"""Hardware-aware execution routing for metagenomic contig processing.

The module converts available host capacity and per-contig sequence complexity
into an explicit execution plan.  Its technical effect is to reserve expensive
matrix/explainability work for ambiguous contigs while preserving a complete,
auditable row for every contig in the canonical ensemble dataset.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import math
import os
from pathlib import Path
import platform
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import psutil

from readfasta import FASTAReader


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ENSEMBLE_INPUT = PROJECT_ROOT / "framework/results/ensemble_dataframe.csv"
ROUTES_OUTPUT = PROJECT_ROOT / "framework/results/execution_routes.csv"
HEAVY_ENSEMBLE_OUTPUT = PROJECT_ROOT / "framework/results/heavy_ensemble_dataframe.csv"
FASTA_INPUT = PROJECT_ROOT / "data/assemblies/metagem_1500/final.contigs.fa"
BIN_COLUMNS = ("metabat_bin", "maxbin_bin", "vamb_bin")
WINDOW_SIZE = 100
KMER_SIZE = 4


@dataclass(frozen=True)
class HardwareCapacity:
    """Snapshot of capacity used to size memory-sensitive work batches.

    `chunk_size` implements the requested bounded allocation rule.  Using
    currently *available* memory, rather than installed memory, prevents a WSL
    guest or a shared workstation from overcommitting RAM during large feature
    matrix and KernelSHAP operations.
    """

    available_ram_gb: float
    total_ram_gb: float
    swap_used_gb: float
    swap_total_gb: float
    cpu_cores: int
    is_wsl: bool
    chunk_size: int


def get_hardware_capacity() -> HardwareCapacity:
    """Read real-time host/guest resources and derive a safe processing chunk.

    The allocation is ``clamp(floor(available_ram_gb / 2) * 250, 50, 2000)``.
    This reserves approximately half of presently available RAM for the OS,
    binners, and file cache, which is particularly important where WSL reports
    a virtualized memory ceiling rather than the host's installed memory.
    """

    memory = psutil.virtual_memory()
    swap = psutil.swap_memory()
    gibibyte = 1024**3
    available_ram_gb = memory.available / gibibyte
    release = platform.release().lower()
    is_wsl = "microsoft" in release or "wsl" in release
    raw_chunk_size = math.floor(available_ram_gb / 2) * 250
    chunk_size = int(np.clip(raw_chunk_size, 50, 2000))

    return HardwareCapacity(
        available_ram_gb=round(available_ram_gb, 3),
        total_ram_gb=round(memory.total / gibibyte, 3),
        swap_used_gb=round(swap.used / gibibyte, 3),
        swap_total_gb=round(swap.total / gibibyte, 3),
        cpu_cores=os.cpu_count() or 1,
        is_wsl=is_wsl,
        chunk_size=chunk_size,
    )


def _shannon_entropy_4mer(sequence: str) -> float:
    """Return normalized Shannon entropy of valid 4-mer frequencies."""

    kmers = [
        sequence[index : index + KMER_SIZE]
        for index in range(len(sequence) - KMER_SIZE + 1)
        if set(sequence[index : index + KMER_SIZE]) <= {"A", "C", "G", "T"}
    ]
    if not kmers:
        return 0.0

    counts = Counter(kmers)
    probabilities = np.fromiter(
        (count / len(kmers) for count in counts.values()),
        dtype=float,
    )
    entropy = -np.sum(probabilities * np.log2(probabilities))
    return float(entropy / math.log2(4**KMER_SIZE))


def sequence_complexity(sequence: str, window_size: int = WINDOW_SIZE) -> tuple[float, float, float]:
    """Calculate GC-window variance, 4-mer entropy, and their composite score.

    GC variance exposes local compositional shifts that make bin membership less
    certain.  Normalized 4-mer entropy captures oligonucleotide diversity.
    Summing the two creates a bounded, sequence-derived complexity signal with
    no hand-tuned absolute biological cutoff.
    """

    normalized = sequence.upper()
    windows = [
        normalized[index : index + window_size]
        for index in range(0, len(normalized), window_size)
        if normalized[index : index + window_size]
    ]
    if not windows:
        return float("nan"), float("nan"), float("nan")

    gc_fractions = np.asarray(
        [
            (window.count("G") + window.count("C")) / len(window)
            for window in windows
        ],
        dtype=float,
    )
    gc_variance = float(np.var(gc_fractions))
    kmer_entropy = _shannon_entropy_4mer(normalized)
    return gc_variance, kmer_entropy, float(gc_variance + kmer_entropy)


def _assignment_metrics(row: pd.Series) -> tuple[float, int, str | None]:
    """Measure base-binner disagreement and identify a valid consensus bin."""

    assignments = [str(row[column]) for column in BIN_COLUMNS if pd.notna(row[column])]
    if not assignments:
        return 1.0, 0, None

    counts = Counter(assignments)
    consensus_bin, consensus_count = counts.most_common(1)[0]
    ambiguity = 1.0 - (consensus_count / len(assignments))
    return float(ambiguity), consensus_count, consensus_bin


def route_contig_execution(
    ensemble: pd.DataFrame,
    sequences: Mapping[str, str],
) -> tuple[pd.DataFrame, float | None]:
    """Tag contigs as direct consensus or heavy adaptive-ensemble work.

    The low-complexity boundary is the observed first quartile of the current
    dataset's composite complexity scores.  This empirical calibration avoids a
    fixed biological threshold and adapts to read length and assembly ecology.
    A fast path requires at least two base-binner assignments: unanimous
    available assignments route immediately, while low-complexity majority
    consensus also bypasses high-cost work.  Missing sequence or assignment
    evidence is conservatively routed to ``HEAVY_ENSEMBLE``.
    """

    missing_columns = set(BIN_COLUMNS) - set(ensemble.columns)
    if missing_columns:
        raise ValueError(f"Missing binner columns for routing: {sorted(missing_columns)}")

    routed = ensemble.copy()
    metrics = [sequence_complexity(sequences[str(contig_id)]) if str(contig_id) in sequences else (np.nan, np.nan, np.nan)
               for contig_id in routed["contig_id"]]
    routed[["gc_window_variance", "kmer_entropy_4", "complexity_score"]] = pd.DataFrame(
        metrics,
        index=routed.index,
    )

    finite_complexity = routed["complexity_score"].dropna()
    complexity_low_threshold = (
        float(finite_complexity.quantile(0.25)) if not finite_complexity.empty else None
    )
    assignment_metrics = routed.apply(_assignment_metrics, axis=1)
    routed[["tool_ambiguity", "tool_consensus_count", "consensus_bin"]] = pd.DataFrame(
        assignment_metrics.tolist(),
        index=routed.index,
    )
    routed["tool_consensus_count"] = routed["tool_consensus_count"].astype(int)

    unanimous = (routed["tool_ambiguity"] == 0.0) & (routed["tool_consensus_count"] >= 2)
    low_complexity_consensus = (
        routed["complexity_score"].lt(complexity_low_threshold)
        if complexity_low_threshold is not None
        else pd.Series(False, index=routed.index)
    ) & (routed["tool_consensus_count"] >= 2)
    routed["execution_route"] = np.where(
        unanimous | low_complexity_consensus,
        "FAST_PATH",
        "HEAVY_ENSEMBLE",
    )
    routed["routing_reason"] = np.select(
        [unanimous, low_complexity_consensus],
        ["unanimous_base_binner_consensus", "low_complexity_majority_consensus"],
        default="ambiguous_or_high_complexity",
    )
    return routed, complexity_low_threshold


def _load_sequences(fasta_path: Path) -> dict[str, str]:
    """Load sequences once so routing does not repeatedly scan the assembly."""

    return {
        str(contig["id"]): str(contig["sequence"])
        for contig in FASTAReader(fasta_path).read_contigs()
    }


def prepare_execution_routes(
    ensemble_path: Path = ENSEMBLE_INPUT,
    fasta_path: Path = FASTA_INPUT,
) -> tuple[HardwareCapacity, pd.DataFrame]:
    """Persist routed canonical and heavy-only datasets for pipeline consumers.

    The canonical CSV retains every contig plus route metadata for downstream
    joins.  The heavy-only CSV is a compact workset intended for operations such
    as SHAP and high-dimensional matrix analysis; excluding fast-path rows from
    it reduces RAM pressure without losing traceability in final CSV products.
    """

    if not ensemble_path.exists():
        raise FileNotFoundError(f"Ensemble input not found: {ensemble_path}")
    if not fasta_path.exists():
        raise FileNotFoundError(f"Assembly FASTA not found: {fasta_path}")

    capacity = get_hardware_capacity()
    ensemble = pd.read_csv(ensemble_path)
    if "contig_id" not in ensemble.columns:
        raise ValueError("Ensemble input must contain a contig_id column.")

    routed, complexity_low_threshold = route_contig_execution(
        ensemble,
        _load_sequences(fasta_path),
    )
    ensemble_path.parent.mkdir(parents=True, exist_ok=True)
    routed.to_csv(ensemble_path, index=False)
    routed.loc[routed["execution_route"] == "HEAVY_ENSEMBLE"].to_csv(
        HEAVY_ENSEMBLE_OUTPUT,
        index=False,
    )
    route_columns: Sequence[str] = (
        "contig_id", "execution_route", "routing_reason", "consensus_bin",
        "tool_ambiguity", "tool_consensus_count", "gc_window_variance",
        "kmer_entropy_4", "complexity_score",
    )
    route_report = routed.loc[:, route_columns].copy()
    route_report["complexity_low_threshold"] = complexity_low_threshold
    route_report["hardware_chunk_size"] = capacity.chunk_size
    route_report.to_csv(ROUTES_OUTPUT, index=False)

    return capacity, routed


def format_startup_telemetry(capacity: HardwareCapacity, routed: pd.DataFrame) -> str:
    """Format a concise, operationally useful routing and capacity report."""

    fast_count = int((routed["execution_route"] == "FAST_PATH").sum())
    heavy_count = int((routed["execution_route"] == "HEAVY_ENSEMBLE").sum())
    environment = "WSL" if capacity.is_wsl else "native"
    return (
        "[RESOURCE-ENGINE] "
        f"RAM available: {capacity.available_ram_gb:.2f} GB / {capacity.total_ram_gb:.2f} GB | "
        f"CPU cores: {capacity.cpu_cores} | Swap used: {capacity.swap_used_gb:.2f} GB | "
        f"Environment: {environment}\n"
        f"[RESOURCE-ENGINE] Routed {fast_count} contigs to FAST_PATH, "
        f"{heavy_count} to HEAVY_ENSEMBLE | Batch Size: {capacity.chunk_size}"
    )


def hardware_metadata(capacity: HardwareCapacity) -> dict[str, float | int | bool]:
    """Expose the immutable capacity snapshot for structured telemetry clients."""

    return asdict(capacity)
