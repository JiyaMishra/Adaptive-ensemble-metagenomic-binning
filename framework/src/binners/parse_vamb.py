from pathlib import Path
import pandas as pd


def parse_vamb(vamb_file, output_file):
    vamb_path = Path(vamb_file)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not vamb_path.exists():
        pd.DataFrame(columns=["Contig", "Bin"]).to_csv(output_path, index=False)
        print(f"Saved VAMB assignments (empty): {output_file}")
        return

    df = pd.read_csv(vamb_file, sep="\t")
    print("Columns detected:")
    print(df.columns.tolist())

    contig_col = next((c for c in ["contigname", "contigName", "contig_name", "contig"] if c in df.columns), None)
    cluster_col = next((c for c in ["clustername", "cluster_name", "cluster", "bin"] if c in df.columns), None)

    if contig_col is None or cluster_col is None:
        raise ValueError(f"Unexpected VAMB output format: {list(df.columns)}")

    result = pd.DataFrame()
    result["Contig"] = df[contig_col]
    result["Bin"] = "vamb_bin_" + df[cluster_col].astype(str)

    result.to_csv(output_file, index=False)
    print(f"Saved VAMB assignments: {output_file}")


if __name__ == "__main__":
    PROJECT_ROOT = Path(__file__).resolve().parents[3]

    candidate_files = [
        PROJECT_ROOT / "framework/results/vamb/vae_clusters_unsplit.tsv",
        PROJECT_ROOT / "framework/results/vamb/clusters.tsv",
        PROJECT_ROOT / "framework/results/vamb/vae_clusters.tsv",
    ]
    vamb_file = next((p for p in candidate_files if p.exists()), candidate_files[0])
    output = PROJECT_ROOT / "framework/results/vamb/vamb_assignments.csv"

    parse_vamb(vamb_file, output)