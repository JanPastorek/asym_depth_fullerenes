import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-glob", default="results/metrics_*.csv")
    parser.add_argument("--output-dir", default="results")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = sorted(Path().glob(args.input_glob))
    if not paths:
        raise SystemExit("No input files found.")

    frames = [pd.read_csv(path) for path in paths]
    data = pd.concat(frames, ignore_index=True)
    data = data.dropna(subset=["hex_neighbor_index", "asymmetric_depth"])

    pearson = data["hex_neighbor_index"].corr(
        data["asymmetric_depth"], method="pearson")
    spearman = data["hex_neighbor_index"].corr(
        data["asymmetric_depth"], method="spearman")

    summary = {
        "n": int(len(data)),
        "pearson": float(pearson),
        "spearman": float(spearman),
    }

    combined_path = output_dir / "metrics_all.csv"
    summary_path = output_dir / "correlation_summary.json"

    data.to_csv(combined_path, index=False)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
