#!/usr/bin/env python
"""
Phase 5: candidate analysis and interpretability.

Example:
    python experiments/03_candidate_analysis.py \
        --frontier results/phase4_pareto_full/pareto_frontier.csv \
        --output results/phase5_candidate_analysis
"""

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.analysis import (
    baseline_descriptor_rows,
    load_candidate_rows,
    mutation_enrichment,
    position_frequency_matrix,
    summarize_candidates,
    write_dict_rows,
)


def current_git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def write_position_heatmap(path: Path, frequency_rows: list[dict]) -> bool:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    positions = sorted({row["position"] for row in frequency_rows})
    residues = sorted({row["residue"] for row in frequency_rows})
    matrix = [
        [
            next(
                (
                    row["frequency"]
                    for row in frequency_rows
                    if row["position"] == position and row["residue"] == residue
                ),
                0.0,
            )
            for position in positions
        ]
        for residue in residues
    ]

    plt.figure(figsize=(7, max(3, len(residues) * 0.35)), dpi=150)
    plt.imshow(matrix, aspect="auto", vmin=0, vmax=1)
    plt.colorbar(label="Frequency")
    plt.xticks(range(len(positions)), positions)
    plt.yticks(range(len(residues)), residues)
    plt.xlabel("Position")
    plt.ylabel("Residue")
    plt.title("Pareto Frontier Position Frequencies")
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    return True


def write_tradeoff_plot(path: Path, summary_rows: list[dict]) -> bool:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(6, 5), dpi=150)
    for row in summary_rows:
        plt.scatter(row["penetration_objective"], row["functional_objective"], s=48)
        plt.text(
            row["penetration_objective"] + 0.003,
            row["functional_objective"],
            row["sequence"],
            fontsize=8,
        )
    plt.xlabel("Penetration objective")
    plt.ylabel("Functional-preservation objective")
    plt.title("Pareto Candidate Trade-Offs")
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    return True


def main(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    candidates = load_candidate_rows(args.frontier)
    candidate_summary = summarize_candidates(candidates)
    frequency_rows = position_frequency_matrix([row.sequence for row in candidates])
    mutation_rows = mutation_enrichment([row.sequence for row in candidates])
    baseline_rows = baseline_descriptor_rows(args.palmitoylated_smiles)

    write_dict_rows(output_dir / "candidate_summary.csv", candidate_summary)
    write_dict_rows(output_dir / "position_frequency_matrix.csv", frequency_rows)
    write_dict_rows(output_dir / "mutation_enrichment.csv", mutation_rows)
    write_dict_rows(output_dir / "baseline_comparison.csv", baseline_rows)

    heatmap_written = write_position_heatmap(
        output_dir / "figures" / "position_frequency_heatmap.png",
        frequency_rows,
    )
    tradeoff_written = write_tradeoff_plot(
        output_dir / "figures" / "candidate_tradeoffs.png",
        candidate_summary,
    )

    metadata = {
        "frontier": args.frontier,
        "palmitoylated_smiles": args.palmitoylated_smiles,
        "candidate_count": len(candidates),
        "mutation_count": len(mutation_rows),
        "baseline_count": len(baseline_rows),
        "heatmap_written": heatmap_written,
        "tradeoff_written": tradeoff_written,
        "git_commit": current_git_commit(),
        "python": sys.version,
        "platform": platform.platform(),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))

    print("Candidate analysis complete")
    print(f"Candidates analyzed: {len(candidates)}")
    print(f"Mutation rows: {len(mutation_rows)}")
    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 5 candidate analysis")
    parser.add_argument("--frontier", default="results/phase4_pareto_full/pareto_frontier.csv")
    parser.add_argument("--output", default="results/phase5_candidate_analysis")
    parser.add_argument(
        "--palmitoylated-smiles",
        default="data/molecules/matrixyl_palmitoylated.smi",
    )
    main(parser.parse_args())
