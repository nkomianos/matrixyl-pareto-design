#!/usr/bin/env python
"""
Phase 6: exploratory conformer validation.

Example:
    python experiments/04_structural_validation.py \
        --candidates results/phase5_candidate_analysis/candidate_summary.csv \
        --output results/phase6_structural_validation
"""

import argparse
import csv
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.analysis import load_candidate_rows, write_dict_rows
from src.chemistry import read_smiles
from src.structure import summarize_sequence_conformers, summarize_smiles_conformers


def current_git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def load_priority_labels(path: str | Path) -> dict[str, str]:
    with Path(path).open() as handle:
        return {
            row["sequence"]: row.get("priority_label", "")
            for row in csv.DictReader(handle)
        }


def summary_to_row(summary, *, sequence: str, priority_label: str = "", penetration_score=None) -> dict:
    return {
        "name": summary.name,
        "sequence": sequence,
        "source": summary.source,
        "priority_label": priority_label,
        "penetration_score": penetration_score,
        "conformer_count": summary.conformer_count,
        "min_radius_of_gyration": summary.min_radius_of_gyration,
        "mean_radius_of_gyration": summary.mean_radius_of_gyration,
        "max_radius_of_gyration": summary.max_radius_of_gyration,
        "std_radius_of_gyration": summary.std_radius_of_gyration,
        "min_energy": summary.min_energy,
        "mean_energy": summary.mean_energy,
    }


def write_compactness_plot(path: Path, rows: list[dict]) -> bool:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    plot_rows = [row for row in rows if row["penetration_score"] not in (None, "")]
    plt.figure(figsize=(7, 5), dpi=150)
    for row in plot_rows:
        plt.scatter(row["mean_radius_of_gyration"], row["penetration_score"], s=48)
        plt.text(
            row["mean_radius_of_gyration"] + 0.02,
            row["penetration_score"],
            row["sequence"],
            fontsize=8,
        )
    plt.xlabel("Mean radius of gyration")
    plt.ylabel("Penetration score")
    plt.title("Conformer Compactness vs Penetration Score")
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    return True


def main(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    candidate_rows = load_candidate_rows(args.candidates)
    priority_labels = load_priority_labels(args.candidates)
    selected = candidate_rows[: args.top_k]

    summary_rows = []
    conformer_rows = []
    for candidate in selected:
        summary = summarize_sequence_conformers(
            candidate.sequence,
            name=candidate.sequence,
            num_conformers=args.num_conformers,
            seed=args.seed,
        )
        summary_rows.append(
            summary_to_row(
                summary,
                sequence=candidate.sequence,
                priority_label=priority_labels.get(candidate.sequence, ""),
                penetration_score=candidate.penetration_objective,
            )
        )
        for metric in summary.metrics:
            conformer_rows.append(
                {
                    "name": summary.name,
                    "sequence": candidate.sequence,
                    **metric.to_dict(),
                }
            )

    if "KTTKS" not in {row.sequence for row in selected}:
        core = summarize_sequence_conformers(
            "KTTKS",
            name="Matrixyl core",
            num_conformers=args.num_conformers,
            seed=args.seed,
        )
        summary_rows.append(summary_to_row(core, sequence="KTTKS", priority_label="function_anchor"))

    smiles, name = read_smiles(args.palmitoylated_smiles)
    pal = summarize_smiles_conformers(
        smiles,
        name=name,
        source="pubchem_smiles",
        num_conformers=max(2, min(args.num_conformers, args.palmitoylated_conformers)),
        seed=args.seed,
    )
    summary_rows.append(summary_to_row(pal, sequence="KTTKS", priority_label="lipidated_baseline"))

    write_dict_rows(output_dir / "conformer_summary.csv", summary_rows)
    write_dict_rows(output_dir / "conformer_metrics.csv", conformer_rows)
    figure_written = write_compactness_plot(
        output_dir / "figures" / "compactness_vs_penetration.png",
        summary_rows,
    )

    metadata = {
        "candidates": args.candidates,
        "palmitoylated_smiles": args.palmitoylated_smiles,
        "top_k": args.top_k,
        "num_conformers": args.num_conformers,
        "palmitoylated_conformers": args.palmitoylated_conformers,
        "seed": args.seed,
        "summary_count": len(summary_rows),
        "conformer_metric_count": len(conformer_rows),
        "figure_written": figure_written,
        "git_commit": current_git_commit(),
        "python": sys.version,
        "platform": platform.platform(),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))

    print("Structural validation complete")
    print(f"Summaries written: {len(summary_rows)}")
    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 6 structural validation")
    parser.add_argument("--candidates", default="results/phase5_candidate_analysis/candidate_summary.csv")
    parser.add_argument("--output", default="results/phase6_structural_validation")
    parser.add_argument("--palmitoylated-smiles", default="data/molecules/matrixyl_palmitoylated.smi")
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--num-conformers", type=int, default=10)
    parser.add_argument("--palmitoylated-conformers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    main(parser.parse_args())
