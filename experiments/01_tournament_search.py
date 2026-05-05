#!/usr/bin/env python
"""
Phase 3: deterministic tournament-search baseline.

Example:
    python experiments/01_tournament_search.py \
        --sequence data/sequences/matrixyl_core.fasta \
        --generations 25 \
        --population 30 \
        --output results/phase3_tournament
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

from src.candidates import CandidateEvaluation, CandidateEvaluator
from src.search_space import SearchSpaceRules
from src.tournament_search import TournamentSearchConfig, TournamentSearchOptimizer


def load_sequence(fasta_path: str) -> str:
    return "".join(
        line.strip()
        for line in Path(fasta_path).read_text().splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()


def current_git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def parse_locked_positions(value: str) -> tuple[int, ...]:
    if not value:
        return ()
    return tuple(int(part.strip()) for part in value.split(",") if part.strip())


def flatten_candidate(candidate: CandidateEvaluation, generation: int | None = None) -> dict:
    penetration = candidate.penetration
    descriptors = penetration.descriptors if penetration else None
    row = {
        "generation": generation,
        "sequence": candidate.sequence,
        "modification": candidate.modification,
        "is_valid": candidate.is_valid,
        "failed_filters": ";".join(candidate.failed_filters),
        "optimization_score": candidate.optimization_score,
        "functional_score": candidate.functional_preservation.score,
        "edit_distance": candidate.functional_preservation.edit_distance,
        "penetration_score": penetration.score if penetration else None,
    }
    if descriptors:
        row.update(
            {
                "formula": descriptors.formula,
                "molecular_weight": descriptors.molecular_weight,
                "tpsa": descriptors.tpsa,
                "logp": descriptors.logp,
                "hbd": descriptors.hbd,
                "hba": descriptors.hba,
                "rotatable_bonds": descriptors.rotatable_bonds,
                "formal_charge": descriptors.formal_charge,
            }
        )
    return row


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_convergence_plot(path: Path, summaries: list[dict]) -> bool:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    generations = [summary["generation"] for summary in summaries]
    best_scores = [summary["best_score"] for summary in summaries]
    mean_scores = [summary["mean_score"] for summary in summaries]

    plt.figure(figsize=(7, 4), dpi=150)
    plt.plot(generations, best_scores, label="Best score")
    plt.plot(generations, mean_scores, label="Mean valid score")
    plt.xlabel("Generation")
    plt.ylabel("Optimization score")
    plt.title("Tournament Search Convergence")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    return True


def main(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_sequence = load_sequence(args.sequence)
    search_space = SearchSpaceRules(
        reference_sequence=reference_sequence,
        min_length=len(reference_sequence),
        max_length=len(reference_sequence),
        max_edit_distance=args.max_edit_distance,
        locked_positions=parse_locked_positions(args.locked_positions),
    )
    config = TournamentSearchConfig(
        population_size=args.population,
        generations=args.generations,
        tournament_size=args.tournament_size,
        mutation_rate=args.mutation_rate,
        elite_count=args.elite_count,
        seed=args.seed,
    )
    optimizer = TournamentSearchOptimizer(
        config=config,
        evaluator=CandidateEvaluator(search_space=search_space),
    )

    result = optimizer.run()

    summaries = [summary.to_dict() for summary in result.summaries]
    evaluated_rows = [
        flatten_candidate(candidate, generation=generation_index)
        for generation_index, generation in enumerate(result.evaluations_by_generation)
        for candidate in generation
    ]
    top_rows = [
        flatten_candidate(candidate)
        for candidate in result.final_population[: args.top_k]
    ]
    random_rows = [
        flatten_candidate(candidate)
        for candidate in result.random_baseline[: args.top_k]
    ]

    write_csv(output_dir / "evaluated_candidates.csv", evaluated_rows)
    write_csv(output_dir / "top_candidates.csv", top_rows)
    write_csv(output_dir / "random_baseline_top_candidates.csv", random_rows)
    write_csv(output_dir / "convergence.csv", summaries)

    figure_written = write_convergence_plot(output_dir / "figures" / "convergence.png", summaries)

    metadata = {
        "reference_sequence": reference_sequence,
        "config": result.config.__dict__,
        "search_space": search_space.to_dict() if hasattr(search_space, "to_dict") else {
            "reference_sequence": search_space.reference_sequence,
            "min_length": search_space.min_length,
            "max_length": search_space.max_length,
            "max_edit_distance": search_space.max_edit_distance,
            "locked_positions": search_space.locked_positions,
        },
        "git_commit": current_git_commit(),
        "python": sys.version,
        "platform": platform.platform(),
        "figure_written": figure_written,
        "best": result.best.to_dict(),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))
    (output_dir / "config.json").write_text(json.dumps(result.config.__dict__, indent=2))

    print("Tournament search complete")
    print(f"Best sequence: {result.best.sequence}")
    print(f"Best optimization score: {result.best.optimization_score:.4f}")
    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 3 tournament-search baseline")
    parser.add_argument("--sequence", default="data/sequences/matrixyl_core.fasta")
    parser.add_argument("--output", default="results/phase3_tournament")
    parser.add_argument("--population", type=int, default=30)
    parser.add_argument("--generations", type=int, default=25)
    parser.add_argument("--tournament-size", type=int, default=3)
    parser.add_argument("--mutation-rate", type=float, default=0.2)
    parser.add_argument("--elite-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-edit-distance", type=int, default=2)
    parser.add_argument("--locked-positions", default="")
    parser.add_argument("--top-k", type=int, default=10)
    main(parser.parse_args())
