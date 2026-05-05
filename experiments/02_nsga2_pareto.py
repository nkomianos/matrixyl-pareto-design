#!/usr/bin/env python
"""
Phase 4: Pareto optimization for Matrixyl-family peptides.

Example:
    python experiments/02_nsga2_pareto.py \
        --sequence data/sequences/matrixyl_core.fasta \
        --population 50 \
        --generations 50 \
        --output results/phase4_pareto
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

from src.candidates import CandidateEvaluator
from src.pareto_search import ParetoCandidate, ParetoSearchConfig, ParetoSearchOptimizer
from src.search_space import SearchSpaceRules


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


def flatten_candidate(candidate: ParetoCandidate, generation: int | None = None) -> dict:
    evaluation = candidate.evaluation
    penetration = evaluation.penetration
    descriptors = penetration.descriptors if penetration else None
    row = {
        "generation": generation,
        "sequence": evaluation.sequence,
        "is_valid": evaluation.is_valid,
        "failed_filters": ";".join(evaluation.failed_filters),
        "pareto_rank": candidate.rank,
        "crowding_distance": candidate.crowding_distance,
        "penetration_objective": candidate.objectives.penetration,
        "functional_objective": candidate.objectives.functional_preservation,
        "synthesis_objective": candidate.objectives.synthesis_feasibility,
        "scalar_optimization_score": evaluation.optimization_score,
        "edit_distance": evaluation.functional_preservation.edit_distance,
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


def write_frontier_plot(path: Path, frontier_rows: list[dict], population_rows: list[dict]) -> bool:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(6, 5), dpi=150)
    plt.scatter(
        [float(row["penetration_objective"]) for row in population_rows],
        [float(row["functional_objective"]) for row in population_rows],
        s=24,
        alpha=0.35,
        label="Final population",
    )
    plt.scatter(
        [float(row["penetration_objective"]) for row in frontier_rows],
        [float(row["functional_objective"]) for row in frontier_rows],
        s=44,
        label="Pareto frontier",
    )
    plt.xlabel("Penetration objective")
    plt.ylabel("Functional-preservation objective")
    plt.title("Pareto Frontier")
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
    config = ParetoSearchConfig(
        population_size=args.population,
        generations=args.generations,
        tournament_size=args.tournament_size,
        mutation_rate=args.mutation_rate,
        seed=args.seed,
    )
    optimizer = ParetoSearchOptimizer(
        config=config,
        evaluator=CandidateEvaluator(search_space=search_space),
    )
    result = optimizer.run()

    frontier_rows = [flatten_candidate(candidate) for candidate in result.frontier]
    final_population_rows = [flatten_candidate(candidate) for candidate in result.final_population]
    evaluated_rows = [
        flatten_candidate(candidate, generation=generation_index)
        for generation_index, generation in enumerate(result.evaluations_by_generation)
        for candidate in generation
    ]
    summaries = [summary.to_dict() for summary in result.summaries]

    write_csv(output_dir / "pareto_frontier.csv", frontier_rows)
    write_csv(output_dir / "final_population.csv", final_population_rows)
    write_csv(output_dir / "evaluated_candidates.csv", evaluated_rows)
    write_csv(output_dir / "convergence.csv", summaries)
    figure_written = write_frontier_plot(
        output_dir / "figures" / "pareto_frontier.png",
        frontier_rows,
        final_population_rows,
    )

    metadata = {
        "reference_sequence": reference_sequence,
        "config": config.__dict__,
        "search_space": {
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
        "frontier_size": len(result.frontier),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))
    (output_dir / "config.json").write_text(json.dumps(config.__dict__, indent=2))

    print("Pareto optimization complete")
    print(f"Frontier size: {len(result.frontier)}")
    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 4 Pareto optimization")
    parser.add_argument("--sequence", default="data/sequences/matrixyl_core.fasta")
    parser.add_argument("--output", default="results/phase4_pareto")
    parser.add_argument("--population", type=int, default=50)
    parser.add_argument("--generations", type=int, default=50)
    parser.add_argument("--tournament-size", type=int, default=2)
    parser.add_argument("--mutation-rate", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-edit-distance", type=int, default=2)
    parser.add_argument("--locked-positions", default="")
    main(parser.parse_args())
