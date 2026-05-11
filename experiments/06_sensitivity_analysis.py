#!/usr/bin/env python3
"""
Sensitivity Analysis: Vary TPSA/MW/LogP penalty weights and re-rank candidates.
Demonstrates robustness of Pareto frontier to parameter variations.
"""

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from datetime import datetime

from src.candidates import CandidateEvaluator, SearchSpaceRules
from src.chemistry import descriptors_from_sequence


def run_sensitivity_analysis(
    frontier_csv: Path,
    output_dir: Path,
    weight_variations: dict,
    seed: int = 42,
) -> dict:
    """
    Load frontier candidates and re-evaluate with varied penalty weights.

    Args:
        frontier_csv: Path to pareto_frontier.csv
        output_dir: Output directory
        weight_variations: Dict of parameter name -> [list of multipliers]
        seed: Random seed

    Returns:
        Sensitivity analysis results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load frontier
    frontier_df = pd.read_csv(frontier_csv)
    sequences = frontier_df['sequence'].tolist()

    # Baseline evaluator (default weights, default Matrixyl search space)
    baseline_eval = CandidateEvaluator()

    # Store results for each parameter sweep
    sensitivity_results = {}

    # Weight ranges to test (±30% from default assumed weights)
    weight_scenarios = {
        "tpsa_penalty_baseline": 1.0,
        "tpsa_penalty_minus30": 0.7,
        "tpsa_penalty_plus30": 1.3,
        "mw_penalty_baseline": 1.0,
        "mw_penalty_minus30": 0.7,
        "mw_penalty_plus30": 1.3,
        "logp_penalty_baseline": 1.0,
        "logp_penalty_minus30": 0.7,
        "logp_penalty_plus30": 1.3,
    }

    # For each scenario, create a modified evaluator by re-scoring
    all_scenario_results = []

    for scenario_name, scenario_config in [
        ("baseline", {"tpsa_mult": 1.0, "mw_mult": 1.0, "logp_mult": 1.0}),
        ("tpsa_minus30", {"tpsa_mult": 0.7, "mw_mult": 1.0, "logp_mult": 1.0}),
        ("tpsa_plus30", {"tpsa_mult": 1.3, "mw_mult": 1.0, "logp_mult": 1.0}),
        ("mw_minus30", {"tpsa_mult": 1.0, "mw_mult": 0.7, "logp_mult": 1.0}),
        ("mw_plus30", {"tpsa_mult": 1.0, "mw_mult": 1.3, "logp_mult": 1.0}),
        ("logp_minus30", {"tpsa_mult": 1.0, "mw_mult": 1.0, "logp_mult": 0.7}),
        ("logp_plus30", {"tpsa_mult": 1.0, "mw_mult": 1.0, "logp_mult": 1.3}),
    ]:
        scenario_results = []

        for seq in sequences:
            result = baseline_eval.evaluate_sequence(seq)
            if result.is_valid and result.penetration is not None:
                desc = result.penetration.descriptors

                # Re-score with modified weights by perturbing the per-descriptor
                # penalty contributions and recombining multiplicatively.
                penalties = result.penetration.penalties
                tpsa_pen = penalties.get("tpsa", 0.0) * scenario_config["tpsa_mult"]
                mw_pen = penalties.get("molecular_weight", 0.0) * scenario_config["mw_mult"]
                logp_pen = penalties.get("logp", 0.0) * scenario_config["logp_mult"]
                other_pen = sum(
                    v for k, v in penalties.items()
                    if k not in {"tpsa", "molecular_weight", "logp"}
                )
                # Clamp each perturbed penalty to [0, 1] before recombining
                tpsa_pen = min(1.0, max(0.0, tpsa_pen))
                mw_pen = min(1.0, max(0.0, mw_pen))
                logp_pen = min(1.0, max(0.0, logp_pen))
                # Recombine as the original penetration score does: product of (1 - penalty)
                adjusted_penetration = max(
                    0.01,
                    (1.0 - tpsa_pen) * (1.0 - mw_pen) * (1.0 - logp_pen) * (1.0 - other_pen)
                )

                scenario_results.append({
                    "sequence": seq,
                    "scenario": scenario_name,
                    "penetration_score_original": result.penetration.score,
                    "penetration_score_adjusted": adjusted_penetration,
                    "functional_score": result.functional_preservation.score,
                    "molecular_weight": desc.molecular_weight,
                    "tpsa": desc.tpsa,
                    "logp": desc.logp,
                })

        all_scenario_results.extend(scenario_results)

    # Convert to DataFrame
    sensitivity_df = pd.DataFrame(all_scenario_results)

    # Save sensitivity results
    sensitivity_csv = output_dir / "sensitivity_analysis.csv"
    sensitivity_df.to_csv(sensitivity_csv, index=False)

    # Analyze rank stability
    rank_stability = []
    baseline_ranking = sensitivity_df[sensitivity_df['scenario'] == 'baseline'].sort_values(
        'penetration_score_adjusted', ascending=False
    )['sequence'].reset_index(drop=True).to_dict()
    baseline_ranking = {seq: rank for rank, seq in baseline_ranking.items()}

    for scenario in sensitivity_df['scenario'].unique():
        if scenario == "baseline":
            continue
        scenario_df = sensitivity_df[sensitivity_df['scenario'] == scenario].sort_values(
            'penetration_score_adjusted', ascending=False
        )
        scenario_ranking = scenario_df['sequence'].reset_index(drop=True).to_dict()
        scenario_ranking = {seq: rank for rank, seq in scenario_ranking.items()}

        rank_changes = []
        for seq in baseline_ranking:
            baseline_rank = baseline_ranking.get(seq)
            scenario_rank = scenario_ranking.get(seq)
            if baseline_rank is not None and scenario_rank is not None:
                rank_changes.append({
                    "sequence": seq,
                    "baseline_rank": baseline_rank,
                    "scenario_rank": scenario_rank,
                    "rank_change": abs(baseline_rank - scenario_rank),
                })

        rank_stability.append({
            "scenario": scenario,
            "mean_rank_change": np.mean([r['rank_change'] for r in rank_changes]) if rank_changes else 0,
            "max_rank_change": max([r['rank_change'] for r in rank_changes]) if rank_changes else 0,
            "top_candidate_unchanged": baseline_ranking[sequences[0]] == scenario_ranking.get(sequences[0], -1),
        })

    rank_stability_df = pd.DataFrame(rank_stability)
    rank_stability_csv = output_dir / "rank_stability.csv"
    rank_stability_df.to_csv(rank_stability_csv, index=False)

    # Summary statistics
    summary = {
        "num_candidates_analyzed": len(sequences),
        "num_scenarios": len(sensitivity_df['scenario'].unique()),
        "mean_rank_change_across_scenarios": float(rank_stability_df['mean_rank_change'].mean()),
        "max_rank_change_observed": float(rank_stability_df['max_rank_change'].max()),
        "top_candidate_stability": float(rank_stability_df['top_candidate_unchanged'].sum() / len(rank_stability_df)),
        "conclusion": (
            "Frontier ranking is robust to ±30% parameter variations. "
            "Core candidates maintain top positions across all scenarios."
            if rank_stability_df['max_rank_change'].max() <= 3
            else "Frontier ranking shows sensitivity to parameter choices."
        ),
    }

    # Write summary
    metadata = {
        "timestamp": datetime.now().isoformat(),
        "frontier_input": str(frontier_csv),
        "seed": seed,
        "summary": summary,
    }

    metadata_json = output_dir / "run_metadata.json"
    with open(metadata_json, "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"Sensitivity analysis complete")
    print(f"Results: {sensitivity_csv}")
    print(f"Rank stability: {rank_stability_csv}")
    print(f"Summary: {metadata_json}")
    print(f"\nKey finding: {summary['conclusion']}")

    return summary


def main(args):
    summary = run_sensitivity_analysis(
        frontier_csv=Path(args.frontier),
        output_dir=Path(args.output),
        weight_variations={},
        seed=args.seed,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sensitivity analysis for Pareto frontier",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--frontier",
        type=str,
        default="results/phase2_pareto/pareto_frontier.csv",
        help="Path to pareto_frontier.csv",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results/sensitivity_analysis",
        help="Output directory",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed",
    )

    main(parser.parse_args())
