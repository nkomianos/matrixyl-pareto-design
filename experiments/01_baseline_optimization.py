#!/usr/bin/env python
"""
Phase 1: Baseline single-objective optimization using NSGA-II.

This experiment:
1. Loads the baseline Matrixyl sequence
2. Initializes oracle and constraints
3. Runs NSGA-II to optimize for penetration and binding
4. Outputs top candidates to CSV

Usage:
    python 01_baseline_optimization.py \
        --sequence ../data/sequences/baseline_matrixyl.fasta \
        --generations 50 \
        --population 30 \
        --output results/phase1
"""

import argparse
import sys
import os
import json
import numpy as np
import pandas as pd
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.oracle import EnsembleOracle
from src.constraints import PhysicochemicalCalculator
from src.evolutionary_algorithm import EvolutionaryOptimizer


def load_sequence(fasta_path: str) -> str:
    """Load sequence from FASTA file."""
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence.upper()


def main(args):
    print("=" * 60)
    print("PHASE 1: Baseline Optimization")
    print("=" * 60)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load baseline sequence
    print(f"\n[1/4] Loading sequence from {args.sequence}...")
    baseline_seq = load_sequence(args.sequence)
    print(f"Baseline sequence: {baseline_seq}")
    print(f"Length: {len(baseline_seq)}")

    # Initialize oracle and constraints
    print(f"\n[2/4] Initializing oracle and constraints...")
    print("(This may take a few minutes on first run to download model weights)")
    try:
        oracle = EnsembleOracle(
            models=['facebook/esm2_t33_650M_UR50D'],
            device='cuda' if args.device == 'cuda' else 'cpu',
            batch_size=args.batch_size
        )
    except Exception as e:
        print(f"WARNING: Could not load ensemble oracle: {e}")
        print("Falling back to dummy oracle for testing...")
        oracle = None

    constraints = PhysicochemicalCalculator()

    # Compute baseline properties
    print(f"\n[3/4] Computing baseline properties...")
    baseline_props = constraints.compute_amino_acid_properties(baseline_seq)
    baseline_penetration = constraints.compute_penetration_score(baseline_seq)

    print(f"  Molecular Weight: {baseline_props['mw']:.1f} Da")
    print(f"  H-bond Donors: {baseline_props['hbd']}")
    print(f"  H-bond Acceptors: {baseline_props['hba']}")
    print(f"  LogP: {baseline_props['logp']:.2f}")
    print(f"  Flexibility: {baseline_props['flexibility']:.2f}")
    print(f"  Penetration Score: {baseline_penetration:.3f}")

    if oracle is not None:
        try:
            binding_scores, binding_unc = oracle.score_binding_affinity([baseline_seq])
            print(f"  Binding Affinity: {binding_scores[0]:.3f} ± {binding_unc[0]:.3f}")
        except Exception as e:
            print(f"  Binding Affinity: (error computing) {e}")

    # Run optimization
    print(f"\n[4/4] Running NSGA-II optimization...")
    print(f"  Population: {args.population}")
    print(f"  Generations: {args.generations}")
    print(f"  Objectives: penetration, binding")

    if oracle is None:
        print("\nSkipping optimization (oracle unavailable). Test completed.")
        return

    optimizer = EvolutionaryOptimizer(oracle, constraints, baseline_seq)

    try:
        pareto_sequences, pareto_objectives = optimizer.optimize_nsga2(
            objectives=['penetration', 'binding'],
            population_size=args.population,
            max_generations=args.generations
        )

        # Compute properties for all Pareto sequences
        results = []
        for i, seq in enumerate(pareto_sequences):
            props = constraints.compute_amino_acid_properties(seq)
            penetration = constraints.compute_penetration_score(seq)

            binding_scores, binding_unc = oracle.score_binding_affinity([seq])

            results.append({
                'rank': i + 1,
                'sequence': seq,
                'length': len(seq),
                'penetration_score': penetration,
                'binding_score': binding_scores[0],
                'binding_uncertainty': binding_unc[0],
                'mw': props['mw'],
                'hbd': props['hbd'],
                'hba': props['hba'],
                'logp': props['logp'],
                'charge': props['charge'],
                'flexibility': props['flexibility'],
            })

        # Save results
        results_df = pd.DataFrame(results)
        csv_path = output_dir / 'pareto_frontier.csv'
        results_df.to_csv(csv_path, index=False)
        print(f"\n✓ Results saved to {csv_path}")

        # Save config
        config = {
            'baseline_sequence': baseline_seq,
            'population_size': args.population,
            'max_generations': args.generations,
            'objectives': ['penetration', 'binding'],
            'n_pareto_solutions': len(pareto_sequences),
        }
        config_path = output_dir / 'config.json'
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)

        print(f"\nSummary:")
        print(f"  Pareto-optimal solutions: {len(pareto_sequences)}")
        print(f"  Best penetration: {results_df['penetration_score'].max():.3f}")
        print(f"  Best binding: {results_df['binding_score'].max():.3f}")
        print(f"\nTop 5 candidates:")
        print(results_df[['sequence', 'penetration_score', 'binding_score']].head(5).to_string(index=False))

    except Exception as e:
        print(f"\nERROR during optimization: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Phase 1: Baseline optimization')
    parser.add_argument('--sequence', default='../data/sequences/baseline_matrixyl.fasta',
                        help='Path to baseline sequence (FASTA)')
    parser.add_argument('--generations', type=int, default=50,
                        help='Number of generations')
    parser.add_argument('--population', type=int, default=30,
                        help='Population size')
    parser.add_argument('--output', default='../results/phase1',
                        help='Output directory')
    parser.add_argument('--device', default='cuda',
                        help='Compute device (cuda/cpu)')
    parser.add_argument('--batch-size', type=int, default=16,
                        help='Batch size for oracle inference')

    args = parser.parse_args()
    main(args)
