#!/usr/bin/env python3
"""Generate publication-ready figures for the paper.

Reads phase-1 through phase-3 result CSVs from ``results/`` and writes the
five manuscript figures directly into ``manuscript/figures/`` so the LaTeX
build always sees the same artifacts the script produced. All paths are
anchored to the repository root, so the script can be invoked from any CWD.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

sns.set_theme(style="whitegrid", palette="husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
OUTPUT_DIR = REPO_ROOT / "manuscript" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def fig_pareto_frontier():
    """Figure 1: Pareto frontier - penetration vs functional preservation."""
    df = pd.read_csv(RESULTS_DIR / "phase2_pareto" / "pareto_frontier.csv")

    fig, ax = plt.subplots(figsize=(8, 6))

    # All evaluated candidates in background
    all_df = pd.read_csv(RESULTS_DIR / "phase2_pareto" / "evaluated_candidates.csv")
    ax.scatter(all_df['penetration_objective'], all_df['functional_objective'],
              alpha=0.15, s=20, c='gray', label='Evaluated candidates')

    # Pareto frontier in foreground
    scatter = ax.scatter(df['penetration_objective'], df['functional_objective'],
                        s=150, c=df['edit_distance'], cmap='viridis',
                        edgecolors='black', linewidth=1.5, alpha=0.8,
                        label='Pareto frontier')

    # Annotate key candidates
    for idx, row in df.iterrows():
        if row['sequence'] in ['PTTPS', 'KTTKS', 'KTTPS', 'KTTPP']:
            ax.annotate(row['sequence'],
                       xy=(row['penetration_objective'], row['functional_objective']),
                       xytext=(5, 5), textcoords='offset points', fontsize=9,
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    ax.set_xlabel('Penetration Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Functional Preservation Score', fontsize=12, fontweight='bold')
    ax.set_title('Pareto Frontier: Multi-Objective Optimization', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)

    cbar = plt.colorbar(scatter, ax=ax, label='Edit Distance from KTTKS')
    ax.legend(loc='lower left', fontsize=10)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "01_pareto_frontier.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Figure 1: Pareto frontier")


def fig_convergence_comparison():
    """Figure 2: Convergence curves - Tournament GA vs NSGA-II."""
    conv_ga = pd.read_csv(RESULTS_DIR / "phase1_tournament" / "convergence.csv")
    conv_pareto = pd.read_csv(RESULTS_DIR / "phase2_pareto" / "convergence.csv")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Tournament GA convergence
    ax1.plot(conv_ga['generation'], conv_ga['best_score'],
            marker='o', linewidth=2, markersize=4, label='Best fitness')
    ax1.fill_between(conv_ga['generation'],
                    conv_ga['best_score'], conv_ga['mean_score'],
                    alpha=0.3, label='Population range')
    ax1.set_xlabel('Generation', fontsize=11)
    ax1.set_ylabel('Optimization Score', fontsize=11)
    ax1.set_title('Tournament Search (Single-Objective)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # NSGA-II convergence
    if 'frontier_size' in conv_pareto.columns:
        ax2.plot(conv_pareto['generation'], conv_pareto['frontier_size'],
                marker='s', linewidth=2, markersize=4, label='Frontier size')
        ax2.set_ylabel('Pareto Frontier Size', fontsize=11)
    else:
        # Fallback: plot best penetration over generations
        all_df = pd.read_csv(RESULTS_DIR / "phase2_pareto" / "evaluated_candidates.csv")
        best_by_gen = all_df.groupby('generation')['penetration_objective'].max()
        ax2.plot(best_by_gen.index, best_by_gen.values,
                marker='s', linewidth=2, markersize=4, label='Best penetration')
        ax2.set_ylabel('Best Penetration Score', fontsize=11)

    ax2.set_xlabel('Generation', fontsize=11)
    ax2.set_title('NSGA-II Multi-Objective Optimization', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "02_convergence.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Figure 2: Convergence curves")


def fig_mutation_analysis():
    """Figure 3: Mutation enrichment heatmap."""
    df = pd.read_csv(RESULTS_DIR / "phase3_analysis" / "mutation_enrichment.csv")
    frontier = pd.read_csv(RESULTS_DIR / "phase2_pareto" / "pareto_frontier.csv")

    fig, ax = plt.subplots(figsize=(10, 6))

    # Create position × mutation heatmap
    positions = df['position'].unique()
    mutations = df['mutation'].unique()

    heatmap_data = np.zeros((len(positions), len(mutations)))
    for i, pos in enumerate(sorted(positions)):
        for j, mut in enumerate(sorted(mutations)):
            freq = df[(df['position'] == pos) & (df['mutation'] == mut)]['frequency_among_candidates'].values
            heatmap_data[i, j] = freq[0] if len(freq) > 0 else 0

    im = ax.imshow(heatmap_data, cmap='YlOrRd', aspect='auto')
    ax.set_xticks(range(len(sorted(mutations))))
    ax.set_yticks(range(len(sorted(positions))))
    ax.set_xticklabels(sorted(mutations), rotation=45, ha='right')
    ax.set_yticklabels(sorted(positions))
    ax.set_xlabel('Mutation', fontsize=12, fontweight='bold')
    ax.set_ylabel('Position', fontsize=12, fontweight='bold')
    ax.set_title('Mutation Enrichment on Pareto Frontier', fontsize=13, fontweight='bold')

    # Add text annotations
    for i in range(len(sorted(positions))):
        for j in range(len(sorted(mutations))):
            text = ax.text(j, i, f'{heatmap_data[i, j]:.1%}',
                          ha="center", va="center", color="black", fontsize=9)

    cbar = plt.colorbar(im, ax=ax, label='Frequency')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "03_mutations.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Figure 3: Mutation enrichment")


# Note: a per-descriptor bar-panel comparison was retired from the manuscript;
# the same information is presented in Table 2 (baseline comparison) and the
# baseline_comparison.csv data file under results/phase3_analysis/.


def fig_search_space_visualization():
    """Figure 4: Search space visualization - all candidates colored by edit distance."""
    all_df = pd.read_csv(RESULTS_DIR / "phase2_pareto" / "evaluated_candidates.csv")
    frontier = pd.read_csv(RESULTS_DIR / "phase2_pareto" / "pareto_frontier.csv")

    fig, ax = plt.subplots(figsize=(10, 8))

    # All candidates
    scatter1 = ax.scatter(all_df['penetration_objective'], all_df['molecular_weight'],
                         c=all_df['edit_distance'], cmap='viridis',
                         s=30, alpha=0.5, label='Explored candidates')

    # Frontier
    scatter2 = ax.scatter(frontier['penetration_objective'], frontier['molecular_weight'],
                         c=frontier['edit_distance'], cmap='viridis',
                         s=200, edgecolors='red', linewidth=2, alpha=0.9,
                         label='Pareto frontier', marker='*')

    # Baseline
    baseline_mw = 563.65  # KTTKS MW
    baseline_pen = 0.395  # KTTKS penetration
    ax.scatter([baseline_pen], [baseline_mw], s=300, marker='D',
              color='gold', edgecolors='black', linewidth=2,
              label='Baseline (KTTKS)', zorder=10)

    ax.set_xlabel('Penetration Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Molecular Weight (Da)', fontsize=12, fontweight='bold')
    ax.set_title('Search Space Exploration: Edit Distance Effect', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(scatter1, ax=ax, label='Edit Distance')
    ax.legend(loc='best', fontsize=10)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "04_search_space.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Figure 4: Search space visualization")


def fig_edit_distance_distribution():
    """Figure 5: Edit distance distribution in frontier vs full space."""
    all_df = pd.read_csv(RESULTS_DIR / "phase2_pareto" / "evaluated_candidates.csv")
    frontier = pd.read_csv(RESULTS_DIR / "phase2_pareto" / "pareto_frontier.csv")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Full space distribution
    ed_counts_all = all_df['edit_distance'].value_counts().sort_index()
    ax1.bar(ed_counts_all.index, ed_counts_all.values, color='steelblue', alpha=0.7, edgecolor='black')
    ax1.set_xlabel('Edit Distance', fontsize=11)
    ax1.set_ylabel('Count', fontsize=11)
    ax1.set_title('Edit Distance Distribution (All Evaluated)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')

    # Frontier distribution
    ed_counts_frontier = frontier['edit_distance'].value_counts().sort_index()
    ax2.bar(ed_counts_frontier.index, ed_counts_frontier.values, color='coral', alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Edit Distance', fontsize=11)
    ax2.set_ylabel('Count', fontsize=11)
    ax2.set_title('Edit Distance Distribution (Pareto Frontier)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "05_edit_distance.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Figure 5: Edit distance distribution")


def main():
    """Generate all figures."""
    print("Generating publication figures...")

    try:
        fig_pareto_frontier()
        fig_convergence_comparison()
        fig_mutation_analysis()
        fig_search_space_visualization()
        fig_edit_distance_distribution()

        print(f"\n✅ All figures saved to {OUTPUT_DIR}/")
    except Exception as e:
        print(f"❌ Error generating figures: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
