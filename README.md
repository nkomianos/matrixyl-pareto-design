# Computational Evolution of Collagen-Stimulating Peptides for Enhanced Epidermal Permeability

## Overview

This project applies evolutionary algorithms and protein language models to optimize cosmetic peptides for enhanced skin penetration while maintaining collagen-binding activity. We target the critical challenge in peptide therapeutics: crossing the stratum corneum barrier without losing functional structure.

### The Problem
Cosmetic peptides (e.g., Matrixyl/palmitoyl pentapeptide-4) work in vitro but struggle with in vivo skin penetration. The barrier properties that make skin protective also prevent functional molecules from entering. Current approaches rely on chemical permeation enhancers or delivery vehicles, which add cost and complexity.

### The Solution
**Computational peptide optimization** using:
- Evolutionary algorithms (CMA-ES, NSGA-II multi-objective)
- Protein language models as structural/functional oracles
- Physicochemical constraints tied to transdermal delivery literature
- Iterative refinement balancing penetration vs. collagen-binding function

## Target Peptide

**Starting sequence**: Palmitoyl pentapeptide-4 (Matrixyl)
- Current MW: ~573 Da (slightly above ideal threshold)
- Known collagen-binding motif: well-characterized in literature
- Commercial baseline: ~10+ years in formulations, proven safe

## Optimization Constraints

| Property | Target | Rationale |
|---|---|---|
| TPSA | < 140 Ų | Literature: threshold for skin permeability |
| Molecular Weight | < 500 Da | Passive diffusion dominates; Matrixyl is 573 Da |
| LogP | 1–3 | Balance: hydrophobic enough to cross membrane, hydrophilic enough to dissolve |
| H-bond Donors/Acceptors | Minimize in core motif | Reduce interaction with stratum corneum polysaccharides |
| Flexibility (Gyration Radius) | Minimize | Lower conformational entropy = easier penetration |
| Collagen Binding | Preserve | Non-negotiable: functional constraint |

## Proposed Improvements to Baseline Concept

### 1. **Multi-Objective Optimization (NSGA-II)**
- Optimize Pareto frontier between penetration and binding, not single weighted sum
- Avoid getting trapped in local optima that sacrifice binding for penetration
- Enables trade-off analysis: "which solutions maintain 80% binding while gaining 50% penetration?"

### 2. **Ensemble Protein Language Models + Uncertainty**
- Use ESM-2, ProtBERT, and OmegaFold in weighted ensemble
- Quantify prediction uncertainty; prioritize high-confidence candidates
- Reduces bias from single model training data

### 3. **Specialized Skin Permeability Models**
- Integrate data-driven models predicting transdermal flux from TPSA + MW + LogP
- Layer on top of PLM predictions for double-validation
- Consider human skin vs. rat skin extrapolation factors

### 4. **Active Learning Loop**
- Don't optimize blindly; identify most informative sequences
- Allocate computational budget to high-uncertainty, high-impact designs
- Reduces total number of evaluations needed

### 5. **Hierarchical Evolutionary Strategy**
- Phase 1: Evolve for penetration properties (unconstrained)
- Phase 2: Constrain best penetrators, refine binding affinity
- Prevents premature convergence to locally optimal solutions

### 6. **Structure-Activity Relationship (SAR) Analysis**
- Evolutionary footprinting: which positions/mutations drive improvements?
- Alanine scanning on top candidates
- Guide interpretability for medicinal chemistry

### 7. **Molecular Dynamics Validation**
- Run 100+ ns MD simulations on top 5 candidates
- Measure RMSD stability, radius of gyration, H-bond networks
- Validate computational predictions before experimental testing

### 8. **Experimental Validation Pipeline**
- **In vitro**: Franz cell diffusion assays on human skin
- **Positive controls**: known skin-penetrating peptides
- **Negative controls**: original Matrixyl
- **Metrics**: cumulative penetration, lag time, flux rate

## Project Structure

```
Bio-paper/
├── README.md                          # This file
├── CLAUDE.md                          # Technical documentation
├── data/
│   ├── sequences/
│   │   ├── baseline_matrixyl.fasta   # Original sequence
│   │   └── optimized_candidates.csv  # Top evolved sequences
│   └── physicochemical/
│       └── training_data.csv          # TPSA/MW/LogP/permeability
├── src/
│   ├── oracle.py                      # Protein language model wrapper
│   ├── evolutionary_algorithm.py      # NSGA-II / CMA-ES implementation
│   ├── constraints.py                 # Physicochemical property calculators
│   ├── active_learning.py             # Uncertainty-driven selection
│   └── validation.py                  # MD simulations, metrics
├── experiments/
│   ├── 01_baseline_optimization.py    # Initial CMA-ES run
│   ├── 02_pareto_evolution.py         # NSGA-II multi-objective
│   ├── 03_uncertainty_sampling.py     # Active learning
│   └── 04_validation.py               # MD + metrics
├── results/
│   └── (output logs, candidate sequences, metrics)
└── scripts/
    ├── download_plm_weights.py        # Setup ESM-2, ProtBERT
    └── generate_hpc_job.py            # Cluster submission templates
```

## Experiment Plan

### Phase 1: Baseline Single-Objective Optimization (Week 1–2)
- CMA-ES on penetration score (TPSA + MW + LogP + gyration)
- Constraint: must maintain ≥80% collagen binding (via PLM)
- Population size: 50, generations: 200
- **Deliverable**: top 20 candidates, penetration vs. binding scatter plot

### Phase 2: Multi-Objective Pareto Optimization (Week 3–4)
- NSGA-II optimizing (penetration, binding, synthesis feasibility)
- 3-objective optimization
- **Deliverable**: Pareto frontier visualization, trade-off analysis

### Phase 3: Active Learning & Refinement (Week 5–6)
- Sample candidates from Pareto frontier using uncertainty
- Ensemble predictions vs. individual model
- Iterative refinement
- **Deliverable**: uncertainty-weighted candidate list

### Phase 4: Validation (Week 7–8)
- Top 5 candidates: MD simulations (100+ ns each)
- SAR analysis: alanine scanning, position importance
- Comparison with Matrixyl baseline
- **Deliverable**: validation report, figures for paper

### Phase 5: Experimental Design (Week 9–10)
- Design in vitro Franz cell assay protocol
- Select 2–3 top candidates for synthesis
- **Deliverable**: experimental methods, candidate justification

## Success Metrics

1. **Computational**:
   - Achieve ≥2–3× improvement in skin penetration score vs. Matrixyl
   - Maintain ≥85% collagen binding activity
   - All candidates pass MD stability checks

2. **Publication**:
   - Novel optimization framework (ensemble + multi-objective + active learning)
   - Clear SAR interpretation
   - Reproducible pipeline (open-source code)

3. **Commercial**:
   - Identified candidates feasible to synthesize
   - Patent landscape analysis (no conflicts)
   - Cost-benefit vs. formulation complexity

## Dependencies

- **Structural Biology**: BioPython, py3Dmol, ProDy
- **ML/Optimization**: ESM-2, scikit-optimize, DEAP (evolutionary algorithms), pymoo (NSGA-II)
- **Chemistry**: RDKit (TPSA, LogP, MW), ProLIF (binding predictions)
- **MD Simulations**: GROMACS or OpenMM + MDTraj
- **Data**: UniProt, PubChem, literature datasets on skin permeability

## References

### Key Literature
1. Matrixyl mechanism: Blanes-Mira et al., *Peptides* (2002)
2. Skin permeability: Potts & Guy, *Pharm. Res.* (1992); TPSA threshold
3. Peptide optimization: Leman et al., *Chem. Rev.* (2023) on protein design
4. Active learning: Osuna et al., *ML for Drug Discovery* (2024)
5. Multi-objective optimization: Deb et al., NSGA-II paper (2002)

## Authors & Contact

This project is designed for publication in *Nature Computational Science* or *ACS Synthetic Biology* (target audience: computational biologists + cosmetic chemists).

---

**Status**: Planning phase  
**Last Updated**: 2026-05-05
