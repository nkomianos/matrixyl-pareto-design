# Computational Optimization Results: Matrixyl Analogs for Enhanced Epidermal Permeability

> **Reference note.** This is an archival working document retained for development tracking. The published findings, with the figures and tables that accompany them, are in the manuscript at `manuscript/main.tex`. If a number disagrees between this file and the manuscript, the manuscript is authoritative.

**Pipeline status:** all phases complete (Tournament GA + Pareto + Analysis + Structural Validation + Experimental Design + Sensitivity Analysis).
**Initial run date:** 2026-05-05; sensitivity analysis re-run 2026-05-11.

---

## Executive Summary

A constrained evolutionary optimization workflow was used to discover Matrixyl-family (`KTTKS`) peptide analogs with improved topical-delivery descriptors while preserving functional motifs. **Real computational results** show:

- **9-member Pareto frontier** separating permeability vs. motif-preservation trade-offs
- **Top candidate**: `PTTPS` (+74% penetration descriptor vs KTTKS, 67% functional preservation)
- **Conservative candidate**: `KTTPS` (+38% penetration descriptor vs KTTKS, 83% functional preservation)
- **Dominant mutation**: K4P (lysine→proline at position 4, present in 56% of frontier)
- **Synthesis panel**: 3 candidates selected for experimental validation

---

## Phase 1: Tournament Search (Single-Objective GA)

**Algorithm**: Mutation-only tournament selection with elitism  
**Parameters**: 100 population, 100 generations, tournament size 3, seed 42  
**Search Space**: Fixed-length `KTTKS` with max edit distance 2 (3,706 valid sequences)

**Results**:
- Best sequence: `PTTPS`
- Best optimization score: **0.4578**
- Penetration score: **0.6866**
- Functional score: **0.6667**
- Edit distance: 2

**Descriptor Profile for PTTPS**:
```
Molecular Weight:     501.54 Da   (↓ -10.4% vs KTTKS)
TPSA:                 217.63 Ų   (↓ -25.6% vs KTTKS)
LogP:                -3.98       (better hydrophobicity balance)
H-bond donors:        8          (↓ -27.3% vs KTTKS)
H-bond acceptors:     9          (↓ -18.2% vs KTTKS)
Rotatable bonds:      11         (↓ -45.0% vs KTTKS)
Formal charge:        0
```

---

## Phase 2: Pareto Frontier (Multi-Objective Optimization)

**Algorithm**: NSGA-II with non-dominated sorting and crowding distance  
**Parameters**: 100 population, 100 generations, seed 42  
**Objectives**: Penetration score, Functional preservation, Synthesis feasibility  

**Pareto Frontier (9 Candidates)**:

| Sequence | Penetration | Function | Synthesis | MW (Da) | TPSA (Ų) | Edit Dist |
|----------|-------------|----------|-----------|---------|----------|-----------|
| **PTTPS** | **0.687** | 0.667 | 1.0 | **501.5** | **217.6** | 2 |
| **KTTPP** | 0.678 | 0.670 | 1.0 | 542.6 | 228.6 | 2 |
| KTTPG | 0.635 | 0.678 | 1.0 | 502.6 | 237.4 | 2 |
| KTTPA | 0.635 | 0.687 | 1.0 | 516.6 | 237.4 | 2 |
| KTTKP | 0.504 | 0.836 | 1.0 | 573.7 | 263.4 | 1 |
| **KTTPS** | 0.546 | **0.833** | 1.0 | 532.6 | 257.6 | 1 |
| KTTKG | 0.435 | 0.845 | 1.0 | 533.6 | 272.2 | 1 |
| KTTKA | 0.434 | 0.854 | 1.0 | 547.7 | 272.2 | 1 |
| **KTTKS** | 0.395 | **1.000** | 1.0 | 563.7 | 292.5 | 0 |

**Key Insight**: The frontier clearly separates three strategy zones:
- **High-penetration trade-offs** (`PTTPS`, `KTTPP`): 67-68% penetration, 67% function
- **Balanced conservative** (`KTTPS`): 55% penetration, 83% function
- **Function anchor** (`KTTKS`): 40% penetration, 100% function (no mutations)

---

## Phase 3: Mutation Analysis

**Dominant Mutations** (Pareto frontier, 9 candidates):

| Mutation | Position | Count | Frequency | Effect |
|----------|----------|-------|-----------|--------|
| **K4P** | 4 (Lysine→Proline) | 5 | **55.6%** | Removes positive charge; reduces TPSA & H-bonds |
| S5A | 5 (Serine→Alanine) | 2 | 22.2% | Reduces H-donors |
| S5G | 5 (Serine→Glycine) | 2 | 22.2% | Increases flexibility |
| S5P | 5 (Serine→Proline) | 2 | 22.2% | Rigidifies C-terminus |
| K1P | 1 (Lysine→Proline) | 1 | 11.1% | Only in `PTTPS` (high-penetration extreme) |

**Interpretation**: K4P is the most robust permeability-improving mutation. Position 5 substitutions are variable and context-dependent.

---

## Phase 4: Structural Validation (RDKit Conformers)

**Method**: ETKDG conformer ensemble (8 conformers per sequence, 100-step MMFF optimization)  
**Key Metric**: Radius of gyration (RoG) - measure of conformational compactness

**Compactness Results**:

| Sequence | Role | Mean RoG (Å) | vs KTTKS | Notes |
|----------|------|-------------|---------|-------|
| **PTTPS** | High-penetration | **5.05** | -6.9% | **Most compact** |
| **KTTPS** | Balanced | **5.08** | -6.3% | More rigid backbone |
| KTTPP | Backup | 5.25 | -3.3% | Moderate compactness |
| KTTKS | Baseline | 5.42 | -- | Reference |
| Pal-KTTKS | Lipidated | 8.38 | +54.5% | **Large spatial extent** |

**Interpretation**: K4P replacement (→ PTTPS) reduces conformational size, potentially improving diffusion. Pal-KTTKS remains spatially large despite improved LogP, confirming that lipidation creates a different permeability problem.

---

## Phase 5: Experimental Validation Design

**Selected Synthesis Panel** (3 candidates):

### 1. PTTPS (High-Penetration Trade-off)
- **Rationale**: Highest predicted permeability; tests aggressive motif editing
- **Strengths**: Best descriptor profile; highest RoG improvement
- **Risk**: Modest functional preservation (67%); may lose collagen signaling in vitro
- **Use**: Proof-of-concept for penetration-optimized design

### 2. KTTPS (Balanced Conservative)
- **Rationale**: Single-mutation analog; strong motif preservation (83%)
- **Strengths**: Minimal synthetic cost; robust compactness; likely active if tested
- **Risk**: More modest permeability gain (+38% vs KTTKS) than PTTPS (+74%)
- **Use**: Conservative follow-up; likely to show both permeability AND functional activity

### 3. KTTPP (Backup High-Penetration)
- **Rationale**: Alternative high-penetration strategy; strong trade-off value
- **Strengths**: Similar penetration to PTTPS; more function preservation (67% vs 67%)
- **Risk**: Dual prolines may affect solubility
- **Use**: Backup if PTTPS synthesis/stability issues arise

**Controls**:
- Vehicle blank
- KTTKS (Matrixyl core)
- Pal-KTTKS (commercial lipidated Matrixyl)
- GPKGDPGA (collagen-like sequence control)

---

## Baseline Comparison

| Peptide | MW (Da) | TPSA (Ų) | LogP | HBD | Penetration | Function |
|---------|---------|----------|------|-----|-------------|----------|
| KTTKS (core) | 563.65 | 292.45 | -4.65 | 11 | **0.395** | 1.0 |
| Pal-KTTKS | 802.07 | 295.53 | +0.99 | 11 | 0.508 | 1.0 |
| **PTTPS** | **501.54** | **217.63** | -3.98 | **8** | **0.687** | 0.667 |
| **KTTPS** | 532.60 | 257.64 | -4.27 | 9 | **0.546** | **0.833** |

**Key Findings**:
- PTTPS achieves **74% higher penetration score** than KTTKS baseline
- KTTPS achieves **38% penetration gain** with minimal function loss
- Pal-KTTKS shows that lipidation improves LogP but cannot overcome MW/TPSA penalties alone
- The frontier reveals penetration-preservation trade-offs cannot be avoided, only managed

---

## Experimental Next Steps

### Franz Diffusion Cell Assay (In Vitro Permeation)
- Candidates: PTTPS, KTTPS, KTTPP + controls
- Skin: Excised human cadaver or reconstructed skin (RhE)
- Receptor: Saline ± formulation enhancers (DMSO, propylene glycol)
- Quantification: LC-MS/MS or validated HPLC
- Duration: 24-48 hours; sample at 2, 4, 8, 24, 48h

### Collagen/Procollagen Assay (Functional Preservation)
- Cell type: Human dermal fibroblasts (HDF)
- Readout: Procollagen Type I secretion or soluble collagen I/III by Sirius Red
- Doses: 0.1, 1, 10 µM (solubility-permitting)
- Positive control: Ascorbic acid or TGF-β
- Expected outcome: KTTPS should show activity; PTTPS activity questionable

### Safety Gating
- Fibroblast/keratinocyte viability (MTT, resazurin)
- LDH release (membrane integrity)
- Reconstructed human epidermis (RhE) irritation if candidates show activity

---

## Computational Methods Summary

### Permeability Scoring
- **Descriptors**: RDKit molecular weight, TPSA, LogP, H-bond donors/acceptors, rotatable bonds
- **Penalty function**: Distance from literature-informed optima (MW < 500, TPSA < 140, LogP 1-3)
- **Score range**: [0, 1] where 1 = optimal permeability

### Functional Preservation Scoring
- **Components**: 
  - Sequence identity (KTTKS anchor)
  - Edit distance (max 2)
  - BLOSUM62 conservative-substitution similarity
  - Length consistency
- **Score range**: [0, 1] where 1 = identical to KTTKS

### Synthesis Feasibility
- **Score**: 1.0 if canonical amino acids only
- **Cost proxy**: Not yet integrated; all candidates feasible for standard peptide synthesis

### Optimization Algorithms
- **Phase 1**: Tournament selection GA (mutation-only, elitism)
- **Phase 2**: NSGA-II (non-dominated sorting, crowding distance)
- **Validation**: Exhaustive enumeration of 3,706-candidate space (confirmed perfect Pareto recovery)

---

## Limitations and Caveats

1. **Descriptor thresholds** (MW < 500, TPSA < 140) come from small-molecule transdermal literature and may not perfectly apply to peptides
2. **Functional-preservation score is sequence-based** (motif similarity, edit distance, BLOSUM62) — adding a learned scoring layer such as a protein language model is plausible future work but is outside the scope of this study
3. **Conformer ensembles** from RDKit are exploratory; real 3D dynamics require MD simulations
4. **No explicit solubility modeling**; some candidates may be insoluble in standard formulations (e.g., PTTPS with multiple prolines)
5. **Palmitoylation effect** cannot be captured by sequence alone; Pal-KTTKS represents a different chemical species

---

## Manuscript figures

The five publication figures are at `manuscript/figures/` and are regenerated from the per-phase CSVs by `scripts/generate_figures.py`:

1. **Pareto frontier** (penetration vs. functional preservation)
2. **Convergence curves** (Tournament GA and NSGA-II)
3. **Mutation enrichment heatmap** (position-frequency)
4. **Search space visualization** (full 3,706-candidate space)
5. **Edit-distance distribution** (frontier vs. full space)

---

## Reproducibility

- **Random seeds**: Fixed (seed=42 for GA, seed=42 for conformers)
- **Environment**: Python 3.10, RDKit 2023.09, NumPy 1.24.3
- **Data checksums**: Available in run metadata files
- **Code version**: https://github.com/nkomianos/matrixyl-pareto-design (branch `main`)

---

## Next Immediate Actions

1. **Synthesis quotes**: Contact CRO for cost estimates for PTTPS, KTTPS, KTTPP
2. **Formulation development**: Test solubility in PBS, ethanol, DMSO, propylene glycol
3. **Positive controls**: Order ascorbic acid, L-carnitine, or other known collagen-active peptides
4. **Skin source**: Confirm availability of fresh cadaver or RhE tissue
5. **Draft manuscript**: Intro + Methods written; Results ready for integration

---

**Data generated**: deterministic evolutionary search + RDKit descriptor calculations on commodity CPU hardware
**Status**: ready for experimental validation
