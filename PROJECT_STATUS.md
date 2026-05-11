# Project Status: Computational Peptide Design for Topical Delivery

**Date**: 2026-05-11
**Status**: computational framework complete; preprint draft pending bioRxiv submission
**Repository**: https://github.com/nkomianos/matrixyl-pareto-design (branch `main`)

---

## What this project built

A **reproducible computational optimization pipeline** for discovering cosmetic peptide analogs with improved topical-delivery descriptors, starting from Matrixyl (KTTKS) baseline.

### Pipeline Phases (All Complete)

| Phase | Component | Algorithm | Output |
|-------|-----------|-----------|--------|
| **1** | Single-objective optimization | Tournament GA | Best: PTTPS (0.4578) |
| **2** | Multi-objective optimization | NSGA-II | 9-member Pareto frontier |
| **3** | Mutation analysis | Enrichment statistics | K4P dominant (56%) |
| **4** | Structural validation | RDKit conformers | Compactness metrics |
| **5** | Experimental design | Panel selection | 3 candidates + controls |
| **6** | Sensitivity analysis | Parameter sweeps (±30 %) | Frontier stable under 5 of 6 perturbations |
| **ED3** | Extended search space | Pareto (edit distance 3) | 20 candidates, diminishing returns |

---

## Key results

The Pareto frontier, K4P design driver, sensitivity-analysis findings, and structural compactness numbers are documented in the manuscript (`manuscript/main.tex`) and in [`COMPUTATIONAL_RESULTS.md`](COMPUTATIONAL_RESULTS.md). Three candidates are recommended for experimental validation: KTTPS (single-mutation, balanced), PTTPS (aggressive), and KTTPP (backup). This file is kept for project-status tracking only and should not be used as a primary results reference.

---

## Deliverables

### Documentation
- **[COMPUTATIONAL_RESULTS.md](COMPUTATIONAL_RESULTS.md)** — Archival working document of technical results
- **[PROJECT_STATUS.md](PROJECT_STATUS.md)** — this file (project tracking only)

### Code & Scripts
- **experiments/06_sensitivity_analysis.py** — Parameter robustness analysis
- **scripts/generate_figures.py** — figure-generation script for the manuscript
- **src/** — Core modules (all 48 unit tests passing)
  - `chemistry.py` — RDKit descriptor calculations
  - `constraints.py` — Penetration & function scoring
  - `candidates.py` — Sequence evaluation framework
  - `tournament_search.py` — GA implementation
  - `pareto_search.py` — NSGA-II implementation
  - `analysis.py` — Mutation enrichment & statistics
  - `structure.py` — Conformer generation

### Data & Results
```
results/
├── phase1_tournament/          # GA single-objective (best: PTTPS)
├── phase2_pareto/              # NSGA-II frontier (9 candidates)
├── phase3_analysis/            # Mutation enrichment analysis
├── phase4_structure/           # Conformer ensemble metrics
├── phase5_experimental/        # Synthesis panel + controls
├── figures/                    # 5 manuscript figures (300 dpi PNG)
│   ├── 01_pareto_frontier.png
│   ├── 02_convergence.png
│   ├── 03_mutations.png
│   ├── 04_search_space.png
│   └── 05_edit_distance.png
```

### CSV Outputs
- **pareto_frontier.csv** — 9 candidates with all descriptors & scores
- **mutation_enrichment.csv** — Position-by-mutation frequency matrix
- **synthesis_candidates.csv** — Top 3 with rationales
- **baseline_comparison.csv** — KTTKS, Pal-KTTKS, evolved analogs
- **candidate_summary.csv** — Detailed candidate profiles
- **convergence.csv** — Generation-by-generation optimization progress

---

## Ready for experimental validation

### Synthesis Panel
| Candidate | Role | Penetration | Function | MW | Rationale |
|-----------|------|-------------|----------|-----|-----------|
| **PTTPS** | High-penetration | 0.687 | 0.667 | 501.5 | Test aggressive optimization |
| **KTTPS** | Conservative | 0.546 | 0.833 | 532.6 | **Primary candidate** (balanced) |
| **KTTPP** | Backup | 0.678 | 0.670 | 542.6 | Alternative if PTTPS problematic |

### Controls
- Vehicle blank
- KTTKS (Matrixyl core)
- Pal-KTTKS (commercial lipidated Matrixyl)
- GPKGDPGA (collagen-like control)

### Proposed Assays
1. **Franz diffusion cell** — In vitro permeation across excised skin or RhE
2. **Fibroblast collagen assay** — Procollagen I/III secretion (if budget allows)
3. **Safety gating** — MTT viability, LDH release

---

## Summary numbers

| Metric | Value | Notes |
|--------|-------|-------|
| Search space size | 3,706 | Edit distance ≤ 2 from KTTKS |
| Pareto frontier | 9 candidates | 100% validated (exhaustive comparison) |
| Dominant mutation | K4P | Present in 55.6% of frontier |
| Best penetration gain | +74% | PTTPS vs KTTKS baseline |
| Recommended gain | +38% | KTTPS vs KTTKS (balanced risk-return) |
| Function preservation (best) | 100% | KTTKS (unmodified) |
| Function preservation (rec.) | 83% | KTTPS (1 mutation) |
| Computational time | ~20 min | Phases 1-5 on a single CPU core |
| Unit test coverage | 48 tests | 100% pass rate |
| Code quality | Reproducible | Fixed random seeds, full metadata logging |

---

## Pre-submission self-check

The deterministic phases of the pipeline (Phases 1--5 and the Phase 6 sensitivity sweep) are complete and committed. The Pareto frontier has been validated against an exhaustive enumeration of the 3,706-candidate search space. The 48-test unit-test suite passes from a clean clone. Code, data, and results are licensed under MIT; the preprint is distributed under CC BY 4.0. The Zenodo snapshot for the release accompanying the preprint is at [10.5281/zenodo.20126787](https://doi.org/10.5281/zenodo.20126787).

---

## Next steps

### Immediate (This Week)
1. **Get synthesis quotes** from peptide CRO
   - PTTPS, KTTPS, KTTPP (2-3 week turnaround typical)
   - Controls: KTTKS, Pal-KTTKS, GPKGDPGA

2. **Confirm lab resources**
   - [ ] Franz cell access or RhE tissue availability
   - [ ] LC-MS/MS or HPLC capability
   - [ ] Fibroblast cell line ready
   - [ ] Budget for controls

### Short Term (1-2 Weeks)
3. **Manuscript polish**
   - The full manuscript is at `manuscript/main.tex` and compiles cleanly via `tectonic main.tex`
   - 5 figures are at `manuscript/figures/`
   - Tables (Pareto frontier, mutations, RoG, sensitivity, validation) are in `manuscript/main.tex`

4. **Choose target journal**
   - Recommended: *Journal of Chemical Information and Modeling* (impact ~3.8)
   - Alternative: *Molecular Informatics*, *Computational Biology and Chemistry*
   - Timeline: Submission in 2-3 weeks, review ~3 months

### Medium Term (1-2 Months)
5. **Experimental validation**
   - Receive synthesis
   - Run Franz permeation assays
   - Collect preliminary fibroblast data
   - Use results as "future work" section in paper

---

## How to use these deliverables

**For writing the paper:**
1. The manuscript source lives at `manuscript/main.tex` and is the canonical reference
2. Background and supporting numbers are in `COMPUTATIONAL_RESULTS.md`
3. Figures are at `manuscript/figures/`

**For colleagues/reviewers:**
- Point to COMPUTATIONAL_RESULTS.md for full technical description
- Reference GitHub branch for reproducibility (all code + results)
- Share PDF of generated figures

**For experimental planning:**
- Use synthesis_candidates.csv for CRO quote request
- Reference assay_matrix.csv for lab protocol planning
- Share protocol_draft.md with collaborators

---

## Scientific limitations

The full discussion of limitations lives in the manuscript (§Limitations and Future Work). The principal items are: descriptor thresholds (MW < 500, TPSA < 140) are soft penalties inherited from small-molecule transdermal literature; the functional-preservation score is a sequence-based proxy (motif similarity, BLOSUM62), not a learned biological oracle; conformer ensembles from RDKit ETKDG are exploratory and not substitutes for molecular dynamics; skin-specific protein binding, enzymatic degradation, and formulation effects are not modelled; all claims here are computational and require wet-lab validation.

---

## Summary

The pipeline is complete and reproducible. The 9-member Pareto frontier matches the ground-truth Pareto set obtained by exhaustive enumeration; the K4P substitution recurs across all robust perturbations of the descriptor weights; three candidates (PTTPS, KTTPS, KTTPP) are proposed for experimental validation. All numerical findings are in `manuscript/main.tex`; this document is internal project tracking and is not intended as a results reference.

---

**Computational work**: commodity CPU hardware (no GPU required).
**Tests**: 48 unit tests, fully reproducible.
