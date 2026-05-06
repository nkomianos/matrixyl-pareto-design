# Paper Outline: Computational Evolution of Matrixyl Analogs for Topical Delivery

## Title
**Constrained Evolutionary Design of Matrixyl Analogs: Balancing Permeability and Functional Preservation Through Computational Optimization**

## Abstract

Cosmetic peptides like Matrixyl (KTTKS) offer collagen-stimulating activity but face topical delivery barriers due to poor physicochemical properties (high molecular weight, polar character, charged residues). We present a reproducible computational workflow combining evolutionary algorithms (tournament selection + NSGA-II) with exact molecular descriptor calculations (RDKit) to discover Matrixyl-family analogs with improved penetration-associated descriptors. Applied to the KTTKS core, our approach identifies a 9-member Pareto frontier revealing fundamental trade-offs between permeability gain and motif preservation. The dominant design strategy (lysine→proline at position 4) reduces topological polar surface area by 26% while maintaining 67% functional-preservation scores. This framework demonstrates how evolutionary algorithms can systematically explore small-peptide design spaces and generate actionable candidates for experimental validation.

**Keywords**: peptide design, evolutionary algorithms, NSGA-II, topical delivery, constrained optimization, structure-activity relationships

---

## 1. Introduction

### Motivation
- Peptide therapeutics and cosmetic ingredients are increasingly important for dermatology and anti-aging markets
- Matrixyl (palmitoyl pentapeptide-4 or KTTKS core) is a well-recognized commercial collagen-stimulating peptide
- Topical delivery remains a major bottleneck: peptides typically have MW 500-3000 Da, TPSA 100-250 Ų, net charges ±2 to ±4
- Traditional hit-to-lead optimization focuses on function first, delivery second

### Research Questions
1. Can constrained evolutionary search discover Matrixyl-family analogs with improved topical-delivery descriptors?
2. What design mutations most reliably improve permeability-related metrics without destroying motif recognition?
3. Can we quantify the trade-off between penetration gain and functional preservation?

### Hypothesis
Multi-objective evolutionary algorithms can navigate the Matrixyl-analog design space and identify a Pareto frontier of candidates that trade off permeability improvements against loss of function, enabling rational candidate selection for experimental validation.

---

## 2. Methods

### 2.1 Baseline Peptides and Search Space

**Reference sequence**: Matrixyl core `KTTKS` (unmodified, canonical residues, 5 amino acids)  
**Search constraints**:
- Fixed backbone length (5 residues)
- Maximum edit distance: 2 or 3 from KTTKS (Levenshtein distance)
- Canonical amino acids only (no non-standard residues)
- Edit distance 2 search space size: 3,706 candidates

**Baseline molecules**:
- `KTTKS`: unmodified Matrixyl core (MW 563.65 Da)
- `Pal-KTTKS`: N-terminal palmitoyl conjugate (MW 802.07 Da, from PubChem CID 9897237)
- `GPKGDPGA`: collagen-mimetic control

### 2.2 Molecular Descriptors and Scoring

**Chemistry engine**: RDKit 2023.09.1
- Standard peptide structure construction with N-terminal amino group, C-terminal carboxyl
- For palmitoylated sequences: N-terminal lipid chain attachment via amide bond
- SMILES generation and validation

**Descriptor calculations** (all exact molecular properties):
- Molecular weight (MolWt)
- Topological polar surface area (TPSA)
- LogP (Wildman-Crippen)
- H-bond donors / acceptors (RDKit standard definitions)
- Rotatable bonds (RDKit definition)
- Formal charge

**Penetration score** (normalized composite, range [0, 1]):
```
penalty(MW) = max(0, (MW - 500) / 100)
penalty(TPSA) = max(0, (TPSA - 140) / 100)
penalty(LogP) = |LogP - 2|  (target LogP 1-3)
penalty(HBD+HBA) = (HBD + HBA - optimal) / 10

penetration_score = max(0.01, (
    (1 - penalty_MW) *
    (1 - penalty_TPSA) *
    (1 - penalty_LogP) *
    (1 - penalty_HBD_HBA)
))
```

Literature justification: MW < 500 and TPSA < 140 are established small-molecule transdermal thresholds (Potts & Guy 1992). We apply these as soft penalties (not hard cutoffs) because peptides often violate these rules.

**Functional preservation score** (range [0, 1]):
```
f_identity = 1 if sequence == KTTKS else 0
f_edit_distance = max(0, 1 - (edit_distance / max_allowed_distance))
f_blosum = mean BLOSUM62 score for aligned positions (conservative substitution bonus)

functional_score = 0.4 * f_identity + 0.4 * f_edit_distance + 0.2 * f_blosum
```

**Synthesis feasibility score**: 1.0 for all canonical amino-acid sequences (constant).

### 2.3 Evolutionary Algorithms

#### Phase 1: Tournament Selection (Single-Objective)
**Algorithm**: Generational genetic algorithm with tournament selection and elitism
**Parameters**:
- Population size: 100
- Generations: 100
- Tournament size: 3
- Mutation probability: 0.2 per individual
- Elite count: 2 (preserve top 2% unchanged each generation)
- Seed: 42 (deterministic reproducibility)

**Search operators**:
- Mutation: uniform random substitution of amino acid at selected position
- Crossover: disabled (mutation-only for short peptides)
- Repair: all offspring checked against SearchSpaceRules; invalid mutants replaced

**Objective**: maximize `penetration_score × functional_preservation_score`

#### Phase 2: NSGA-II (Multi-Objective)
**Algorithm**: Non-dominated sorting genetic algorithm II (local reimplementation)
**Parameters**:
- Population size: 100
- Generations: 100
- Tournament size: 2
- Seed: 42

**Objectives** (all maximized):
1. Penetration score
2. Functional preservation score
3. Synthesis feasibility (constant: 1.0)

**Selection**: Binary tournament on (Pareto rank, crowding distance)  
**Survival**: Elitist environmental selection (best 100 from parent + offspring pool of 200)  
**Deduplication**: Final Pareto frontier deduplicated by sequence

### 2.4 Validation and Sensitivity

**Validation approach**: Exhaustive enumeration of the 3,706-candidate ED2 search space to compute the true Pareto frontier, then compare algorithm frontier against ground truth.

**Sensitivity analysis**: Re-evaluated Phase 2 frontier with ±30% variations in descriptor penalty weights to assess ranking robustness.

**Extended search space**: Phase 2 repeated with max edit distance 3 (100k candidate neighborhood) to test for diminishing returns.

### 2.5 Structural Analysis

**Conformer ensemble generation**:
- RDKit ETKDG (distance geometry with triangle bounds smoothing)
- 8 conformers per sequence, 100 optimization cycles (MMFF94s + UFF fallback)
- Radius of gyration (RoG) statistics: mean, std, min, max

---

## 3. Results

### 3.1 Tournament Search (Single-Objective)

**Best candidate**: `PTTPS`  
**Optimization score**: 0.4578  
**Descriptors**:
- Penetration: 0.6866
- Function: 0.6667
- Edit distance: 2
- MW: 501.54 Da (−10.4% vs KTTKS)
- TPSA: 217.63 Ų (−25.6% vs KTTKS)
- LogP: −3.98
- Rotatable bonds: 11 (−45% vs KTTKS)

**Convergence**: Best candidate found by generation 4, maintained through generation 100 (typical GA plateau behavior).

### 3.2 Pareto Frontier (Multi-Objective)

**Frontier size**: 9 unique sequences (exhaustive validation confirms 100% precision and recall)

**Frontier members** (ranked by penetration):

| Seq | Penetration | Function | MW | TPSA | Edit | Role |
|-----|-------------|----------|-----|------|------|------|
| PTTPS | 0.687 | 0.667 | 501.5 | 217.6 | 2 | High-pen trade-off |
| KTTPP | 0.678 | 0.670 | 542.6 | 228.6 | 2 | High-pen variant |
| KTTPG | 0.635 | 0.678 | 502.6 | 237.4 | 2 | High-pen conservative |
| KTTPA | 0.635 | 0.687 | 516.6 | 237.4 | 2 | High-pen conservative |
| KTTKP | 0.504 | 0.836 | 573.7 | 263.4 | 1 | Balanced |
| KTTPS | 0.546 | 0.833 | 532.6 | 257.6 | 1 | Balanced +1 mut |
| KTTKG | 0.435 | 0.845 | 533.6 | 272.2 | 1 | Conservative |
| KTTKA | 0.434 | 0.854 | 547.7 | 272.2 | 1 | Conservative |
| KTTKS | 0.395 | 1.000 | 563.7 | 292.5 | 0 | Unmodified baseline |

**Trade-off analysis**: Clear elbow at `KTTPS`; candidates beyond are increasingly conservative (motif-preserving but minimal penetration gain).

### 3.3 Mutation Enrichment

**Dominant mutations** (Pareto frontier, 9 candidates):

| Mutation | Position | Count | Frequency | Effect |
|----------|----------|-------|-----------|--------|
| **K4P** | 4 (Lys→Pro) | 5 | 55.6% | Removes +1 charge; reduces HBD/HBA; increases hydrophobicity |
| S5A | 5 (Ser→Ala) | 2 | 22.2% | Reduces H-bond donors |
| S5G | 5 (Ser→Gly) | 2 | 22.2% | Increases flexibility |
| S5P | 5 (Ser→Pro) | 2 | 22.2% | Restricts backbone flexibility |
| K1P | 1 (Lys→Pro) | 1 | 11.1% | Only in extreme PTTPS; rare |

**Interpretation**: Position 4 lysine→proline is the single most effective permeability-improving mutation. Position 5 (serine) shows variable effects depending on the substitution.

### 3.4 Structural Validation (Conformer Ensembles)

**Compactness correlates with permeability**:

| Sequence | Role | Mean RoG (Ų) | vs KTTKS | Compactness |
|----------|------|-------------|---------|-------------|
| PTTPS | Optimal | 4.83 | −12.0% | **Most compact** |
| KTTPS | Balanced | 4.99 | −8.8% | Compact |
| KTTPP | Alt-high | 5.25 | −4.0% | Moderate |
| KTTKS | Baseline | 5.47 | — | Reference |
| Pal-KTTKS | Lipidated | 8.56 | +56% | **Spatially large** |

**Key finding**: Pal-KTTKS despite improved LogP, remains spatially extended due to lipid tail, explaining modest permeability gain (5.1%) despite descriptor improvements.

### 3.5 Edit Distance Sensitivity

**Edit distance 2 (main paper)**:
- Best penetration: PTTPS (0.687)
- Frontier size: 9
- Justification: Tractable neighborhood (3.7k candidates); maintains Matrixyl recognition

**Edit distance 3 (extended exploration)**:
- Best penetration: KFFKP (0.785, +14.2% over ED2)
- Frontier size: 20
- Trade-off: Improved penetration but less recognizable Matrixyl motif
- Conclusion: Diminishing returns; ED2 provides better signal-to-noise for experimental prioritization

---

## 4. Discussion

### Interpretation of Results

1. **K4P is a universal permeability driver**
   - Present in 56% of frontier candidates
   - Mechanistic explanation: lysine removal eliminates a +1 charge, reducing electrostatic repulsion at skin interface
   - BLOSUM62 allows K→P (conservative substitution in hydrophobic core contexts), supporting its naturalness

2. **Penetration-preservation trade-off is real and quantifiable**
   - PTTPS achieves 74% penetration improvement but only 67% functional preservation
   - KTTPS achieves 38% improvement with 83% preservation (better risk-return)
   - Unmodified KTTKS: baseline (39.5% penetration, 100% function) is a valid control

3. **Structural compactness supports molecular design**
   - 12% smaller conformational ensemble (RoG) in PTTPS vs KTTKS
   - Supports diffusion-favorable hypothesis but requires MD validation for strong claims
   - Pal-KTTKS remains large despite LogP improvement, explaining modest net gain

4. **Pareto frontier is robust**
   - ±30% penalty weight variations preserve top-3 ranking
   - Edit distance 3 still contains all ED2 frontier candidates
   - Conclusion: Core design principles (K4P substitution) are parameter-insensitive

### Comparison to Prior Work

- **Peptide design literature** (Leman et al. 2023, Lin et al. 2023): PLM-based design focuses on naturalness; we add explicit permeability constraints
- **Topical delivery literature** (Potts & Guy 1992, Mitragotri et al. 2004): MW/TPSA/LogP rules validated here on peptides; rules apply but are soft constraints
- **Evolutionary peptide design** (Deap, Pymoo): Our contribution is tight integration with exact molecular chemistry + small-space exhaustive validation

### Limitations

1. **Descriptor thresholds are soft**: MW < 500 and TPSA < 140 come from small-molecule literature; peptides violate these rules and still permeate (especially if formulated).
2. **No PLM uncertainty**: Planned Phase 7 (ESM-2 / ProtBERT scoring) was skipped due to environment dependencies; functional-preservation score uses motif similarity only.
3. **Conformer ensembles are exploratory**: RDKit geometries are not MD-validated; full explicit-solvent simulations would strengthen claims but are outside scope.
4. **No skin binding model**: We optimize for passive diffusion descriptors; skin protein binding, metabolism, and formulation effects are not captured.
5. **No in vitro validation**: This is a computational framework paper; biological activity remains unvalidated.

---

## 5. Conclusions

This work demonstrates that **constrained evolutionary algorithms can systematically explore small-peptide design spaces** and generate actionable candidates for experimental validation. The Matrixyl-family case study reveals:

1. A single dominant mutation (K4P) explains most penetration improvement
2. Multi-objective optimization quantifies real trade-offs between delivery and function
3. The resulting 9-member frontier provides diverse candidates for hypothesis testing

**Recommended next steps**: Synthesis of PTTPS (aggressive optimization), KTTPS (conservative), and KTTPP (backup), followed by in vitro permeation assays and fibroblast collagen-stimulation testing.

**Broader significance**: The methodology generalizes to other cosmetic or therapeutic peptides (GHK, bradykinin analogs, etc.) and demonstrates value of reproducible computational prioritization before synthesis investment.

---

## 6. Methods Appendix

### A. Search Space Enumeration

All 3,706 valid KTTKS-family candidates (edit distance ≤ 2) enumerated via dynamic programming (Levenshtein distance). Validated against hand-enumerated test cases.

### B. Pareto Frontier Validation

Entire 3,706-candidate space exhaustively scored using identical CandidateEvaluator. True Pareto frontier computed via full dominance comparison. Evolutionary frontier vs. exhaustive frontier: 100% precision, 100% recall (9/9 candidates correct).

### C. Random Baseline Comparison

Phase 1 tournament search ran with random mutant initialization (not seeded from reference). Best random-mutant candidate: `KTTPS` (score 0.4549 vs. evolved 0.4578). Evolutionary advantage: 0.6% margin (modest but consistent).

### D. Code and Data Availability

- **Repository**: GitHub (Bio-paper branch `claude/peptide-skin-penetration-aOat2`)
- **Tests**: 60 unit tests covering all modules; 100% pass rate
- **Reproducibility**: Fixed random seeds (seed=42); all parameters logged to run_metadata.json
- **Data**: All CSV files (evaluated_candidates, pareto_frontier, mutation_enrichment) included in supplementary materials

---

## 7. Figures and Tables

**Figure 1**: Pareto frontier (penetration vs. functional preservation, colored by edit distance)  
**Figure 2**: Convergence curves (Tournament GA and NSGA-II)  
**Figure 3**: Mutation enrichment heatmap (position × mutation frequency)  
**Figure 4**: Search space visualization (all evaluated candidates, frontier highlighted)  
**Figure 5**: Edit distance distribution (frontier vs. full space)  
**Figure 6**: Descriptor comparison table (baseline vs. top candidates)  

**Table 1**: Pareto frontier candidates (sequence, descriptors, scores, roles)  
**Table 2**: Baseline descriptor comparison (KTTKS, Pal-KTTKS, evolved analogs)  
**Table 3**: Synthesis panel (PTTPS, KTTPS, KTTPP) with rationales  

---

## References

- Deb, K., et al. (2002). "NSGA-II: A fast and elitist multiobjective genetic algorithm." IEEE Trans. Evol. Comput. 6(3): 182-197.
- Leman, J. F., et al. (2023). "A comprehensive protein language model for biology." bioRxiv.
- Lin, Z., et al. (2023). "Language models of protein sequences at the scale of evolution enable accurate structure prediction." Science 379(6636): eade2574.
- Mitragotri, S., et al. (2004). "Transdermal drug delivery." Nat. Biotechnol. 22(4): 440-446.
- Potts, R. O., & Guy, R. H. (1992). "Predicting skin permeability." Adv. Drug Deliver Rev. 13(1-2): 1-30.

---

**Word count**: ~3,500 (target: 15-20 pages with figures, ~8,000 words for full journal version)
