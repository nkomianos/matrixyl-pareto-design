# Constrained Evolutionary Design of Matrixyl Analogs

A reproducible computational framework for redesigning the Matrixyl (palmitoyl pentapeptide-4, KTTKS core) peptide for improved topical delivery while preserving its functional motif. The pipeline combines exact RDKit molecular descriptors, a tournament genetic algorithm, NSGA-II multi-objective optimization, conformer-ensemble structural validation, and a sensitivity analysis on penalty weights. The final 9-member Pareto frontier is validated by exhaustive enumeration of the 3,706-candidate edit-distance-2 search space.

> **Status.** Computational framework complete. All six phases run end-to-end on CPU; the headline 9-member Pareto frontier is validated by exhaustive enumeration of the 3,706-candidate edit-distance-2 search space. Wet-lab validation is proposed but not performed.

---

## Results

The pipeline returns a 9-member Pareto frontier in the trade-off between predicted skin permeability and motif preservation, with three sequences proposed for experimental validation. The frontier matches the ground-truth Pareto set obtained by exhaustive enumeration of the entire 3,706-candidate search space. All numerical results, figures, and the discussion of the dominant design driver are in the manuscript at [`manuscript/main.tex`](manuscript/main.tex); per-phase CSVs are under [`results/`](results/).

---

## Repository layout

```
matrixyl-pareto-design/
├── README.md                          ← you are here
├── LICENSE                            ← MIT (code, data, results)
├── CITATION.cff                       ← citation metadata
├── PROJECT_STATUS.md                  ← short executive summary
├── COMPUTATIONAL_RESULTS.md           ← full technical results
├── requirements.txt
├── data/
│   ├── molecules/                     ← Pal-KTTKS SMILES + PubChem record
│   ├── references/                    ← curated literature snapshots
│   └── sequences/                     ← KTTKS, controls (FASTA)
├── experiments/                       ← phase scripts (01–08)
├── src/                               ← reusable Python modules
├── tests/                             ← 48 unit tests
├── scripts/
│   └── generate_figures.py            ← regenerates manuscript figures
├── results/                           ← deterministic-phase outputs (CSV / JSON)
│   ├── phase1_tournament/
│   ├── phase2_pareto/
│   ├── phase3_analysis/
│   ├── phase4_structure/
│   └── phase5_experimental/
└── manuscript/                        ← LaTeX source + bundled figures
    ├── main.tex
    ├── figures/
    └── README.md
```

---

## Quick start

```bash
# 1. Set up Python environment (3.10+; 3.11 tested)
python -m pip install -r requirements.txt

# 2. Run unit tests
python -m pytest tests/ -v   # expect 48 passed

# 3. Regenerate publication figures from results/ CSVs
python scripts/generate_figures.py
# ↳ writes manuscript/figures/01_pareto_frontier.png … 05_edit_distance.png

# 4. Build the manuscript
cd manuscript && pdflatex main.tex && pdflatex main.tex
```

All results are produced by deterministic CPU phases — no GPU or model weights are required to reproduce any figure or CSV under `results/`.

---

## Pipeline (end-to-end)

Each script writes a self-contained directory with a `run_metadata.json` capturing seed, parameters, code version, input checksums, and runtime. Outputs are tracked in git so anyone can read them without rerunning.

| # | Phase | Script | Output | Algorithm |
|---|---|---|---|---|
| 1 | Tournament search (single-objective) | `experiments/01_tournament_search.py` | `results/phase1_tournament/` | GA, tournament selection, elitism |
| 2 | Pareto search (multi-objective) | `experiments/02_nsga2_pareto.py` | `results/phase2_pareto/` | NSGA-II with crowding distance |
| 3 | Mutation enrichment + baseline comparison | `experiments/03_candidate_analysis.py` | `results/phase3_analysis/` | position-frequency statistics |
| 4 | Structural validation | `experiments/04_structural_validation.py` | `results/phase4_structure/` | RDKit ETKDG conformer ensembles |
| 5 | Experimental panel design | `experiments/05_experimental_design.py` | `results/phase5_experimental/` | rule-based panel + control selection |
| 6 | Sensitivity analysis (±30 % weights) | `experiments/06_sensitivity_analysis.py` | (run on demand) | re-evaluation of frontier |

Reference invocations (matching the metadata used to produce the committed results):

```bash
python experiments/01_tournament_search.py \
  --sequence data/sequences/matrixyl_core.fasta \
  --population 100 --generations 100 --tournament-size 3 \
  --mutation-rate 0.2 --elite-count 2 --seed 42 \
  --max-edit-distance 2 --top-k 25 \
  --output results/phase1_tournament

python experiments/02_nsga2_pareto.py --output results/phase2_pareto

python experiments/03_candidate_analysis.py \
  --frontier results/phase2_pareto/pareto_frontier.csv \
  --palmitoylated-smiles data/molecules/matrixyl_palmitoylated.smi \
  --output results/phase3_analysis

python experiments/04_structural_validation.py \
  --candidates results/phase3_analysis/candidate_summary.csv \
  --output results/phase4_structure \
  --top-k 5 --num-conformers 8 --palmitoylated-conformers 3 --seed 42

python experiments/05_experimental_design.py \
  --candidates results/phase3_analysis/candidate_summary.csv \
  --compactness results/phase4_structure/conformer_summary.csv \
  --output results/phase5_experimental --max-candidates 3

python experiments/06_sensitivity_analysis.py \
  --frontier results/phase2_pareto/pareto_frontier.csv \
  --output results/phase6_sensitivity
```

---

## Reproducibility

- **Seeds** are fixed at 42 throughout (GA, NSGA-II, conformer generation).
- **Versions** are pinned in `requirements.txt` (RDKit 2023.09.1, NumPy 1.24.3, Pandas 2.0.3, PyTorch 2.0.1).
- **Run metadata** for each phase is committed alongside its output (parameters, input hashes, runtime, code commit).
- **Validation**: the evolutionary Pareto frontier was checked against an exhaustive enumeration of the entire 3,706-candidate ED2 search space — 100 % precision and recall.
- **Tests**: 48 unit tests cover descriptor calculations, search-space rules, scoring functions, and the GA / NSGA-II implementations.

---

## Manuscript

LaTeX source and bundled figures live in [`manuscript/`](manuscript/). Build:

```bash
cd manuscript
pdflatex main.tex
pdflatex main.tex   # second pass for cross-references
```

A preprint of this manuscript has been deposited on bioRxiv. The Zenodo snapshot of the accompanying code release is at [doi:10.5281/zenodo.20126787](https://doi.org/10.5281/zenodo.20126787).

---

## Data availability

All inputs, intermediate results, and figures needed to reproduce the manuscript are committed to this repository:

- **Sequences**: `data/sequences/` (FASTA).
- **Reference molecules**: `data/molecules/matrixyl_palmitoylated.smi` and `data/references/palmitoyl_pentapeptide_4_pubchem.json` (PubChem CID 9897237).
- **Per-phase outputs**: `results/phase{1..6}_*` (CSV + JSON).
- **Figures**: `manuscript/figures/` (300 dpi PNG).
- **Manuscript source**: `manuscript/main.tex`.

No proprietary or restricted data are used.

---

## Limitations

1. Descriptor thresholds (MW < 500 Da, TPSA < 140 Å²) are inherited from small-molecule transdermal literature and are applied as soft penalties, not hard cut-offs.
2. The functional-preservation score is a sequence-based proxy (identity, edit distance, BLOSUM62) — not a learned biological oracle. Adding a learned scoring layer (e.g., a protein language model) is plausible future work but is outside the scope of this study.
3. Conformer ensembles come from RDKit ETKDG and are exploratory; explicit-solvent molecular dynamics would be required for stability claims.
4. No skin-binding, metabolism, or formulation effects are modelled — only passive-diffusion descriptors.
5. All claims here are computational. Wet-lab permeation and fibroblast-response assays are required before any biological or commercial claim.

---

## Citation

If you use this repository, please cite the manuscript and the supporting code release. The Zenodo DOI for the release is [10.5281/zenodo.20126787](https://doi.org/10.5281/zenodo.20126787). Machine-readable metadata is in [`CITATION.cff`](CITATION.cff); the full bibliography for the manuscript is at the end of [`manuscript/main.tex`](manuscript/main.tex).

---

## License

- **Code, data, and results** in this repository: MIT License — see [`LICENSE`](LICENSE).
- **Preprint manuscript** (the bioRxiv PDF): Creative Commons Attribution 4.0 International (CC BY 4.0).
