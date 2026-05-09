# CLAUDE.md — developer quickstart for the Bio-paper repository

This file is read by AI coding assistants (and humans) at the start of a session.
It describes what this repository actually contains today; for narrative results
and the bioRxiv-bound manuscript, read [`README.md`](README.md) and the docs
listed below.

## What this repository is

A reproducible computational pipeline that redesigns the Matrixyl peptide
(palmitoyl pentapeptide-4, KTTKS core) for improved topical delivery while
preserving the functional motif. Methodology: tournament GA, NSGA-II,
RDKit-based descriptor scoring, conformer-ensemble structural validation, and
a sensitivity analysis on penalty weights. The 9-member Pareto frontier is
validated against an exhaustive enumeration of the 3,706-candidate edit-
distance-2 search space.

The pipeline is **deterministic** and **CPU-only**. There are no learned
oracles, no protein language models, and no GPU dependencies in the published
results. References to PLMs in the manuscript are limited to academic context
and future work.

## Source-of-truth docs

- [`README.md`](README.md) — entry point, headline results, quickstart, full
  pipeline command list, license info.
- [`PROJECT_STATUS.md`](PROJECT_STATUS.md) — short executive summary.
- [`COMPUTATIONAL_RESULTS.md`](COMPUTATIONAL_RESULTS.md) — full technical
  results.
- [`PAPER_DRAFT_OUTLINE.md`](PAPER_DRAFT_OUTLINE.md) — manuscript outline.
- [`manuscript/main.tex`](manuscript/main.tex) — LaTeX source for the bioRxiv
  preprint.

## Module map (what's actually in `src/`)

| Module | Purpose |
|---|---|
| `chemistry.py` | RDKit-backed molecular descriptor calculation (MW, TPSA, LogP, HBD/HBA, rotatable bonds, formal charge) for peptide sequences |
| `constraints.py` | Soft-penalty penetration scoring composed from the chemistry descriptors |
| `function_scores.py` | Functional-preservation score (identity, edit distance, BLOSUM62) |
| `search_space.py` | Fixed-length max-edit-distance search rules and exhaustive enumeration |
| `candidates.py` | `CandidateEvaluator` glue that scores a sequence against all constraints |
| `tournament_search.py` | Phase 1: single-objective tournament GA with elitism |
| `pareto_search.py` | Phase 2: NSGA-II with non-dominated sorting and crowding distance |
| `analysis.py` | Phase 3: mutation-enrichment statistics and baseline comparison |
| `structure.py` | Phase 4: RDKit ETKDG conformer ensembles + radius-of-gyration |
| `experimental_design.py` | Phase 5: rule-based synthesis-panel and control selection |

## Pipeline (six phases)

```
experiments/01_tournament_search.py     → results/phase1_tournament/
experiments/02_nsga2_pareto.py          → results/phase2_pareto/
experiments/03_candidate_analysis.py    → results/phase3_analysis/
experiments/04_structural_validation.py → results/phase4_structure/
experiments/05_experimental_design.py   → results/phase5_experimental/
experiments/06_sensitivity_analysis.py  → results/phase6_sensitivity/
```

Each script writes a `run_metadata.json` capturing seed, parameters, input
hashes, and runtime.

Reference invocations are listed in [`README.md`](README.md#pipeline-end-to-end).

## Key invariants

- Random seed `42` is fixed throughout (GA, NSGA-II, conformer generation).
- All published results live under `results/` and are committed.
- All five publication figures live under `manuscript/figures/`. They are
  produced by `scripts/generate_figures.py`, which writes there directly so
  the LaTeX build always sees the artifacts the script produced.
- Path handling in `scripts/generate_figures.py` is anchored to the repository
  root via `Path(__file__).resolve().parent.parent`, so it works from any CWD.

## Running tests and figures

```bash
python -m pytest tests/ -q     # 48 unit tests, all green on CPU
python scripts/generate_figures.py
```

## Common tasks

### Add a new physicochemical constraint

1. Implement the calculator in `src/constraints.py` (return float in `[0, 1]`).
2. Compose it into `PenetrationScore` so it participates in the existing
   product penalty.
3. Add a unit test in `tests/test_penetration_scoring.py`.

### Tune the GA / NSGA-II

Adjust hyperparameters in the experiment scripts (`experiments/01_*.py`,
`experiments/02_*.py`) — they are CLI flags, not source-level constants.

### Add a new figure

1. Add a function to `scripts/generate_figures.py` that reads from
   `RESULTS_DIR / "phaseN_*"` and writes to `OUTPUT_DIR / "NN_name.png"`.
2. Add the matching `\includegraphics{figures/NN_name.png}` reference in
   `manuscript/main.tex`.

## Status

- 48 unit tests pass.
- All six phases reproduce the committed results from a clean clone
  (`pip install -r requirements.txt && pytest`).
- Manuscript (`manuscript/main.tex`) is bioRxiv-format with line numbers,
  ORCIDs, full author/competing-interest/funding/data-availability statements,
  and `[Affiliation]` placeholders for the two authors to fill in before
  submission.

## What this repository does **not** contain

To set expectations for any AI assistant or contributor:

- No learned scoring layer (no PLM, no oracle, no surrogate model).
- No GPU code.
- No molecular-dynamics simulations.
- No closed-loop / autonomous design components.
- No experimental data — wet-lab work is proposed but not executed.

If you find references to any of the above in the codebase or docs, they are
stale and should be removed.
