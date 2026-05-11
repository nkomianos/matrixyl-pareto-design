# Manuscript

LaTeX source for the paper *Constrained Evolutionary Design of Matrixyl Analogs: Balancing Permeability and Functional Preservation Through Computational Optimization*.

## Files

- `main.tex` — full manuscript
- `figures/` — publication-quality figures (300 dpi PNG)
  - `01_pareto_frontier.png` — Pareto frontier scatter (penetration vs. function)
  - `02_convergence.png` — GA and NSGA-II convergence curves
  - `03_mutations.png` — mutation enrichment heatmap
  - `04_search_space.png` — search-space visualization with frontier highlighted
  - `05_edit_distance.png` — edit-distance distribution

## Build

Local:

```bash
cd manuscript
pdflatex main.tex
pdflatex main.tex   # second pass for cross-references
```

Overleaf: upload the `manuscript/` folder as a new project; set `main.tex` as the root document.

## Regenerating figures

Figures are produced from CSVs under `../results/` by:

```bash
python ../scripts/generate_figures.py
```

The script writes directly into `manuscript/figures/`, so no manual copying is required.

## Source of truth

Code, unit tests, and full results: <https://github.com/nkomianos/matrixyl-pareto-design> (branch `main`).
