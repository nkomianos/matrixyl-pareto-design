# Matrixyl Peptide Optimization - Overleaf Project

Complete LaTeX manuscript ready for journal submission.

## Files

- `main.tex` — Complete manuscript (760 lines)
- `figures/` — Publication-quality figures (300 dpi PNG)
  - `01_pareto_frontier.png` — Pareto frontier scatter
  - `02_convergence.png` — Convergence curves
  - `03_mutations.png` — Mutation enrichment heatmap
  - `04_search_space.png` — Search space visualization
  - `05_edit_distance.png` — Edit distance distribution

## To Compile

1. Upload to Overleaf or use pdflatex locally:
   ```bash
   pdflatex main.tex
   pdflatex main.tex  # Run twice for references
   ```

2. Output: `main.pdf`

## Document Structure

- Title & Abstract with keywords
- Introduction (motivation, prior work, research questions)
- Methods (algorithms, descriptor calculations, validation)
- Results (6 phases: GA, Pareto frontier, mutations, RoG, edit distance, sensitivity)
- Discussion (interpretation, comparison to literature, limitations)
- Conclusions & Next Steps
- Appendices (search space enumeration, code availability)
- References (15+ citations)

## Submission Checklist

- [x] All figures embedded and compiling
- [x] All tables formatted
- [x] Citations complete
- [x] No LaTeX errors
- [x] Ready for journal submission

## Journal Recommendations

- Journal of Chemical Information and Modeling (JCIM)
- Molecular Informatics
- Computational Biology and Chemistry

All supporting data and code available on GitHub:
https://github.com/nkomianos/Bio-paper (branch `claude/peptide-skin-penetration-aOat2`)
