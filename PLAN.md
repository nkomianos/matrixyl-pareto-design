# PLAN: Computational Evolution of Collagen-Stimulating Peptides for Enhanced Epidermal Permeability

## 1. Executive Summary

This repository is an early scaffold for a computational peptide-design paper. The intended thesis is strong: start from a recognizable cosmetic peptide family, evolve analogs under skin-penetration constraints, and preserve collagen-stimulating activity well enough to justify topical formulation relevance.

The current codebase captures the broad idea, but it is not yet a scientifically defensible pipeline. It contains:

- A `README.md` and `CLAUDE.md` that describe the target concept and intended modules.
- A baseline FASTA in `data/sequences/baseline_matrixyl.fasta`.
- Three source modules: `src/constraints.py`, `src/oracle.py`, and `src/evolutionary_algorithm.py`.
- One experiment entrypoint: `experiments/01_baseline_optimization.py`.
- No tests, no generated results, no active-learning module, no validation module, and no paper artifacts yet.

The immediate priority is to turn the scaffold into a reproducible computational study with corrected scientific assumptions, exact chemistry descriptors, a working evolutionary loop, transparent baselines, and figures that can support a research-paper narrative.

## 2. Scientific Correction Required Before Coding Forward

The repository currently treats `GPKGDPGA` as Matrixyl / palmitoyl pentapeptide-4. That is likely incorrect for the planned paper.

Matrixyl, commonly described as palmitoyl pentapeptide-4, is generally represented as `Pal-KTTKS`:

- Peptide core: `KTTKS`
- Lipidated cosmetic ingredient: palmitoyl-KTTKS
- Approximate molecular weight: about 563 Da for the unmodified peptide core and about 802 Da for the palmitoylated molecule
- Common CAS for palmitoyl pentapeptide-4: `214047-00-4`

The current FASTA contains:

```text
>Matrixyl|Palmitoyl-Pentapeptide-4|CAS-823801-91-6
GPKGDPGA
```

This looks more like a collagen-like motif fragment than Matrixyl. We need to decide whether the paper is optimizing:

1. `KTTKS` / `Pal-KTTKS` analogs as Matrixyl-family collagen-stimulating matrikines.
2. `GPKGDP` / collagen-mimetic peptides as collagen-binding or collagen-interacting motifs.
3. A broader benchmark where Matrixyl is the commercial comparator and collagen-like motifs are a second design family.

Recommended direction: use Matrixyl-family `KTTKS` as the primary commercial baseline because cosmetic brands will recognize it, and treat `GPKGDP` only as a separate collagen-like control if the paper needs one.

## 3. Repository Audit

### Current Files

`README.md`

- Strong high-level narrative and useful phase outline.
- Mentions modules and scripts that do not exist yet.
- Contains a Matrixyl molecular-weight statement that should be corrected or made precise by distinguishing unmodified `KTTKS` from `Pal-KTTKS`.
- Uses the TPSA threshold but currently includes a corrupted unit character. Use `A^2` or `angstrom^2` in ASCII.

`CLAUDE.md`

- Provides useful project context, architectural intent, and expected workflows.
- Describes RDKit descriptor calculation, PLM ensemble, NSGA-II, CMA-ES, active learning, and MD validation.
- Much of it is aspirational and not yet implemented.

`data/sequences/baseline_matrixyl.fasta`

- Needs scientific correction.
- Should be replaced or supplemented with explicit baselines:
  - `matrixyl_core.fasta`: `KTTKS`
  - `matrixyl_palmitoylated.smi` or structured metadata for `Pal-KTTKS`
  - optional `collagen_like_control.fasta`: `GPKGDPGA`

`src/constraints.py`

- Implements a heuristic sequence-property calculator.
- Does not use RDKit yet.
- Molecular-weight calculation is not chemically reliable. It uses residue masses and then adds `(n - 1) * water`, which overestimates peptide mass if the table already represents residue masses.
- TPSA is approximated from H-bond counts, and the penetration score uses a different TPSA proxy than `compute_tpsa`.
- LogP is an average Kyte-Doolittle hydrophobicity value, not a molecular LogP.
- H-bond donor and acceptor counts are residue heuristics, not molecule-level descriptor counts.
- Does not support palmitoylation, terminal capping, noncanonical residues, salt state, pH-specific charge, or SMILES/SDF input.

`src/oracle.py`

- Loads ESM-2 and ProtBERT-style transformer models.
- Scores "binding affinity" as embedding cosine similarity to a short reference sequence (`GPKGDP`).
- This is not a validated collagen-binding or collagen-stimulation oracle.
- ProtBERT typically expects spaced amino-acid tokens, so the current tokenizer path may be unreliable.
- Does not cache embeddings, so optimization will be slow and expensive.
- Does not define the biological target clearly: collagen binding, collagen stimulation, fibroblast response, or sequence naturalness.

`src/evolutionary_algorithm.py`

- Contains an NSGA-II wrapper but no implemented CMA-ES path.
- Uses positive penetration and binding scores as objective values even though `pymoo` minimizes by default. Unless objectives are negated, NSGA-II will prefer worse scores.
- Uses `A` as a padding character and then removes every `A`, which prevents alanine from being represented as a valid evolved residue.
- Does not seed the initial population with the reference sequence.
- Does not define custom sampling, crossover, mutation, or repair operators for peptide strings.
- Does not enforce motif preservation or sequence-family constraints.
- Allows length changes conceptually, but the fixed-length integer encoding is not robust.

`experiments/01_baseline_optimization.py`

- The header says "single-objective optimization" but the implementation runs NSGA-II with penetration and binding.
- Falls back to `oracle = None` if model loading fails and then skips optimization.
- Imports through `src`, which currently makes even lightweight modules fail if optional optimizer dependencies are missing.
- Writes CSV and config outputs but does not generate plots, logs, environment metadata, dataset checksums, or paper-ready artifacts.

`requirements.txt`

- Lists the right families of dependencies, but the current environment does not have `pymoo` installed.
- Pins may need refresh and platform validation, especially for RDKit on macOS.
- There is no lockfile, no `pyproject.toml`, and no environment setup guide beyond `requirements.txt`.

### Runtime Finding

In the current environment, importing `src.constraints` through the package fails because `src/__init__.py` eagerly imports `src.evolutionary_algorithm`, which imports `pymoo`, and `pymoo` is not installed. This means simple descriptor work is blocked by optional optimizer dependencies.

Using a direct file import for `src/constraints.py`, the current heuristic calculator reports for `GPKGDPGA`:

- Length: 8
- Molecular weight: about 805.4 Da by the current heuristic
- TPSA estimate: 270
- LogP proxy: -1.25
- Penetration score: about 0.39

These numbers should not be used in the paper until the chemistry engine is replaced with exact or at least well-documented molecular descriptor calculations.

## 4. Core Research Question

Can a constrained evolutionary search discover Matrixyl-family peptide analogs with improved topical-delivery descriptors while preserving a computational proxy for collagen-stimulating function?

The paper should avoid overclaiming direct biological efficacy from PLM scores. A defensible framing is:

- Primary computational objective: improve physicochemical descriptors linked to epidermal permeability.
- Functional preservation objective: preserve the matrikine-like signal motif and sequence plausibility, with optional docking or learned surrogate models if a credible target or dataset is available.
- Validation objective: nominate candidates for synthesis and in vitro testing rather than claim confirmed biological activity.

## 5. Hypotheses

H1: Constrained sequence optimization can produce `KTTKS`-family analogs with lower TPSA, fewer H-bond donors/acceptors, and improved LogP-like behavior relative to the unmodified Matrixyl core.

H2: Multi-objective search will produce a Pareto frontier showing a measurable trade-off between permeability-like descriptors and motif/function preservation.

H3: Explicit chemistry-aware descriptors will rank candidates differently than amino-acid heuristic descriptors, demonstrating why exact molecular representation matters for cosmetic peptide optimization.

H4: Lipidation state changes the optimization problem enough that `KTTKS` and `Pal-KTTKS` must be analyzed as separate baselines.

## 6. Target Optimization Constraints

Use raw descriptors and normalized scores. Do not hide the raw values.

Minimum descriptor set:

- Molecular weight: target below 500 Da for passive diffusion, with a separate discussion for lipidated Matrixyl at about 802 Da.
- TPSA: target below 140 A^2.
- LogP or cLogP: target 1 to 3 for topical/transdermal balance.
- H-bond donors and acceptors: minimize, while preserving required motif chemistry.
- Net charge: near neutral is generally preferred for passive diffusion; report at pH 5.5 and pH 7.4.
- Rotatable bonds: lower is generally preferable for permeability.
- Radius of gyration or compactness: estimate from conformer ensembles, not only sequence length.
- Synthetic feasibility: avoid difficult or costly residues/modifications unless justified.

Important nuance: peptides often violate small-molecule rules. The paper should present these constraints as design pressures, not absolute laws.

## 7. Proposed Architecture

### Data Layer

Create explicit baseline and metadata files:

- `data/sequences/matrixyl_core.fasta`
- `data/sequences/collagen_like_control.fasta`
- `data/molecules/matrixyl_palmitoylated.smi`
- `data/baselines/baseline_descriptors.csv`
- `data/literature/literature_claims.yaml`

Each baseline should include:

- Name
- Sequence
- Modification state
- CAS, if known
- Literature source
- Descriptor values
- Role in experiments

### Chemistry Layer

Replace sequence-only descriptor heuristics with a molecule-aware calculator:

- Convert canonical peptide sequences to molecular structures with explicit termini.
- Support N-terminal palmitoylation and optional acetylation/amidation.
- Use RDKit for molecular weight, TPSA, HBD, HBA, LogP, rotatable bonds, formal charge, and molecular formula.
- Generate conformers for approximate radius of gyration.
- Keep the current heuristic calculator only as a fast prefilter, clearly labeled as heuristic.

Suggested module:

- `src/chemistry.py`: molecule construction, modifications, RDKit descriptors.
- `src/constraints.py`: scoring functions built on descriptor records.

### Functional Oracle Layer

The current PLM embedding-similarity approach is acceptable as a placeholder but not as a final "binding affinity" model.

Recommended oracle hierarchy:

1. Motif preservation score: exact, interpretable constraint for `KTTKS`-family analogs.
2. PLM naturalness score: pseudo-likelihood or masked-residue plausibility from ESM-2.
3. Similarity-to-baseline score: embedding similarity to Matrixyl core, clearly described as a sequence plausibility proxy.
4. Optional biological surrogate: train or calibrate against literature/in vitro collagen-stimulation data if a dataset can be assembled.
5. Optional docking/MD: only if a credible receptor or interaction model is defined.

Suggested module:

- `src/oracle.py`: abstract oracle interface and PLM implementation.
- `src/function_scores.py`: motif preservation, baseline similarity, naturalness, uncertainty.

### Optimization Layer

Start with a simple, robust genetic algorithm before adding CMA-ES:

- Population initialized from the reference sequence plus conservative mutants.
- Mutation operators:
  - conservative substitution
  - hydrophobicity-directed substitution
  - donor/acceptor-reducing substitution
  - terminal modification toggles
  - short insertion/deletion only when length changes are part of the hypothesis
- Crossover: optional for such short peptides; mutation-only evolution may be easier to interpret.
- Constraints:
  - locked or partially locked motif positions
  - maximum edit distance from baseline
  - valid amino acids or allowed noncanonical residues
  - descriptor hard filters
- Selection:
  - tournament selection for single weighted objective
  - NSGA-II for Pareto optimization

Fix before using NSGA-II:

- Negate objectives for `pymoo` minimization.
- Replace `A` padding with an explicit gap token outside the amino-acid alphabet.
- Add custom sampling and mutation operators.
- Return raw and normalized objectives consistently.
- Save every evaluated candidate for auditability.

CMA-ES can be added later if optimizing continuous latent vectors, descriptor weights, or amino-acid embedding coordinates. It is not the first choice for direct discrete sequence edits.

### Experiment Layer

Create four reproducible experiment scripts:

- `experiments/00_baseline_descriptors.py`
- `experiments/01_tournament_search.py`
- `experiments/02_nsga2_pareto.py`
- `experiments/03_candidate_analysis.py`

Later:

- `experiments/04_md_validation.py`
- `experiments/05_experimental_design.py`

Each experiment should write:

- `config.yaml`
- `evaluated_candidates.csv`
- `top_candidates.csv`
- `pareto_frontier.csv`, where applicable
- `figures/*.png` and `figures/*.svg`
- `run_metadata.json` with timestamp, git commit, Python version, package versions, and random seed

### Analysis Layer

Add reusable plotting and analysis helpers:

- `src/analysis.py`
- `src/plotting.py`

Required figures for the paper:

- Figure 1: Conceptual workflow from Matrixyl baseline to evolved candidates.
- Figure 2: Baseline descriptor comparison for `KTTKS`, `Pal-KTTKS`, and candidate families.
- Figure 3: Optimization trajectory across generations.
- Figure 4: Pareto frontier: permeability score vs functional-preservation score.
- Figure 5: Sequence logo or mutation heatmap for top candidates.
- Figure 6: Candidate shortlist table with raw descriptors and rationale.
- Optional Figure 7: MD or conformer compactness comparison.

## 8. Implementation Roadmap

### Phase 0: Scientific Reset and Repository Hygiene

Goal: make the project scientifically coherent and runnable.

Tasks:

- Correct baseline files and distinguish `KTTKS`, `Pal-KTTKS`, and `GPKGDPGA`.
- Update `README.md` language so it does not claim implemented modules that do not exist.
- Move aspirational content into this plan or mark it as future work.
- Add `pyproject.toml` or environment setup instructions.
- Make `src/__init__.py` avoid eager importing optional heavy modules.
- Add a minimal `tests/` directory.
- Add a reproducibility section and expected commands.

Acceptance criteria:

- `python -c "from src.constraints import PhysicochemicalCalculator"` works without optimizer dependencies installed.
- Baseline sequence identities and molecular weights are documented.
- `pytest` runs at least the descriptor unit tests.

### Phase 1: Exact Baseline Descriptor Engine

Goal: generate trustworthy baseline chemistry.

Tasks:

- Implement RDKit-backed descriptor calculation.
- Add support for terminal states and palmitoylation.
- Validate descriptor outputs against PubChem or supplier data for `Pal-KTTKS`.
- Write unit tests for molecular weight, TPSA, HBD, HBA, LogP, charge, and rotatable bonds.
- Generate `results/baselines/baseline_descriptors.csv`.

Acceptance criteria:

- Descriptor outputs match known references within documented tolerance.
- Raw descriptor values are saved for each baseline.
- Heuristic descriptor code is no longer used for final reported results.

### Phase 2: Search Space and Scoring Definition

Goal: define what evolution is allowed to change.

Status: complete for the first evolutionary baseline. Exact descriptor-based permeability scoring, search-space validation, functional-preservation scoring, and the shared candidate schema are implemented and tested.

Tasks:

- Decide allowed sequence lengths and whether positions in `KTTKS` are locked, soft-locked, or mutable. Completed via `SearchSpaceRules`.
- Define allowed amino acids and optional modifications. Completed for canonical fixed-length sequence candidates; chemical modifications remain represented through molecule files such as `Pal-KTTKS`.
- Define permeability score from raw descriptors. Completed in `src/constraints.py` via RDKit descriptor scoring.
- Define functional-preservation score from motif preservation, edit distance, PLM similarity, and PLM uncertainty. Completed for motif preservation, edit distance, length consistency, and BLOSUM62-style conservative substitution. PLM scoring remains future work.
- Add hard filters for invalid or scientifically irrelevant candidates. Completed via `SearchSpaceRules.validate()`.
- Create a candidate schema shared across all experiments. Completed via `CandidateEvaluator` and `CandidateEvaluation`.

Acceptance criteria:

- Every candidate has sequence, modification state, raw descriptors, normalized scores, and validity flags.
- The score function can explain why a candidate passed or failed.
- Baseline candidates are evaluated through the same path as evolved candidates.

### Phase 3: Working Evolutionary Baseline

Goal: produce the first real candidate set.

Status: complete for the first deterministic baseline. The Phase 3 baseline uses deterministic tournament selection with elitism, mutation-only fixed-length peptide edits, and `CandidateEvaluator` as the only scoring path.

Tasks:

- Implement tournament-selection GA. Completed in `src/tournament_search.py`.
- Seed the initial population from Matrixyl core and conservative mutants. Completed with deterministic valid mutant initialization.
- Add deterministic random seeds. Completed through `TournamentSearchConfig.seed`.
- Save every generation. Completed in `TournamentSearchResult.evaluations_by_generation` and `experiments/01_tournament_search.py`.
- Add progress metrics and convergence plots. Completed with generation summaries, convergence CSV, and optional matplotlib convergence plot.
- Compare against a random-mutant baseline. Completed through `TournamentSearchResult.random_baseline` and `random_baseline_top_candidates.csv`.

Acceptance criteria:

- Running `experiments/01_tournament_search.py` produces candidate CSVs and figures.
- Best candidates improve permeability score over baseline without violating function constraints.
- Random baseline comparison is included.

Implementation notes:

- Use tournament size as the selection-pressure control, with deterministic winner selection inside each sampled tournament.
- Preserve the top 1 to 5 percent of each generation through elitism; for small populations, preserve at least one elite.
- Start with mutation-only sequence evolution because the `KTTKS` search space is short and interpretability matters more than recombination at this stage.
- Keep sequence length fixed for the first baseline and rely on `SearchSpaceRules` for hard validity checks.
- Use fixed random seeds in tests and experiment configs.

Completed:

- Added `src/tournament_search.py` with `TournamentSearchConfig`, `TournamentSearchOptimizer`, `GenerationSummary`, and `TournamentSearchResult`.
- Added deterministic population initialization from the Matrixyl core plus valid mutants.
- Added deterministic tournament selection, elitism, mutation-only reproduction, per-generation summaries, and random-baseline reporting.
- Added `experiments/01_tournament_search.py` to write evaluated candidates, top candidates, random-baseline candidates, convergence data, config, run metadata, and an optional convergence plot.
- Added unit tests for deterministic initialization, reproducible runs, valid final populations, and random-baseline reporting.

Verification:

- Focused tournament tests passed.
- Smoke command passed:
  `PYTHONDONTWRITEBYTECODE=1 python experiments/01_tournament_search.py --sequence data/sequences/matrixyl_core.fasta --population 8 --generations 2 --seed 123 --output results/phase3_smoke`
- Smoke result best sequence: `KTPAS`; best optimization score: `0.4233`.
- Full suite passed with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests`.
- Final suite status: 27 tests passed.
- No linter diagnostics were reported for `src/tournament_search.py`, `experiments/01_tournament_search.py`, `tests/test_tournament_search.py`, or `PLAN.md`.

Remaining work before Phase 4:

- Run a real Phase 3 experiment with larger population and generation settings.
- Review the candidate distribution and decide whether the hard max edit distance of 2 is too restrictive or appropriate for the paper's first result.
- Add plotting/analysis helpers if the generated CSVs need richer publication figures.

### 2026-05-05: Phase 3 Extensive Experiment

Research and sizing:

- Fixed-length `KTTKS` with max edit distance 2 gives a valid canonical neighborhood of 3,706 candidates.
- Genetic-algorithm parameter guidance supports population sizes in the 50-100 range for small combinatorial problems, with elitism and fixed seeds for reproducibility.
- Extensive run settings: population 100, generations 100, tournament size 3, mutation rate 0.2, elite count 2, seed 42.

Experiment gate:

- Run the full unit suite before the experiment.
- Write outputs to `results/phase3_tournament_extensive`.
- Inspect best candidate, convergence, unique candidate count, edit-distance distribution, descriptor ranges, and random-baseline comparison.
- Run the full unit suite again after the experiment.
- Decide whether max edit distance 2 is sufficient before moving to Phase 4.

Completed:

- Preflight suite passed: 27 tests.
- Ran:
  `PYTHONDONTWRITEBYTECODE=1 python experiments/01_tournament_search.py --sequence data/sequences/matrixyl_core.fasta --population 100 --generations 100 --tournament-size 3 --mutation-rate 0.2 --elite-count 2 --seed 42 --max-edit-distance 2 --top-k 25 --output results/phase3_tournament_extensive`
- Output files were written under `results/phase3_tournament_extensive`.
- The convergence figure was written successfully.
- Post-run suite passed: 27 tests.
- No linter diagnostics were reported for the phase files.

Results:

- Best evolved sequence: `PTTPS`.
- Best optimization score: `0.4578`.
- Best functional-preservation score: `0.6667`.
- Best penetration score: `0.6866`.
- Best edit distance from `KTTKS`: 2.
- Descriptor profile for `PTTPS`: MW `501.54`, TPSA `217.63`, RDKit MolLogP `-3.98`, HBD `8`, HBA `9`, rotatable bonds `11`, formal charge `0`.
- The best sequence first appeared by generation 4 and remained the best through generation 100.
- The run evaluated 10,100 rows and covered 2,438 unique sequences out of the 3,706-candidate max-edit-2 neighborhood.
- Unique evaluated edit-distance distribution: edit 0 = 1, edit 1 = 94, edit 2 = 2,343.
- Best random-baseline candidate: `KTTPS` with optimization score `0.4549`.
- Exhaustive enumeration of all 3,706 valid candidates confirmed `PTTPS` is the global best under the current scalar scoring contract.

Interpretation:

- The tournament search is working, but the neighborhood is small enough that exhaustive enumeration is feasible and useful as a control.
- The current scalar objective rewards lower polarity and H-bond burden but still leaves the best candidate just above the MW target and far above the TPSA target.
- Max edit distance 2 is appropriate for the first conservative baseline, but Phase 4 should use Pareto analysis instead of a single scalar so we can inspect the trade-off between permeability gain and motif preservation.
- Before paper claims, we should compare `PTTPS` against `KTTPS`, `KTTPP`, `KPTPS`, and `KTPPS` in a Pareto/frontier view rather than declare a single winner.

### Phase 4: Pareto Optimization

Goal: produce a paper-ready trade-off analysis.

Status: implemented and smoke-tested. A dependency-light NSGA-II-style optimizer now replaces the old `pymoo` scaffold for this phase, so objective direction, sequence validity, and mutation behavior are explicit and unit-tested.

Research notes:

- NSGA-II uses non-dominated sorting, crowding distance, binary tournament parent selection, and elitist survival from the combined parent plus offspring population.
- All Phase 4 objectives are maximized: permeability score, functional-preservation score, and synthesis-feasibility score.
- Pareto frontiers are the right framing because peptide design involves competing objectives; a scalar "winner" can hide loss of motif preservation or synthesis risk.

Tasks:

- Fix and harden the NSGA-II implementation. Completed via `src/pareto_search.py`.
- Use correct minimization signs. Completed by using explicit maximization in local Pareto utilities.
- Implement custom peptide sampling, repair, mutation, and decoding. Completed by reusing fixed-length mutation and `SearchSpaceRules`.
- Track penetration, functional preservation, synthesis feasibility, and uncertainty. Completed for penetration, function, and synthesis; uncertainty remains future PLM work.
- Plot Pareto frontiers and label baseline positions. Completed for frontier plotting; explicit baseline labels remain a plotting enhancement.

Acceptance criteria:

- `experiments/02_nsga2_pareto.py` produces a stable Pareto frontier.
- The frontier shows meaningful trade-offs.
- Top candidates can be traced back to exact mutations and descriptor changes.

Implementation gate:

- Add unit-tested Pareto utilities: objective extraction, dominance, non-dominated sorting, crowding distance, and survival selection.
- Add a Pareto optimizer that uses `CandidateEvaluator` and fixed-length mutation only.
- Add a Phase 4 experiment script that writes `pareto_frontier.csv`, `final_population.csv`, `evaluated_candidates.csv`, `convergence.csv`, metadata, and a penetration-vs-function plot.
- Run focused tests after the module is added, then run the full test suite and a smoke Pareto experiment.

Completed:

- Added `src/pareto_search.py` with objective extraction, synthesis-feasibility scoring, dominance, non-dominated sorting, crowding distance, rank assignment, survival selection, and `ParetoSearchOptimizer`.
- Added `experiments/02_nsga2_pareto.py` to write `pareto_frontier.csv`, `final_population.csv`, `evaluated_candidates.csv`, `convergence.csv`, metadata, config, and `figures/pareto_frontier.png`.
- Added `tests/test_pareto_search.py` for synthesis score, objective extraction, maximized dominance, non-dominated sorting, crowding distance, and reproducible optimizer behavior.
- Exported Pareto APIs from `src/__init__.py`.

Verification:

- Focused Pareto tests passed.
- Full suite passed with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests`.
- Final suite status after Phase 4 implementation: 33 tests passed.
- No linter diagnostics were reported for Pareto modules, tests, experiment script, or `PLAN.md`.
- Smoke command passed:
  `PYTHONDONTWRITEBYTECODE=1 python experiments/02_nsga2_pareto.py --sequence data/sequences/matrixyl_core.fasta --population 20 --generations 5 --seed 123 --output results/phase4_pareto_smoke`
- Smoke frontier size: 6.
- Smoke frontier examples included `KTTKS` as the function-preserving reference and `KTAPS` as a higher-penetration / lower-function trade-off candidate.

Next Phase 4 step:

- Run a full Pareto experiment with population 100, generations 100, seed 42, max edit distance 2, and compare the resulting frontier against the exhaustive Phase 3 scalar optimum `PTTPS`.
- Add baseline labels to the plot if the frontier figure will be used directly in the paper.

### 2026-05-05: Phase 4 Full-Scale Pareto Test

Research and validation notes:

- Standard Pareto validation uses convergence and spread indicators such as hypervolume, generational distance, inverse generational distance, and frontier coverage.
- Because the current fixed-length max-edit-2 search space is only 3,706 candidates, exhaustive enumeration is a stronger validation control than approximate distance metrics.
- The full Phase 4 test does not need PLM/model stubs. It uses deterministic local objectives only: RDKit penetration score, functional-preservation score, and synthesis-feasibility score.

Experiment gate:

- Run the full unit suite before the experiment.
- Run `experiments/02_nsga2_pareto.py` with population 100, generations 100, seed 42, max edit distance 2.
- Write outputs to `results/phase4_pareto_full`.
- Enumerate all 3,706 valid candidates and compute the true Pareto frontier.
- Compare the evolutionary frontier against the exhaustive frontier for coverage and missing sequences.
- Inspect whether `PTTPS`, the Phase 3 scalar optimum, appears on the true and evolved Pareto frontiers.
- Run the full unit suite again after the experiment.

Robustness fix during full run:

- The first full Pareto run exposed duplicate sequences in survivor selection, inflating the reported frontier size.
- Fixed `src/pareto_search.py` to deduplicate candidates by sequence before Pareto operations and to inject valid immigrant mutants when diversity collapses.
- Added a regression test to ensure final populations and frontiers contain unique sequences.
- Regression and full suite passed after the fix.

Completed:

- Preflight suite passed: 33 tests.
- Ran:
  `PYTHONDONTWRITEBYTECODE=1 python experiments/02_nsga2_pareto.py --sequence data/sequences/matrixyl_core.fasta --population 100 --generations 100 --tournament-size 2 --mutation-rate 0.2 --seed 42 --max-edit-distance 2 --output results/phase4_pareto_full`
- Output files were written under `results/phase4_pareto_full`.
- The Pareto frontier figure was written successfully.
- Wrote exhaustive validation outputs:
  - `results/phase4_pareto_full/exhaustive_pareto_frontier.csv`
  - `results/phase4_pareto_full/frontier_comparison.json`
- Post-run suite passed: 34 tests.
- No linter diagnostics were reported for the phase files.

Results:

- Exhaustive candidate count: 3,706.
- True exhaustive Pareto frontier size: 9.
- Evolved Pareto frontier size: 9.
- Evolved frontier precision against exhaustive frontier: 1.0.
- Exhaustive frontier recall: 1.0.
- Missing true frontier sequences: none.
- Extra evolved frontier sequences: none.
- Phase 3 scalar optimum `PTTPS` appears on both the evolved and exhaustive Pareto frontiers.
- Final frontier sequences:
  - `PTTPS`
  - `KTTPP`
  - `KTTKS`
  - `KTTPS`
  - `KTTPA`
  - `KTTKA`
  - `KTTKP`
  - `KTTKG`
  - `KTTPG`

Interpretation:

- The Phase 4 evolutionary Pareto optimizer recovered the exact true frontier for the current fixed-length max-edit-2 search space.
- `KTTKS` remains the function-preserving anchor, while `PTTPS` is the strongest high-penetration trade-off under current deterministic objectives.
- The frontier provides a better paper story than the Phase 3 scalar winner alone because it separates permeability gain from motif preservation.
- The next work should move to Phase 5 analysis: sequence logos/mutation enrichment, candidate table, and paper-ready figures comparing Matrixyl core, Pal-KTTKS, the Phase 3 scalar optimum, and the Phase 4 Pareto frontier.

### Phase 5: Candidate Analysis and Interpretability

Goal: turn optimized outputs into scientific claims.

Status: implemented and run on the exact Phase 4 Pareto frontier. Phase 5 converts the frontier into mutation-enrichment tables, position-frequency summaries, baseline comparisons, candidate priority labels, and paper-facing figures.

Research notes:

- Sequence logos and position-frequency matrices are standard tools for visualizing protein and peptide motifs.
- Enrichment/depletion views are useful when optimized variants share a reference motif and the goal is to explain which positions changed.
- Pareto-front candidate prioritization should present raw descriptors and trade-off objectives, not just a single scalar rank.

Tasks:

- Create sequence logos for top candidates. Completed as position-frequency tables and a lightweight heatmap.
- Run mutation enrichment analysis by position. Completed.
- Perform in silico alanine or conservative-substitution scans.
- Compare candidates against Matrixyl core, palmitoylated Matrixyl, random mutants, and collagen-like controls. Completed for Matrixyl core and Pal-KTTKS; random/control comparisons remain a follow-up analysis.
- Flag candidates likely to be too polar, too hydrophobic, too large, too charged, or hard to synthesize. Completed via candidate warning flags.

Acceptance criteria:

- Candidate selection rationale is reproducible.
- The top 5 to 10 candidates have a concise explanation grounded in raw descriptors.
- The paper can explain which mutations drove predicted permeability gains.

Implementation gate:

- Add analysis helpers for loading candidate CSVs, position-frequency matrices, mutation-enrichment rows, descriptor warning flags, and candidate priority labels.
- Add unit tests for every helper.
- Add `experiments/03_candidate_analysis.py` to read `results/phase4_pareto_full/pareto_frontier.csv` and write analysis artifacts.
- Generate tables and figures under `results/phase5_candidate_analysis`.
- Run full tests before and after the analysis script.

Completed:

- Added `src/analysis.py` with helpers for loading candidate rows, position-frequency matrices, mutation enrichment, descriptor warning flags, priority labels, candidate summaries, baseline descriptor rows, and CSV writing.
- Added `tests/test_analysis.py` covering analysis helper behavior.
- Added `experiments/03_candidate_analysis.py`.
- Ran:
  `PYTHONDONTWRITEBYTECODE=1 python experiments/03_candidate_analysis.py --frontier results/phase4_pareto_full/pareto_frontier.csv --output results/phase5_candidate_analysis --palmitoylated-smiles data/molecules/matrixyl_palmitoylated.smi`
- Generated:
  - `results/phase5_candidate_analysis/candidate_summary.csv`
  - `results/phase5_candidate_analysis/position_frequency_matrix.csv`
  - `results/phase5_candidate_analysis/mutation_enrichment.csv`
  - `results/phase5_candidate_analysis/baseline_comparison.csv`
  - `results/phase5_candidate_analysis/figures/position_frequency_heatmap.png`
  - `results/phase5_candidate_analysis/figures/candidate_tradeoffs.png`
  - `results/phase5_candidate_analysis/run_metadata.json`

Verification:

- Focused analysis tests passed.
- Full suite passed with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests`.
- Final suite status after Phase 5 implementation: 40 tests passed.
- No linter diagnostics were reported for `src/analysis.py`, `tests/test_analysis.py`, `experiments/03_candidate_analysis.py`, or `PLAN.md`.

Results:

- Candidates analyzed: 9.
- Mutation enrichment rows: 5.
- Dominant frontier mutation: `K4P`, present in 5 of 9 candidates.
- Position 5 substitutions split across `S5A`, `S5G`, and `S5P`, each present in 2 of 9 candidates.
- `K1P` appears only in `PTTPS`, the high-penetration trade-off candidate.
- Priority labels:
  - High-penetration trade-offs: `PTTPS`, `KTTPP`.
  - Balanced candidates: `KTTPS`, `KTTKP`.
  - Function anchor: `KTTKS`.
- Baseline comparison:
  - Matrixyl core `KTTKS`: penetration score `0.3949`, functional score `1.0`.
  - Pal-KTTKS: penetration score `0.5077`, functional score `1.0`.
- Pal-KTTKS improves the descriptor score relative to unmodified `KTTKS` mainly through LogP, while retaining very high MW and rotatable-bond burden.

Interpretation:

- The current frontier suggests that replacing the fourth-position lysine with proline is the strongest recurring move for lowering descriptor penalties while retaining a Matrixyl-like short motif.
- `PTTPS` is a clear high-penetration computational trade-off, but it has only moderate functional-preservation score and should not be presented as a biological winner without validation.
- `KTTPS` and `KTTKP` are more conservative follow-up candidates because they preserve more of the Matrixyl motif while improving descriptor score.
- The next step is either candidate selection for structural/conformer validation or adding random/control comparison plots before drafting result figures.

### Phase 6: Optional Structural and MD Validation

Goal: add credibility without overbuilding.

Status: completed for exploratory conformer validation. Full MD remains optional future work.

Research notes:

- RDKit ETKDG is a standard conformer-generation method and ETKDGv3 is commonly used for more complex flexible molecules.
- Radius of gyration is a useful exploratory compactness descriptor because it approximates molecular size and diffusion-relevant spatial extent.
- Flexible peptide conformer predictions are limited by force-field and sampling accuracy, so these metrics should support exploratory ranking only.

Tasks:

- Generate conformer ensembles for top candidates. Completed with RDKit ETKDG.
- Estimate radius of gyration and compactness. Completed.
- If needed, run short explicit-solvent MD or membrane-interface simulations.
- Avoid claiming receptor-level function unless a credible biological target is defined.

Acceptance criteria:

- Structural validation supports stability/compactness claims.
- MD outputs are treated as exploratory unless run length and setup are publication-grade.

Implementation gate:

- Add a conformer module that builds RDKit molecules from sequence or SMILES, embeds multiple conformers, optimizes with MMFF where possible and UFF fallback, and reports radius of gyration statistics.
- Add unit tests for conformer generation, compactness metrics, and deterministic seeds.
- Add `experiments/04_structural_validation.py` to validate frontier candidates plus Matrixyl core and Pal-KTTKS.
- Run full tests before and after the structural validation script.

Completed:

- Added `src/structure.py` with RDKit conformer ensemble generation, ETKDGv3 embedding, MMFF optimization with UFF fallback, and radius-of-gyration summary statistics.
- Added `tests/test_structure.py` covering deterministic sequence conformer summaries, invalid conformer counts, and Pal-KTTKS SMILES conformers.
- Added `experiments/04_structural_validation.py`.
- Ran:
  `PYTHONDONTWRITEBYTECODE=1 python experiments/04_structural_validation.py --candidates results/phase5_candidate_analysis/candidate_summary.csv --output results/phase6_structural_validation --top-k 5 --num-conformers 8 --palmitoylated-conformers 3 --seed 42`
- Generated:
  - `results/phase6_structural_validation/conformer_summary.csv`
  - `results/phase6_structural_validation/conformer_metrics.csv`
  - `results/phase6_structural_validation/run_metadata.json`
  - `results/phase6_structural_validation/figures/compactness_vs_penetration.png`

Results:

- `PTTPS` had the most compact mean conformer ensemble among the evaluated candidates: mean radius of gyration 4.83 vs. 5.47 for the unmodified Matrixyl core.
- `KTTPS`, the balanced one-mutation candidate, also improved compactness relative to `KTTKS`: mean radius of gyration 4.99 vs. 5.47.
- Pal-KTTKS was much larger in 3D extent, with mean radius of gyration 8.56, consistent with the lipid tail increasing spatial size and flexibility despite improving LogP.
- These results support using `PTTPS` as a high-penetration structural trade-off and `KTTPS` as the more conservative candidate for follow-up validation.

Verification:

- Focused structural tests passed with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -p 'test_structure.py'`.
- Full suite passed before and after the structural run with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests`.
- Final suite status: 43 tests passed.
- No linter diagnostics were reported for `src/structure.py`, `tests/test_structure.py`, `experiments/04_structural_validation.py`, `src/__init__.py`, or `PLAN.md`.

Limitations:

- RDKit conformer ensembles are exploratory and should not be framed as MD-level stability evidence.
- Full explicit-solvent or membrane-interface MD remains optional and should only be run if the paper needs stronger claims about conformational stability or membrane interaction.

### Phase 7: Experimental Validation Design

Goal: prepare the wet-lab bridge.

Status: completed for experimental validation planning.

Research notes:

- Franz diffusion cells are the standard in vitro format for estimating topical permeation across excised skin, reconstructed skin, or synthetic membranes.
- Peptide permeation studies need method-specific quantification, ideally LC-MS/MS or validated HPLC, plus solubility and receptor-medium stability checks.
- Collagen-stimulating claims should be tested in human dermal fibroblasts with procollagen type I, collagen I/III, or related secretion assays rather than inferred from sequence similarity alone.
- Safety gating should include fibroblast and keratinocyte viability, LDH or membrane-integrity readouts, and reconstructed human epidermis irritation before any strong cosmetic claim.

Tasks:

- Select 2 to 3 synthesis candidates plus controls. Completed.
- Propose Franz diffusion cell protocol. Completed.
- Include vehicle/formulation assumptions. Completed.
- Include fibroblast collagen I/III or procollagen assay design. Completed.
- Include cytotoxicity and irritation controls. Completed.

Acceptance criteria:

- The paper has a credible "future experimental validation" section.
- Candidate selection is realistic for synthesis and topical testing.

Implementation gate:

- Add a validation-design helper module that joins Phase 5 candidate scores with Phase 6 compactness and selects a small synthesis panel.
- Generate structured outputs: synthesis candidate table, assay matrix, experimental protocol draft, and run metadata.
- Add focused tests for candidate selection, controls, assay matrix completeness, and protocol content.
- Run full tests before and after generating the Phase 7 design artifacts.

Completed:

- Added `src/experimental_design.py` with helpers for joining candidate and compactness outputs, selecting a synthesis panel, defining controls, defining the assay matrix, and drafting protocol text.
- Added `tests/test_experimental_design.py` covering candidate selection, invalid panel size, required controls, assay coverage, and protocol guardrails.
- Added `experiments/05_experimental_design.py`.
- Ran:
  `PYTHONDONTWRITEBYTECODE=1 python experiments/05_experimental_design.py --candidates results/phase5_candidate_analysis/candidate_summary.csv --compactness results/phase6_structural_validation/conformer_summary.csv --output results/phase7_experimental_design --max-candidates 3`
- Generated:
  - `results/phase7_experimental_design/synthesis_candidates.csv`
  - `results/phase7_experimental_design/experimental_controls.csv`
  - `results/phase7_experimental_design/assay_matrix.csv`
  - `results/phase7_experimental_design/protocol_draft.md`
  - `results/phase7_experimental_design/run_metadata.json`

Selected synthesis panel:

- `PTTPS`: high-penetration trade-off, highest predicted permeability score, mean conformer radius of gyration 4.83.
- `KTTPS`: balanced conservative one-mutation analog, stronger motif preservation than the high-penetration trade-offs, mean conformer radius of gyration 4.99.
- `KTTPP`: backup high-penetration analog with strong scalar trade-off value, mean conformer radius of gyration 5.31.

Controls:

- Vehicle blank.
- `KTTKS` Matrixyl core.
- Pal-KTTKS commercial/lipidated Matrixyl comparator.
- `GPKGDPGA` collagen-like sequence control.
- Positive collagen assay control such as ascorbic acid or TGF-beta, assay-dependent.
- Cytotoxicity positive control selected by the testing lab.

Assay design:

- Franz diffusion permeation with cumulative receptor recovery, flux, skin retention, and mass balance.
- Vehicle and receptor compatibility checks for solubility, recovery, and stability.
- Human dermal fibroblast collagen/procollagen assay using non-toxic concentrations.
- Fibroblast and keratinocyte cytotoxicity using MTT/resazurin plus LDH release.
- Reconstructed human epidermis irritation for finalist safety positioning.

Verification:

- Focused Phase 7 tests passed with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -p 'test_experimental_design.py'`.
- Full suite passed before and after design generation with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests`.
- Final suite status: 48 tests passed.
- No linter diagnostics were reported for `src/experimental_design.py`, `tests/test_experimental_design.py`, `experiments/05_experimental_design.py`, `src/__init__.py`, or `PLAN.md`.

Limitations:

- Phase 7 outputs are a proposed validation plan, not experimental evidence.
- The paper should clearly separate computational candidate prioritization from wet-lab confirmation of permeation, collagen/procollagen activity, and safety.

### Phase 8: PLM Oracle and GPU Run Preparation

Goal: make the repository ready for a real GPU-backed PLM scoring run.

Status: completed for GPU-run preparation. The actual PLM scoring run is ready but not launched on this non-CUDA machine.

Research notes:

- ProtBERT expects uppercase amino acids with spaces between residues and uncommon amino acids such as `U`, `Z`, `O`, and `B` replaced with `X`.
- ESM-2 Hugging Face models can consume protein sequence strings directly; batching should still be explicit and conservative because padding and memory use dominate GPU throughput.
- The previous oracle's `GPKGDP` reference is not appropriate for Matrixyl-family functional preservation; Phase 8 should default to `KTTKS` and allow explicit reference override.
- The fixed Phase 4 search space has only 3,706 valid max-edit-2 candidates, so a cached full-space PLM scoring run is feasible on a GPU once the oracle is hardened.

Tasks:

- Replace eager PLM imports with lazy dependency loading so tests and descriptor code still work without GPU dependencies loaded. Completed.
- Add model-aware sequence formatting for ESM-2 and ProtBERT. Completed.
- Add deterministic embedding caching to support resumable GPU scoring. Completed.
- Add a CPU-safe preflight script that reports CUDA availability, model configuration, candidate counts, and output paths without loading all weights unless requested. Completed.
- Add a GPU-run scoring script that reads candidate sequences, scores them against the `KTTKS` reference, and writes PLM scores/uncertainties. Completed.
- Add unit tests for oracle prep without downloading or loading real models. Completed.

Acceptance criteria:

- Full unit suite passes without loading Hugging Face weights.
- `src.oracle.EnsembleOracle` is ready to load real models when run on a GPU host.
- A preflight command can be run locally to verify candidate inputs and expected GPU-run configuration.
- The actual GPU command is documented and produces a reproducible output path.

Completed:

- Rebuilt `src/oracle.py` so heavy `torch` and `transformers` imports are lazy and only occur when `EnsembleOracle` is instantiated.
- Changed the default PLM reference sequence from the incorrect `GPKGDP` proxy to the Matrixyl-family core `KTTKS`.
- Added `ProteinModelSpec`, model-aware input formatting, ProtBERT spaced-token formatting, uncommon amino-acid replacement for PLM inputs, cosine similarity scoring, and `EmbeddingCache`.
- Added `src/plm_pipeline.py` for CPU-safe torch/CUDA detection, candidate manifest generation, fixed-length search-space enumeration, and preflight report construction.
- Added `tests/test_oracle_prep.py` and `tests/test_plm_pipeline.py`.
- Added `experiments/06_plm_gpu_preflight.py` and `experiments/07_score_plm_oracle.py`.

Audit hardening:

- Reviewed the implementation up to the GPU step for import safety, scoring correctness, candidate enumeration, output consistency, and GPU execution readiness.
- Updated PLM embedding pooling to exclude tokenizer special tokens from the mean-pooled sequence embedding.
- Updated `experiments/07_score_plm_oracle.py` so models are scored sequentially by default, reducing GPU memory pressure versus loading ESM-2 and ProtBERT together. `--load-models-together` remains available for larger GPUs.
- Added a runtime invariant that rejects model batches if the number of returned embeddings does not match the number of input sequences.
- Added tests covering special-token pooling and embedding batch-size mismatch failures.

Prepared GPU inputs:

- Ran:
  `PYTHONDONTWRITEBYTECODE=1 python experiments/06_plm_gpu_preflight.py --enumerate-search-space --output results/phase8_plm_preflight --cache-dir results/phase8_plm_cache --batch-size 16 --reference-sequence KTTKS`
- Generated:
  - `results/phase8_plm_preflight/candidate_manifest.csv`
  - `results/phase8_plm_preflight/preflight_report.json`
- Candidate manifest contains all 3,706 fixed-length max-edit-2 Matrixyl-family sequences.
- Preflight detected `torch` but no CUDA device on the current machine, so the heavy model run was not launched locally.

GPU command to run on a CUDA host:

```bash
PYTHONDONTWRITEBYTECODE=1 python experiments/07_score_plm_oracle.py \
  --candidates results/phase8_plm_preflight/candidate_manifest.csv \
  --output results/phase8_plm_preflight/gpu_scores \
  --cache-dir results/phase8_plm_cache \
  --batch-size 16 \
  --reference-sequence KTTKS \
  --require-cuda \
  --models facebook/esm2_t33_650M_UR50D \
  --models rostlab/prot_bert
```

Expected GPU outputs:

- `results/phase8_plm_preflight/gpu_scores/plm_scores.csv`
- `results/phase8_plm_preflight/gpu_scores/run_metadata.json`
- Per-model embedding cache files under `results/phase8_plm_cache`.

Verification:

- Focused oracle prep tests passed with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -p 'test_oracle_prep.py'`.
- Focused PLM pipeline tests passed with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -p 'test_plm_pipeline.py'`.
- Full suite passed before and after preflight with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests`.
- Final audit suite status: 60 tests passed.
- No linter diagnostics were reported for `src/oracle.py`, `src/plm_pipeline.py`, `tests/test_oracle_prep.py`, `tests/test_plm_pipeline.py`, `experiments/06_plm_gpu_preflight.py`, `experiments/07_score_plm_oracle.py`, `src/__init__.py`, or `PLAN.md`.

GPU readiness checklist:

- Candidate manifest has 3,706 rows and 3,706 unique sequences.
- Required anchor and shortlisted sequences are present: `KTTKS`, `PTTPS`, `KTTPS`, and `KTTPP`.
- Preflight report sequence count matches the manifest.
- Preflight command includes `--require-cuda`, so accidental CPU execution of the heavy run is blocked.
- Current local machine has `torch` available but no CUDA device; run the scoring command on a CUDA host.

Next after GPU scores:

- Join `plm_scores.csv` back to the Phase 4/5 candidate tables.
- Re-rank the Pareto frontier using descriptor score, motif score, PLM similarity, and PLM uncertainty.
- Decide whether the synthesis panel remains `PTTPS`, `KTTPS`, `KTTPP` or changes after real-model evidence.

## 9. Paper Plan

### Working Title

Computational Evolution of Collagen-Stimulating Peptides for Enhanced Epidermal Permeability

### Abstract Claim

A constrained evolutionary design workflow can identify Matrixyl-family analogs with improved permeability-associated physicochemical profiles while preserving motif-level similarity to the collagen-stimulating parent peptide.

### Introduction

Key points:

- Cosmetic peptides are commercially successful but delivery-limited.
- The stratum corneum imposes strong constraints on polarity, size, charge, and flexibility.
- Matrixyl is a recognizable baseline for collagen-stimulating peptide design.
- Prior peptide-design workflows often optimize function first and delivery second.
- This paper frames delivery as a first-class optimization objective.

### Methods

Sections:

- Baseline peptide and modification states.
- Molecular descriptor calculation.
- Permeability score construction.
- Functional-preservation oracle.
- Evolutionary optimization algorithms.
- Candidate filtering and ranking.
- Statistical and random-baseline comparisons.
- Optional conformer or MD validation.

### Results

Planned results:

- Baseline descriptor audit shows why unmodified and palmitoylated Matrixyl pose different delivery problems.
- Tournament search improves descriptor score over random mutants.
- NSGA-II reveals a Pareto frontier between permeability and functional preservation.
- Mutation analysis identifies positions or residue classes responsible for improved delivery descriptors.
- Candidate shortlist balances permeability, motif preservation, and synthesis feasibility.

### Discussion

Key discussion points:

- The workflow is useful because topical peptides need formulation-aware design.
- PLM and motif scores are computational proxies, not proof of biological activity.
- Peptide permeability rules are imperfect, especially for lipidated and formulated ingredients.
- The approach is most valuable as a candidate-prioritization engine before synthesis.

### Limitations

Include explicitly:

- Descriptor thresholds come from small-molecule and transdermal literature and may not fully transfer to peptides.
- PLM embedding similarity is not a direct collagen-stimulation assay.
- Palmitoylation changes skin partitioning and cannot be represented by sequence alone.
- Experimental permeability and fibroblast-response assays remain necessary.

## 10. Immediate Next Steps

Recommended order:

1. Correct baseline sequence data and README claims.
2. Add tests and fix package imports so lightweight modules work without heavy optional dependencies.
3. Implement exact RDKit descriptor calculation for `KTTKS` and `Pal-KTTKS`.
4. Define the candidate schema and scoring contract.
5. Implement tournament-search baseline.
6. Fix NSGA-II objective signs and sequence encoding.
7. Generate first baseline descriptor table and first optimization run.
8. Draft paper outline with placeholders for figures.

## 11. Engineering Risks

High-priority risks:

- Wrong baseline sequence will invalidate the study if not fixed.
- Heuristic descriptors are not strong enough for publication claims.
- Positive objectives in `pymoo` currently optimize in the wrong direction.
- Package import coupling blocks lightweight usage.
- PLM similarity to `GPKGDP` is not a credible Matrixyl functional oracle.

Medium-priority risks:

- RDKit peptide construction with modifications may require careful validation.
- Very short peptides can make PLM embeddings noisy or biologically ambiguous.
- Optimizing for passive diffusion may fight the known reason Matrixyl is palmitoylated.
- Candidate search may collapse to trivial small hydrophobic sequences unless motif constraints are explicit.

## 12. Definition of Done for a First Preprint-Ready Computational Study

The repository is ready to support a first computational manuscript when:

- Baselines are scientifically correct and cited.
- Descriptor calculations are exact, tested, and reproducible.
- Search algorithms run from clean setup commands.
- Random and baseline controls are included.
- Candidate rankings include raw descriptors, normalized objectives, and function-preservation scores.
- Figures are generated by scripts, not hand-made.
- Limitations are explicit and avoid overclaiming biological activity.
- The top candidates are plausible for synthesis and experimental follow-up.

## 13. Implementation Log

### 2026-05-05: Phase 0 and Phase 1 Foundation

Completed:

- Corrected the Matrixyl baseline from the collagen-like `GPKGDPGA` motif to the Matrixyl peptide core `KTTKS`.
- Preserved `GPKGDPGA` as `data/sequences/collagen_like_control.fasta`.
- Added `data/sequences/README.md` to distinguish sequence-only `KTTKS` from lipidated `Pal-KTTKS`.
- Made package-root imports lazy for optional heavy modules so lightweight descriptor code does not require optimizer dependencies.
- Added RDKit-backed descriptor calculation in `src/chemistry.py`.
- Added `data/molecules/matrixyl_palmitoylated.smi` and PubChem reference targets for palmitoyl pentapeptide-4.
- Added unit tests for baseline sequences, lightweight imports, and RDKit descriptors.

Verification:

- `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests` passed with 7 tests.
- RDKit was installed from the project pin, and NumPy was aligned to the project pin to avoid RDKit ABI warnings.

### 2026-05-05: Phase 2 Scoring Contract

Research notes:

- Permeability scoring should use raw RDKit descriptors as design pressures, not as a direct biological permeability claim.
- Literature and descriptor APIs support using molecular weight, LogP, TPSA, H-bond donors/acceptors, charge, and rotatable bonds as the first scoring layer.
- Peptide and lipidated-peptide caveats must remain explicit because small-molecule rules such as MW < 500 Da and TPSA < 140 A^2 are not absolute for formulated topical peptides.

Completed:

- Added an exact descriptor-based `PenetrationScore` contract in `src/constraints.py`.
- Updated `compute_penetration_score`, `compute_molecular_weight`, `compute_tpsa`, and `compute_logp` to use RDKit descriptors for canonical peptide sequences.
- Kept the old residue-table score as `compute_heuristic_penetration_score` for prefiltering only.
- Added tests proving exact scoring uses RDKit descriptors, exposes component penalties, and evaluates `KTTKS` and `Pal-KTTKS` through the same scoring path.

Verification:

- `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests` passed with 12 tests.
- No linter diagnostics were reported for `src/constraints.py`, `src/__init__.py`, `tests/test_penetration_scoring.py`, or `PLAN.md`.

Remaining Phase 2 work:

- Define sequence search-space rules around `KTTKS`: locked positions, soft-lock positions, allowed substitutions, length bounds, and optional modifications.
- Add functional-preservation scoring before running evolution so optimization cannot collapse to trivial permeable but nonfunctional sequences.
- Add a shared candidate schema that combines raw descriptors, normalized scores, hard-filter flags, and functional-preservation metadata.

### 2026-05-05: Phase 2 Search Space and Functional Preservation - In Progress

Research notes:

- `KTTKS` is the Matrixyl-family matrikine core; palmitoylation should be represented as an N-terminal chemical modification, not as a FASTA character.
- Initial evolution should stay conservative: fixed-length `KTTKS`-family sequences, limited edit distance, canonical amino acids, and optional locked positions.
- BLOSUM62-style conservative-substitution scoring is an appropriate transparent proxy for motif preservation until experimental labels or a validated collagen-stimulation model exist.

Implementation gate:

- Add a reusable search-space validator that reports hard-filter failures instead of silently rejecting sequences.
- Add functional-preservation scoring with identity, edit distance, length consistency, and conservative substitution components.
- Add a shared candidate-evaluation schema combining search validity, RDKit permeability score, functional-preservation score, and an optimization-ready scalar.
- Run unit tests after each module is added and update this log with results.

Completed:

- Added `src/search_space.py` with fixed-length Matrixyl-family defaults, canonical amino-acid validation, max edit-distance filtering, optional locked positions, and explicit failure reasons.
- Added `src/function_scores.py` with `FunctionalPreservationScorer`, combining sequence identity, edit distance, length consistency, and normalized BLOSUM62 conservative-substitution similarity.
- Added `src/candidates.py` with `CandidateEvaluator`, which evaluates every sequence through search validity, functional-preservation scoring, exact RDKit permeability scoring, and a simple optimization scalar.
- Exported the new Phase 2 APIs from `src/__init__.py`.
- Added unit tests for search-space validation, functional-preservation scoring, and candidate evaluation.

Verification:

- Focused tests passed after each module was added.
- Full suite passed with `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests`.
- Final Phase 2 suite status: 23 tests passed.
- No linter diagnostics were reported for new Phase 2 modules, tests, or `PLAN.md`.

Phase 3 entry point:

- The first evolutionary baseline should use `CandidateEvaluator.evaluate_sequence()` as its only scoring path.
- Start with fixed-length `KTTKS`-family sequences and max edit distance 2 from the reference.
- Treat `optimization_score = penetration_score * functional_preservation_score` as the scalar objective for the tournament-search baseline.
- Keep PLM oracle integration out of the first baseline run until the deterministic evolutionary loop and random-mutant comparison are tested.

