# Computational Evolution of Collagen-Stimulating Peptides for Enhanced Epidermal Permeability

This repository supports a computational bio-research paper on Matrixyl-family peptide optimization for topical delivery. The current implementation builds a reproducible pipeline from corrected baselines through RDKit descriptor scoring, evolutionary search, Pareto analysis, structural/conformer checks, experimental validation planning, and GPU-ready protein language model (PLM) preflight.

The actual GPU PLM scoring run has not been executed in this local environment because no CUDA device is available. The repository is prepared for that run.

## Scientific Framing

The primary baseline is the Matrixyl / palmitoyl pentapeptide-4 peptide core:

- Matrixyl core: `KTTKS`
- Lipidated commercial comparator: `Pal-KTTKS`, represented by SMILES in `data/molecules/matrixyl_palmitoylated.smi`
- Collagen-like control retained separately: `GPKGDPGA`

The project optimizes permeability-associated descriptors while preserving motif-level similarity to the Matrixyl core. The current functional score is a computational preservation proxy, not a direct collagen-stimulation assay.

## Current Status

Completed:

- Corrected baseline sequence data from the original collagen-like fragment to `KTTKS`.
- Added RDKit-backed descriptor calculation for peptide FASTA and SMILES inputs.
- Added exact penetration scoring using molecular weight, TPSA, LogP, hydrogen bonding, rotatable bonds, and formal charge.
- Added fixed-length Matrixyl-family search-space validation with max edit distance 2.
- Added transparent functional-preservation scoring using identity, edit distance, length consistency, and BLOSUM62-style substitution similarity.
- Added deterministic tournament search and dependency-light Pareto optimization.
- Added candidate analysis, mutation enrichment, baseline comparison, and figure generation.
- Added RDKit conformer ensemble validation with radius-of-gyration compactness summaries.
- Added experimental validation design artifacts for synthesis candidates, controls, Franz diffusion, fibroblast collagen/procollagen, cytotoxicity, and irritation assays.
- Hardened the PLM oracle for a real GPU run with ESM-2 and ProtBERT, batching, model-aware input formatting, sequential model loading, special-token-safe pooling, and embedding caching.
- Added 60 unit tests covering the implemented pipeline.

Pending:

- Run the real GPU-backed PLM scoring job.
- Join PLM scores back into the candidate tables and re-rank the Pareto frontier.
- Decide whether the synthesis panel changes after PLM evidence.
- Execute wet-lab validation if moving beyond a computational manuscript.

## Optimization Targets

The descriptor score treats these as design pressures rather than absolute biological laws:

- TPSA: target below 140 A^2.
- Molecular weight: target below 500 Da for passive diffusion, with separate interpretation for lipidated `Pal-KTTKS`.
- LogP: target 1 to 3 for topical/transdermal balance.
- Hydrogen-bond donors and acceptors: minimized while preserving peptide function.
- Formal charge: near neutral preferred.
- Rotatable bonds and radius of gyration: lower values preferred as compactness/flexibility proxies.
- Functional preservation: motif similarity to `KTTKS`, plus later PLM similarity once GPU scoring is run.

## Repository Layout

```text
Bio-paper/
├── PLAN.md
├── README.md
├── data/
│   ├── molecules/
│   │   └── matrixyl_palmitoylated.smi
│   ├── references/
│   │   └── palmitoyl_pentapeptide_4_pubchem.json
│   └── sequences/
│       ├── baseline_matrixyl.fasta
│       ├── collagen_like_control.fasta
│       ├── matrixyl_core.fasta
│       └── README.md
├── experiments/
│   ├── 01_tournament_search.py
│   ├── 02_nsga2_pareto.py
│   ├── 03_candidate_analysis.py
│   ├── 04_structural_validation.py
│   ├── 05_experimental_design.py
│   ├── 06_plm_gpu_preflight.py
│   └── 07_score_plm_oracle.py
├── src/
│   ├── analysis.py
│   ├── candidates.py
│   ├── chemistry.py
│   ├── constraints.py
│   ├── experimental_design.py
│   ├── function_scores.py
│   ├── oracle.py
│   ├── pareto_search.py
│   ├── plm_pipeline.py
│   ├── search_space.py
│   ├── structure.py
│   └── tournament_search.py
└── tests/
```

`results/` is generated output and is not part of the source tree.

## Setup

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

For the GPU PLM run, install a CUDA-compatible PyTorch build for the target machine. First-time Hugging Face model downloads require network access unless model snapshots are already cached and `--local-files-only` is used.

## Verification

Run the full unit suite:

```bash
python -m pytest tests/ -v
```

Current local audit status:

- 60 tests passing.
- CPU-safe PLM preflight passes.
- Local machine has `torch` installed but no CUDA device.

## Pipeline Commands

### Phase 3: Tournament Search

```bash
python experiments/01_tournament_search.py \
  --sequence data/sequences/matrixyl_core.fasta \
  --population 100 \
  --generations 100 \
  --tournament-size 3 \
  --mutation-rate 0.2 \
  --elite-count 2 \
  --seed 42 \
  --max-edit-distance 2 \
  --top-k 25 \
  --output results/phase3_tournament_extensive
```

### Phase 4: Pareto Search

```bash
python experiments/02_nsga2_pareto.py \
  --output results/phase4_pareto_full
```

### Phase 5: Candidate Analysis

```bash
python experiments/03_candidate_analysis.py \
  --frontier results/phase4_pareto_full/pareto_frontier.csv \
  --output results/phase5_candidate_analysis \
  --palmitoylated-smiles data/molecules/matrixyl_palmitoylated.smi
```

### Phase 6: Structural Validation

```bash
python experiments/04_structural_validation.py \
  --candidates results/phase5_candidate_analysis/candidate_summary.csv \
  --output results/phase6_structural_validation \
  --top-k 5 \
  --num-conformers 8 \
  --palmitoylated-conformers 3 \
  --seed 42
```

### Phase 7: Experimental Design

```bash
python experiments/05_experimental_design.py \
  --candidates results/phase5_candidate_analysis/candidate_summary.csv \
  --compactness results/phase6_structural_validation/conformer_summary.csv \
  --output results/phase7_experimental_design \
  --max-candidates 3
```

### Phase 8: GPU Preflight

This prepares the full fixed-length, max-edit-2 Matrixyl candidate manifest without loading model weights:

```bash
python experiments/06_plm_gpu_preflight.py \
  --enumerate-search-space \
  --output results/phase8_plm_preflight \
  --cache-dir results/phase8_plm_cache \
  --batch-size 16 \
  --reference-sequence KTTKS
```

Expected generated files:

- `results/phase8_plm_preflight/candidate_manifest.csv`
- `results/phase8_plm_preflight/preflight_report.json`

The manifest should contain 3,706 unique candidate sequences.

### Phase 8: Real GPU PLM Scoring

Run this on a CUDA host:

```bash
python experiments/07_score_plm_oracle.py \
  --candidates results/phase8_plm_preflight/candidate_manifest.csv \
  --output results/phase8_plm_preflight/gpu_scores \
  --cache-dir results/phase8_plm_cache \
  --batch-size 16 \
  --reference-sequence KTTKS \
  --require-cuda \
  --models facebook/esm2_t33_650M_UR50D \
  --models rostlab/prot_bert
```

Expected generated files:

- `results/phase8_plm_preflight/gpu_scores/plm_scores.csv`
- `results/phase8_plm_preflight/gpu_scores/run_metadata.json`
- Per-model embedding cache files in `results/phase8_plm_cache`

By default, the GPU scorer loads models sequentially to reduce VRAM pressure. Use `--load-models-together` only on a host with enough GPU memory.

## Current Candidate Shortlist

The deterministic pipeline currently nominates:

- `PTTPS`: high-penetration trade-off.
- `KTTPS`: balanced conservative one-mutation analog.
- `KTTPP`: backup high-penetration analog.

This shortlist should be revisited after PLM scoring is joined back into the Phase 4/5 candidate tables.

## Experimental Validation Plan

The generated wet-lab planning artifacts include:

- Synthesis candidates and controls.
- Franz diffusion permeation design.
- Vehicle and receptor compatibility checks.
- Human dermal fibroblast collagen/procollagen assay.
- Fibroblast and keratinocyte cytotoxicity.
- Reconstructed human epidermis irritation testing.

These are proposed validation designs, not experimental results.

## Important Limitations

- Descriptor thresholds come from small-molecule and transdermal literature and may not fully transfer to short peptides.
- PLM embedding similarity is a proxy for sequence-level preservation, not direct collagen stimulation.
- Palmitoylation changes partitioning, formulation behavior, and 3D extent; it cannot be represented by sequence alone.
- RDKit conformer ensembles are exploratory and should not be framed as MD-level stability evidence.
- Experimental permeability and fibroblast-response assays are required before biological or commercial claims.

## References to Cite

- Deb et al., NSGA-II: A fast and elitist multiobjective genetic algorithm (2002).
- Lin et al., ESM-2 language models of protein sequences at evolutionary scale (2023).
- Potts and Guy, Predicting skin permeability (1992).
- Mitragotri et al., Transdermal drug delivery (2004).
- Blanes-Mira et al., Matrixyl / palmitoyl pentapeptide activity literature (2002).

---

Status: GPU-ready computational pipeline; real PLM scoring pending on CUDA hardware.  
Last updated: 2026-05-05
