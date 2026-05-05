# CLAUDE.md: Technical Guidance for Peptide Optimization Project

## Quick Start for Future Sessions

This document guides you through the codebase structure, architectural decisions, and common tasks.

### Project Goals
1. Optimize Matrixyl (palmitoyl pentapeptide-4) for skin penetration
2. Maintain collagen-binding function
3. Use evolutionary algorithms + protein language models
4. Validate computationally and experimentally

### Key Decisions Made Upfront

#### Why NSGA-II + CMA-ES?
- **CMA-ES**: Excellent for single-objective optimization with continuous variables (physicochemical properties). Fast convergence.
- **NSGA-II**: Better for multi-objective (penetration vs. binding). Maintains diverse Pareto frontier.
- **Hybrid strategy**: Run CMA-ES first (fast exploration), then NSGA-II for refinement.

#### Why Ensemble of PLMs?
- ESM-2, ProtBERT, OmegaFold have different training data and biases
- Voting ensemble reduces outliers
- Uncertainty (std dev across models) flags sequences needing validation

#### Physicochemical Constraints
All calculations via RDKit after converting peptide sequence to SMILES:
```python
# Example constraint calculation
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

mol = Chem.MolFromSmiles(smiles)
tpsa = Descriptors.TPSA(mol)
mw = Descriptors.MolWt(mol)
logp = Crippen.MolLogP(mol)
```

**Note**: SMILES generation from peptide sequences is lossy (doesn't capture all stereochemistry). We handle this by using PLM embeddings as primary oracle; physicochemical properties are secondary constraints.

#### Why Not Just Machine Learning?
Training a single regression model on limited peptide data → overfitting. PLMs pre-trained on millions of proteins provide much stronger generalization. Evolutionary algorithms let us search the space intelligently without needing labeled training data.

---

## Architecture Overview

### Core Modules

#### `src/oracle.py`
Wraps protein language models (ESM-2, ProtBERT, OmegaFold).
```python
class EnsembleOracle:
    def __init__(self, models=['esm2', 'protbert', 'omegafold']):
        # Load pretrained weights
        
    def score_binding_affinity(self, sequence: str) -> float:
        # Returns (score, uncertainty) tuple
        # score: ensemble mean across models
        # uncertainty: std dev (used for active learning)
        
    def validate_structure(self, sequence: str) -> bool:
        # Check secondary structure predictions
        # Reject if predicted as intrinsically disordered
```

**Weights**: Download via `scripts/download_plm_weights.py` (ESM-2 is ~1.3 GB per model).

#### `src/evolutionary_algorithm.py`
NSGA-II implementation using pymoo. Multi-objective fitness function:
```python
def fitness_function(sequence):
    penetration_score = compute_penetration_score(sequence)
    binding_score = oracle.score_binding_affinity(sequence)
    synthesis_cost = estimate_synthesis_complexity(sequence)
    
    return [penetration_score, binding_score, -synthesis_cost]
    # NSGA-II will maximize all three
```

**Mutation operators**:
- Position substitution (50%): random amino acid swap
- Insertion/deletion (20%): add or remove residue
- Structure-aware mutation (30%): swap to similar physicochemical residues

#### `src/constraints.py`
Physicochemical property calculators. All return float in [0, 1] for optimization:
```python
def penetration_score(sequence):
    tpsa = compute_tpsa(sequence)
    mw = compute_molecular_weight(sequence)
    logp = compute_logp(sequence)
    flexibility = compute_gyration_radius(sequence)
    
    # Penalty functions (minimize distance from optimal range)
    tpsa_penalty = max(0, (tpsa - 140) / 100)
    mw_penalty = max(0, (mw - 500) / 100)
    logp_penalty = min(abs(logp - 1.5), abs(logp - 3))
    
    return (1 - tpsa_penalty) * (1 - mw_penalty) * (1 - logp_penalty) * ...
```

#### `src/active_learning.py`
Identifies most informative sequences from Pareto frontier:
```python
def select_next_batch(pareto_frontier, budget=10):
    # For each sequence, compute uncertainty (std dev across ensemble)
    uncertainties = [oracle.score_binding_affinity(seq)[1] for seq in pareto_frontier]
    
    # Select sequences with highest uncertainty + high penetration
    # This targets the "risky" candidates that could be wins or losses
    return frontier[argsort(uncertainty)[-budget:]]
```

#### `src/validation.py`
MD simulation wrapper and metrics:
```python
def run_md_simulation(sequence, timesteps=100000):
    # Generate 3D structure (via OmegaFold or AlphaFold2)
    # Run GROMACS/OpenMM
    # Compute RMSD, Rg, H-bond count over time
    # Return stability metrics
    
def compare_to_baseline(optimized_seq, matrixyl_seq):
    # Side-by-side metrics comparison
    # Statistical significance testing
```

---

## Running Experiments

### Setup (One-time)
```bash
# Download model weights
python scripts/download_plm_weights.py --models esm2 protbert

# Install dependencies
pip install -r requirements.txt

# Verify setup
python -c "from src.oracle import EnsembleOracle; o = EnsembleOracle(); print('OK')"
```

### Phase 1: Baseline CMA-ES
```bash
cd experiments/
python 01_baseline_optimization.py \
    --sequence data/sequences/baseline_matrixyl.fasta \
    --generations 200 \
    --population 50 \
    --output results/phase1_cmaes
```
**Expected runtime**: 2–4 hours on GPU (8 GB VRAM).
**Outputs**: `pareto_frontier.csv`, `convergence_plot.png`

### Phase 2: NSGA-II Multi-Objective
```bash
python 02_pareto_evolution.py \
    --init_population results/phase1_cmaes/top_20.csv \
    --generations 150 \
    --output results/phase2_nsga2
```

### Phase 3: Active Learning
```bash
python 03_uncertainty_sampling.py \
    --frontier results/phase2_nsga2/pareto_frontier.csv \
    --budget 50 \
    --output results/phase3_active_learning
```

### Phase 4: Validation
```bash
python 04_validation.py \
    --candidates results/phase3_active_learning/top_5.csv \
    --baseline data/sequences/baseline_matrixyl.fasta \
    --output results/phase4_validation
```

---

## Common Tasks

### Adding a New Constraint
1. Implement calculator in `src/constraints.py`
2. Add to fitness function in `src/evolutionary_algorithm.py`
3. Test with dummy sequence first:
```python
test_seq = "MKVLWAALLVTNVLSAAPK"
print(f"New constraint: {new_constraint(test_seq)}")
```

### Debugging a Poor Convergence
- Check oracle predictions: are they stable across sequences?
  ```python
  for seq in random_sequences:
      score, unc = oracle.score_binding_affinity(seq)
      print(f"{seq}: {score:.3f} ± {unc:.3f}")
  ```
- Verify physicochemical constraints aren't too restrictive (compare ranges to literature)
- Increase population size or mutation rate in `evolutionary_algorithm.py`

### Computing Cost
- **ESM-2 inference**: ~50 ms per sequence on V100
- **NSGA-II (100 gen, pop 50)**: ~5,000 sequences evaluated = 250 seconds ≈ 4 min
- **Full pipeline**: ~30 hours CPU-equivalent on single GPU

---

## Validation & Output Standards

### Sequence Representation
- **Input**: FASTA format (canonical amino acids only)
- **Output**: CSV with columns `[sequence, penetration_score, binding_score, uncertainty, tpsa, mw, logp, flexibility]`

### Plot Standards
All figures should be publication-ready (300 dpi, colorblind-friendly palettes):
- Pareto frontier: scatter plot (penetration vs. binding), color-code by uncertainty
- Convergence: line plot with mean ± std over generations
- SAR: bar plot of residue importance across top candidates

### Reproducibility Checklist
Before committing results:
- [ ] Random seed fixed (`np.random.seed(42)`)
- [ ] All hyperparameters logged to `results/config.yaml`
- [ ] Code version (`git describe --tags`)
- [ ] Dataset checksums (MD5 of input sequences)
- [ ] Time-to-completion logged

---

## Literature & References

### Key Papers to Cite
1. **Multi-objective optimization**: Deb et al., "NSGA-II: A fast and elitist multiobjective genetic algorithm" (2002)
2. **Protein language models**: Lin et al., "Language models of protein sequences at the scale of evolution enable accurate structure prediction" (2023, ESM-2)
3. **Skin permeability**: Potts & Guy, "Predicting skin permeability" (1992); Mitragotri et al., "Transdermal drug delivery" (2004)
4. **Peptide design**: Leman et al., "A comprehensive protein language model" (2023)

### Datasets to Consider
- **UniProt**: background distribution of amino acids
- **PDBbind**: collagen-binding motif reference structures
- **DrugBank/PubChem**: MW/LogP/TPSA distributions for known drugs
- **Skin penetration literature**: QSAR models from Potts, Mitragotri, Suhnel

---

## Troubleshooting

### GPU Out of Memory
- Reduce batch size in `oracle.py` (default: 32)
- Use `torch.cuda.empty_cache()` between iterations

### Model Weights Won't Download
- Check internet connection and storage space (~5 GB total)
- Manually download from Hugging Face: `huggingface-cli download facebook/esm2_t33_650M_UR50D`

### Divergent Evolution
- Reduce mutation rate (in `evolutionary_algorithm.py`: `mutation_prob` from 0.5 to 0.3)
- Increase population size (pop 50 → 100)
- Add elitism constraint (keep top 10% unchanged each generation)

---

## Next Steps After Initial Results

1. **Experimental validation**: Submit top 3 candidates for synthesis
2. **Patent search**: Check for existing peptide patents conflicting with optimized sequences
3. **Cost analysis**: Synthesize top 1 candidate, measure cost per gram vs. original Matrixyl
4. **Publication**: Write up as 15–20 page paper with computational methods + results sections
5. **Open-source release**: Clean up code, add unit tests, publish on GitHub

---

**Document Status**: v1.0 (Initial planning)  
**Last Updated**: 2026-05-05  
**Maintainer**: (You)
