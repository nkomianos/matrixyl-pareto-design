# Sequence Baselines

This directory separates the commercial Matrixyl-family baseline from collagen-like control motifs.

- `baseline_matrixyl.fasta` and `matrixyl_core.fasta` contain the peptide core `KTTKS`, commonly associated with Matrixyl / palmitoyl pentapeptide-4 before adding the N-terminal palmitoyl group.
- `collagen_like_control.fasta` contains `GPKGDPGA`, the collagen-like motif previously stored as the Matrixyl baseline.

The commercial ingredient is lipidated (`Pal-KTTKS`). Sequence-only FASTA files cannot encode the palmitoyl chain, so molecular descriptor work should represent the lipidated form with a molecule format such as SMILES or an RDKit-built structure.
