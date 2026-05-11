# Sequence Baselines

This directory separates the commercial Matrixyl-family baseline from collagen-like control motifs.

- `matrixyl_core.fasta` contains the peptide core `KTTKS`, the bioactive cargo of palmitoyl pentapeptide-4 before N-terminal lipidation.
- `collagen_like_control.fasta` contains `GPKGDPGA`, a collagen-like sequence used as a motif control.

The commercial ingredient is lipidated (`Pal-KTTKS`). Sequence-only FASTA files cannot encode the palmitoyl chain, so molecular descriptor work should represent the lipidated form with a molecule format such as SMILES or an RDKit-built structure.
