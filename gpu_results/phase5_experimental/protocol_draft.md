# Phase 7 Experimental Validation Design

## Synthesis Panel

- PTTPS: high_penetration_tradeoff.
  Rationale: Highest predicted permeability score; tests whether aggressive motif editing translates to better delivery.
  Scores: penetration 0.687, function 0.667, mean Rg 5.05.
- KTTPS: balanced_conservative.
  Rationale: One-mutation candidate with stronger motif preservation and improved compactness relative to KTTKS.
  Scores: penetration 0.546, function 0.833, mean Rg 5.08.
- KTTPP: backup_high_penetration.
  Rationale: Backup analog with strong computed trade-off value if synthesis, solubility, or assay behavior limits the primary candidates.
  Scores: penetration 0.678, function 0.670, mean Rg 5.25.

## Core Experimental Flow

1. Confirm synthesis feasibility, purity, identity, aqueous/formulation solubility, and short-term stability for each candidate.
2. Run Franz diffusion permeation with KTTKS, Pal-KTTKS, vehicle blank, and selected analogs; quantify peptide by validated LC-MS/MS or HPLC.
3. Measure fibroblast collagen/procollagen response only at concentrations shown to be non-toxic.
4. Gate cosmetic safety with fibroblast/keratinocyte viability, LDH release, and reconstructed human epidermis irritation testing.

## Vehicle Assumptions

- Start with a simple aqueous or hydroalcoholic screening vehicle compatible with peptide stability and analytical recovery.
- Treat final cream or serum formulation as a second-stage experiment because excipients can dominate skin delivery.
- Maintain sink conditions in the receptor phase and verify mass balance at the end of each Franz diffusion run.

## Reporting Guidance

- Report permeation, skin retention, cytotoxicity, and collagen/procollagen effects separately.
- Do not claim retained biological activity from permeation alone; require fibroblast assay evidence.
- Frame RDKit conformer metrics as exploratory candidate-selection support, not MD-grade stability validation.
