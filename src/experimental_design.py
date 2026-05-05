"""
Experimental validation design helpers for the peptide optimization study.

The goal is to produce reproducible wet-lab planning artifacts from the
computational candidate tables, not to replace protocol optimization by a lab.
"""

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional

from .analysis import CandidateRow


@dataclass(frozen=True)
class CompactnessRow:
    sequence: str
    name: str
    source: str
    mean_radius_of_gyration: float
    min_radius_of_gyration: float
    max_radius_of_gyration: float

    @classmethod
    def from_csv_row(cls, row: dict) -> "CompactnessRow":
        return cls(
            sequence=row["sequence"],
            name=row["name"],
            source=row["source"],
            mean_radius_of_gyration=float(row["mean_radius_of_gyration"]),
            min_radius_of_gyration=float(row["min_radius_of_gyration"]),
            max_radius_of_gyration=float(row["max_radius_of_gyration"]),
        )


@dataclass(frozen=True)
class SynthesisCandidate:
    sequence: str
    role: str
    rationale: str
    penetration_score: float
    functional_score: float
    scalar_score: float
    edit_distance: int
    molecular_weight: float
    tpsa: float
    logp: float
    mean_radius_of_gyration: Optional[float] = None

    def to_dict(self) -> dict:
        return asdict(self)


def load_compactness_rows(path: str | Path) -> dict[str, CompactnessRow]:
    with Path(path).open() as handle:
        rows = [CompactnessRow.from_csv_row(row) for row in csv.DictReader(handle)]
    return {row.sequence: row for row in rows if row.source == "rdkit_fasta"}


def _candidate_to_synthesis_row(
    candidate: CandidateRow,
    *,
    role: str,
    rationale: str,
    compactness: Optional[CompactnessRow],
) -> SynthesisCandidate:
    return SynthesisCandidate(
        sequence=candidate.sequence,
        role=role,
        rationale=rationale,
        penetration_score=candidate.penetration_objective,
        functional_score=candidate.functional_objective,
        scalar_score=candidate.scalar_optimization_score,
        edit_distance=candidate.edit_distance,
        molecular_weight=candidate.molecular_weight,
        tpsa=candidate.tpsa,
        logp=candidate.logp,
        mean_radius_of_gyration=(
            compactness.mean_radius_of_gyration if compactness is not None else None
        ),
    )


def select_synthesis_candidates(
    candidates: Iterable[CandidateRow],
    compactness_rows: dict[str, CompactnessRow],
    *,
    max_candidates: int = 3,
) -> list[SynthesisCandidate]:
    candidate_list = list(candidates)
    if max_candidates < 2:
        raise ValueError("max_candidates must be at least 2")

    selected: list[SynthesisCandidate] = []
    selected_sequences: set[str] = set()

    def add(candidate: CandidateRow, role: str, rationale: str) -> None:
        if candidate.sequence in selected_sequences or len(selected) >= max_candidates:
            return
        selected.append(
            _candidate_to_synthesis_row(
                candidate,
                role=role,
                rationale=rationale,
                compactness=compactness_rows.get(candidate.sequence),
            )
        )
        selected_sequences.add(candidate.sequence)

    high_penetration = max(candidate_list, key=lambda row: row.penetration_objective)
    add(
        high_penetration,
        "high_penetration_tradeoff",
        "Highest predicted permeability score; tests whether aggressive motif editing translates to better delivery.",
    )

    balanced_pool = [row for row in candidate_list if row.edit_distance <= 1]
    balanced = max(
        balanced_pool,
        key=lambda row: (row.scalar_optimization_score, row.functional_objective),
    )
    add(
        balanced,
        "balanced_conservative",
        "One-mutation candidate with stronger motif preservation and improved compactness relative to KTTKS.",
    )

    backup_pool = [row for row in candidate_list if row.sequence not in selected_sequences]
    if backup_pool and len(selected) < max_candidates:
        backup = max(
            backup_pool,
            key=lambda row: (row.scalar_optimization_score, row.penetration_objective),
        )
        add(
            backup,
            "backup_high_penetration",
            "Backup analog with strong computed trade-off value if synthesis, solubility, or assay behavior limits the primary candidates.",
        )

    return selected


def experimental_controls() -> list[dict]:
    return [
        {
            "control": "vehicle_blank",
            "sample": "Formulation vehicle without peptide",
            "purpose": "Baseline permeation, cytotoxicity, and collagen assay background.",
        },
        {
            "control": "matrixyl_core",
            "sample": "KTTKS",
            "purpose": "Unmodified peptide-core comparator for sequence-level claims.",
        },
        {
            "control": "pal_kttks",
            "sample": "Pal-KTTKS",
            "purpose": "Commercially recognizable lipidated Matrixyl comparator.",
        },
        {
            "control": "collagen_like_sequence",
            "sample": "GPKGDPGA",
            "purpose": "Non-Matrixyl collagen-like sequence control from the original repository baseline.",
        },
        {
            "control": "positive_collagen_control",
            "sample": "Ascorbic acid or TGF-beta, assay-dependent",
            "purpose": "Confirms fibroblast collagen/procollagen assay responsiveness.",
        },
        {
            "control": "cytotoxicity_control",
            "sample": "Lab-standard membrane-disruptive positive control",
            "purpose": "Confirms viability and LDH assays detect cell damage.",
        },
    ]


def assay_matrix() -> list[dict]:
    return [
        {
            "assay": "Franz diffusion permeation",
            "purpose": "Measure peptide delivery across skin or validated membrane model.",
            "samples": "Selected analogs, KTTKS, Pal-KTTKS, vehicle blank",
            "readout": "Cumulative receptor amount, flux, skin retention, mass balance",
            "method_notes": "Use validated LC-MS/MS or HPLC; confirm receptor solubility and peptide stability before the run.",
            "success_signal": "Analog exceeds KTTKS receptor recovery or skin retention without unacceptable mass-balance loss.",
        },
        {
            "assay": "Vehicle and receptor compatibility",
            "purpose": "Avoid false permeation differences from precipitation or degradation.",
            "samples": "All synthesis candidates and controls",
            "readout": "Solubility, recovery, short-term stability",
            "method_notes": "Screen aqueous buffer, PBS plus solubilizer if needed, and final cosmetic vehicle.",
            "success_signal": "Stable quantification over the planned sampling window.",
        },
        {
            "assay": "Human dermal fibroblast collagen/procollagen",
            "purpose": "Test whether analogs preserve collagen-stimulating activity.",
            "samples": "Non-toxic concentration range for each candidate and controls",
            "readout": "Procollagen type I, collagen I/III, or related secretion by ELISA or immunoassay",
            "method_notes": "Normalize to viable cell number or protein content; include positive collagen-control treatment.",
            "success_signal": "Candidate matches or exceeds KTTKS/Pal-KTTKS functional signal at non-toxic doses.",
        },
        {
            "assay": "Fibroblast and keratinocyte cytotoxicity",
            "purpose": "Gate functional interpretation on cell health.",
            "samples": "Dose range bracketing efficacy concentrations",
            "readout": "MTT/resazurin viability and LDH release",
            "method_notes": "Run before or alongside collagen assays to set maximum non-toxic concentration.",
            "success_signal": "At least 80 percent viability and no substantial LDH elevation at test concentration.",
        },
        {
            "assay": "Reconstructed human epidermis irritation",
            "purpose": "Support cosmetic safety positioning before animal or human testing.",
            "samples": "Finalists in intended vehicle plus vehicle blank",
            "readout": "Tissue viability, morphology, optional inflammatory markers",
            "method_notes": "Use a validated reconstructed epidermis model and lab-standard irritant control.",
            "success_signal": "Non-irritant classification at intended topical exposure.",
        },
    ]


def draft_protocol_markdown(
    synthesis_candidates: Iterable[SynthesisCandidate],
    *,
    title: str = "Phase 7 Experimental Validation Design",
) -> str:
    candidates = list(synthesis_candidates)
    lines = [
        f"# {title}",
        "",
        "## Synthesis Panel",
        "",
    ]
    for candidate in candidates:
        rg_text = (
            f"{candidate.mean_radius_of_gyration:.2f}"
            if candidate.mean_radius_of_gyration is not None
            else "not measured"
        )
        lines.extend(
            [
                f"- {candidate.sequence}: {candidate.role}.",
                f"  Rationale: {candidate.rationale}",
                f"  Scores: penetration {candidate.penetration_score:.3f}, function {candidate.functional_score:.3f}, mean Rg {rg_text}.",
            ]
        )

    lines.extend(
        [
            "",
            "## Core Experimental Flow",
            "",
            "1. Confirm synthesis feasibility, purity, identity, aqueous/formulation solubility, and short-term stability for each candidate.",
            "2. Run Franz diffusion permeation with KTTKS, Pal-KTTKS, vehicle blank, and selected analogs; quantify peptide by validated LC-MS/MS or HPLC.",
            "3. Measure fibroblast collagen/procollagen response only at concentrations shown to be non-toxic.",
            "4. Gate cosmetic safety with fibroblast/keratinocyte viability, LDH release, and reconstructed human epidermis irritation testing.",
            "",
            "## Vehicle Assumptions",
            "",
            "- Start with a simple aqueous or hydroalcoholic screening vehicle compatible with peptide stability and analytical recovery.",
            "- Treat final cream or serum formulation as a second-stage experiment because excipients can dominate skin delivery.",
            "- Maintain sink conditions in the receptor phase and verify mass balance at the end of each Franz diffusion run.",
            "",
            "## Reporting Guidance",
            "",
            "- Report permeation, skin retention, cytotoxicity, and collagen/procollagen effects separately.",
            "- Do not claim retained biological activity from permeation alone; require fibroblast assay evidence.",
            "- Frame RDKit conformer metrics as exploratory candidate-selection support, not MD-grade stability validation.",
        ]
    )
    return "\n".join(lines) + "\n"
