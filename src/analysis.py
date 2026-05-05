"""
Candidate analysis and interpretability helpers.
"""

import csv
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from .chemistry import descriptors_from_sequence, descriptors_from_smiles, read_smiles
from .constraints import PhysicochemicalCalculator
from .function_scores import FunctionalPreservationScorer


REFERENCE_SEQUENCE = "KTTKS"


@dataclass(frozen=True)
class CandidateRow:
    sequence: str
    penetration_objective: float
    functional_objective: float
    synthesis_objective: float
    scalar_optimization_score: float
    edit_distance: int
    molecular_weight: float
    tpsa: float
    logp: float
    hbd: int
    hba: int
    rotatable_bonds: int
    formal_charge: int

    @classmethod
    def from_csv_row(cls, row: dict) -> "CandidateRow":
        return cls(
            sequence=row["sequence"],
            penetration_objective=float(row["penetration_objective"]),
            functional_objective=float(row["functional_objective"]),
            synthesis_objective=float(row["synthesis_objective"]),
            scalar_optimization_score=float(row["scalar_optimization_score"]),
            edit_distance=int(row["edit_distance"]),
            molecular_weight=float(row["molecular_weight"]),
            tpsa=float(row["tpsa"]),
            logp=float(row["logp"]),
            hbd=int(row["hbd"]),
            hba=int(row["hba"]),
            rotatable_bonds=int(row["rotatable_bonds"]),
            formal_charge=int(row["formal_charge"]),
        )

    def to_dict(self) -> dict:
        return asdict(self)


def load_candidate_rows(path: str | Path) -> list[CandidateRow]:
    with Path(path).open() as handle:
        return [CandidateRow.from_csv_row(row) for row in csv.DictReader(handle)]


def write_dict_rows(path: str | Path, rows: list[dict]) -> None:
    if not rows:
        return
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def position_frequency_matrix(
    sequences: Iterable[str],
    reference_sequence: str = REFERENCE_SEQUENCE,
) -> list[dict]:
    sequences = list(sequences)
    if not sequences:
        return []
    length = len(reference_sequence)
    if any(len(sequence) != length for sequence in sequences):
        raise ValueError("All sequences must match the reference length")

    rows = []
    for position in range(length):
        counts = Counter(sequence[position] for sequence in sequences)
        total = sum(counts.values())
        for residue, count in sorted(counts.items()):
            rows.append(
                {
                    "position": position + 1,
                    "reference_residue": reference_sequence[position],
                    "residue": residue,
                    "count": count,
                    "frequency": count / total,
                    "is_reference": residue == reference_sequence[position],
                }
            )
    return rows


def mutation_enrichment(
    sequences: Iterable[str],
    reference_sequence: str = REFERENCE_SEQUENCE,
) -> list[dict]:
    sequences = list(sequences)
    if not sequences:
        return []
    length = len(reference_sequence)
    if any(len(sequence) != length for sequence in sequences):
        raise ValueError("All sequences must match the reference length")

    rows = []
    for position in range(length):
        reference_residue = reference_sequence[position]
        mutated = [
            sequence[position]
            for sequence in sequences
            if sequence[position] != reference_residue
        ]
        counts = Counter(mutated)
        for residue, count in sorted(counts.items()):
            rows.append(
                {
                    "position": position + 1,
                    "reference_residue": reference_residue,
                    "mutant_residue": residue,
                    "mutation": f"{reference_residue}{position + 1}{residue}",
                    "count": count,
                    "frequency_among_candidates": count / len(sequences),
                }
            )
    return sorted(rows, key=lambda row: (-row["count"], row["position"], row["mutant_residue"]))


def descriptor_warning_flags(row: CandidateRow) -> tuple[str, ...]:
    flags = []
    if row.molecular_weight > 500:
        flags.append("mw_above_500")
    if row.tpsa > 140:
        flags.append("tpsa_above_140")
    if row.logp < 1:
        flags.append("logp_below_1")
    elif row.logp > 3:
        flags.append("logp_above_3")
    if row.hbd > 5:
        flags.append("hbd_above_5")
    if row.hba > 10:
        flags.append("hba_above_10")
    if row.rotatable_bonds > 10:
        flags.append("rotatable_bonds_above_10")
    if abs(row.formal_charge) > 1:
        flags.append("charge_outside_minus1_to_1")
    return tuple(flags)


def priority_label(row: CandidateRow) -> str:
    flags = descriptor_warning_flags(row)
    if row.functional_objective >= 0.8 and row.penetration_objective >= 0.5:
        return "balanced"
    if row.penetration_objective >= 0.65:
        return "high_penetration_tradeoff"
    if row.functional_objective >= 0.9:
        return "function_anchor"
    if len(flags) <= 3:
        return "descriptor_followup"
    return "exploratory"


def summarize_candidates(rows: list[CandidateRow]) -> list[dict]:
    return [
        {
            **row.to_dict(),
            "priority_label": priority_label(row),
            "warning_flags": ";".join(descriptor_warning_flags(row)),
        }
        for row in sorted(
            rows,
            key=lambda item: (
                item.penetration_objective,
                item.functional_objective,
                item.scalar_optimization_score,
            ),
            reverse=True,
        )
    ]


def baseline_descriptor_rows(palmitoylated_smiles_path: str | Path) -> list[dict]:
    calculator = PhysicochemicalCalculator()
    function_scorer = FunctionalPreservationScorer()
    core_descriptors = descriptors_from_sequence(REFERENCE_SEQUENCE, name="Matrixyl core")
    smiles, name = read_smiles(palmitoylated_smiles_path)
    pal_descriptors = descriptors_from_smiles(smiles, name=name, source="pubchem_smiles")

    rows = []
    for label, descriptors, sequence in (
        ("Matrixyl core", core_descriptors, REFERENCE_SEQUENCE),
        ("Pal-KTTKS", pal_descriptors, REFERENCE_SEQUENCE),
    ):
        score = calculator.score_descriptors(descriptors)
        function_score = function_scorer.score(sequence).score
        rows.append(
            {
                "name": label,
                "sequence": sequence,
                "formula": descriptors.formula,
                "molecular_weight": descriptors.molecular_weight,
                "tpsa": descriptors.tpsa,
                "logp": descriptors.logp,
                "hbd": descriptors.hbd,
                "hba": descriptors.hba,
                "rotatable_bonds": descriptors.rotatable_bonds,
                "formal_charge": descriptors.formal_charge,
                "penetration_score": score.score,
                "functional_score": function_score,
            }
        )
    return rows
