"""
Functional-preservation proxies for Matrixyl-family candidates.

These scores do not claim measured collagen stimulation. They keep early
evolution close to the KTTKS matrikine core until experimental labels or a
validated functional oracle are available.
"""

from dataclasses import asdict, dataclass
from typing import Dict, Optional

from .search_space import DEFAULT_REFERENCE_SEQUENCE, levenshtein_distance, normalize_sequence


AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"

# BLOSUM62 values for the 20 canonical amino acids, ordered by AMINO_ACIDS.
BLOSUM62_ROWS = {
    "A": "4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0",
    "R": "-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3",
    "N": "-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3",
    "D": "-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3",
    "C": "0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1",
    "Q": "-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2",
    "E": "-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2",
    "G": "0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3",
    "H": "-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3",
    "I": "-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3",
    "L": "-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1",
    "K": "-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2",
    "M": "-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1",
    "F": "-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1",
    "P": "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2",
    "S": "1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2",
    "T": "0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0",
    "W": "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3",
    "Y": "-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1",
    "V": "0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4",
}

BLOSUM62 = {
    aa: {
        other: int(value)
        for other, value in zip(AMINO_ACIDS, BLOSUM62_ROWS[aa].split())
    }
    for aa in AMINO_ACIDS
}


@dataclass(frozen=True)
class FunctionalPreservationScore:
    sequence: str
    reference_sequence: str
    score: float
    components: Dict[str, float]
    edit_distance: int

    def to_dict(self) -> dict:
        return asdict(self)


class FunctionalPreservationScorer:
    DEFAULT_WEIGHTS = {
        "identity": 0.45,
        "substitution": 0.35,
        "edit": 0.15,
        "length": 0.05,
    }

    def __init__(
        self,
        reference_sequence: str = DEFAULT_REFERENCE_SEQUENCE,
        weights: Optional[Dict[str, float]] = None,
    ):
        self.reference_sequence = normalize_sequence(reference_sequence)
        self.weights = dict(weights or self.DEFAULT_WEIGHTS)
        total = sum(self.weights.values())
        if not self.reference_sequence:
            raise ValueError("reference_sequence cannot be empty")
        if total <= 0:
            raise ValueError("functional preservation weights must sum to a positive value")
        self.weights = {key: value / total for key, value in self.weights.items()}

    @staticmethod
    def normalized_blosum_similarity(reference_aa: str, candidate_aa: str) -> float:
        if reference_aa == candidate_aa:
            return 1.0
        raw_score = BLOSUM62[reference_aa][candidate_aa]
        self_score = BLOSUM62[reference_aa][reference_aa]
        normalized = (raw_score + 4) / (self_score + 4)
        return min(max(normalized, 0.0), 1.0)

    def score(self, sequence: str) -> FunctionalPreservationScore:
        candidate = normalize_sequence(sequence)
        reference = self.reference_sequence
        compared_length = min(len(candidate), len(reference))

        if compared_length:
            matches = sum(
                1
                for index in range(compared_length)
                if candidate[index] == reference[index]
            )
            identity = matches / len(reference)
            substitution = sum(
                self.normalized_blosum_similarity(reference[index], candidate[index])
                for index in range(compared_length)
            ) / len(reference)
        else:
            identity = 0.0
            substitution = 0.0

        edit_distance = levenshtein_distance(candidate, reference)
        edit_component = 1.0 - min(edit_distance / max(len(candidate), len(reference), 1), 1.0)
        length_component = 1.0 - min(abs(len(candidate) - len(reference)) / len(reference), 1.0)
        components = {
            "identity": identity,
            "substitution": substitution,
            "edit": edit_component,
            "length": length_component,
        }
        score = sum(self.weights[key] * components[key] for key in self.weights)
        return FunctionalPreservationScore(
            sequence=candidate,
            reference_sequence=reference,
            score=min(max(score, 0.0), 1.0),
            components=components,
            edit_distance=edit_distance,
        )
