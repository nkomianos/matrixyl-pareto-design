"""
Search-space rules for Matrixyl-family peptide candidates.

These checks are deliberately cheap and explicit so the evolutionary algorithm
can reject candidates before running heavier descriptor or oracle scoring.
"""

from dataclasses import asdict, dataclass
from typing import Iterable, Optional

from .chemistry import CANONICAL_AMINO_ACIDS


DEFAULT_REFERENCE_SEQUENCE = "KTTKS"


def normalize_sequence(sequence: str) -> str:
    return sequence.strip().upper()


def levenshtein_distance(left: str, right: str) -> int:
    """Compute edit distance without external dependencies."""
    if left == right:
        return 0
    if not left:
        return len(right)
    if not right:
        return len(left)

    previous = list(range(len(right) + 1))
    for i, left_char in enumerate(left, start=1):
        current = [i]
        for j, right_char in enumerate(right, start=1):
            insertion = current[j - 1] + 1
            deletion = previous[j] + 1
            substitution = previous[j - 1] + (left_char != right_char)
            current.append(min(insertion, deletion, substitution))
        previous = current
    return previous[-1]


@dataclass(frozen=True)
class SearchSpaceValidation:
    sequence: str
    is_valid: bool
    failures: tuple[str, ...]
    edit_distance: int
    length: int

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class SearchSpaceRules:
    reference_sequence: str = DEFAULT_REFERENCE_SEQUENCE
    min_length: int = len(DEFAULT_REFERENCE_SEQUENCE)
    max_length: int = len(DEFAULT_REFERENCE_SEQUENCE)
    allowed_amino_acids: frozenset[str] = CANONICAL_AMINO_ACIDS
    max_edit_distance: int = 2
    locked_positions: tuple[int, ...] = ()

    def __post_init__(self):
        reference = normalize_sequence(self.reference_sequence)
        if not reference:
            raise ValueError("reference_sequence cannot be empty")
        invalid = sorted(set(reference) - set(self.allowed_amino_acids))
        if invalid:
            raise ValueError(f"reference_sequence contains invalid amino acids: {invalid}")
        if self.min_length < 1 or self.max_length < self.min_length:
            raise ValueError("length bounds are invalid")
        invalid_locks = [
            position
            for position in self.locked_positions
            if position < 0 or position >= len(reference)
        ]
        if invalid_locks:
            raise ValueError(f"locked_positions outside reference sequence: {invalid_locks}")

        object.__setattr__(self, "reference_sequence", reference)
        object.__setattr__(self, "locked_positions", tuple(sorted(set(self.locked_positions))))

    @classmethod
    def matrixyl_default(cls, locked_positions: Optional[Iterable[int]] = None) -> "SearchSpaceRules":
        """
        Conservative default for early experiments.

        Keep the KTTKS core length fixed, allow canonical amino acids, and cap
        candidates at two edits from the reference. Locked positions can be
        supplied by experiments when testing stricter motif preservation.
        """
        return cls(locked_positions=tuple(locked_positions or ()))

    def validate(self, sequence: str) -> SearchSpaceValidation:
        normalized = normalize_sequence(sequence)
        failures = []

        if len(normalized) < self.min_length:
            failures.append("too_short")
        if len(normalized) > self.max_length:
            failures.append("too_long")

        invalid_aas = sorted(set(normalized) - set(self.allowed_amino_acids))
        if invalid_aas:
            failures.append("invalid_amino_acid")

        edit_distance = levenshtein_distance(normalized, self.reference_sequence)
        if edit_distance > self.max_edit_distance:
            failures.append("too_many_edits")

        if len(normalized) == len(self.reference_sequence):
            for position in self.locked_positions:
                if normalized[position] != self.reference_sequence[position]:
                    failures.append(f"locked_position_{position}")
        elif self.locked_positions:
            failures.append("locked_positions_require_reference_length")

        deduped_failures = tuple(dict.fromkeys(failures))
        return SearchSpaceValidation(
            sequence=normalized,
            is_valid=not deduped_failures,
            failures=deduped_failures,
            edit_distance=edit_distance,
            length=len(normalized),
        )

    def allowed_residues_at(self, position: int) -> tuple[str, ...]:
        if position < 0 or position >= len(self.reference_sequence):
            raise IndexError("position outside reference sequence")
        if position in self.locked_positions:
            return (self.reference_sequence[position],)
        return tuple(sorted(self.allowed_amino_acids))
