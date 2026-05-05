"""
Shared candidate evaluation schema for optimization experiments.
"""

from dataclasses import asdict, dataclass
from typing import Optional

from .constraints import PenetrationScore, PhysicochemicalCalculator
from .function_scores import FunctionalPreservationScore, FunctionalPreservationScorer
from .search_space import SearchSpaceRules, SearchSpaceValidation, normalize_sequence


@dataclass(frozen=True)
class CandidateEvaluation:
    sequence: str
    modification: str
    search_validation: SearchSpaceValidation
    functional_preservation: FunctionalPreservationScore
    penetration: Optional[PenetrationScore]
    is_valid: bool
    failed_filters: tuple[str, ...]
    optimization_score: float

    def to_dict(self) -> dict:
        data = asdict(self)
        data["search_validation"] = self.search_validation.to_dict()
        data["functional_preservation"] = self.functional_preservation.to_dict()
        data["penetration"] = self.penetration.to_dict() if self.penetration else None
        return data


class CandidateEvaluator:
    """
    Evaluate canonical peptide candidates through the Phase 2 scoring contract.

    The optimization score is intentionally simple for Phase 3 baselines:
    permeability score multiplied by functional-preservation score, with hard
    filters setting invalid candidates to zero.
    """

    def __init__(
        self,
        search_space: Optional[SearchSpaceRules] = None,
        physicochemical_calculator: Optional[PhysicochemicalCalculator] = None,
        functional_scorer: Optional[FunctionalPreservationScorer] = None,
    ):
        self.search_space = search_space or SearchSpaceRules.matrixyl_default()
        self.physicochemical_calculator = physicochemical_calculator or PhysicochemicalCalculator()
        self.functional_scorer = functional_scorer or FunctionalPreservationScorer(
            reference_sequence=self.search_space.reference_sequence
        )

    def evaluate_sequence(
        self,
        sequence: str,
        *,
        modification: str = "unmodified",
        name: Optional[str] = None,
    ) -> CandidateEvaluation:
        normalized = normalize_sequence(sequence)
        search_validation = self.search_space.validate(normalized)
        functional_preservation = self.functional_scorer.score(normalized)
        penetration = None
        failed_filters = list(search_validation.failures)

        if search_validation.is_valid:
            try:
                penetration = self.physicochemical_calculator.score_sequence_exact(
                    normalized,
                    name=name or normalized,
                )
            except ValueError:
                failed_filters.append("descriptor_parse_failed")

        is_valid = not failed_filters and penetration is not None
        optimization_score = (
            penetration.score * functional_preservation.score
            if is_valid and penetration is not None
            else 0.0
        )

        return CandidateEvaluation(
            sequence=normalized,
            modification=modification,
            search_validation=search_validation,
            functional_preservation=functional_preservation,
            penetration=penetration,
            is_valid=is_valid,
            failed_filters=tuple(dict.fromkeys(failed_filters)),
            optimization_score=optimization_score,
        )
