"""
Dependency-light Pareto optimization for Matrixyl-family candidates.

This module follows the NSGA-II mechanics we need for the paper: explicit
maximization objectives, non-dominated sorting, crowding distance, binary
tournament parent selection, and elitist survival.
"""

from dataclasses import asdict, dataclass
import math
import random
from statistics import mean
from typing import Iterable, Optional

from .candidates import CandidateEvaluation, CandidateEvaluator
from .search_space import SearchSpaceRules
from .tournament_search import TournamentSearchConfig, TournamentSearchOptimizer


DIFFICULT_SYNTHESIS_RESIDUES = frozenset("CMW")


def synthesis_feasibility_score(sequence: str) -> float:
    """Simple synthesis-risk proxy for short canonical peptide candidates."""
    if not sequence:
        return 0.0
    difficult = sum(1 for residue in sequence if residue in DIFFICULT_SYNTHESIS_RESIDUES)
    return 1.0 - (difficult / len(sequence))


@dataclass(frozen=True)
class ParetoObjectives:
    penetration: float
    functional_preservation: float
    synthesis_feasibility: float

    def values(self) -> tuple[float, float, float]:
        return (
            self.penetration,
            self.functional_preservation,
            self.synthesis_feasibility,
        )

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class ParetoCandidate:
    evaluation: CandidateEvaluation
    objectives: ParetoObjectives
    rank: int = -1
    crowding_distance: float = 0.0

    @property
    def sequence(self) -> str:
        return self.evaluation.sequence

    def with_rank_and_crowding(self, rank: int, crowding_distance: float) -> "ParetoCandidate":
        return ParetoCandidate(
            evaluation=self.evaluation,
            objectives=self.objectives,
            rank=rank,
            crowding_distance=crowding_distance,
        )

    def to_dict(self) -> dict:
        data = self.evaluation.to_dict()
        data["objectives"] = self.objectives.to_dict()
        data["pareto_rank"] = self.rank
        data["crowding_distance"] = self.crowding_distance
        return data


@dataclass(frozen=True)
class ParetoSearchConfig:
    population_size: int = 50
    generations: int = 50
    tournament_size: int = 2
    mutation_rate: float = 0.2
    seed: int = 42
    max_initialization_attempts: int = 5000
    max_mutation_attempts: int = 200

    def __post_init__(self):
        if self.population_size < 2:
            raise ValueError("population_size must be at least 2")
        if self.generations < 0:
            raise ValueError("generations must be non-negative")
        if self.tournament_size < 1:
            raise ValueError("tournament_size must be at least 1")
        if not 0.0 <= self.mutation_rate <= 1.0:
            raise ValueError("mutation_rate must be in [0, 1]")


@dataclass(frozen=True)
class ParetoGenerationSummary:
    generation: int
    frontier_size: int
    best_penetration: float
    best_functional_preservation: float
    best_synthesis_feasibility: float
    mean_penetration: float
    mean_functional_preservation: float
    unique_count: int

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class ParetoSearchResult:
    config: ParetoSearchConfig
    frontier: tuple[ParetoCandidate, ...]
    final_population: tuple[ParetoCandidate, ...]
    summaries: tuple[ParetoGenerationSummary, ...]
    evaluations_by_generation: tuple[tuple[ParetoCandidate, ...], ...]

    def to_dict(self) -> dict:
        return {
            "config": asdict(self.config),
            "frontier": [candidate.to_dict() for candidate in self.frontier],
            "final_population": [candidate.to_dict() for candidate in self.final_population],
            "summaries": [summary.to_dict() for summary in self.summaries],
            "evaluations_by_generation": [
                [candidate.to_dict() for candidate in generation]
                for generation in self.evaluations_by_generation
            ],
        }


def objectives_from_evaluation(evaluation: CandidateEvaluation) -> ParetoObjectives:
    penetration = evaluation.penetration.score if evaluation.penetration else 0.0
    return ParetoObjectives(
        penetration=penetration,
        functional_preservation=evaluation.functional_preservation.score,
        synthesis_feasibility=synthesis_feasibility_score(evaluation.sequence),
    )


def dominates(left: ParetoCandidate, right: ParetoCandidate) -> bool:
    left_values = left.objectives.values()
    right_values = right.objectives.values()
    return all(l >= r for l, r in zip(left_values, right_values)) and any(
        l > r for l, r in zip(left_values, right_values)
    )


def non_dominated_sort(candidates: Iterable[ParetoCandidate]) -> list[list[ParetoCandidate]]:
    population = unique_by_sequence(candidates)
    dominated_by_count = {candidate.sequence: 0 for candidate in population}
    dominates_map = {candidate.sequence: [] for candidate in population}
    fronts: list[list[ParetoCandidate]] = [[]]

    for candidate in population:
        for other in population:
            if candidate.sequence == other.sequence:
                continue
            if dominates(candidate, other):
                dominates_map[candidate.sequence].append(other)
            elif dominates(other, candidate):
                dominated_by_count[candidate.sequence] += 1
        if dominated_by_count[candidate.sequence] == 0:
            fronts[0].append(candidate)

    current_front_index = 0
    while current_front_index < len(fronts) and fronts[current_front_index]:
        next_front = []
        for candidate in fronts[current_front_index]:
            for dominated in dominates_map[candidate.sequence]:
                dominated_by_count[dominated.sequence] -= 1
                if dominated_by_count[dominated.sequence] == 0:
                    next_front.append(dominated)
        if next_front:
            fronts.append(next_front)
        current_front_index += 1

    return fronts


def unique_by_sequence(candidates: Iterable[ParetoCandidate]) -> list[ParetoCandidate]:
    """Deduplicate candidates by sequence before Pareto operations."""
    unique = {}
    for candidate in candidates:
        existing = unique.get(candidate.sequence)
        if existing is None or candidate.objectives.values() > existing.objectives.values():
            unique[candidate.sequence] = candidate
    return list(unique.values())


def assign_crowding_distance(front: list[ParetoCandidate]) -> list[ParetoCandidate]:
    if not front:
        return []
    if len(front) <= 2:
        return [
            candidate.with_rank_and_crowding(candidate.rank, math.inf)
            for candidate in front
        ]

    distances = {candidate.sequence: 0.0 for candidate in front}
    objective_getters = (
        lambda candidate: candidate.objectives.penetration,
        lambda candidate: candidate.objectives.functional_preservation,
        lambda candidate: candidate.objectives.synthesis_feasibility,
    )

    for getter in objective_getters:
        sorted_front = sorted(front, key=getter)
        distances[sorted_front[0].sequence] = math.inf
        distances[sorted_front[-1].sequence] = math.inf
        min_value = getter(sorted_front[0])
        max_value = getter(sorted_front[-1])
        if max_value == min_value:
            continue
        for index in range(1, len(sorted_front) - 1):
            if math.isinf(distances[sorted_front[index].sequence]):
                continue
            previous_value = getter(sorted_front[index - 1])
            next_value = getter(sorted_front[index + 1])
            distances[sorted_front[index].sequence] += (next_value - previous_value) / (
                max_value - min_value
            )

    return [
        candidate.with_rank_and_crowding(candidate.rank, distances[candidate.sequence])
        for candidate in front
    ]


def rank_population(candidates: Iterable[ParetoCandidate]) -> list[ParetoCandidate]:
    ranked = []
    for rank, front in enumerate(non_dominated_sort(unique_by_sequence(candidates))):
        ranked_front = [
            candidate.with_rank_and_crowding(rank, candidate.crowding_distance)
            for candidate in front
        ]
        ranked.extend(assign_crowding_distance(ranked_front))
    return sorted(
        ranked,
        key=lambda candidate: (
            -candidate.rank,
            candidate.crowding_distance,
            candidate.objectives.penetration,
            candidate.sequence,
        ),
        reverse=True,
    )


class ParetoSearchOptimizer:
    """NSGA-II-style optimizer over valid fixed-length peptide sequences."""

    def __init__(
        self,
        evaluator: Optional[CandidateEvaluator] = None,
        config: Optional[ParetoSearchConfig] = None,
    ):
        self.evaluator = evaluator or CandidateEvaluator()
        self.config = config or ParetoSearchConfig()
        self.search_space: SearchSpaceRules = self.evaluator.search_space
        self.rng = random.Random(self.config.seed)
        self._mutator = TournamentSearchOptimizer(
            evaluator=self.evaluator,
            config=TournamentSearchConfig(
                population_size=max(self.config.population_size, 2),
                generations=0,
                tournament_size=2,
                mutation_rate=self.config.mutation_rate,
                elite_count=1,
                seed=self.config.seed,
                max_initialization_attempts=self.config.max_initialization_attempts,
                max_mutation_attempts=self.config.max_mutation_attempts,
            ),
        )
        self._mutator.rng = self.rng

    def evaluate_sequence(self, sequence: str) -> ParetoCandidate:
        evaluation = self.evaluator.evaluate_sequence(sequence)
        return ParetoCandidate(
            evaluation=evaluation,
            objectives=objectives_from_evaluation(evaluation),
        )

    def evaluate_population(self, sequences: Iterable[str]) -> list[ParetoCandidate]:
        seen = set()
        candidates = []
        for sequence in sequences:
            if sequence in seen:
                continue
            seen.add(sequence)
            candidate = self.evaluate_sequence(sequence)
            if candidate.evaluation.is_valid:
                candidates.append(candidate)
        return rank_population(candidates)

    def initialize_population(self) -> list[str]:
        return self._mutator.initialize_population()

    def _select_parent(self, population: list[ParetoCandidate]) -> ParetoCandidate:
        contestants = self.rng.sample(population, min(self.config.tournament_size, len(population)))
        return sorted(
            contestants,
            key=lambda candidate: (
                -candidate.rank,
                candidate.crowding_distance,
                candidate.objectives.penetration,
                candidate.sequence,
            ),
            reverse=True,
        )[0]

    def _make_offspring(self, population: list[ParetoCandidate]) -> list[str]:
        offspring = []
        seen = set()
        attempts = 0
        while len(offspring) < self.config.population_size:
            attempts += 1
            if attempts > self.config.max_initialization_attempts:
                raise RuntimeError("Could not create enough Pareto offspring")
            parent = self._select_parent(population)
            child = self._mutator.mutate_sequence(parent.sequence, force_mutation=True)
            if child not in seen and self.search_space.validate(child).is_valid:
                offspring.append(child)
                seen.add(child)
        return offspring

    def _fill_with_immigrants(self, population: list[ParetoCandidate]) -> list[ParetoCandidate]:
        filled = list(population)
        seen = {candidate.sequence for candidate in filled}
        attempts = 0

        while len(filled) < self.config.population_size:
            attempts += 1
            if attempts > self.config.max_initialization_attempts:
                raise RuntimeError("Could not maintain a full diverse Pareto population")
            seed_sequence = self.search_space.reference_sequence
            if filled:
                seed_sequence = self.rng.choice(filled).sequence
            immigrant = self._mutator.mutate_sequence(seed_sequence, force_mutation=True)
            if immigrant in seen or not self.search_space.validate(immigrant).is_valid:
                continue
            candidate = self.evaluate_sequence(immigrant)
            if candidate.evaluation.is_valid:
                filled.append(candidate)
                seen.add(candidate.sequence)

        return rank_population(filled)

    def _survival_select(self, candidates: list[ParetoCandidate]) -> list[ParetoCandidate]:
        ranked = rank_population(unique_by_sequence(candidates))
        selected = []
        fronts = non_dominated_sort(ranked)
        for rank, front in enumerate(fronts):
            ranked_front = [
                candidate.with_rank_and_crowding(rank, candidate.crowding_distance)
                for candidate in front
            ]
            crowded_front = assign_crowding_distance(ranked_front)
            if len(selected) + len(crowded_front) <= self.config.population_size:
                selected.extend(crowded_front)
            else:
                remaining = self.config.population_size - len(selected)
                selected.extend(
                    sorted(
                        crowded_front,
                        key=lambda candidate: (
                            candidate.crowding_distance,
                            candidate.objectives.penetration,
                            candidate.sequence,
                        ),
                        reverse=True,
                    )[:remaining]
                )
                break
        return self._fill_with_immigrants(rank_population(selected))

    def summarize_generation(
        self,
        generation: int,
        population: list[ParetoCandidate],
    ) -> ParetoGenerationSummary:
        frontier_size = sum(1 for candidate in population if candidate.rank == 0)
        return ParetoGenerationSummary(
            generation=generation,
            frontier_size=frontier_size,
            best_penetration=max(candidate.objectives.penetration for candidate in population),
            best_functional_preservation=max(
                candidate.objectives.functional_preservation for candidate in population
            ),
            best_synthesis_feasibility=max(
                candidate.objectives.synthesis_feasibility for candidate in population
            ),
            mean_penetration=mean(candidate.objectives.penetration for candidate in population),
            mean_functional_preservation=mean(
                candidate.objectives.functional_preservation for candidate in population
            ),
            unique_count=len({candidate.sequence for candidate in population}),
        )

    def run(self) -> ParetoSearchResult:
        population = self.evaluate_population(self.initialize_population())
        evaluations_by_generation = []
        summaries = []

        for generation in range(self.config.generations + 1):
            population = rank_population(population)
            evaluations_by_generation.append(tuple(population))
            summaries.append(self.summarize_generation(generation, population))
            if generation < self.config.generations:
                offspring = self.evaluate_population(self._make_offspring(population))
                population = self._survival_select(population + offspring)

        final_population = rank_population(population)
        frontier = tuple(candidate for candidate in final_population if candidate.rank == 0)
        return ParetoSearchResult(
            config=self.config,
            frontier=frontier,
            final_population=tuple(final_population),
            summaries=tuple(summaries),
            evaluations_by_generation=tuple(evaluations_by_generation),
        )
