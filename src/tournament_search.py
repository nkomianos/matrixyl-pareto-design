"""
Deterministic tournament-search baseline for Matrixyl-family peptide evolution.
"""

from dataclasses import asdict, dataclass
import random
from statistics import mean
from typing import Optional

from .candidates import CandidateEvaluation, CandidateEvaluator
from .search_space import SearchSpaceRules


@dataclass(frozen=True)
class TournamentSearchConfig:
    population_size: int = 30
    generations: int = 25
    tournament_size: int = 3
    mutation_rate: float = 0.2
    elite_count: int = 1
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
        if self.elite_count < 1 or self.elite_count >= self.population_size:
            raise ValueError("elite_count must be at least 1 and less than population_size")


@dataclass(frozen=True)
class GenerationSummary:
    generation: int
    best_sequence: str
    best_score: float
    mean_score: float
    valid_count: int
    unique_count: int

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class TournamentSearchResult:
    config: TournamentSearchConfig
    best: CandidateEvaluation
    final_population: tuple[CandidateEvaluation, ...]
    summaries: tuple[GenerationSummary, ...]
    evaluations_by_generation: tuple[tuple[CandidateEvaluation, ...], ...]
    random_baseline: tuple[CandidateEvaluation, ...]

    def to_dict(self) -> dict:
        return {
            "config": asdict(self.config),
            "best": self.best.to_dict(),
            "final_population": [candidate.to_dict() for candidate in self.final_population],
            "summaries": [summary.to_dict() for summary in self.summaries],
            "evaluations_by_generation": [
                [candidate.to_dict() for candidate in generation]
                for generation in self.evaluations_by_generation
            ],
            "random_baseline": [candidate.to_dict() for candidate in self.random_baseline],
        }


class TournamentSearchOptimizer:
    """Mutation-only tournament selection using the shared candidate evaluator."""

    def __init__(
        self,
        evaluator: Optional[CandidateEvaluator] = None,
        config: Optional[TournamentSearchConfig] = None,
    ):
        self.evaluator = evaluator or CandidateEvaluator()
        self.config = config or TournamentSearchConfig()
        self.rng = random.Random(self.config.seed)
        self.search_space: SearchSpaceRules = self.evaluator.search_space

    @staticmethod
    def _rank(evaluations: list[CandidateEvaluation]) -> list[CandidateEvaluation]:
        return sorted(
            evaluations,
            key=lambda candidate: (candidate.optimization_score, candidate.sequence),
            reverse=True,
        )

    def evaluate_population(self, sequences: list[str]) -> list[CandidateEvaluation]:
        return [
            self.evaluator.evaluate_sequence(sequence)
            for sequence in sequences
        ]

    def summarize_generation(
        self,
        generation: int,
        evaluations: list[CandidateEvaluation],
    ) -> GenerationSummary:
        ranked = self._rank(evaluations)
        valid_scores = [
            candidate.optimization_score
            for candidate in evaluations
            if candidate.is_valid
        ]
        return GenerationSummary(
            generation=generation,
            best_sequence=ranked[0].sequence,
            best_score=ranked[0].optimization_score,
            mean_score=mean(valid_scores) if valid_scores else 0.0,
            valid_count=sum(candidate.is_valid for candidate in evaluations),
            unique_count=len({candidate.sequence for candidate in evaluations}),
        )

    def initialize_population(self) -> list[str]:
        reference = self.search_space.reference_sequence
        population = [reference]
        seen = {reference}
        attempts = 0

        while len(population) < self.config.population_size:
            attempts += 1
            if attempts > self.config.max_initialization_attempts:
                raise RuntimeError("Could not initialize a full valid population")
            candidate = self.mutate_sequence(reference, force_mutation=True)
            if candidate not in seen and self.search_space.validate(candidate).is_valid:
                population.append(candidate)
                seen.add(candidate)

        return population

    def mutate_sequence(self, sequence: str, force_mutation: bool = True) -> str:
        reference = self.search_space.reference_sequence
        positions = list(range(len(reference)))

        for _ in range(self.config.max_mutation_attempts):
            mutated = list(sequence)
            changed = False

            for position in positions:
                if self.rng.random() < self.config.mutation_rate:
                    options = [
                        residue
                        for residue in self.search_space.allowed_residues_at(position)
                        if residue != mutated[position]
                    ]
                    if options:
                        mutated[position] = self.rng.choice(options)
                        changed = True

            if force_mutation and not changed:
                mutable_positions = [
                    position
                    for position in positions
                    if len(self.search_space.allowed_residues_at(position)) > 1
                ]
                if not mutable_positions:
                    return sequence
                position = self.rng.choice(mutable_positions)
                options = [
                    residue
                    for residue in self.search_space.allowed_residues_at(position)
                    if residue != mutated[position]
                ]
                mutated[position] = self.rng.choice(options)

            candidate = "".join(mutated)
            if self.search_space.validate(candidate).is_valid:
                return candidate

        return sequence

    def tournament_select(self, evaluations: list[CandidateEvaluation]) -> CandidateEvaluation:
        tournament_size = min(self.config.tournament_size, len(evaluations))
        contestants = self.rng.sample(evaluations, tournament_size)
        return self._rank(contestants)[0]

    def _next_generation(self, evaluations: list[CandidateEvaluation]) -> list[str]:
        ranked = self._rank(evaluations)
        next_sequences = [candidate.sequence for candidate in ranked[:self.config.elite_count]]
        seen = set(next_sequences)

        while len(next_sequences) < self.config.population_size:
            parent = self.tournament_select(evaluations)
            child = self.mutate_sequence(parent.sequence, force_mutation=True)
            if child not in seen:
                next_sequences.append(child)
                seen.add(child)

        return next_sequences

    def random_baseline(self) -> list[CandidateEvaluation]:
        sequences = self.initialize_population()
        return self._rank(self.evaluate_population(sequences))

    def run(self) -> TournamentSearchResult:
        population = self.initialize_population()
        evaluations_by_generation = []
        summaries = []

        for generation in range(self.config.generations + 1):
            evaluations = self._rank(self.evaluate_population(population))
            evaluations_by_generation.append(tuple(evaluations))
            summaries.append(self.summarize_generation(generation, evaluations))

            if generation < self.config.generations:
                population = self._next_generation(evaluations)

        final_population = evaluations_by_generation[-1]
        all_evaluations = [
            candidate
            for generation in evaluations_by_generation
            for candidate in generation
        ]
        best = self._rank(all_evaluations)[0]

        return TournamentSearchResult(
            config=self.config,
            best=best,
            final_population=final_population,
            summaries=tuple(summaries),
            evaluations_by_generation=tuple(evaluations_by_generation),
            random_baseline=tuple(self.random_baseline()),
        )
