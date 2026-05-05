import math
import unittest

from src.candidates import CandidateEvaluator
from src.pareto_search import (
    ParetoSearchConfig,
    ParetoSearchOptimizer,
    assign_crowding_distance,
    dominates,
    non_dominated_sort,
    objectives_from_evaluation,
    rank_population,
    synthesis_feasibility_score,
    unique_by_sequence,
)


class ParetoSearchTests(unittest.TestCase):
    def setUp(self):
        self.evaluator = CandidateEvaluator()

    def candidate(self, sequence):
        evaluation = self.evaluator.evaluate_sequence(sequence)
        return ParetoSearchOptimizer(self.evaluator).evaluate_sequence(evaluation.sequence)

    def test_synthesis_feasibility_penalizes_difficult_residues(self):
        self.assertEqual(synthesis_feasibility_score("KTTKS"), 1.0)
        self.assertEqual(synthesis_feasibility_score("KTTMS"), 0.8)
        self.assertEqual(synthesis_feasibility_score(""), 0.0)

    def test_objectives_are_extracted_from_candidate_evaluation(self):
        evaluation = self.evaluator.evaluate_sequence("KTTKS")
        objectives = objectives_from_evaluation(evaluation)

        self.assertEqual(objectives.penetration, evaluation.penetration.score)
        self.assertEqual(objectives.functional_preservation, 1.0)
        self.assertEqual(objectives.synthesis_feasibility, 1.0)

    def test_dominance_uses_maximized_objectives(self):
        strong = self.candidate("KTTPS")
        weak = self.candidate("KTTMS")

        self.assertTrue(dominates(strong, weak))
        self.assertFalse(dominates(weak, strong))

    def test_non_dominated_sort_identifies_frontier(self):
        candidates = [self.candidate(sequence) for sequence in ["KTTKS", "KTTPS", "PTTPS"]]
        fronts = non_dominated_sort(candidates)

        frontier_sequences = {candidate.sequence for candidate in fronts[0]}
        self.assertIn("KTTKS", frontier_sequences)
        self.assertIn("PTTPS", frontier_sequences)

    def test_crowding_distance_preserves_extremes(self):
        ranked = rank_population(
            [self.candidate(sequence) for sequence in ["KTTKS", "KTTPS", "PTTPS", "KTTPP"]]
        )
        frontier = [candidate for candidate in ranked if candidate.rank == 0]
        crowded = assign_crowding_distance(frontier)

        self.assertTrue(any(math.isinf(candidate.crowding_distance) for candidate in crowded))

    def test_unique_by_sequence_deduplicates_frontier_inputs(self):
        duplicate = self.candidate("PTTPS")
        unique = unique_by_sequence([duplicate, duplicate, self.candidate("KTTKS")])

        self.assertEqual({candidate.sequence for candidate in unique}, {"PTTPS", "KTTKS"})

    def test_pareto_run_is_reproducible_and_valid(self):
        config = ParetoSearchConfig(
            population_size=12,
            generations=3,
            tournament_size=2,
            mutation_rate=0.35,
            seed=17,
        )

        first = ParetoSearchOptimizer(config=config).run()
        second = ParetoSearchOptimizer(config=config).run()

        self.assertEqual(
            [candidate.sequence for candidate in first.frontier],
            [candidate.sequence for candidate in second.frontier],
        )
        self.assertEqual(len(first.final_population), config.population_size)
        self.assertEqual(
            len({candidate.sequence for candidate in first.final_population}),
            config.population_size,
        )
        self.assertEqual(
            len({candidate.sequence for candidate in first.frontier}),
            len(first.frontier),
        )
        self.assertTrue(all(candidate.evaluation.is_valid for candidate in first.final_population))
        self.assertGreaterEqual(len(first.frontier), 1)
        self.assertEqual(len(first.summaries), config.generations + 1)


if __name__ == "__main__":
    unittest.main()
