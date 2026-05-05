import unittest

from src.tournament_search import TournamentSearchConfig, TournamentSearchOptimizer


class TournamentSearchTests(unittest.TestCase):
    def test_initial_population_is_deterministic_and_valid(self):
        config = TournamentSearchConfig(population_size=8, generations=2, seed=7)
        first = TournamentSearchOptimizer(config=config).initialize_population()
        second = TournamentSearchOptimizer(config=config).initialize_population()

        self.assertEqual(first, second)
        self.assertEqual(first[0], "KTTKS")
        self.assertEqual(len(first), 8)
        self.assertEqual(len(set(first)), 8)

        optimizer = TournamentSearchOptimizer(config=config)
        self.assertTrue(all(optimizer.search_space.validate(sequence).is_valid for sequence in first))

    def test_run_is_reproducible_with_fixed_seed(self):
        config = TournamentSearchConfig(
            population_size=10,
            generations=3,
            tournament_size=3,
            mutation_rate=0.4,
            seed=99,
        )

        first = TournamentSearchOptimizer(config=config).run()
        second = TournamentSearchOptimizer(config=config).run()

        self.assertEqual(first.best.sequence, second.best.sequence)
        self.assertEqual(first.best.optimization_score, second.best.optimization_score)
        self.assertEqual(
            [summary.to_dict() for summary in first.summaries],
            [summary.to_dict() for summary in second.summaries],
        )

    def test_run_keeps_valid_population_and_history(self):
        config = TournamentSearchConfig(population_size=12, generations=4, seed=11)
        result = TournamentSearchOptimizer(config=config).run()

        self.assertEqual(len(result.summaries), config.generations + 1)
        self.assertEqual(len(result.evaluations_by_generation), config.generations + 1)
        self.assertEqual(len(result.final_population), config.population_size)
        self.assertTrue(result.best.is_valid)
        self.assertTrue(all(candidate.is_valid for candidate in result.final_population))

    def test_random_baseline_is_reported(self):
        config = TournamentSearchConfig(population_size=8, generations=1, seed=5)
        result = TournamentSearchOptimizer(config=config).run()

        self.assertEqual(len(result.random_baseline), config.population_size)
        self.assertTrue(all(candidate.is_valid for candidate in result.random_baseline))


if __name__ == "__main__":
    unittest.main()
