import unittest

from src.function_scores import FunctionalPreservationScorer


class FunctionalPreservationTests(unittest.TestCase):
    def test_reference_scores_perfectly(self):
        score = FunctionalPreservationScorer().score("KTTKS")

        self.assertEqual(score.score, 1.0)
        self.assertEqual(score.edit_distance, 0)
        self.assertEqual(score.components["identity"], 1.0)
        self.assertEqual(score.components["substitution"], 1.0)

    def test_conservative_substitution_scores_above_nonconservative(self):
        scorer = FunctionalPreservationScorer()

        conservative = scorer.score("KSTKS")
        nonconservative = scorer.score("KWTKS")

        self.assertGreater(conservative.score, nonconservative.score)
        self.assertEqual(conservative.edit_distance, 1)
        self.assertEqual(nonconservative.edit_distance, 1)

    def test_length_changes_are_penalized(self):
        scorer = FunctionalPreservationScorer()

        same_prefix_extra_residue = scorer.score("KTTKSS")
        exact = scorer.score("KTTKS")

        self.assertLess(same_prefix_extra_residue.score, exact.score)
        self.assertLess(same_prefix_extra_residue.components["length"], 1.0)


if __name__ == "__main__":
    unittest.main()
