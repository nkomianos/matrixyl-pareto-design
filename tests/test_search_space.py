import unittest

from src.search_space import SearchSpaceRules, levenshtein_distance


class SearchSpaceTests(unittest.TestCase):
    def test_levenshtein_distance(self):
        self.assertEqual(levenshtein_distance("KTTKS", "KTTKS"), 0)
        self.assertEqual(levenshtein_distance("KTTKS", "KSTKS"), 1)
        self.assertEqual(levenshtein_distance("KTTKS", "KTTKSS"), 1)

    def test_default_matrixyl_rules_accept_reference_and_close_mutant(self):
        rules = SearchSpaceRules.matrixyl_default()

        self.assertTrue(rules.validate("KTTKS").is_valid)
        self.assertTrue(rules.validate("KSTKS").is_valid)

    def test_default_matrixyl_rules_reject_invalid_candidates(self):
        rules = SearchSpaceRules.matrixyl_default()

        self.assertIn("invalid_amino_acid", rules.validate("KTTKZ").failures)
        self.assertIn("too_long", rules.validate("KTTKSS").failures)
        self.assertIn("too_many_edits", rules.validate("AAAAA").failures)

    def test_locked_positions_are_enforced(self):
        rules = SearchSpaceRules.matrixyl_default(locked_positions=(0, 4))

        self.assertTrue(rules.validate("KTTKS").is_valid)
        self.assertIn("locked_position_0", rules.validate("RTTKS").failures)
        self.assertEqual(rules.allowed_residues_at(0), ("K",))
        self.assertGreater(len(rules.allowed_residues_at(1)), 1)


if __name__ == "__main__":
    unittest.main()
