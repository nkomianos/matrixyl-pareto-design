import unittest

from src.candidates import CandidateEvaluator
from src.search_space import SearchSpaceRules


class CandidateEvaluationTests(unittest.TestCase):
    def test_reference_candidate_uses_shared_schema(self):
        evaluation = CandidateEvaluator().evaluate_sequence("KTTKS", name="Matrixyl core")

        self.assertTrue(evaluation.is_valid)
        self.assertEqual(evaluation.sequence, "KTTKS")
        self.assertEqual(evaluation.modification, "unmodified")
        self.assertEqual(evaluation.failed_filters, ())
        self.assertIsNotNone(evaluation.penetration)
        self.assertEqual(evaluation.functional_preservation.score, 1.0)
        self.assertEqual(
            evaluation.optimization_score,
            evaluation.penetration.score * evaluation.functional_preservation.score,
        )

    def test_invalid_candidate_reports_filters_without_descriptor_scoring(self):
        evaluation = CandidateEvaluator().evaluate_sequence("AAAAA")

        self.assertFalse(evaluation.is_valid)
        self.assertIn("too_many_edits", evaluation.failed_filters)
        self.assertIsNone(evaluation.penetration)
        self.assertEqual(evaluation.optimization_score, 0.0)

    def test_locked_position_policy_flows_into_candidate_evaluation(self):
        search_space = SearchSpaceRules.matrixyl_default(locked_positions=(0,))
        evaluation = CandidateEvaluator(search_space=search_space).evaluate_sequence("RTTKS")

        self.assertFalse(evaluation.is_valid)
        self.assertIn("locked_position_0", evaluation.failed_filters)

    def test_candidate_to_dict_is_flattenable_for_csv_rows(self):
        row = CandidateEvaluator().evaluate_sequence("KSTKS").to_dict()

        self.assertEqual(row["sequence"], "KSTKS")
        self.assertIsInstance(row["search_validation"], dict)
        self.assertIsInstance(row["functional_preservation"], dict)
        self.assertIsInstance(row["penetration"], dict)


if __name__ == "__main__":
    unittest.main()
