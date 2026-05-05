from pathlib import Path
import unittest

from src.analysis import (
    CandidateRow,
    baseline_descriptor_rows,
    descriptor_warning_flags,
    mutation_enrichment,
    position_frequency_matrix,
    priority_label,
    summarize_candidates,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def candidate_row(**overrides):
    values = {
        "sequence": "PTTPS",
        "penetration_objective": 0.68,
        "functional_objective": 0.67,
        "synthesis_objective": 1.0,
        "scalar_optimization_score": 0.45,
        "edit_distance": 2,
        "molecular_weight": 501.5,
        "tpsa": 217.6,
        "logp": -3.9,
        "hbd": 8,
        "hba": 9,
        "rotatable_bonds": 11,
        "formal_charge": 0,
    }
    values.update(overrides)
    return CandidateRow(**values)


class AnalysisTests(unittest.TestCase):
    def test_position_frequency_matrix(self):
        rows = position_frequency_matrix(["KTTKS", "PTTPS"], reference_sequence="KTTKS")

        position_1 = [row for row in rows if row["position"] == 1]
        self.assertEqual({row["residue"] for row in position_1}, {"K", "P"})
        self.assertTrue(all(row["frequency"] == 0.5 for row in position_1))

    def test_mutation_enrichment(self):
        rows = mutation_enrichment(["KTTKS", "PTTPS", "KTTPP"], reference_sequence="KTTKS")
        mutations = {row["mutation"]: row for row in rows}

        self.assertIn("K1P", mutations)
        self.assertIn("K4P", mutations)
        self.assertIn("S5P", mutations)
        self.assertEqual(mutations["K4P"]["count"], 2)

    def test_descriptor_warning_flags(self):
        flags = descriptor_warning_flags(candidate_row())

        self.assertIn("mw_above_500", flags)
        self.assertIn("tpsa_above_140", flags)
        self.assertIn("logp_below_1", flags)
        self.assertIn("hbd_above_5", flags)

    def test_priority_label(self):
        self.assertEqual(priority_label(candidate_row()), "high_penetration_tradeoff")
        self.assertEqual(
            priority_label(
                candidate_row(
                    sequence="KTTPS",
                    penetration_objective=0.55,
                    functional_objective=0.83,
                )
            ),
            "balanced",
        )
        self.assertEqual(
            priority_label(
                candidate_row(
                    sequence="KTTKS",
                    penetration_objective=0.39,
                    functional_objective=1.0,
                )
            ),
            "function_anchor",
        )

    def test_summarize_candidates_sorts_by_penetration_then_function(self):
        rows = summarize_candidates(
            [
                candidate_row(sequence="KTTKS", penetration_objective=0.39, functional_objective=1.0),
                candidate_row(sequence="PTTPS", penetration_objective=0.68, functional_objective=0.67),
            ]
        )

        self.assertEqual(rows[0]["sequence"], "PTTPS")
        self.assertIn("priority_label", rows[0])
        self.assertIn("warning_flags", rows[0])

    def test_baseline_descriptor_rows_include_core_and_palmitoylated(self):
        rows = baseline_descriptor_rows(REPO_ROOT / "data" / "molecules" / "matrixyl_palmitoylated.smi")
        names = {row["name"] for row in rows}

        self.assertEqual(names, {"Matrixyl core", "Pal-KTTKS"})
        pal = next(row for row in rows if row["name"] == "Pal-KTTKS")
        self.assertGreater(pal["molecular_weight"], 800)


if __name__ == "__main__":
    unittest.main()
