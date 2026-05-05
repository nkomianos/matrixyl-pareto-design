import unittest

from src.analysis import CandidateRow
from src.experimental_design import (
    CompactnessRow,
    assay_matrix,
    draft_protocol_markdown,
    experimental_controls,
    select_synthesis_candidates,
)


def candidate(sequence, penetration, functional, scalar, edit_distance):
    return CandidateRow(
        sequence=sequence,
        penetration_objective=penetration,
        functional_objective=functional,
        synthesis_objective=1.0,
        scalar_optimization_score=scalar,
        edit_distance=edit_distance,
        molecular_weight=500.0,
        tpsa=220.0,
        logp=-3.5,
        hbd=8,
        hba=9,
        rotatable_bonds=12,
        formal_charge=0,
    )


class ExperimentalDesignTests(unittest.TestCase):
    def test_select_synthesis_candidates_balances_penetration_and_function(self):
        rows = [
            candidate("PTTPS", 0.68, 0.66, 0.45, 2),
            candidate("KTTPS", 0.55, 0.83, 0.46, 1),
            candidate("KTTPP", 0.67, 0.67, 0.44, 2),
            candidate("KTTKS", 0.39, 1.00, 0.39, 0),
        ]
        compactness = {
            "PTTPS": CompactnessRow("PTTPS", "PTTPS", "rdkit_fasta", 4.83, 4.2, 5.3),
            "KTTPS": CompactnessRow("KTTPS", "KTTPS", "rdkit_fasta", 4.99, 4.5, 5.6),
        }

        selected = select_synthesis_candidates(rows, compactness, max_candidates=3)
        sequences = [row.sequence for row in selected]

        self.assertEqual(sequences[0], "PTTPS")
        self.assertIn("KTTPS", sequences)
        self.assertEqual(len(sequences), len(set(sequences)))
        self.assertEqual(selected[0].mean_radius_of_gyration, 4.83)

    def test_selection_rejects_too_small_panel(self):
        with self.assertRaises(ValueError):
            select_synthesis_candidates([], {}, max_candidates=1)

    def test_controls_include_core_commercial_and_safety_controls(self):
        controls = {row["control"] for row in experimental_controls()}

        self.assertIn("matrixyl_core", controls)
        self.assertIn("pal_kttks", controls)
        self.assertIn("vehicle_blank", controls)
        self.assertIn("cytotoxicity_control", controls)

    def test_assay_matrix_covers_permeation_function_and_safety(self):
        assays = [row["assay"] for row in assay_matrix()]

        self.assertTrue(any("Franz" in assay for assay in assays))
        self.assertTrue(any("collagen" in assay.lower() for assay in assays))
        self.assertTrue(any("cytotoxicity" in assay.lower() for assay in assays))
        self.assertTrue(any("epidermis irritation" in assay.lower() for assay in assays))

    def test_protocol_mentions_quantification_controls_and_limits(self):
        selected = [
            candidate("PTTPS", 0.68, 0.66, 0.45, 2),
            candidate("KTTPS", 0.55, 0.83, 0.46, 1),
        ]
        synthesis_rows = select_synthesis_candidates(selected, {}, max_candidates=2)
        protocol = draft_protocol_markdown(synthesis_rows)

        self.assertIn("LC-MS/MS or HPLC", protocol)
        self.assertIn("Pal-KTTKS", protocol)
        self.assertIn("Do not claim retained biological activity from permeation alone", protocol)
        self.assertIn("reconstructed human epidermis", protocol)


if __name__ == "__main__":
    unittest.main()
