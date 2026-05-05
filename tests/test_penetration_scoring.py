from pathlib import Path
import unittest

from src.chemistry import read_smiles
from src.constraints import PenetrationScore, PhysicochemicalCalculator


REPO_ROOT = Path(__file__).resolve().parents[1]
MOLECULE_DIR = REPO_ROOT / "data" / "molecules"


class PenetrationScoringTests(unittest.TestCase):
    def setUp(self):
        self.calculator = PhysicochemicalCalculator()

    def test_sequence_score_uses_rdkit_descriptors(self):
        score = self.calculator.score_sequence_exact("KTTKS", name="Matrixyl core")

        self.assertIsInstance(score, PenetrationScore)
        self.assertEqual(score.descriptors.name, "Matrixyl core")
        self.assertAlmostEqual(score.descriptors.molecular_weight, 563.653, places=3)
        self.assertAlmostEqual(score.descriptors.tpsa, 292.45, places=2)
        self.assertAlmostEqual(score.descriptors.logp, -4.651, places=3)
        self.assertGreaterEqual(score.score, 0.0)
        self.assertLessEqual(score.score, 1.0)

    def test_legacy_public_descriptor_methods_are_exact(self):
        self.assertAlmostEqual(self.calculator.compute_molecular_weight("KTTKS"), 563.653, places=3)
        self.assertAlmostEqual(self.calculator.compute_tpsa("KTTKS"), 292.45, places=2)
        self.assertAlmostEqual(self.calculator.compute_logp("KTTKS"), -4.651, places=3)
        self.assertEqual(
            self.calculator.compute_penetration_score("KTTKS"),
            self.calculator.score_sequence_exact("KTTKS").score,
        )

    def test_score_exposes_component_penalties(self):
        score = self.calculator.score_sequence_exact("KTTKS")

        self.assertEqual(
            set(score.penalties),
            {
                "tpsa",
                "mw",
                "logp",
                "hbd",
                "hba",
                "rotatable_bonds",
                "formal_charge",
            },
        )
        self.assertEqual(score.penalties["tpsa"], 1.0)
        self.assertEqual(score.penalties["logp"], 1.0)
        self.assertEqual(score.penalties["formal_charge"], 0.0)

    def test_palmitoylation_score_uses_same_contract(self):
        smiles, name = read_smiles(MOLECULE_DIR / "matrixyl_palmitoylated.smi")

        core_score = self.calculator.score_sequence_exact("KTTKS")
        palmitoylated_score = self.calculator.score_smiles_exact(
            smiles,
            name=name,
            source="pubchem_smiles",
        )

        self.assertEqual(palmitoylated_score.descriptors.name, "Pal-KTTKS")
        self.assertAlmostEqual(palmitoylated_score.descriptors.molecular_weight, 802.068, places=3)
        self.assertAlmostEqual(palmitoylated_score.descriptors.logp, 0.988, places=3)
        self.assertGreater(
            palmitoylated_score.score,
            core_score.score,
            "Palmitoylation should improve the LogP component despite MW/TPSA liabilities.",
        )

    def test_heuristic_score_remains_available_for_prefiltering(self):
        heuristic = self.calculator.compute_heuristic_penetration_score("KTTKS")
        exact = self.calculator.compute_penetration_score("KTTKS")

        self.assertGreaterEqual(heuristic, 0.0)
        self.assertLessEqual(heuristic, 1.0)
        self.assertNotEqual(heuristic, exact)


if __name__ == "__main__":
    unittest.main()
