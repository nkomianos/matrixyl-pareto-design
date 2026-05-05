from pathlib import Path
import unittest

from src.chemistry import read_smiles
from src.structure import summarize_sequence_conformers, summarize_smiles_conformers


REPO_ROOT = Path(__file__).resolve().parents[1]


class StructureTests(unittest.TestCase):
    def test_sequence_conformer_summary_is_deterministic(self):
        first = summarize_sequence_conformers("KTTKS", num_conformers=3, seed=7)
        second = summarize_sequence_conformers("KTTKS", num_conformers=3, seed=7)

        self.assertEqual(first.conformer_count, 3)
        self.assertEqual(second.conformer_count, 3)
        self.assertAlmostEqual(first.mean_radius_of_gyration, second.mean_radius_of_gyration)
        self.assertGreater(first.min_radius_of_gyration, 0)
        self.assertLessEqual(first.min_radius_of_gyration, first.max_radius_of_gyration)

    def test_invalid_conformer_count_is_rejected(self):
        with self.assertRaises(ValueError):
            summarize_sequence_conformers("KTTKS", num_conformers=0)

    def test_palmitoylated_smiles_conformers(self):
        smiles, name = read_smiles(REPO_ROOT / "data" / "molecules" / "matrixyl_palmitoylated.smi")
        summary = summarize_smiles_conformers(smiles, name=name, num_conformers=2, seed=11)

        self.assertEqual(summary.name, "Pal-KTTKS")
        self.assertEqual(summary.conformer_count, 2)
        self.assertGreater(summary.mean_radius_of_gyration, 0)
        self.assertTrue(all(metric.force_field in {"MMFF", "UFF", "none"} for metric in summary.metrics))


if __name__ == "__main__":
    unittest.main()
