import json
from pathlib import Path
import unittest

from src.chemistry import (
    descriptors_from_sequence,
    descriptors_from_smiles,
    read_smiles,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
MOLECULE_DIR = REPO_ROOT / "data" / "molecules"
REFERENCE_DIR = REPO_ROOT / "data" / "references"


class ChemistryDescriptorTests(unittest.TestCase):
    def test_matrixyl_core_rdkit_descriptors(self):
        descriptors = descriptors_from_sequence("KTTKS", name="Matrixyl core")

        self.assertEqual(descriptors.formula, "C23H45N7O9")
        self.assertAlmostEqual(descriptors.molecular_weight, 563.653, places=3)
        self.assertAlmostEqual(descriptors.tpsa, 292.45, places=2)
        self.assertEqual(descriptors.hbd, 11)
        self.assertEqual(descriptors.hba, 11)
        self.assertEqual(descriptors.rotatable_bonds, 20)
        self.assertEqual(descriptors.formal_charge, 0)

    def test_palmitoylated_matrixyl_matches_pubchem_reference(self):
        smiles, name = read_smiles(MOLECULE_DIR / "matrixyl_palmitoylated.smi")
        descriptors = descriptors_from_smiles(smiles, name=name, source="pubchem_smiles")
        reference = json.loads(
            (REFERENCE_DIR / "palmitoyl_pentapeptide_4_pubchem.json").read_text()
        )["properties"]

        self.assertEqual(name, "Pal-KTTKS")
        self.assertEqual(descriptors.formula, reference["molecular_formula"])
        self.assertAlmostEqual(
            descriptors.molecular_weight,
            reference["molecular_weight"],
            delta=0.1,
        )
        self.assertAlmostEqual(descriptors.tpsa, reference["tpsa"], delta=1.0)
        self.assertEqual(descriptors.hbd, reference["h_bond_donor_count"])
        self.assertEqual(descriptors.rotatable_bonds, reference["rotatable_bond_count"])

        # RDKit and PubChem use different HBA and logP implementations. Lock the
        # RDKit values here while retaining PubChem values in the reference JSON.
        self.assertEqual(descriptors.hba, 11)
        self.assertAlmostEqual(descriptors.logp, 0.988, places=3)

    def test_invalid_sequence_is_rejected(self):
        with self.assertRaises(ValueError):
            descriptors_from_sequence("KTTKZ")


if __name__ == "__main__":
    unittest.main()
