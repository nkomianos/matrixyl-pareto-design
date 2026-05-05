import csv
from pathlib import Path
import tempfile
import unittest

from src.plm_pipeline import (
    build_preflight_report,
    enumerate_fixed_length_search_space,
    load_unique_sequences_from_csv,
    write_sequence_manifest,
)
from src.search_space import SearchSpaceRules


class PLMPipelineTests(unittest.TestCase):
    def test_enumerates_matrixyl_max_edit_two_space(self):
        sequences = enumerate_fixed_length_search_space(SearchSpaceRules.matrixyl_default())

        self.assertEqual(len(sequences), 3706)
        self.assertIn("KTTKS", sequences)
        self.assertIn("PTTPS", sequences)

    def test_enumeration_respects_locked_positions(self):
        rules = SearchSpaceRules.matrixyl_default(locked_positions=(0,))
        sequences = enumerate_fixed_length_search_space(rules)

        self.assertTrue(all(sequence[0] == "K" for sequence in sequences))

    def test_manifest_round_trip_deduplicates_sequences(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "manifest.csv"
            write_sequence_manifest(path, ["KTTKS", "PTTPS", "KTTKS"])

            loaded = load_unique_sequences_from_csv(path)

        self.assertEqual(loaded, ["KTTKS", "PTTPS"])

    def test_load_unique_sequences_from_existing_csv_shape(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "candidates.csv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sequence", "score"])
                writer.writeheader()
                writer.writerow({"sequence": "KTTKS", "score": "1"})
                writer.writerow({"sequence": "KTTKS", "score": "1"})
                writer.writerow({"sequence": "PTTPS", "score": "2"})

            self.assertEqual(load_unique_sequences_from_csv(path), ["KTTKS", "PTTPS"])

    def test_preflight_report_is_cpu_safe(self):
        report = build_preflight_report(
            sequences=["KTTKS", "PTTPS"],
            output_dir="results/phase8_plm_gpu",
            model_names=["facebook/esm2_t33_650M_UR50D", "rostlab/prot_bert"],
            batch_size=4,
        )

        self.assertEqual(report["sequence_count"], 2)
        self.assertEqual(report["reference_sequence"], "KTTKS")
        self.assertEqual(report["batch_size"], 4)
        self.assertEqual(report["models"][1]["input_format"], "spaced")
        self.assertIn("torch_available", report["torch"])


if __name__ == "__main__":
    unittest.main()
