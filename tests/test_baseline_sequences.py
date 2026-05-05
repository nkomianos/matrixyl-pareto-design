from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
SEQUENCE_DIR = REPO_ROOT / "data" / "sequences"


def read_fasta_sequence(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith(">")
    )


class BaselineSequenceTests(unittest.TestCase):
    def test_matrixyl_baseline_uses_kttks_core(self):
        self.assertEqual(read_fasta_sequence(SEQUENCE_DIR / "baseline_matrixyl.fasta"), "KTTKS")
        self.assertEqual(read_fasta_sequence(SEQUENCE_DIR / "matrixyl_core.fasta"), "KTTKS")

    def test_collagen_like_control_preserves_previous_sequence(self):
        self.assertEqual(
            read_fasta_sequence(SEQUENCE_DIR / "collagen_like_control.fasta"),
            "GPKGDPGA",
        )


if __name__ == "__main__":
    unittest.main()
