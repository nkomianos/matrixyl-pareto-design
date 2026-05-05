"""
Utilities for preparing and running GPU-backed PLM scoring jobs.
"""

import csv
from dataclasses import asdict, dataclass
from itertools import combinations, product
from pathlib import Path
from typing import Iterable, Optional

from .oracle import DEFAULT_MODEL_NAMES, DEFAULT_REFERENCE_SEQUENCE, model_specs_from_names
from .search_space import SearchSpaceRules


@dataclass(frozen=True)
class TorchStatus:
    torch_available: bool
    cuda_available: bool
    device_count: int
    device_names: tuple[str, ...]
    torch_version: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)


def detect_torch_status() -> TorchStatus:
    try:
        import torch
    except ImportError:
        return TorchStatus(False, False, 0, ())

    cuda_available = torch.cuda.is_available()
    device_count = torch.cuda.device_count() if cuda_available else 0
    device_names = (
        tuple(torch.cuda.get_device_name(index) for index in range(device_count))
        if cuda_available
        else ()
    )
    return TorchStatus(
        torch_available=True,
        cuda_available=cuda_available,
        device_count=device_count,
        device_names=device_names,
        torch_version=torch.__version__,
    )


def load_unique_sequences_from_csv(path: str | Path, column: str = "sequence") -> list[str]:
    with Path(path).open() as handle:
        sequences = {
            row[column].strip().upper()
            for row in csv.DictReader(handle)
            if row.get(column, "").strip()
        }
    return sorted(sequences)


def enumerate_fixed_length_search_space(
    rules: Optional[SearchSpaceRules] = None,
) -> list[str]:
    rules = rules or SearchSpaceRules.matrixyl_default()
    reference = rules.reference_sequence
    if rules.min_length != len(reference) or rules.max_length != len(reference):
        raise ValueError("Only fixed-length substitution spaces can be enumerated")

    sequences = {reference}
    mutable_positions = [
        position
        for position in range(len(reference))
        if position not in set(rules.locked_positions)
    ]
    for edit_count in range(1, rules.max_edit_distance + 1):
        for positions in combinations(mutable_positions, edit_count):
            residue_options = [
                [
                    residue
                    for residue in rules.allowed_residues_at(position)
                    if residue != reference[position]
                ]
                for position in positions
            ]
            for replacements in product(*residue_options):
                candidate = list(reference)
                for position, residue in zip(positions, replacements):
                    candidate[position] = residue
                sequence = "".join(candidate)
                if rules.validate(sequence).is_valid:
                    sequences.add(sequence)
    return sorted(sequences)


def write_sequence_manifest(path: str | Path, sequences: Iterable[str]) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence"])
        writer.writeheader()
        for sequence in sorted(set(sequences)):
            writer.writerow({"sequence": sequence})


def build_preflight_report(
    *,
    sequences: Iterable[str],
    output_dir: str | Path,
    model_names: Iterable[str] = DEFAULT_MODEL_NAMES,
    reference_sequence: str = DEFAULT_REFERENCE_SEQUENCE,
    batch_size: int = 16,
    cache_dir: str | Path = "results/phase8_plm_cache",
) -> dict:
    sequence_list = sorted(set(sequences))
    model_specs = model_specs_from_names(tuple(model_names))
    torch_status = detect_torch_status()
    output_path = Path(output_dir)
    return {
        "sequence_count": len(sequence_list),
        "reference_sequence": reference_sequence,
        "batch_size": batch_size,
        "models": [spec.to_dict() for spec in model_specs],
        "cache_dir": str(cache_dir),
        "output_dir": str(output_path),
        "expected_scores_csv": str(output_path / "plm_scores.csv"),
        "expected_metadata_json": str(output_path / "run_metadata.json"),
        "torch": torch_status.to_dict(),
    }
