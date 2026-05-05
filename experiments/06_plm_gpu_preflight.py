#!/usr/bin/env python
"""
Phase 8: CPU-safe preflight for the GPU PLM scoring run.

This script does not load model weights. It prepares the candidate manifest and
reports whether the current machine has torch/CUDA available.
"""

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.oracle import DEFAULT_MODEL_NAMES
from src.plm_pipeline import (
    build_preflight_report,
    enumerate_fixed_length_search_space,
    load_unique_sequences_from_csv,
    write_sequence_manifest,
)
from src.search_space import SearchSpaceRules


def current_git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def load_or_enumerate_sequences(args) -> list[str]:
    if args.enumerate_search_space:
        rules = SearchSpaceRules.matrixyl_default()
        return enumerate_fixed_length_search_space(rules)
    return load_unique_sequences_from_csv(args.candidates)


def main(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    sequences = load_or_enumerate_sequences(args)
    manifest_path = output_dir / "candidate_manifest.csv"
    write_sequence_manifest(manifest_path, sequences)

    model_names = tuple(args.models or DEFAULT_MODEL_NAMES)
    report = build_preflight_report(
        sequences=sequences,
        output_dir=output_dir,
        model_names=model_names,
        reference_sequence=args.reference_sequence,
        batch_size=args.batch_size,
        cache_dir=args.cache_dir,
    )
    report.update(
        {
            "candidate_manifest": str(manifest_path),
            "source_candidates": "enumerated_search_space" if args.enumerate_search_space else args.candidates,
            "expected_scores_csv": str(output_dir / "gpu_scores" / "plm_scores.csv"),
            "expected_metadata_json": str(output_dir / "gpu_scores" / "run_metadata.json"),
            "git_commit": current_git_commit(),
            "python": sys.version,
            "platform": platform.platform(),
            "gpu_command": (
                "PYTHONDONTWRITEBYTECODE=1 python experiments/07_score_plm_oracle.py "
                f"--candidates {manifest_path} "
                f"--output {output_dir / 'gpu_scores'} "
                f"--cache-dir {args.cache_dir} "
                f"--batch-size {args.batch_size} "
                f"--reference-sequence {args.reference_sequence} "
                "--require-cuda "
                + " ".join(f"--models {model}" for model in model_names)
            ),
        }
    )
    (output_dir / "preflight_report.json").write_text(json.dumps(report, indent=2))

    print("PLM GPU preflight complete")
    print(f"Sequences prepared: {report['sequence_count']}")
    print(f"CUDA available: {report['torch']['cuda_available']}")
    print(f"Candidate manifest: {manifest_path}")
    print("GPU command:")
    print(report["gpu_command"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare PLM GPU scoring run")
    parser.add_argument("--candidates", default="results/phase5_candidate_analysis/candidate_summary.csv")
    parser.add_argument("--output", default="results/phase8_plm_preflight")
    parser.add_argument("--cache-dir", default="results/phase8_plm_cache")
    parser.add_argument("--reference-sequence", default="KTTKS")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--models", action="append")
    parser.add_argument(
        "--enumerate-search-space",
        action="store_true",
        help="Prepare all fixed-length max-edit-2 Matrixyl candidates instead of reading input CSV.",
    )
    main(parser.parse_args())
