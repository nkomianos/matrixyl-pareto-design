#!/usr/bin/env python
"""
Phase 8: score candidate sequences with real protein language models.

Run this on the GPU host after `experiments/06_plm_gpu_preflight.py`.
"""

import argparse
import csv
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.oracle import DEFAULT_MODEL_NAMES, EnsembleOracle
from src.plm_pipeline import detect_torch_status, load_unique_sequences_from_csv


def current_git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def write_scores(path: Path, sequences: list[str], scores, uncertainties) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sequence", "plm_similarity_score", "plm_uncertainty"],
        )
        writer.writeheader()
        for sequence, score, uncertainty in zip(sequences, scores, uncertainties):
            writer.writerow(
                {
                    "sequence": sequence,
                    "plm_similarity_score": float(score),
                    "plm_uncertainty": float(uncertainty),
                }
            )


def score_with_sequential_models(args, sequences: list[str], model_names: tuple[str, ...]):
    import numpy as np

    model_scores = []
    devices = []
    for model_name in model_names:
        oracle = EnsembleOracle(
            models=[model_name],
            reference_sequence=args.reference_sequence,
            device=args.device,
            batch_size=args.batch_size,
            cache_dir=args.cache_dir,
            local_files_only=args.local_files_only,
        )
        scores, _ = oracle.score_binding_affinity(sequences)
        model_scores.append(scores)
        devices.append(oracle.device)
        torch_module = oracle.torch
        del oracle
        if torch_module.cuda.is_available():
            torch_module.cuda.empty_cache()

    stacked = np.vstack(model_scores)
    return stacked.mean(axis=0), stacked.std(axis=0), devices


def main(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    sequences = load_unique_sequences_from_csv(args.candidates)
    if args.limit:
        sequences = sequences[: args.limit]

    model_names = tuple(args.models or DEFAULT_MODEL_NAMES)
    torch_status = detect_torch_status()
    if args.require_cuda and not torch_status.cuda_available:
        raise SystemExit("CUDA is required for this run, but torch.cuda.is_available() is false.")

    if args.load_models_together:
        oracle = EnsembleOracle(
            models=list(model_names),
            reference_sequence=args.reference_sequence,
            device=args.device,
            batch_size=args.batch_size,
            cache_dir=args.cache_dir,
            local_files_only=args.local_files_only,
        )
        scores, uncertainties = oracle.score_binding_affinity(sequences)
        devices = [oracle.device]
    else:
        scores, uncertainties, devices = score_with_sequential_models(args, sequences, model_names)

    scores_path = output_dir / "plm_scores.csv"
    write_scores(scores_path, sequences, scores, uncertainties)

    metadata = {
        "candidate_source": args.candidates,
        "sequence_count": len(sequences),
        "reference_sequence": args.reference_sequence,
        "models": list(model_names),
        "device": oracle.device,
        "batch_size": args.batch_size,
        "cache_dir": args.cache_dir,
        "local_files_only": args.local_files_only,
        "load_models_together": args.load_models_together,
        "devices": devices,
        "scores_csv": str(scores_path),
        "torch": torch_status.to_dict(),
        "git_commit": current_git_commit(),
        "python": sys.version,
        "platform": platform.platform(),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))

    print("PLM scoring complete")
    print(f"Sequences scored: {len(sequences)}")
    print(f"Scores written to: {scores_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Score candidates with PLM oracle")
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--output", default="results/phase8_plm_gpu_scores")
    parser.add_argument("--cache-dir", default="results/phase8_plm_cache")
    parser.add_argument("--reference-sequence", default="KTTKS")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--models", action="append")
    parser.add_argument("--device", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--local-files-only", action="store_true")
    parser.add_argument("--require-cuda", action="store_true")
    parser.add_argument(
        "--load-models-together",
        action="store_true",
        help="Load all PLMs into memory at once. By default, models are scored sequentially to reduce GPU memory pressure.",
    )
    main(parser.parse_args())
