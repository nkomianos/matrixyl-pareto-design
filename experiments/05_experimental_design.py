#!/usr/bin/env python
"""
Phase 7: generate experimental validation design artifacts.

Example:
    python experiments/05_experimental_design.py \
        --candidates results/phase5_candidate_analysis/candidate_summary.csv \
        --compactness results/phase6_structural_validation/conformer_summary.csv \
        --output results/phase7_experimental_design
"""

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.analysis import load_candidate_rows, write_dict_rows
from src.experimental_design import (
    assay_matrix,
    draft_protocol_markdown,
    experimental_controls,
    load_compactness_rows,
    select_synthesis_candidates,
)


def current_git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def main(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    candidates = load_candidate_rows(args.candidates)
    compactness = load_compactness_rows(args.compactness)
    synthesis_candidates = select_synthesis_candidates(
        candidates,
        compactness,
        max_candidates=args.max_candidates,
    )
    controls = experimental_controls()
    assays = assay_matrix()
    protocol = draft_protocol_markdown(synthesis_candidates)

    write_dict_rows(
        output_dir / "synthesis_candidates.csv",
        [candidate.to_dict() for candidate in synthesis_candidates],
    )
    write_dict_rows(output_dir / "experimental_controls.csv", controls)
    write_dict_rows(output_dir / "assay_matrix.csv", assays)
    (output_dir / "protocol_draft.md").write_text(protocol)

    metadata = {
        "candidate_source": args.candidates,
        "compactness_source": args.compactness,
        "max_candidates": args.max_candidates,
        "selected_sequences": [candidate.sequence for candidate in synthesis_candidates],
        "control_count": len(controls),
        "assay_count": len(assays),
        "git_commit": current_git_commit(),
        "python": sys.version,
        "platform": platform.platform(),
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))

    print("Experimental validation design complete")
    print("Selected candidates:", ", ".join(metadata["selected_sequences"]))
    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 7 experimental validation design")
    parser.add_argument("--candidates", default="results/phase5_candidate_analysis/candidate_summary.csv")
    parser.add_argument("--compactness", default="results/phase6_structural_validation/conformer_summary.csv")
    parser.add_argument("--output", default="results/phase7_experimental_design")
    parser.add_argument("--max-candidates", type=int, default=3)
    main(parser.parse_args())
