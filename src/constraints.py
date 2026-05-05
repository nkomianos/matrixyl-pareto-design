"""
Physicochemical property calculators for peptide skin penetration.
All properties normalized to [0, 1] for fitness functions.
"""

from dataclasses import asdict, dataclass
import logging
import numpy as np
from typing import Dict

from .chemistry import (
    MolecularDescriptors,
    descriptors_from_sequence,
    descriptors_from_smiles,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PenetrationScore:
    """
    Normalized permeability-oriented score with raw descriptor context.

    The score is a design heuristic, not a measured permeability coefficient.
    Higher is better, and penalties are normalized to [0, 1].
    """

    score: float
    weighted_penalty: float
    penalties: Dict[str, float]
    descriptors: MolecularDescriptors

    def to_dict(self) -> dict:
        data = asdict(self)
        data["descriptors"] = self.descriptors.to_dict()
        return data


class PhysicochemicalCalculator:
    """
    Compute peptide properties relevant to transdermal delivery.
    """

    # Fast residue-table prefilter ranges. Final paper scoring should use
    # EXACT_TARGETS through score_descriptors().
    TARGETS = {
        'tpsa': (80, 140),          # A^2; optimal < 140
        'mw': (200, 500),           # Da; optimal < 500
        'logp': (1.0, 3.0),         # optimal 1-3
        'hbd': (0, 5),              # H-bond donors; minimize
        'hba': (0, 10),             # H-bond acceptors; minimize
    }

    # RDKit descriptor targets used for publication-grade candidate scoring.
    EXACT_TARGETS = {
        'tpsa': (0, 140),            # A^2; lower is generally better for passive diffusion
        'mw': (0, 500),              # Da
        'logp': (1.0, 3.0),          # topical/transdermal balance window
        'hbd': (0, 5),
        'hba': (0, 10),
        'rotatable_bonds': (0, 10),
        'formal_charge': (-1, 1),
    }

    EXACT_WEIGHTS = {
        'tpsa': 0.25,
        'mw': 0.25,
        'logp': 0.20,
        'hbd': 0.10,
        'hba': 0.10,
        'rotatable_bonds': 0.05,
        'formal_charge': 0.05,
    }

    @staticmethod
    def compute_amino_acid_properties(sequence: str) -> Dict[str, float]:
        """
        Compute properties directly from amino acid sequence (no SMILES needed).
        This avoids lossy conversion to SMILES.

        Args:
            sequence: amino acid sequence (string)

        Returns:
            dict with properties: mw, hbd, hba, logp, charge, flexibility
        """
        # Amino acid molecular weights (monoisotopic, excluding water)
        aa_weights = {
            'A': 71.04, 'R': 156.10, 'N': 114.04, 'D': 115.03,
            'C': 103.01, 'Q': 128.06, 'E': 129.04, 'G': 57.02,
            'H': 137.06, 'I': 113.08, 'L': 113.08, 'K': 128.09,
            'M': 131.04, 'F': 147.07, 'P': 97.05, 'S': 87.03,
            'T': 101.05, 'W': 186.08, 'Y': 163.06, 'V': 99.07,
        }

        # H-bond donors (including backbone NH groups)
        aa_hbd = {
            'A': 1, 'R': 4, 'N': 2, 'D': 1, 'C': 1, 'Q': 2, 'E': 1,
            'G': 1, 'H': 2, 'I': 1, 'L': 1, 'K': 2, 'M': 1, 'F': 1,
            'P': 0, 'S': 2, 'T': 2, 'W': 2, 'Y': 2, 'V': 1,
        }

        # H-bond acceptors
        aa_hba = {
            'A': 1, 'R': 2, 'N': 2, 'D': 3, 'C': 1, 'Q': 3, 'E': 3,
            'G': 1, 'H': 2, 'I': 1, 'L': 1, 'K': 1, 'M': 1, 'F': 1,
            'P': 1, 'S': 2, 'T': 2, 'W': 2, 'Y': 2, 'V': 1,
        }

        # Logp per residue (Kyte-Doolittle hydrophobicity, approximate)
        aa_logp = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
        }

        # Charge at pH 7.4 (pKa values)
        aa_charge = {
            'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0, 'Q': 0, 'E': -1,
            'G': 0, 'H': 0.1, 'I': 0, 'L': 0, 'K': 1, 'M': 0, 'F': 0,
            'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0,
        }

        # Flexibility (normalized disorder propensity)
        aa_disorder = {
            'A': 0.5, 'R': 0.8, 'N': 0.8, 'D': 0.9, 'C': 0.3,
            'Q': 0.8, 'E': 0.9, 'G': 0.9, 'H': 0.5, 'I': 0.3,
            'L': 0.3, 'K': 0.8, 'M': 0.4, 'F': 0.3, 'P': 0.9,
            'S': 0.8, 'T': 0.8, 'W': 0.3, 'Y': 0.6, 'V': 0.3,
        }

        n = len(sequence)
        if n == 0:
            return None

        # Sum properties
        total_mw = sum(aa_weights.get(aa, 0) for aa in sequence)
        total_mw += (n - 1) * 18.015  # Add water molecules between residues

        total_hbd = sum(aa_hbd.get(aa, 0) for aa in sequence)
        total_hba = sum(aa_hba.get(aa, 0) for aa in sequence)
        avg_logp = np.mean([aa_logp.get(aa, 0) for aa in sequence])
        total_charge = sum(aa_charge.get(aa, 0) for aa in sequence)
        avg_disorder = np.mean([aa_disorder.get(aa, 0) for aa in sequence])

        return {
            'mw': total_mw,
            'hbd': total_hbd,
            'hba': total_hba,
            'logp': avg_logp,
            'charge': total_charge,
            'flexibility': avg_disorder,
            'length': n,
        }

    @staticmethod
    def penalty_score(value: float, target_range: tuple) -> float:
        """
        Compute penalty for being outside target range.
        Returns 0 if in range, increases quadratically outside.

        Args:
            value: measured value
            target_range: (min, max) tuple

        Returns:
            penalty: 0 (perfect) to 1 (very bad)
        """
        min_val, max_val = target_range
        if min_val <= value <= max_val:
            return 0.0

        if value < min_val:
            penalty = ((min_val - value) / min_val) ** 2
        else:
            penalty = ((value - max_val) / max_val) ** 2

        return min(penalty, 1.0)

    def compute_heuristic_penetration_score(self, sequence: str) -> float:
        """
        Fast residue-table prefilter for skin penetration potential.

        This is retained for rough exploration only. Use
        score_sequence_exact() or score_smiles_exact() for reported results.
        Higher is better. Range [0, 1].

        Based on:
        - TPSA < 140 A^2
        - MW < 500 Da
        - LogP 1-3
        - Minimize H-bond donors/acceptors
        - Flexibility (not too rigid, not too disordered)
        """
        props = self.compute_amino_acid_properties(sequence)
        if props is None:
            return 0.0

        # Individual penalties
        tpsa_penalty = self.penalty_score(props['hba'] * 20, self.TARGETS['tpsa'])  # Approximate TPSA
        mw_penalty = self.penalty_score(props['mw'], self.TARGETS['mw'])
        logp_penalty = self.penalty_score(props['logp'], self.TARGETS['logp'])

        hbd_penalty = min(props['hbd'] / 5.0, 1.0)  # Penalize >5 donors
        hba_penalty = min(props['hba'] / 10.0, 1.0)  # Penalize >10 acceptors

        # Flexibility: target middle ground (0.3-0.7)
        flexibility_penalty = abs(props['flexibility'] - 0.5) / 0.5

        # Combine penalties (weighted)
        total_penalty = (
            0.25 * tpsa_penalty +
            0.25 * mw_penalty +
            0.15 * logp_penalty +
            0.15 * hbd_penalty +
            0.10 * hba_penalty +
            0.10 * flexibility_penalty
        )

        return 1.0 - min(total_penalty, 1.0)

    def score_descriptors(self, descriptors: MolecularDescriptors) -> PenetrationScore:
        """
        Score exact RDKit descriptors against topical-delivery design targets.

        The output intentionally exposes every component penalty so candidates
        can be audited instead of treated as black-box scalar scores.
        """
        raw_values = {
            'tpsa': descriptors.tpsa,
            'mw': descriptors.molecular_weight,
            'logp': descriptors.logp,
            'hbd': descriptors.hbd,
            'hba': descriptors.hba,
            'rotatable_bonds': descriptors.rotatable_bonds,
            'formal_charge': descriptors.formal_charge,
        }
        penalties = {
            key: self.penalty_score(value, self.EXACT_TARGETS[key])
            for key, value in raw_values.items()
        }
        weighted_penalty = sum(
            self.EXACT_WEIGHTS[key] * penalties[key]
            for key in self.EXACT_WEIGHTS
        )
        bounded_penalty = min(max(weighted_penalty, 0.0), 1.0)
        return PenetrationScore(
            score=1.0 - bounded_penalty,
            weighted_penalty=bounded_penalty,
            penalties=penalties,
            descriptors=descriptors,
        )

    def score_sequence_exact(self, sequence: str, name: str = None) -> PenetrationScore:
        """Compute exact RDKit descriptors for a canonical peptide sequence and score them."""
        return self.score_descriptors(descriptors_from_sequence(sequence, name=name))

    def score_smiles_exact(self, smiles: str, name: str, source: str = "smiles") -> PenetrationScore:
        """Compute exact RDKit descriptors for a molecule SMILES and score them."""
        return self.score_descriptors(descriptors_from_smiles(smiles, name=name, source=source))

    def compute_penetration_score(self, sequence: str) -> float:
        """
        Publication-grade penetration score for canonical peptide sequences.

        This method now uses RDKit descriptors through score_sequence_exact().
        Use compute_heuristic_penetration_score() for the legacy residue-table
        approximation.
        """
        return self.score_sequence_exact(sequence).score

    def compute_tpsa(self, sequence: str) -> float:
        """Return RDKit TPSA for an unmodified canonical peptide sequence."""
        return self.score_sequence_exact(sequence).descriptors.tpsa

    def compute_molecular_weight(self, sequence: str) -> float:
        """Return RDKit molecular weight for an unmodified canonical peptide sequence."""
        return self.score_sequence_exact(sequence).descriptors.molecular_weight

    def compute_logp(self, sequence: str) -> float:
        """Return RDKit Crippen MolLogP for an unmodified canonical peptide sequence."""
        return self.score_sequence_exact(sequence).descriptors.logp

    def compute_gyration_radius(self, sequence: str) -> float:
        """
        Estimate radius of gyration from sequence length.
        Simplified: Rg ≈ sqrt(5/3) * L^0.6 / sqrt(12π) (Gaussian chain)

        Returns normalized value [0, 1] where 1 = very flexible/extended
        """
        length = len(sequence)
        # For typical peptides: Rg ~ 0.2-0.5 nm for 5-20 residues
        # Normalize: 0 = very compact, 1 = very extended
        rg_estimate = (0.2 * length) ** 0.6
        return min(rg_estimate / 2.0, 1.0)  # Cap at 1.0
