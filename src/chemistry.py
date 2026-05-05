"""
RDKit-backed molecular descriptors for peptide baselines and candidates.

The sequence-only calculator in constraints.py is a fast heuristic. This module
is the chemistry-grounded path for values that should appear in paper figures.
"""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional


PALMITOYL_KTTKS_SMILES = (
    "CCCCCCCCCCCCCCCC(=O)N[C@@H](CCCCN)C(=O)N[C@@H]([C@@H](C)O)"
    "C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CO)C(=O)O"
)

CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


class RDKitUnavailableError(ImportError):
    """Raised when RDKit-backed descriptor calculation is requested without RDKit."""


@dataclass(frozen=True)
class MolecularDescriptors:
    name: str
    source: str
    formula: str
    molecular_weight: float
    exact_molecular_weight: float
    tpsa: float
    logp: float
    hbd: int
    hba: int
    rotatable_bonds: int
    formal_charge: int
    heavy_atom_count: int
    sequence: Optional[str] = None
    smiles: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)


def _rdkit_modules():
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    except ImportError as exc:
        raise RDKitUnavailableError(
            "RDKit is required for exact molecular descriptors. "
            "Install the project's chemistry dependencies with requirements.txt."
        ) from exc

    return Chem, Crippen, Descriptors, rdMolDescriptors


def molecule_from_sequence(sequence: str):
    """Build an RDKit molecule from an unmodified canonical amino-acid sequence."""
    Chem, _, _, _ = _rdkit_modules()
    normalized = sequence.strip().upper()
    invalid = sorted(set(normalized) - CANONICAL_AMINO_ACIDS)
    if not normalized or invalid:
        raise ValueError(f"Invalid canonical peptide sequence: {sequence!r}")

    mol = Chem.MolFromFASTA(normalized)
    if mol is None:
        raise ValueError(f"RDKit could not parse peptide sequence: {sequence!r}")
    return mol


def molecule_from_smiles(smiles: str):
    """Build an RDKit molecule from a SMILES string."""
    Chem, _, _, _ = _rdkit_modules()
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError("RDKit could not parse SMILES string")
    return mol


def descriptors_from_molecule(
    mol,
    *,
    name: str,
    source: str,
    sequence: Optional[str] = None,
    smiles: Optional[str] = None,
) -> MolecularDescriptors:
    """Compute the descriptor set used by permeability scoring and reports."""
    _, Crippen, Descriptors, rdMolDescriptors = _rdkit_modules()
    return MolecularDescriptors(
        name=name,
        source=source,
        formula=rdMolDescriptors.CalcMolFormula(mol),
        molecular_weight=Descriptors.MolWt(mol),
        exact_molecular_weight=Descriptors.ExactMolWt(mol),
        tpsa=Descriptors.TPSA(mol),
        logp=Crippen.MolLogP(mol),
        hbd=rdMolDescriptors.CalcNumHBD(mol),
        hba=rdMolDescriptors.CalcNumHBA(mol),
        rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol),
        formal_charge=sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        heavy_atom_count=mol.GetNumHeavyAtoms(),
        sequence=sequence,
        smiles=smiles,
    )


def descriptors_from_sequence(sequence: str, name: Optional[str] = None) -> MolecularDescriptors:
    normalized = sequence.strip().upper()
    mol = molecule_from_sequence(normalized)
    return descriptors_from_molecule(
        mol,
        name=name or normalized,
        source="rdkit_fasta",
        sequence=normalized,
    )


def descriptors_from_smiles(smiles: str, name: str, source: str = "smiles") -> MolecularDescriptors:
    normalized = smiles.strip()
    mol = molecule_from_smiles(normalized)
    return descriptors_from_molecule(
        mol,
        name=name,
        source=source,
        smiles=normalized,
    )


def read_smiles(path: str | Path) -> tuple[str, str]:
    """
    Read a simple .smi file as (smiles, name).

    The expected format is one non-comment line with SMILES followed by an
    optional molecule name separated by whitespace.
    """
    for line in Path(path).read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split(maxsplit=1)
        smiles = parts[0]
        name = parts[1] if len(parts) > 1 else Path(path).stem
        return smiles, name

    raise ValueError(f"No SMILES record found in {path}")
