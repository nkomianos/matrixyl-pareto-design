"""
Exploratory conformer validation for peptide candidates.

These metrics are intended for candidate triage, not publication-grade MD.
"""

from dataclasses import asdict, dataclass
from statistics import mean, pstdev
from typing import Optional

from .chemistry import molecule_from_sequence, molecule_from_smiles


@dataclass(frozen=True)
class ConformerMetrics:
    conformer_id: int
    radius_of_gyration: float
    energy: Optional[float]
    force_field: str
    optimization_converged: bool

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class ConformerEnsembleSummary:
    name: str
    source: str
    conformer_count: int
    min_radius_of_gyration: float
    mean_radius_of_gyration: float
    max_radius_of_gyration: float
    std_radius_of_gyration: float
    min_energy: Optional[float]
    mean_energy: Optional[float]
    metrics: tuple[ConformerMetrics, ...]

    def to_dict(self) -> dict:
        data = asdict(self)
        data["metrics"] = [metric.to_dict() for metric in self.metrics]
        return data


def _rdkit_conformer_modules():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
    except ImportError as exc:
        raise ImportError("RDKit is required for conformer validation") from exc
    return Chem, AllChem, rdMolDescriptors


def _optimize_conformer(mol, conformer_id: int, max_iterations: int) -> tuple[Optional[float], str, bool]:
    _, AllChem, _ = _rdkit_conformer_modules()
    try:
        properties = AllChem.MMFFGetMoleculeProperties(mol)
        if properties is not None:
            force_field = AllChem.MMFFGetMoleculeForceField(
                mol,
                properties,
                confId=conformer_id,
            )
            if force_field is not None:
                status = force_field.Minimize(maxIts=max_iterations)
                return force_field.CalcEnergy(), "MMFF", status == 0
    except Exception:
        pass

    try:
        force_field = AllChem.UFFGetMoleculeForceField(mol, confId=conformer_id)
        if force_field is not None:
            status = force_field.Minimize(maxIts=max_iterations)
            return force_field.CalcEnergy(), "UFF", status == 0
    except Exception:
        pass

    return None, "none", False


def summarize_conformer_ensemble(
    mol,
    *,
    name: str,
    source: str,
    num_conformers: int = 20,
    seed: int = 42,
    max_optimization_iterations: int = 200,
) -> ConformerEnsembleSummary:
    Chem, AllChem, rdMolDescriptors = _rdkit_conformer_modules()
    if num_conformers < 1:
        raise ValueError("num_conformers must be at least 1")

    working_mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.25
    conformer_ids = list(
        AllChem.EmbedMultipleConfs(
            working_mol,
            numConfs=num_conformers,
            params=params,
        )
    )
    if not conformer_ids:
        raise ValueError(f"RDKit could not generate conformers for {name}")

    metrics = []
    for conformer_id in conformer_ids:
        energy, force_field, converged = _optimize_conformer(
            working_mol,
            conformer_id,
            max_optimization_iterations,
        )
        rg = rdMolDescriptors.CalcRadiusOfGyration(working_mol, confId=conformer_id)
        metrics.append(
            ConformerMetrics(
                conformer_id=conformer_id,
                radius_of_gyration=rg,
                energy=energy,
                force_field=force_field,
                optimization_converged=converged,
            )
        )

    radii = [metric.radius_of_gyration for metric in metrics]
    energies = [metric.energy for metric in metrics if metric.energy is not None]
    return ConformerEnsembleSummary(
        name=name,
        source=source,
        conformer_count=len(metrics),
        min_radius_of_gyration=min(radii),
        mean_radius_of_gyration=mean(radii),
        max_radius_of_gyration=max(radii),
        std_radius_of_gyration=pstdev(radii) if len(radii) > 1 else 0.0,
        min_energy=min(energies) if energies else None,
        mean_energy=mean(energies) if energies else None,
        metrics=tuple(metrics),
    )


def summarize_sequence_conformers(
    sequence: str,
    *,
    name: Optional[str] = None,
    num_conformers: int = 20,
    seed: int = 42,
) -> ConformerEnsembleSummary:
    normalized = sequence.strip().upper()
    return summarize_conformer_ensemble(
        molecule_from_sequence(normalized),
        name=name or normalized,
        source="rdkit_fasta",
        num_conformers=num_conformers,
        seed=seed,
    )


def summarize_smiles_conformers(
    smiles: str,
    *,
    name: str,
    source: str = "smiles",
    num_conformers: int = 20,
    seed: int = 42,
) -> ConformerEnsembleSummary:
    return summarize_conformer_ensemble(
        molecule_from_smiles(smiles),
        name=name,
        source=source,
        num_conformers=num_conformers,
        seed=seed,
    )
