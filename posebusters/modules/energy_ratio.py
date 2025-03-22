"""Module to check energy of ligand conformations."""

from __future__ import annotations

import logging
from copy import deepcopy
from functools import cache
from math import isfinite

from rdkit import ForceField  # noqa: F401
from rdkit.Chem.inchi import InchiReadWriteError, MolFromInchi
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs, ETKDGv3
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField, UFFOptimizeMoleculeConfs
from rdkit.Chem.rdmolops import AddHs, AssignStereochemistryFrom3D, SanitizeMol

from ..tools.inchi import get_inchi
from ..tools.logging import CaptureLogger

logger = logging.getLogger(__name__)


_warning_prefix = "WARNING: Energy ratio module "
_empty_results = {
    "results": {
        "ensemble_avg_energy": float("nan"),
        "mol_pred_energy": float("nan"),
        "energy_ratio": float("nan"),
        "energy_ratio_passes": float("nan"),
    }
}


def check_energy_ratio(
    mol_pred: Mol,
    threshold_energy_ratio: float = 7.0,
    ensemble_number_conformations: int = 100,
    inchi_strict: bool = False,
    epsilon=1e-10,
) -> dict[str, dict[str, float | bool]]:
    """Check whether the energy of the docked ligand is within user defined range.

    Args:
        mol_pred: Predicted molecule (docked ligand) with exactly one conformer.
        threshold_energy_ratio: Limit above which the energy ratio is deemed to high. Defaults to 7.0.
        ensemble_number_conformations: Number of conformations to generate for the ensemble over which to
            average. Defaults to 100.
        inchi_strict: Whether to treat warnings in the InChI generation as errors. Defaults to False.

    Returns:
        PoseBusters results dictionary.
    """
    mol_pred = deepcopy(mol_pred)

    try:
        assert mol_pred.GetNumConformers() > 0, "Molecule does not have a conformer."
        assert not SanitizeMol(mol_pred, catchErrors=True), "Molecule does not sanitize."
        AddHs(mol_pred, addCoords=True)
    except Exception as e:
        logger.warning(_warning_prefix + "failed because RDKit sanitization failed for molecule: %s", e)
        return _empty_results

    try:
        inchi = get_inchi(mol_pred, inchi_strict=inchi_strict)
    except InchiReadWriteError as e:
        logger.warning(_warning_prefix + "failed because InChI creation failed for molecule: %s", e.args[1])
        return _empty_results
    except Exception as e:
        logger.warning(_warning_prefix + "failed because InChI creation failed for molecule: %s", e)
        return _empty_results

    try:
        observed_energy = get_conf_energy(mol_pred)
    except Exception as e:
        logger.warning(_warning_prefix + "failed to calculate conformation energy for %s: %s", inchi, e)
        observed_energy = float("nan")

    try:
        energies = get_energies(inchi, ensemble_number_conformations)
        mean_energy = sum(energies) / len(energies)
        std_energy = sum((energy - mean_energy) ** 2 for energy in energies) / (len(energies) - 1)
        std_energy = max(epsilon, std_energy)  # clipping
    except Exception as e:
        logger.warning(_warning_prefix + "failed to calculate ensemble conformation energy for %s: %s", inchi, e)
        mean_energy = float("nan")
        std_energy = float("nan")

    if mean_energy == 0:
        logger.warning(_warning_prefix + "calculated average energy of molecule 0 for %s", inchi)
        mean_energy = epsilon  # clipping

    # simple ratio
    ratio = observed_energy / mean_energy
    ratio_passes = ratio <= threshold_energy_ratio if isfinite(ratio) else float("nan")

    # ratio after subtracting mean
    deviation = observed_energy - mean_energy
    relative_deviation = deviation / mean_energy
    relative_deviation_passes = (
        relative_deviation <= threshold_energy_ratio if isfinite(relative_deviation) else float("nan")
    )

    # standard score (ratio after subtracting by population mean and dividing by population std)
    z_value = (observed_energy - mean_energy) / std_energy
    z_value_passes = z_value <= threshold_energy_ratio if isfinite(z_value) else float("nan")

    results = {
        "ensemble_avg_energy": mean_energy,
        "mol_pred_energy": observed_energy,
        "energy_ratio": ratio,
        "relative_deviation": relative_deviation,
        "z_value": z_value,
        "energy_ratio_passes": ratio_passes,
        "relative_deviation_passes": relative_deviation_passes,
        "z_value_passes": z_value_passes,
    }
    return {"results": results}


def get_average_energy(inchi: str, n_confs: int = 50, num_threads: int = 0) -> Mol:
    """Get average energy of an ensemble of molecule conformations."""
    energies = get_energies(inchi, n_confs, num_threads)
    return sum(energies) / len(energies)


@cache
def get_energies(inchi: str, n_confs: int = 50, num_threads: int = 0) -> list[float]:
    """Get energies of an ensemble of molecule conformations."""
    with CaptureLogger():
        mol = MolFromInchi(inchi)
    return new_conformation(mol, n_confs, num_threads)["energies"]


def new_conformation(mol: Mol, n_confs: int = 1, num_threads: int = 0, energy_minimization=True) -> Mol:
    """Generate new conformation(s) for a molecule."""
    assert mol is not None

    etkdg = ETKDGv3()
    etkdg.randomSeed = 42
    etkdg.verbose = False
    etkdg.useRandomCoords = True
    etkdg.numThreads = num_threads

    # prep mol
    mol_etkdg = deepcopy(mol)
    mol_etkdg = AddHs(mol_etkdg, addCoords=True)
    AssignStereochemistryFrom3D(mol_etkdg, replaceExistingTags=True)
    mol_etkdg.RemoveAllConformers()

    # etkdg
    with CaptureLogger():
        cids_etkdg = EmbedMultipleConfs(mol_etkdg, n_confs, etkdg)
        assert len(cids_etkdg) == n_confs, "Failed to generate conformations."

    # energy minimization
    if energy_minimization:
        with CaptureLogger():
            energy_uff = UFFOptimizeMoleculeConfs(mol_etkdg, numThreads=num_threads)
        energy_uff = [v[1] for v in energy_uff]
    else:
        energy_uff = [get_conf_energy(mol_etkdg, conf_id) for conf_id in cids_etkdg]

    return {"mol": mol_etkdg, "energies": energy_uff}


def get_conf_energy(mol: Mol, conf_id: int = -1) -> float:
    """Get energy of a conformation."""
    assert mol is not None

    mol = AddHs(mol, addCoords=True)
    with CaptureLogger():
        uff = UFFGetMoleculeForceField(mol, confId=conf_id)
        e_uff = uff.CalcEnergy()
    return e_uff
