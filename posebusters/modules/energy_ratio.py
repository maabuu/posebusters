"""Module to check energy of ligand conformations."""
from __future__ import annotations

import logging
from copy import deepcopy
from functools import lru_cache

import numpy as np
from rdkit import ForceField  # noqa: F401
from rdkit.Chem.inchi import InchiReadWriteError, MolFromInchi, MolToInchi
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs, ETKDGv3
from rdkit.Chem.rdForceFieldHelpers import (
    UFFGetMoleculeForceField,
    UFFOptimizeMoleculeConfs,
)
from rdkit.Chem.rdmolops import AddHs, AssignStereochemistryFrom3D

from ..tools.logging import CaptureLogger
from ..tools.molecules import assert_sanity

logger = logging.getLogger(__name__)

_empty_results = {
    "results": {
        "ensemble_avg_energy": np.nan,
        "mol_pred_energy": np.nan,
        "energy_ratio": np.nan,
        "energy_ratio_passes": np.nan,
    }
}


def check_energy_ratio(
    mol_pred: Mol,
    threshold_energy_ratio: float = 7.0,
    ensemble_number_conformations: int = 100,
    inchi_strict: bool = False,
):
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
        mol_pred = assert_sanity(mol_pred)
        AddHs(mol_pred, addCoords=True)
    except Exception as e:
        logger.warning(f"Failed to prepare molecule: {e}")
        return _empty_results

    try:
        inchi = get_inchi(mol_pred, inchi_strict=inchi_strict)
    except InchiReadWriteError as e:
        logger.warning(f"Molecule does not sanitize: {e.args[1]}")
        return _empty_results
    except Exception as e:
        logger.warning(f"Molecule does not sanitize: {e}")
        return _empty_results

    try:
        conf_energy = get_conf_energy(mol_pred)
    except Exception as e:
        logger.warning(f"Failed to calculate conformation energy for {inchi}: {e}")
        conf_energy = np.nan

    try:
        avg_energy = float(get_average_energy(inchi, ensemble_number_conformations))
    except Exception as e:
        logger.warning(f"Failed to calculate ensemble conformation energy for {inchi}: {e}")
        avg_energy = np.nan

    pred_factor = conf_energy / avg_energy
    ratio_passes = pred_factor <= threshold_energy_ratio

    results = {
        "ensemble_avg_energy": avg_energy,
        "mol_pred_energy": conf_energy,
        "energy_ratio": pred_factor,
        "energy_ratio_passes": ratio_passes,
    }
    return {"results": results}


def get_inchi(mol: Mol, inchi_strict: bool = False) -> str:
    """Prepare and check InChI of a molecule."""
    with CaptureLogger() as log:
        inchi = MolToInchi(mol, treatWarningAsError=inchi_strict)
        # check inchi because inchi generation does not raise an error if the inchi is invalid
        if MolFromInchi(inchi, sanitize=True) is None:
            raise Exception(log["ERROR"].split("ERROR: ")[-1])
    return inchi


@lru_cache(maxsize=None)
def get_average_energy(inchi: str, n_confs: int = 50, num_threads: int = 0) -> Mol:
    """Get average energy of an ensemble of molecule conformations."""
    with CaptureLogger():
        mol = MolFromInchi(inchi)
    energies = new_conformation(mol, n_confs, num_threads)["energies"]
    return sum(energies) / len(energies)


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
