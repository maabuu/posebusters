"""Module to check RMSD between docked and crystal ligand."""

from __future__ import annotations

import logging
from copy import deepcopy

import numpy as np
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolAlign import CalcRMS, GetBestRMS
from rdkit.Chem.rdmolops import RemoveHs, RemoveStereochemistry
from rdkit.rdBase import LogToPythonLogger

from ..tools.logging import CaptureLogger
from ..tools.molecules import neutralize_atoms, remove_all_charges_and_hydrogens

LogToPythonLogger()
logger = logging.getLogger(__name__)

tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
tautomer_enumerator.SetMaxTautomers(100000)
tautomer_enumerator.SetMaxTransforms(100000)
tautomer_enumerator.SetReassignStereo(True)
tautomer_enumerator.SetRemoveBondStereo(True)
tautomer_enumerator.SetRemoveSp3Stereo(True)


def check_rmsd(
    mol_pred: Mol, mol_true: Mol, rmsd_threshold: float = 2.0, heavy_only: bool = True, choose_by: str = "rmsd"
) -> dict[str, dict[str, bool | float]]:
    """Calculate RMSD and related metrics between predicted molecule and closest ground truth molecule.

    Args:
        mol_pred: Predicted molecule (docked ligand) with exactly one conformer.
        mol_true: Ground truth molecule (crystal ligand) with at least one conformer. If multiple conformers are
            present, the lowest RMSD will be reported.
        rmsd_threshold: Threshold in angstrom for reporting whether RMSD is within threshold. Defaults to 2.0.
        heavy_only: Whether to only consider heavy atoms for RMSD calculation. Defaults to True.
        choose_by: Metric to choose which mol_true conformation to compare to. Defaults to "rmsd".

    Returns:
        PoseBusters results dictionary.
    """
    assert isinstance(mol_true, Mol), "Ground truth molecule is missing."
    assert isinstance(mol_pred, Mol), "Predicted molecule is missing."
    num_conf = mol_true.GetNumConformers()
    assert num_conf > 0, "Crystal ligand needs at least one conformer."
    assert mol_pred.GetNumConformers() == 1, "Docked ligand should only have one conformer."

    rmsds = [robust_rmsd(mol_true, mol_pred, conf_id_probe=i, heavy_only=heavy_only) for i in range(num_conf)]
    kabsch_rmsds = [
        robust_rmsd(mol_true, mol_pred, conf_id_probe=i, kabsch=True, heavy_only=heavy_only) for i in range(num_conf)
    ]
    intercentroids = [
        intercentroid(mol_true, mol_pred, conf_id_probe=i, heavy_only=heavy_only) for i in range(num_conf)
    ]

    if choose_by == "rmsd":
        i = np.argmin(np.nan_to_num(rmsds, nan=np.inf))
    elif choose_by == "kabsch_rmsd":
        i = np.argmin(np.nan_to_num(kabsch_rmsds, nan=np.inf))
    elif choose_by == "centroid_distance":
        i = np.argmin(np.nan_to_num(intercentroids, nan=np.inf))
    else:
        raise ValueError(f"Invalid value {choose_by} for choose_by. Use 'rmsd', 'kabsch_rmsd', 'centroid_distance'.")
    rmsd = rmsds[i]
    kabsch_rmsd = kabsch_rmsds[i]
    centroid_dist = intercentroids[i]

    rmsd_within_threshold = rmsd <= rmsd_threshold

    results = {
        "rmsd": rmsd,
        "kabsch_rmsd": kabsch_rmsd,
        "centroid_distance": centroid_dist,
        "rmsd_within_threshold": rmsd_within_threshold,
    }
    return {"results": results}


def robust_rmsd(  # noqa: PLR0913
    mol_probe: Mol,
    mol_ref: Mol,
    conf_id_probe: int = -1,
    conf_id_ref: int = -1,
    drop_stereo: bool = False,
    heavy_only: bool = True,
    kabsch: bool = False,
    symmetrizeConjugatedTerminalGroups=True,
    **params,
) -> float:
    """RMSD calculation for isomers."""
    mol_probe = deepcopy(mol_probe)
    mol_ref = deepcopy(mol_ref)  # copy mols because rdkit RMSD calculation aligns mols

    if drop_stereo:
        RemoveStereochemistry(mol_probe)
        RemoveStereochemistry(mol_ref)

    if heavy_only:
        mol_probe = RemoveHs(mol_probe, sanitize=False)
        mol_ref = RemoveHs(mol_ref, sanitize=False)

    # combine parameters
    params = dict(symmetrizeConjugatedTerminalGroups=symmetrizeConjugatedTerminalGroups, kabsch=kabsch, **params)

    # calculate RMSD
    rmsd = _call_rdkit_rmsd(mol_probe, mol_ref, conf_id_probe, conf_id_ref, **params)
    if not np.isnan(rmsd):
        return rmsd

    # try again ignoring charges or tautomers
    rmsd = _rmsd_ignoring_charges_and_tautomers(mol_ref, mol_probe, conf_id_ref, conf_id_probe, params)
    if not np.isnan(rmsd):
        return rmsd

    # try assigning the bond orders of one molecule to the other
    try:
        mol_probe_new_bonds = AssignBondOrdersFromTemplate(mol_ref, mol_probe)
    except ValueError:
        return np.nan
    rmsd = _call_rdkit_rmsd(mol_probe_new_bonds, mol_ref, conf_id_probe, conf_id_ref, **params)
    if not np.isnan(rmsd):
        return rmsd

    return np.nan


def _rmsd_ignoring_charges_and_tautomers(
    mol_ref: Mol, mol_probe: Mol, conf_id_ref: int, conf_id_probe: int, params: dict
) -> float:
    # try again but remove charges and hydrogens
    mol_ref_uncharged = remove_all_charges_and_hydrogens(mol_ref)
    mol_probe_uncharged = remove_all_charges_and_hydrogens(mol_probe)
    rmsd = _call_rdkit_rmsd(mol_probe_uncharged, mol_ref_uncharged, conf_id_probe, conf_id_ref, **params)
    if not np.isnan(rmsd):
        return rmsd

    # try again but neutralize atoms
    mol_ref_neutralized = neutralize_atoms(mol_ref)
    mol_probe_neutralized = neutralize_atoms(mol_probe)
    rmsd = _call_rdkit_rmsd(mol_probe_neutralized, mol_ref_neutralized, conf_id_probe, conf_id_ref, **params)
    if not np.isnan(rmsd):
        return rmsd

    # try again but on canonical tautomers
    mol_ref_canonical = tautomer_enumerator.Canonicalize(mol_ref)
    mol_probe_canonical = tautomer_enumerator.Canonicalize(mol_probe)
    rmsd = _call_rdkit_rmsd(mol_probe_canonical, mol_ref_canonical, conf_id_probe, conf_id_ref, **params)
    if not np.isnan(rmsd):
        return rmsd

    # try again but after neutralizing atoms
    mol_ref_neutral_canonical = tautomer_enumerator.Canonicalize(neutralize_atoms(mol_ref))
    mol_probe_neutral_canonical = tautomer_enumerator.Canonicalize(neutralize_atoms(mol_probe))
    rmsd = _call_rdkit_rmsd(
        mol_probe_neutral_canonical, mol_ref_neutral_canonical, conf_id_probe, conf_id_ref, **params
    )
    if not np.isnan(rmsd):
        return rmsd

    return rmsd


def _call_rdkit_rmsd(mol_probe: Mol, mol_ref: Mol, conf_id_probe: int, conf_id_ref: int, **params):
    try:
        with CaptureLogger():
            return _rmsd(mol_probe, mol_ref, conf_id_probe, conf_id_ref, **params)
    except RuntimeError:
        pass
    except ValueError:
        pass

    return np.nan


def _rmsd(mol_probe: Mol, mol_ref: Mol, conf_id_probe: int, conf_id_ref: int, kabsch: bool = False, **params):
    if kabsch is True:
        return GetBestRMS(prbMol=mol_probe, refMol=mol_ref, prbId=conf_id_probe, refId=conf_id_ref, **params)
    return CalcRMS(prbMol=mol_probe, refMol=mol_ref, prbId=conf_id_probe, refId=conf_id_ref, **params)


def intercentroid(
    mol_probe: Mol, mol_ref: Mol, conf_id_probe: int = -1, conf_id_ref: int = -1, heavy_only: bool = True
) -> float:
    """Distance between centroids of two molecules."""

    centroid_probe = get_centroid(mol_probe, heavy_only, conf_id_probe)
    centroid_ref = get_centroid(mol_ref, heavy_only, conf_id_ref)
    return float(np.linalg.norm(centroid_probe - centroid_ref))


def get_centroid(mol: Mol, heavy_only: bool = True, conf_id: int = -1) -> np.ndarray:
    """Get centroid of molecule."""
    pos = mol.GetConformer(conf_id).GetPositions()
    if heavy_only:
        pos = pos[[atom.GetAtomicNum() != 1 for atom in mol.GetAtoms()], :]
    return pos.mean(axis=0)
