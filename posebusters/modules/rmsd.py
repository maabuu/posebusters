"""Module to check RMSD between docked and crystal ligand."""
from __future__ import annotations

import logging
from copy import deepcopy

import numpy as np
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolAlign import CalcRMS, GetBestRMS
from rdkit.Chem.rdmolops import AddHs, RemoveHs
from rdkit.Chem import RemoveStereochemistry
from rdkit.rdBase import LogToPythonLogger

from ..tools.logging import CaptureLogger
from ..tools.molecules import remove_all_charges_and_hydrogens

LogToPythonLogger()
logger = logging.getLogger(__name__)


def check_rmsd(mol_pred: Mol, mol_true: Mol, rmsd_threshold: float = 2.0) -> dict[str, dict[str, bool | float]]:
    """Calculate RMSD and related metrics between predicted molecule and closest ground truth molecule.

    Args:
        mol_pred: Predicted molecule (docked ligand) with exactly one conformer.
        mol_true: Ground truth molecule (crystal ligand) with at least one conformer. If multiple conformers are
            present, the lowest RMSD will be reported.
        rmsd_threshold: Threshold in angstrom for reporting whether RMSD is within threshold. Defaults to 2.0.

    Returns:
        PoseBusters results dictionary.
    """
    assert isinstance(mol_true, Mol), "Ground truth molecule is missing."
    assert isinstance(mol_pred, Mol), "Predicted molecule is missing."
    num_conf = mol_true.GetNumConformers()
    assert num_conf > 0, "Crystal ligand needs at least one conformer."
    assert mol_pred.GetNumConformers() == 1, "Docked ligand should only have one conformer."

    rmsds = [robust_rmsd(mol_true, mol_pred, conf_id_probe=i) for i in range(num_conf)]
    kabsch_rmsd = [robust_rmsd(mol_true, mol_pred, conf_id_probe=i, kabsch=True) for i in range(num_conf)]

    i = np.argmin(np.nan_to_num(rmsds, nan=np.inf))
    rmsd = rmsds[i]
    kabsch_rmsd = kabsch_rmsd[i]

    rmsd_within_threshold = rmsd <= rmsd_threshold

    results = {"rmsd": rmsd, "kabsch_rmsd": kabsch_rmsd, "rmsd_within_threshold": rmsd_within_threshold}
    return {"results": results}


def robust_rmsd(
    mol_probe: Mol,
    mol_ref: Mol,
    conf_id_probe: int = -1,
    conf_id_ref: int = -1,
    drop_stereo: bool = True,
    heavy_only: bool = True,
    kabsch: bool = False,
    symmetrizeConjugatedTerminalGroups=True,
    **params,
) -> float:
    """RMSD calculation that handles errors and returns NaN if RMSD cannot be calculated."""
    mol_probe = deepcopy(mol_probe)
    mol_ref = deepcopy(mol_ref)

    # if drop_stereo:
    #     RemoveStereochemistry(mol_probe)
    #     RemoveStereochemistry(mol_ref)
    # if heavy_only:
    #     mol_probe = remove_all_charges_and_hydrogens(mol_probe)
    #     mol_ref = remove_all_charges_and_hydrogens(mol_ref)

    if heavy_only:
        mol_probe = RemoveHs(mol_probe)
        mol_ref = RemoveHs(mol_ref)

    params = dict(symmetrizeConjugatedTerminalGroups=symmetrizeConjugatedTerminalGroups, **params)

    try:
        return _call_rdkit_rmsd(mol_probe, mol_ref, conf_id_probe, conf_id_ref, params, kabsch)
    except RuntimeError as error:
        pass
    except ValueError as error:
        pass

    # try again but on canonical tautomers
    mol_ref = rdMolStandardize.TautomerEnumerator().Canonicalize(mol_ref)
    mol_probe = rdMolStandardize.TautomerEnumerator().Canonicalize(mol_probe)

    try:
        return _call_rdkit_rmsd(mol_probe, mol_ref, conf_id_probe, conf_id_ref, params, kabsch)
    except RuntimeError as error:
        logger.info(f"Could not calculate RMSD due to {error}")
        return np.nan
    except ValueError as error:
        logger.info(f"Could not calculate RMSD due to {error}")
        return np.nan


def _call_rdkit_rmsd(mol_probe: Mol, mol_ref: Mol, conf_id_probe: int, conf_id_ref: int, params, kabsch: bool = False):
    if kabsch is True:
        with CaptureLogger():
            return GetBestRMS(prbMol=mol_probe, refMol=mol_ref, prbId=conf_id_probe, refId=conf_id_ref, **params)
    else:
        with CaptureLogger():
            return CalcRMS(prbMol=mol_probe, refMol=mol_ref, prbId=conf_id_probe, refId=conf_id_ref, **params)
