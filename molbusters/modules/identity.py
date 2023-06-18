"""Module to check identity of docked and crystal ligand."""
from __future__ import annotations

import logging
from copy import deepcopy
from re import findall
from typing import Any

from rdkit.Chem.inchi import MolFromInchi, MolToInchi
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
from rdkit.Chem.rdchem import AtomValenceException, Mol
from rdkit.Chem.rdmolops import (
    AssignStereochemistryFrom3D,
    RemoveHs,
    RemoveStereochemistry,
)
from rdkit.rdBase import LogToPythonLogger

from ..tools.logging import CaptureLogger
from ..tools.molecules import (
    add_stereo_hydrogens,
    assert_sanity,
    neutralize_atoms,
    remove_isotopic_info,
)

LogToPythonLogger()
logger = logging.getLogger(__name__)


def check_identity(mol_pred: Mol, mol_true: Mol, inchi_options: str = "") -> dict[str, Any]:
    """Check two molecules are identical (docking relevant identity).

    Args:
        mol_pred: Predicted molecule (docked ligand).
        mol_true:Ground truth molecule (crystal ligand) with a conformer.
        inchi_options: String of options to pass to the InChI module. Defaults to "".

    Returns:
        MolBusters results dictionary.
    """
    # generate inchis
    inchi_crystal = standardize_and_get_inchi(mol_true, options=inchi_options)
    inchi_docked = standardize_and_get_inchi(mol_pred, options=inchi_options)

    # check inchis are valid
    inchi_crystal_valid = _is_valid_inchi(inchi_crystal)
    inchi_docked_valid = _is_valid_inchi(inchi_docked)

    # compare inchis
    if inchi_crystal_valid and inchi_docked_valid:
        inchi_comparison = _compare_inchis(inchi_docked, inchi_crystal)
    else:
        inchi_comparison = {}

    results = {
        "inchi_crystal_valid": inchi_crystal_valid,
        "inchi_docked_valid": inchi_docked_valid,
        "inchi_crystal": inchi_crystal,
        "inchi_docked": inchi_docked,
        **inchi_comparison,
    }

    return {"results": results}


# isotopic layer /i (for isotopes) not present if no isotopic information in molecule
# fixed hydrogen layer /f (for tautomers) not part of standard inchi string
# reconnected layer /r (for metal ions) not part of standard inchi string
standard_layers = ["=", "/", "/c", "/h", "/q", "/p", "/t", "/b", "/m", "/s"]
layer_names = {
    "=": "inchi_version",
    "/": "formula",
    "/c": "connections",
    "/h": "hydrogens",
    "/q": "net_charge",
    "/p": "protons",
    "/t": "tetrahedral",
    "/b": "double_bonds",
    "/m": "stereo_relative_inchi_enatiomer",
    "/s": "stereo_absolute",
    "/i": "isotopic",
}
stereo_layers = ["tetrahedral", "double_bonds", "stereo_relative_inchi_enatiomer", "stereo_absolute"]


def standardize_and_get_inchi(mol: Mol, options: str = "", log_level=None, warnings_as_errors=False) -> str:
    """Return InChI after standardising molecule and inferring stereo from coordinates."""
    mol = deepcopy(mol)
    mol = assert_sanity(mol)

    # standardize molecule and assign stereo from 3D coordinates
    mol = remove_isotopic_info(mol)
    RemoveStereochemistry(mol)
    mol = RemoveHs(mol)
    try:
        mol = neutralize_atoms(mol)
    except AtomValenceException:
        logger.warning("Failed to neutralize molecule. Using uncharger. InChI check might fail.")
        mol = Uncharger().uncharge(mol)
    mol = add_stereo_hydrogens(mol)
    AssignStereochemistryFrom3D(mol, replaceExistingTags=True)

    with CaptureLogger():
        inchi = MolToInchi(mol, options=options, logLevel=log_level, treatWarningAsError=warnings_as_errors)

    return inchi


def _compare_inchis(inchi_true: str, inchi_pred: str, layers: list[str] = standard_layers) -> dict[str, bool]:
    results = {}

    # fast when overall InChI is the same
    if inchi_true == inchi_pred:
        results["inchi_overall"] = True
        for layer in layers:
            results[layer_names[layer]] = True
            results["stereo"] = True
        return results

    # otherwise comparison by layer
    results["inchi_overall"] = False
    layers_true = _split_inchi(inchi_true)
    layers_pred = _split_inchi(inchi_pred)
    assert "/" in layers_true, "Molecular formula layer missing from InChI string"
    for layer in layers:
        name = layer_names[layer]
        if (layer not in layers_true) or (layer not in layers_pred):
            results[name] = (layer not in layers_true) and (layer not in layers_pred)
        else:
            results[name] = layers_true[layer] == layers_pred[layer]
    results["stereo"] = all(results[name] for name in stereo_layers if name in results)  # combine stereo layers

    return results


def _split_inchi(inchi: str) -> dict[str, str]:
    """Split the standard InChI string without isotopic info into its layers."""
    if not inchi.startswith("InChI="):
        raise ValueError("InChI string must start with 'InChI='")

    # inchi always InChi=1S/...formula.../...layer.../...layer.../...layer...
    version = inchi[6:].split(r"/", 1)[0]
    layers = findall(r"(?=.*)(\/[a-z]{0,1})(.*?)(?=\/.*|$)", inchi[6:])

    inchi_parts = {"=": version}
    for prefix, layer in layers:
        # standard inchi strings without isotopic info have each layer no more than once
        assert prefix not in inchi_parts, f"Layer {prefix} more than once!"
        inchi_parts[prefix] = layer

    return inchi_parts


def _is_valid_inchi(inchi: str) -> bool:
    """Check that InChI can be parsed and sanitization does not fail."""
    try:
        mol = MolFromInchi(inchi)
        assert_sanity(mol)
        assert mol is not None
        return True
    except Exception:
        return False
