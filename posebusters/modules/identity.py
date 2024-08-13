"""Module to check identity of docked and crystal ligand."""

from __future__ import annotations

import logging
from typing import Any

from rdkit.Chem.rdchem import Mol
from rdkit.rdBase import LogToPythonLogger

from ..tools.inchi import is_valid_inchi, split_inchi, standardize_and_get_inchi

LogToPythonLogger()
logger = logging.getLogger(__name__)


def check_identity(mol_pred: Mol, mol_true: Mol, inchi_options: str = "") -> dict[str, Any]:
    """Check two molecules are identical (docking relevant identity).

    Args:
        mol_pred: Predicted molecule (docked ligand).
        mol_true:Ground truth molecule (crystal ligand) with a conformer.
        inchi_options: String of options to pass to the InChI module. Defaults to "".

    Returns:
        PoseBusters results dictionary.
    """
    # generate inchis
    inchi_crystal = standardize_and_get_inchi(mol_true, options=inchi_options)
    inchi_docked = standardize_and_get_inchi(mol_pred, options=inchi_options)

    # check inchis are valid
    inchi_crystal_valid = is_valid_inchi(inchi_crystal)
    inchi_docked_valid = is_valid_inchi(inchi_docked)

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
    "/b": "stereo_dbond",  # double bond (Z/E) stereochemistry
    "/t": "stereo_sp3",  # tetrahderal stereochemistry
    "/m": "stereo_sp3_inverted",
    "/s": "stereo_type",
    "/i": "isotopic",
}
stereo_all_layers = ["stereo_dbond", "stereo_sp3", "stereo_sp3_inverted", "stereo_type"]
stereo_tetrahedral_layers = ["stereo_sp3", "stereo_sp3_inverted", "stereo_type"]


def _compare_inchis(inchi_true: str, inchi_pred: str, layers: list[str] = standard_layers) -> dict[str, bool]:
    results = {}

    # fast return when overall InChI is the same
    if inchi_true == inchi_pred:
        results["inchi_overall"] = True
        for layer in layers:
            results[layer_names[layer]] = True
            results["stereo_tetrahedral"] = True
            results["stereo"] = True
        return results

    # otherwise comparison by layer
    results["inchi_overall"] = False
    layers_true = split_inchi(inchi_true)
    layers_pred = split_inchi(inchi_pred)
    assert "/" in layers_true, "Molecular formula layer missing from InChI string"
    for layer in layers:
        name = layer_names[layer]
        if (layer not in layers_true) or (layer not in layers_pred):
            results[name] = (layer not in layers_true) and (layer not in layers_pred)
        else:
            results[name] = layers_true[layer] == layers_pred[layer]

    # combine stereo layers (pass if not present)
    results["stereo_tetrahedral"] = all(results.get(name, True) for name in stereo_tetrahedral_layers)
    results["stereo"] = all(results.get(name, True) for name in stereo_all_layers)

    return results
