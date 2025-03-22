from __future__ import annotations

import logging
import re
from copy import deepcopy
from re import findall

from rdkit.Chem.inchi import MolFromInchi, MolToInchi
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
from rdkit.Chem.rdchem import AtomValenceException, Mol
from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D, RemoveHs, RemoveStereochemistry, SanitizeMol

from .logging import CaptureLogger
from .molecules import add_stereo_hydrogens, neutralize_atoms, remove_isotopic_info

logger = logging.getLogger(__name__)


def get_inchi(mol: Mol, inchi_strict: bool = False) -> str:
    """Get inchi of a molecule."""

    with CaptureLogger() as log:
        SanitizeMol(mol)
        inchi = MolToInchi(mol, treatWarningAsError=inchi_strict)
        # check inchi because inchi generation does not raise an error if the inchi is invalid
        if MolFromInchi(inchi, sanitize=True) is None:
            warnings = re.split(r"\[.+?\] WARNING: ", log.get("WARNING", ""))[1:]
            errors = re.split(r"\[.+?\] ERROR: ", log.get("ERROR", ""))[1:]
            message = f"InChI ERRORs: {'; '.join(errors)} WARNINGs: {'; '.join(warnings)}"
            raise Exception(message)
    return inchi


def standardize_and_get_inchi(mol: Mol, options: str = "", log_level=None, warnings_as_errors=False) -> str:
    """Return InChI after standardising molecule and inferring stereo from coordinates."""

    mol = deepcopy(mol)

    if flags := SanitizeMol(mol, catchErrors=True):
        logger.debug("Cannot get InChI because molecule doesn't sanitize due to %s.", flags)
        return ""

    # standardise molecule
    mol = remove_isotopic_info(mol)

    # assign stereo from 3D coordinates only if 3D coordinates are present
    has_pose = mol.GetNumConformers() > 0
    if has_pose:
        RemoveStereochemistry(mol)

    mol = RemoveHs(mol)
    try:
        mol = neutralize_atoms(mol)
    except AtomValenceException:
        logger.warning("Failed to neutralize molecule. Using uncharger. InChI check might fail.")
        mol = Uncharger().uncharge(mol)
    mol = add_stereo_hydrogens(mol)

    if has_pose:
        AssignStereochemistryFrom3D(mol, replaceExistingTags=True)

    with CaptureLogger():
        inchi = MolToInchi(mol, options=options, logLevel=log_level, treatWarningAsError=warnings_as_errors)

    return inchi


def is_valid_inchi(inchi: str) -> bool:
    """Check that InChI can be parsed and sanitization does not fail."""
    try:
        mol = MolFromInchi(inchi)
        assert mol is not None, "Molecule is None."
        assert not SanitizeMol(mol, catchErrors=True), "Molecule does not sanitize."
        return True
    except Exception:
        return False


def split_inchi(inchi: str) -> dict[str, str]:
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
