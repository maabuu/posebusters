"""Module to check chemistry of molecules."""
import pandas as pd
from rdkit import RDLogger
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import DetectChemistryProblems, GetMolFrags, SanitizeFlags


def check_chemistry(mol_pred: Mol) -> Mol:
    """Check chemical sanity of molecule.

    Args:
        mol_pred: Molecule.

    Returns:
        PoseBusters results dictionary.
    """
    assert isinstance(mol_pred, Mol)

    RDLogger.DisableLog("rdApp.*")
    errors = DetectChemistryProblems(mol_pred, sanitizeOps=SanitizeFlags.SANITIZE_ALL)
    RDLogger.EnableLog("rdApp.*")

    messages = [error.Message() for error in errors]
    atom_indices = [e.GetAtomIndices() if hasattr(e, "GetAtomIndices") else e.GetAtomIdx() for e in errors]
    types = [error.GetType() for error in errors]

    num_frags = len(GetMolFrags(mol_pred, asMols=False, sanitizeFrags=False))

    results = {
        "passes_valence_checks": "AtomValenceException" not in types,
        "passes_kekulization": "AtomKekulizeException" not in types,
        "passes_rdkit_sanity_checks": len(errors) == 0,
        "all_atoms_connected": num_frags <= 1,
    }
    details = pd.DataFrame(dict(messages=messages, atom_indices=atom_indices, types=types))

    return {"results": results, "details": details}
