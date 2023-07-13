from __future__ import annotations

from rdkit.Chem.rdmolfiles import MolFromSmiles

from posebusters.modules.sanity import check_chemistry


def test_check_rmsd():
    mol_pass = MolFromSmiles("C=C", sanitize=True)
    assert check_chemistry(mol_pass)["results"]["passes_rdkit_sanity_checks"] is True

    mol_pass = MolFromSmiles("c1ccccc1", sanitize=True)
    assert check_chemistry(mol_pass)["results"]["passes_rdkit_sanity_checks"] is True

    mol_pass = MolFromSmiles("CCO", sanitize=True)
    assert check_chemistry(mol_pass)["results"]["passes_rdkit_sanity_checks"] is True

    mol_fail = MolFromSmiles("F=F", sanitize=False)
    assert check_chemistry(mol_fail)["results"]["passes_rdkit_sanity_checks"] is False
    assert check_chemistry(mol_fail)["results"]["passes_valence_checks"] is False
    assert check_chemistry(mol_fail)["results"]["passes_kekulization"] is True

    mol_fail = MolFromSmiles("c1ccccc1c", sanitize=False)
    assert check_chemistry(mol_fail)["results"]["passes_rdkit_sanity_checks"] is False
    assert check_chemistry(mol_fail)["results"]["passes_valence_checks"] is True
    assert check_chemistry(mol_fail)["results"]["passes_kekulization"] is False
