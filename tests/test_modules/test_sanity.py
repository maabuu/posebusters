from __future__ import annotations

from rdkit.Chem.rdmolfiles import MolFromSmiles

from posebusters.modules.sanity import check_chemistry_using_inchi, check_chemistry_using_rdkit


def test_check_rdkit_sanity():
    mol_pass = MolFromSmiles("C=C", sanitize=True)
    assert check_chemistry_using_rdkit(mol_pass)["results"]["passes_rdkit_sanity_checks"] is True

    mol_pass = MolFromSmiles("c1ccccc1", sanitize=True)
    assert check_chemistry_using_rdkit(mol_pass)["results"]["passes_rdkit_sanity_checks"] is True

    mol_pass = MolFromSmiles("CCO", sanitize=True)
    assert check_chemistry_using_rdkit(mol_pass)["results"]["passes_rdkit_sanity_checks"] is True

    mol_fail = MolFromSmiles("F=F", sanitize=False)
    assert check_chemistry_using_rdkit(mol_fail)["results"]["passes_rdkit_sanity_checks"] is False
    assert check_chemistry_using_rdkit(mol_fail)["results"]["passes_valence_checks"] is False
    assert check_chemistry_using_rdkit(mol_fail)["results"]["passes_kekulization"] is True

    mol_fail = MolFromSmiles("c1ccccc1c", sanitize=False)
    assert check_chemistry_using_rdkit(mol_fail)["results"]["passes_rdkit_sanity_checks"] is False
    assert check_chemistry_using_rdkit(mol_fail)["results"]["passes_valence_checks"] is True
    assert check_chemistry_using_rdkit(mol_fail)["results"]["passes_kekulization"] is False


def test_check_chemistry_using_inchi_1afs(mol_pred_1afs_gen87):
    info = check_chemistry_using_inchi(mol_pred_1afs_gen87)

    assert info["results"]["inchi_convertible"] is False


def test_check_chemistry_using_inchi_1afs_94(mol_pred_1afs_gen94):
    info = check_chemistry_using_inchi(mol_pred_1afs_gen94)

    assert info["results"]["inchi_convertible"] is True


def test_check_chemistry_using_inchi_1jn2_3(mol_pred_1jn2_gen3):
    info = check_chemistry_using_inchi(mol_pred_1jn2_gen3)

    assert info["results"]["inchi_convertible"] is False


def test_check_chemistry_using_inchi_1jn2_62(mol_pred_1jn2_gen62):
    info = check_chemistry_using_inchi(mol_pred_1jn2_gen62)

    assert info["results"]["inchi_convertible"] is True
