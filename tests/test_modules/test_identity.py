from __future__ import annotations

from rdkit.Chem.rdmolfiles import MolFromSmiles

from posebusters.modules.identity import (
    check_identity,
    standardize_and_get_inchi,
)
from posebusters.tools.inchi import is_valid_inchi, split_inchi

mol_ethane = MolFromSmiles("CC")
mol_ethene = MolFromSmiles("C=C")
mol_ethyne = MolFromSmiles("C#C")

mol_benzene = MolFromSmiles("C1=CC=CC=C1")
mol_cyclohexane = MolFromSmiles("C1CCCCC1")

mol_ethylamine = MolFromSmiles("CCN")
mol_ethylaminium = MolFromSmiles("CC[NH3+]")
mol_ethanaminide = MolFromSmiles("CC[NH-]")

mol_sp3_anticlockwise = MolFromSmiles("N[C@H](C)C(=O)O")
mol_sp3_clockwise = MolFromSmiles("N[C@@H](C)C(=O)O")

mol_thiamine = MolFromSmiles("CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO")


def test_check_identity_positive():
    out = check_identity(mol_benzene, mol_benzene)
    assert out["results"]["inchi_overall"] is True


def test_check_identity_negative():
    out = check_identity(mol_benzene, mol_cyclohexane)
    assert out["results"]["inchi_overall"] is False
    assert out["results"]["formula"] is False
    assert out["results"]["connections"] is True


def test_check_identity_tetrahedral():
    out = check_identity(mol_sp3_clockwise, mol_sp3_clockwise)
    assert out["results"]["inchi_overall"] is True
    assert out["results"]["stereo_type"] is True
    assert out["results"]["stereo_sp3"] is True
    assert out["results"]["stereo_sp3_inverted"] is True
    assert out["results"]["stereo_tetrahedral"] is True

    out = check_identity(mol_sp3_clockwise, mol_sp3_anticlockwise)
    assert out["results"]["inchi_overall"] is False
    assert out["results"]["stereo_type"] is True
    assert out["results"]["stereo_tetrahedral"] is False


def test_standardize_and_get_inchi():
    assert standardize_and_get_inchi(mol_ethane) == "InChI=1S/C2H6/c1-2/h1-2H3"

    # bond order matters
    assert standardize_and_get_inchi(mol_ethane) != standardize_and_get_inchi(mol_ethene)
    assert standardize_and_get_inchi(mol_benzene) != standardize_and_get_inchi(mol_cyclohexane)

    # protonation state to be normalized
    assert standardize_and_get_inchi(mol_ethylamine) == standardize_and_get_inchi(mol_ethylaminium)
    assert standardize_and_get_inchi(mol_ethylamine) == standardize_and_get_inchi(mol_ethanaminide)


def test_split_inchi():
    inchi_ethane = standardize_and_get_inchi(mol_ethane)
    split_inchi_ethane = split_inchi(inchi_ethane)
    assert split_inchi_ethane["/"] == "C2H6"
    assert split_inchi_ethane["/c"] == "1-2"
    assert split_inchi_ethane["/h"] == "1-2H3"


def test_is_valid_inchi():
    phenol_inchi = "InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H"
    assert is_valid_inchi(phenol_inchi) is True
    assert is_valid_inchi("1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H") is False
    assert is_valid_inchi("InChI=1S/c7-6-4-2-1-3-5-6/h1-5,7H") is False
    assert is_valid_inchi("InChI=1S/C6H6O/c7-6-4-1-3-5-6/h1-5,7H") is False
    assert is_valid_inchi("InChI=1S/C6H6O/c7-6-4-1-3-5-6") is False
