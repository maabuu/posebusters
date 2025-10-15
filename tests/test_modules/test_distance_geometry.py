from __future__ import annotations

from rdkit.Chem.rdmolfiles import MolFromPDBFile

from posebusters.modules.distance_geometry import _two_bonds_to_angle, check_geometry


def test_check_geometry(mol_pm2):
    out = check_geometry(mol_pm2)
    assert out["results"]["bond_lengths_within_bounds"] is True
    assert out["results"]["bond_angles_within_bounds"] is True
    assert out["results"]["no_internal_clash"] is True


def test_two_bonds_to_angle():
    assert _two_bonds_to_angle((0, 1), (0, 1)) is None
    assert _two_bonds_to_angle((0, 1), (1, 0)) is None
    assert _two_bonds_to_angle((0, 1), (1, 2)) == (0, 1, 2)
    assert _two_bonds_to_angle((0, 1), (2, 1)) == (0, 1, 2)
    assert _two_bonds_to_angle((1, 0), (1, 2)) == (0, 1, 2)
    assert _two_bonds_to_angle((1, 0), (2, 1)) == (0, 1, 2)
    assert _two_bonds_to_angle((0, 1), (2, 3)) is None


def test_molecule_issue_80():
    """Test sanitizable molecule from issue 80 that causes a RunTime error in GetMoleculeBoundsMatrix."""

    mol = MolFromPDBFile("tests/conftest/mol_issue_80.pdb", removeHs=False, sanitize=True)
    assert mol is not None
    # should not raise an error
    check_geometry(mol)
