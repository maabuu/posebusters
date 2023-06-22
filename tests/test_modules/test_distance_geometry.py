from __future__ import annotations

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
