from __future__ import annotations

from posebusters.modules.flatness import check_flatness

rings = {
    "aromatic_5_membered_rings_sp2": "[ar5^2]1[ar5^2][ar5^2][ar5^2][ar5^2]1",
    "aromatic_6_membered_rings_sp2": "[ar6^2]1[ar6^2][ar6^2][ar6^2][ar6^2][ar6^2]1",
}
bonds = {"trigonal_planar_double_bonds": "[C;X3;^2](*)(*)=[C;X3;^2](*)(*)"}


def test_check_flatness(mol_pm2, mol_cgb):
    # no rings
    out = check_flatness(mol_cgb)
    assert out["results"]["flatness_passes"] is True

    # rings
    out = check_flatness(mol_pm2)
    assert out["results"]["flatness_passes"] is True


def test_check_flatness_6yr2_t1c(mol_pred_6yr2_t1c, mol_true_6yr2_t1c):
    # pass
    out = check_flatness(mol_true_6yr2_t1c, threshold_flatness=0.5)
    assert out["results"]["flatness_passes"] is True

    # fail
    out = check_flatness(mol_pred_6yr2_t1c)
    assert out["results"]["flatness_passes"] is False

    # fail rings
    out = check_flatness(mol_pred_6yr2_t1c, flat_systems=rings)
    assert out["results"]["flatness_passes"] is False

    # pass bonds
    out = check_flatness(mol_pred_6yr2_t1c, flat_systems=bonds, threshold_flatness=0.5)
    assert out["results"]["flatness_passes"] is False


def test_check_flatness_7ecr_sin(mol_lig_7ecr_sin):
    # pass because no rings or double bonds present
    out = check_flatness(mol_lig_7ecr_sin)
    assert out["results"]["flatness_passes"] is True
