from __future__ import annotations

from posebusters.modules.flatness import check_flatness


def test_check_flatness(mol_pm2, mol_cgb):
    # no rings
    out = check_flatness(mol_cgb)
    assert out["results"]["flatness_passes"] is True

    # rings
    out = check_flatness(mol_pm2)
    assert out["results"]["flatness_passes"] is True


def test_check_flatness_6yr2_t1c(mol_pred_6yr2_t1c, mol_true_6yr2_t1c):
    # pass
    out = check_flatness(mol_true_6yr2_t1c)
    assert out["results"]["flatness_passes"] is True

    # fail
    out = check_flatness(mol_pred_6yr2_t1c)
    assert out["results"]["flatness_passes"] is False
