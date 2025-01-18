from __future__ import annotations

import pytest

from posebusters.modules.rmsd import check_rmsd


def test_check_rmsd(mol_rq3_x00, mol_rq3_x01, mol_rq3_x10):
    out = check_rmsd(mol_rq3_x00, mol_rq3_x00)
    assert out["results"]["rmsd"] == pytest.approx(0.0, abs=1e-6)
    assert out["results"]["rmsd_within_threshold"] is True

    out = check_rmsd(mol_rq3_x00, mol_rq3_x01)
    assert out["results"]["rmsd"] == pytest.approx(1.0, abs=1e-6)

    out = check_rmsd(mol_rq3_x00, mol_rq3_x10)
    assert out["results"]["rmsd"] == pytest.approx(10.0, abs=1e-6)


def test_check_rmsd_1mmv_3ar(mol_true_1mmv_3ar, mol_pred_1mmv_3ar):
    out = check_rmsd(mol_true_1mmv_3ar, mol_pred_1mmv_3ar)
    assert out["results"]["rmsd"] < 10.0
    assert out["results"]["rmsd_within_threshold"] is False


def test_check_rmsd_1of6_dty(mol_true_1of6_dty, mol_pred_1of6_dty):
    out = check_rmsd(mol_true_1of6_dty, mol_pred_1of6_dty)
    assert out["results"]["rmsd"] < 2.0
    assert out["results"]["rmsd_within_threshold"] is True


def test_check_rmsd_1q1g_mti(mol_true_1q1g_mti, mol_pred_1q1g_mti):
    out = check_rmsd(mol_true_1q1g_mti, mol_pred_1q1g_mti)
    assert out["results"]["rmsd"] < 50.0
    assert out["results"]["rmsd_within_threshold"] is False


def test_check_rmsd_6yr2_t1c(mol_true_6yr2_t1c, mol_pred_6yr2_t1c):
    out = check_rmsd(mol_true_6yr2_t1c, mol_pred_6yr2_t1c)
    assert out["results"]["rmsd"] < 50.0
    assert out["results"]["rmsd_within_threshold"] is False


def test_check_rmsd_tautomeric(mol_pred_1g9v, mol_true_1g9v):
    out = check_rmsd(mol_pred_1g9v, mol_true_1g9v)
    assert out["results"]["rmsd"] == pytest.approx(6.97287, abs=1e-4)


def test_check_rmsd_multiple_ground_truth(mol_one_true_1w1p, mol_true_1w1p):
    out = check_rmsd(mol_one_true_1w1p, mol_true_1w1p)
    assert out["results"]["rmsd"] == pytest.approx(0.0, abs=1e-6)
    assert out["results"]["rmsd_within_threshold"] is True


def test_check_rmsd_bond_orders(mol_3wrb_gde_true, mol_3wrb_gde_pred):
    out = check_rmsd(mol_3wrb_gde_true, mol_3wrb_gde_true)
    assert out["results"]["rmsd"] == pytest.approx(0.0, abs=1e-6)

    out = check_rmsd(mol_3wrb_gde_pred, mol_3wrb_gde_pred)
    assert out["results"]["rmsd"] == pytest.approx(0.0, abs=1e-6)

    out = check_rmsd(mol_3wrb_gde_pred, mol_3wrb_gde_true)
    assert out["results"]["rmsd"] == pytest.approx(1.091668, abs=1e-4)

    out = check_rmsd(mol_3wrb_gde_true, mol_3wrb_gde_pred)
    assert out["results"]["rmsd"] == pytest.approx(1.091668, abs=1e-4)
