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


def test_check_rmsd_bug(mol_true_1mmv, mol_pred_1mmv):
    out = check_rmsd(mol_true_1mmv, mol_pred_1mmv)
    assert out["results"]["rmsd"] < 10.0
    assert out["results"]["rmsd_within_threshold"] is False


def test_check_rmsd_tautomeric(mol_pred_1g9v, mol_true_1g9v):
    out = check_rmsd(mol_pred_1g9v, mol_true_1g9v)
    assert out["results"]["rmsd"] == pytest.approx(6.97287, abs=1e-4)


def test_check_rmsd_multiple_ground_truth(mol_one_true_1w1p, mol_true_1w1p):
    out = check_rmsd(mol_one_true_1w1p, mol_true_1w1p)
    assert out["results"]["rmsd"] == pytest.approx(0.0, abs=1e-6)
    assert out["results"]["rmsd_within_threshold"] is True
