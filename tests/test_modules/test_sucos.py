from __future__ import annotations

import pytest

from posebusters.modules.sucos import check_sucos


def test_check_sucos(mol_rq3_x00, mol_rq3_x01, mol_rq3_x10):
    out = check_sucos(mol_rq3_x00, mol_rq3_x00)
    assert out["results"]["sucos"] == pytest.approx(1.0, abs=1e-6)

    out = check_sucos(mol_rq3_x00, mol_rq3_x01)
    assert out["results"]["sucos"] == pytest.approx(0.5, abs=0.1)

    out = check_sucos(mol_rq3_x00, mol_rq3_x10)
    assert out["results"]["sucos"] == pytest.approx(0.0, abs=1e-6)


def test_check_sucos_multiple_ground_truth(mol_one_true_1w1p, mol_true_1w1p):
    out = check_sucos(mol_one_true_1w1p, mol_true_1w1p)
    assert out["results"]["sucos"] == pytest.approx(1.0, abs=1e-6)
