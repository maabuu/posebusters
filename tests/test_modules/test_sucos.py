from __future__ import annotations

import pytest
from rdkit.Chem.rdmolops import AddHs, RemoveHs

from posebusters.modules.sucos import check_sucos


def test_check_sucos(mol_rq3_x00, mol_rq3_x01, mol_rq3_x10, mol_small_7brv_f5r_7wb6_f5r, mol_large_7brv_f5r_7wb6_f5r):
    out = check_sucos(mol_rq3_x00, mol_rq3_x00)
    assert out["results"]["sucos"] == pytest.approx(1.0, abs=1e-6)

    out = check_sucos(mol_rq3_x00, mol_rq3_x01)
    assert out["results"]["sucos"] == pytest.approx(0.5, abs=0.1)

    out = check_sucos(mol_rq3_x00, mol_rq3_x10)
    assert out["results"]["sucos"] == pytest.approx(0.0, abs=1e-6)


def test_check_sucos_hydrogens(mol_small_7brv_f5r_7wb6_f5r, mol_large_7brv_f5r_7wb6_f5r):
    out = check_sucos(
        AddHs(mol_large_7brv_f5r_7wb6_f5r, addCoords=True), AddHs(mol_small_7brv_f5r_7wb6_f5r, addCoords=True)
    )
    assert out["results"]["sucos"] <= 1.0

    out = check_sucos(RemoveHs(mol_large_7brv_f5r_7wb6_f5r), RemoveHs(mol_small_7brv_f5r_7wb6_f5r))
    assert out["results"]["sucos"] <= 1.0

    out = check_sucos(AddHs(mol_large_7brv_f5r_7wb6_f5r, addCoords=True), RemoveHs(mol_small_7brv_f5r_7wb6_f5r))
    assert out["results"]["sucos"] <= 1.0

    out = check_sucos(RemoveHs(mol_small_7brv_f5r_7wb6_f5r), AddHs(mol_large_7brv_f5r_7wb6_f5r, addCoords=True))
    assert out["results"]["sucos"] <= 1.0


def test_check_sucos_multiple_ground_truth(mol_one_true_1w1p, mol_true_1w1p):
    out = check_sucos(mol_one_true_1w1p, mol_true_1w1p)
    assert out["results"]["sucos"] == pytest.approx(1.0, abs=1e-6)


def test_get_sucos_score_065(mol_065, mol_065_left, mol_065_right):
    out = check_sucos(mol_065, mol_065)
    assert out["results"]["sucos"] == pytest.approx(1.0, abs=1e-6)

    out = check_sucos(mol_065_left, mol_065_right)
    assert out["results"]["sucos"] == pytest.approx(0.0, abs=1e-6)


def test_get_sucos_score_TMO(mol_TMO):
    # this molecule has no features
    out = check_sucos(mol_TMO, mol_TMO)
    assert out["results"]["sucos"] == 1.0


def test_get_sucos_score_2YU_HQT(mol_2YU, mol_HQT):
    out = check_sucos(mol_2YU, mol_HQT)
    assert out["results"]["sucos"] <= 1.0
