from __future__ import annotations

import math

import numpy as np

from posebusters.modules.energy_ratio import check_energy_ratio


def test_check_energy_ratio(mol_pm2):
    # molecule has valid conformation
    out = check_energy_ratio(mol_pm2)
    assert math.isfinite(out["results"]["mol_pred_energy"])
    assert math.isfinite(out["results"]["ensemble_avg_energy"])
    assert math.isfinite(out["results"]["energy_ratio"])
    assert out["results"]["energy_ratio_passes"] is True


def test_check_energy_ratio_14gs_0(mol_pred_14gs_gen0):
    # molecule is chemically insane, but conformation looks alright
    out = check_energy_ratio(mol_pred_14gs_gen0)
    assert math.isfinite(out["results"]["mol_pred_energy"])
    assert math.isfinite(out["results"]["ensemble_avg_energy"])
    assert math.isfinite(out["results"]["energy_ratio"])
    assert out["results"]["energy_ratio_passes"] is True


def test_check_energy_ratio_1afs_87(mol_pred_1afs_gen87):
    # nan because molecule cannot be converted to valid InChI
    out = check_energy_ratio(mol_pred_1afs_gen87)
    assert math.isnan(out["results"]["mol_pred_energy"])
    assert math.isnan(out["results"]["ensemble_avg_energy"])
    assert math.isnan(out["results"]["energy_ratio"])
    assert math.isnan(out["results"]["energy_ratio_passes"])


def test_check_energy_ratio_1afs_94(mol_pred_1afs_gen94):
    # molecule is chemically insane, but conformation looks alright
    out = check_energy_ratio(mol_pred_1afs_gen94)
    assert math.isfinite(out["results"]["mol_pred_energy"])
    assert math.isnan(out["results"]["ensemble_avg_energy"])
    assert math.isnan(out["results"]["energy_ratio"])
    assert math.isnan(out["results"]["energy_ratio_passes"])


def test_check_energy_ratio_1jn2_3(mol_pred_1jn2_gen3):
    # nan because molecule cannot be converted to valid InChI
    out = check_energy_ratio(mol_pred_1jn2_gen3)
    assert math.isnan(out["results"]["mol_pred_energy"])
    assert math.isnan(out["results"]["ensemble_avg_energy"])
    assert math.isnan(out["results"]["energy_ratio"])
    assert math.isnan(out["results"]["energy_ratio_passes"])


def test_check_energy_ratio_1jn2_62(mol_pred_1jn2_gen62):
    # molecule is chemically insane, but conformation looks alright
    out = check_energy_ratio(mol_pred_1jn2_gen62)
    assert math.isfinite(out["results"]["mol_pred_energy"])
    assert math.isfinite(out["results"]["ensemble_avg_energy"])
    assert math.isfinite(out["results"]["energy_ratio"])
    assert out["results"]["energy_ratio_passes"] is True


def test_check_energy_ratio_approximate_consistency(mol_1a30_clash_2, mol_1a30_clash_3):
    energy_2 = check_energy_ratio(mol_1a30_clash_2)["results"]["mol_pred_energy"]
    energy_3 = check_energy_ratio(mol_1a30_clash_3)["results"]["mol_pred_energy"]

    # check numbers are approximately equal
    assert np.isclose(energy_2, energy_3)


def test_check_energy_ratio_disconnected_atoms(mol_disconnnected_atoms):
    # no bonds, just disconnected atoms -> energy is 0
    out = check_energy_ratio(mol_disconnnected_atoms)
    assert math.isfinite(out["results"]["mol_pred_energy"])
    assert math.isfinite(out["results"]["ensemble_avg_energy"])
    assert math.isfinite(out["results"]["energy_ratio"])
    assert out["results"]["energy_ratio_passes"] is True
