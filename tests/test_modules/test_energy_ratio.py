from __future__ import annotations

import numpy as np

from posebusters.modules.energy_ratio import check_energy_ratio


def test_check_energy_ratio(mol_pm2):
    out = check_energy_ratio(mol_pm2)
    assert out["results"]["energy_ratio_passes"] is True


def test_check_energy_ratio_approximate_consistency(mol_1a30_clash_2, mol_1a30_clash_3):
    energy_2 = check_energy_ratio(mol_1a30_clash_2)["results"]["mol_pred_energy"]
    energy_3 = check_energy_ratio(mol_1a30_clash_3)["results"]["mol_pred_energy"]

    # check numbers are approximately equal
    assert np.isclose(energy_2, energy_3)
