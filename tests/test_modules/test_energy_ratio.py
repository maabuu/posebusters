from __future__ import annotations

from molbuster.modules.energy_ratio import check_energy_ratio


def test_check_energy_ratio(mol_pm2):
    out = check_energy_ratio(mol_pm2)
    assert out["results"]["energy_ratio_passes"] is True
