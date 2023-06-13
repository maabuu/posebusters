from __future__ import annotations

from molbusters.modules.flatness import check_flatness


def test_check_flatness(mol_pm2, mol_cgb):
    # no rings
    out = check_flatness(mol_cgb)
    assert out["results"]["flatness_passes"] is True

    # rings
    out = check_flatness(mol_pm2)
    assert out["results"]["flatness_passes"] is True
