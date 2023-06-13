from __future__ import annotations

import pytest

from molbusters.modules.rmsd import check_rmsd


def test_check_rmsd(mol_rq3_x00, mol_rq3_x01, mol_rq3_x10):
    out = check_rmsd(mol_rq3_x00, mol_rq3_x01)
    out["results"]["rmsd"] == pytest.approx(0.0, abs=1e-6)
    out["results"]["rmsd_within_threshold"] is True

    out = check_rmsd(mol_rq3_x00, mol_rq3_x01)
    out["results"]["rmsd"] == pytest.approx(1.0, abs=1e-6)

    out = check_rmsd(mol_rq3_x00, mol_rq3_x10)
    out["results"]["rmsd"] == pytest.approx(10.0, abs=1e-6)
