from __future__ import annotations

from posebusters.modules.volume_overlap import check_volume_overlap


def test_check_volume_overlap(mol_rq3_x00, mol_rq3_x01, mol_rq3_x10):
    out = check_volume_overlap(mol_rq3_x00, mol_rq3_x00)
    assert out["results"]["volume_overlap"] == 1.0
    assert out["results"]["no_volume_clash"] is False

    out = check_volume_overlap(mol_rq3_x00, mol_rq3_x01)
    assert out["results"]["no_volume_clash"] is False

    out = check_volume_overlap(mol_rq3_x00, mol_rq3_x10)
    assert out["results"]["volume_overlap"] == 0.0
    assert out["results"]["no_volume_clash"] is True


def test_check_volume_overlap_without_hydrogens(mol_rq3_x00, mol_rq3_x01, mol_rq3_x10):
    out = check_volume_overlap(mol_rq3_x00, mol_rq3_x00, ignore_hydrogens=True)
    assert out["results"]["volume_overlap"] == 1.0
    assert out["results"]["no_volume_clash"] is False

    out = check_volume_overlap(mol_rq3_x00, mol_rq3_x01, ignore_hydrogens=True)
    assert out["results"]["no_volume_clash"] is False

    out = check_volume_overlap(mol_rq3_x00, mol_rq3_x10, ignore_hydrogens=True)
    assert out["results"]["volume_overlap"] == 0.0
    assert out["results"]["no_volume_clash"] is True
