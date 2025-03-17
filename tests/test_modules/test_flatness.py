from __future__ import annotations

from posebusters.modules.flatness import check_flatness, flat_bonds, flat_rings, nonflat


def test_check_flatness(mol_pm2, mol_cgb):
    # no rings
    out = check_flatness(mol_cgb)
    assert out["results"]["flatness_passes"] is True

    # rings
    out = check_flatness(mol_pm2)
    assert out["results"]["flatness_passes"] is True


def test_check_flatness_6yr2_t1c(mol_pred_6yr2_t1c, mol_true_6yr2_t1c):
    # pass
    out = check_flatness(mol_true_6yr2_t1c, threshold_flatness=0.5)
    assert out["results"]["flatness_passes"] is True

    # fail
    out = check_flatness(mol_pred_6yr2_t1c)
    assert out["results"]["flatness_passes"] is False

    # fail rings
    out = check_flatness(mol_pred_6yr2_t1c, flat_systems=flat_rings)
    assert out["results"]["flatness_passes"] is False

    # pass bonds
    out = check_flatness(mol_pred_6yr2_t1c, flat_systems=flat_bonds, threshold_flatness=0.5)
    assert out["results"]["flatness_passes"] is False


def test_check_flatness_7ecr_sin(mol_lig_7ecr_sin):
    # pass because no rings or double bonds present
    out = check_flatness(mol_lig_7ecr_sin)
    assert out["results"]["flatness_passes"] is True


def test_check_nonflatness_synthetic(mols_flat_etkdgv3, mols_nonflat_etkdgv3):
    for mol in mols_flat_etkdgv3:
        # these are flat and should be ignored by the patterns
        out = check_flatness(mol, flat_systems=nonflat, check_nonflat=True)
        assert out["results"]["flatness_passes"]

    for mol in mols_nonflat_etkdgv3:
        # these should be non-flat
        out = check_flatness(mol, flat_systems=nonflat, check_nonflat=True)
        assert out["results"]["flatness_passes"]


def test_check_nonflatness_awj_ideal(mol_awj_ideal):
    # RDKit does not recognize the 5-membered ring as aromatic (2 < 4n + 2 for all n)
    out = check_flatness(mol_awj_ideal)
    assert out["results"]["flatness_passes"] is True
    out = check_flatness(mol_awj_ideal, flat_systems=nonflat, check_nonflat=True, threshold_flatness=0.05)
    assert out["results"]["flatness_passes"] is True


def test_check_nonflatness_o2u_ideal(mol_o2u_ideal):
    out = check_flatness(mol_o2u_ideal)
    assert out["results"]["flatness_passes"] is True
    out = check_flatness(mol_o2u_ideal, flat_systems=nonflat, check_nonflat=True, threshold_flatness=0.05)
    assert out["results"]["flatness_passes"] is True


def test_check_nonflatness_pip_ideal(mol_pip_ideal):
    out = check_flatness(mol_pip_ideal)
    assert out["results"]["flatness_passes"] is True
    out = check_flatness(mol_pip_ideal, flat_systems=nonflat, check_nonflat=True, threshold_flatness=0.05)
    assert out["results"]["flatness_passes"] is True


def test_check_nonflatness_pip_wrong(mol_pip_wrong):
    out = check_flatness(mol_pip_wrong)
    assert out["results"]["flatness_passes"] is True
    out = check_flatness(mol_pip_wrong, flat_systems=nonflat, check_nonflat=True, threshold_flatness=0.05)
    assert out["results"]["flatness_passes"] is False
