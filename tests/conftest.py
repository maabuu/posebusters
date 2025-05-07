from __future__ import annotations

import pytest
from rdkit import Chem
from rdkit.Chem.rdDistGeom import EmbedMolecule, srETKDGv3
from rdkit.Chem.rdmolfiles import MolFromMol2File, MolFromMolFile, MolFromPDBFile, SDMolSupplier


@pytest.fixture
def mol_cgb():
    return MolFromMolFile("tests/conftest/mol_CGB.sdf")


@pytest.fixture
def mol_pm2():
    return MolFromMolFile("tests/conftest/mol_PM2.sdf")


@pytest.fixture
def mol_rq3_x00():
    return MolFromMolFile("tests/conftest/mol_RQ3_x00.sdf")


@pytest.fixture
def mol_rq3_x01():
    return MolFromMolFile("tests/conftest/mol_RQ3_x01.sdf")


@pytest.fixture
def mol_rq3_x10():
    return MolFromMolFile("tests/conftest/mol_RQ3_x10.sdf")


@pytest.fixture
def ligand_2bm2():
    return MolFromMolFile("tests/conftest/mol_2bm2_ligand.sdf")


@pytest.fixture
def protein_2bm2():
    return MolFromPDBFile(
        "tests/conftest/mol_2bm2_protein_one_ligand_removed.pdb", sanitize=False, removeHs=False, proximityBonding=False
    )


@pytest.fixture
def mol_1a30_ligand():
    return MolFromMol2File("tests/conftest/mol_1a30_ligand.mol2", sanitize=False)


@pytest.fixture
def mol_1a30_clash_2():
    return MolFromMolFile("tests/conftest/mol_1a30_clash_2.sdf", sanitize=False)


@pytest.fixture
def mol_1a30_clash_3():
    return MolFromMolFile("tests/conftest/mol_1a30_clash_3.sdf", sanitize=False)


@pytest.fixture
def mol_pred_14gs_gen0():
    return MolFromMolFile("tests/conftest/mol_14gs_0_out.sdf", sanitize=True)


@pytest.fixture
def mol_pred_1afs_gen87():
    return MolFromMolFile("tests/conftest/mol_1afs_87_out.sdf", sanitize=True)


@pytest.fixture
def mol_pred_1afs_gen94():
    return MolFromMolFile("tests/conftest/mol_1afs_94_out.sdf", sanitize=True)


@pytest.fixture
def mol_pred_1jn2_gen3():
    return MolFromMolFile("tests/conftest/mol_1jn2_3_out.sdf", sanitize=True)


@pytest.fixture
def mol_pred_1jn2_gen62():
    return MolFromMolFile("tests/conftest/mol_1jn2_62_out.sdf", sanitize=True)


@pytest.fixture
def mol_calcium():
    return MolFromMolFile("tests/conftest/mol_calcium.sdf", sanitize=True)


@pytest.fixture
def mol_cholesterol():
    return MolFromMolFile("tests/conftest/mol_cholesterol.sdf", sanitize=True)


@pytest.fixture
def mol_true_1g9v():
    path = "tests/conftest/1G9V_RQ3/1G9V_RQ3_ligands.sdf"
    with SDMolSupplier(str(path), sanitize=True, removeHs=True, strictParsing=True) as supplier:
        mol = next(supplier)
        while mol is None and supplier.atEnd() is False:
            mol = next(supplier)
        for mol_next in supplier:
            if mol_next is not None:
                mol.AddConformer(mol_next.GetConformer(), assignId=True)
        return mol


@pytest.fixture
def mol_one_true_1g9v():
    return MolFromMolFile("tests/conftest/1G9V_RQ3/1G9V_RQ3_ligand.sdf", sanitize=True)


@pytest.fixture
def mol_pred_1g9v():
    return MolFromMolFile("tests/conftest/1G9V_RQ3/1G9V_RQ3_gold_redock.sdf", sanitize=True)


@pytest.fixture
def mol_true_1w1p():
    path = "tests/conftest/1W1P_GIO/1W1P_GIO_ligands.sdf"
    with SDMolSupplier(str(path), sanitize=True, removeHs=True, strictParsing=True) as supplier:
        mol = next(supplier)
        while mol is None and supplier.atEnd() is False:
            mol = next(supplier)
        for mol_next in supplier:
            if mol_next is not None:
                mol.AddConformer(mol_next.GetConformer(), assignId=True)
        return mol


@pytest.fixture
def mol_one_true_1w1p():
    return MolFromMolFile("tests/conftest/1W1P_GIO/1W1P_GIO_ligand.sdf", sanitize=True)


@pytest.fixture
def mol_true_1mmv_3ar():
    return MolFromMolFile("tests/conftest/1MMV_3AR/1MMV_3AR_crystal.sdf", sanitize=True)


@pytest.fixture
def mol_pred_1mmv_3ar():
    return MolFromMolFile("tests/conftest/1MMV_3AR/1MMV_3AR_tankbind.sdf", sanitize=True)


@pytest.fixture
def mol_true_1of6_dty():
    return MolFromMolFile("tests/conftest/1OF6_DTY/1OF6_DTY_crystal.sdf", sanitize=True)


@pytest.fixture
def mol_pred_1of6_dty():
    return MolFromMolFile("tests/conftest/1OF6_DTY/1OF6_DTY_vina.sdf", sanitize=True)


@pytest.fixture
def mol_true_1q1g_mti():
    return MolFromMolFile("tests/conftest/1Q1G_MTI/1Q1G_MTI_crystal.sdf", sanitize=True)


@pytest.fixture
def mol_pred_1q1g_mti():
    return MolFromMolFile("tests/conftest/1Q1G_MTI/1Q1G_MTI_tankbind.sdf", sanitize=True)


@pytest.fixture
def mol_true_6yr2_t1c():
    return MolFromMolFile("tests/conftest/6YR2_T1C/6YR2_T1C_crystal.sdf", sanitize=True)


@pytest.fixture
def mol_pred_6yr2_t1c():
    return MolFromMolFile("tests/conftest/6YR2_T1C/6YR2_T1C_tankbind.sdf", sanitize=True)


@pytest.fixture
def mol_lig_7ecr_sin():
    return MolFromMolFile("tests/conftest/7ECR_SIN/7ECR_SIN_ligand.sdf", sanitize=True)


@pytest.fixture
def mol_cond_7ecr_sin():
    return MolFromPDBFile("tests/conftest/7ECR_SIN/7ECR_SIN_protein.pdb", sanitize=False, proximityBonding=False)


@pytest.fixture
def mol_lig_7cnq_g8x():
    return MolFromMolFile("tests/conftest/7CNQ_G8X/7CNQ_G8X_ligand.sdf", sanitize=True)


@pytest.fixture
def mol_cond_7cnq_g8x():
    return MolFromPDBFile("tests/conftest/7CNQ_G8X/7CNQ_G8X_protein.pdb", sanitize=False, proximityBonding=False)


@pytest.fixture
def mol_lig_7ztl_bcn():
    return MolFromMolFile("tests/conftest/7ZTL_BCN/7ZTL_BCN_ligand.sdf", sanitize=True)


@pytest.fixture
def mol_cond_7ztl_bcn():
    return MolFromPDBFile("tests/conftest/7ZTL_BCN/7ZTL_BCN_protein.pdb", sanitize=False, proximityBonding=False)


@pytest.fixture
def mol_disconnnected_atoms():
    return MolFromMolFile("tests/conftest/mol_disconnected_atoms.sdf", sanitize=True)


@pytest.fixture
def mol_small_7brv_f5r_7wb6_f5r():
    return MolFromMolFile("tests/conftest/7BRV_F5R_7WB6_F5R/7BRV_F5R_7WB6_F5R_smaller_ligand.sdf", removeHs=True)


@pytest.fixture
def mol_large_7brv_f5r_7wb6_f5r():
    return MolFromMolFile("tests/conftest/7BRV_F5R_7WB6_F5R/7BRV_F5R_7WB6_F5R_larger_ligand.sdf", removeHs=True)


@pytest.fixture
def mol_065():
    return MolFromMolFile("tests/conftest/mol_065_ideal.sdf")


@pytest.fixture
def mol_065_left():
    return MolFromMolFile("tests/conftest/mol_065_left.sdf")


@pytest.fixture
def mol_065_right():
    return MolFromMolFile("tests/conftest/mol_065_right.sdf")


@pytest.fixture
def mol_TMO():
    return MolFromMolFile("tests/conftest/mol_TMO.sdf")


@pytest.fixture
def mol_2YU():
    return MolFromMolFile("tests/conftest/mol_2YU.mol")


@pytest.fixture
def mol_HQT():
    return MolFromMolFile("tests/conftest/mol_HQT.mol")


@pytest.fixture
def mol_3wrb_gde_true():
    return MolFromPDBFile("tests/conftest/mol_3WRB_1_GDE_0_ligand_true.pdb")


@pytest.fixture
def mol_3wrb_gde_pred():
    return MolFromMolFile("tests/conftest/mol_3WRB_1_GDE_0_ligand_pred.sdf")


def embed_mol(smi: str) -> Chem.Mol:
    hmol = Chem.AddHs(Chem.MolFromSmiles(smi))
    ps = srETKDGv3()
    ps.randomState = 42
    _ = EmbedMolecule(hmol, params=ps)
    return hmol


@pytest.fixture()
def mols_flat_etkdgv3():
    # non-aromatic flat
    smis_flat = ["C1=CC2=CC=CC2=C1", "C1=C[CH]C=C1", "N1C=CNC=C1", "C1C=CN=CN1"]
    # aromatic flat
    smis_flat.extend(["c1ccccc1", "c1cnc[nH]1"])
    return [embed_mol(smi) for smi in smis_flat]


@pytest.fixture()
def mols_nonflat_etkdgv3():
    # just nonflat
    smis_flat = ["C1CC(C)NC1", "N1C=CCC=C1", "C1C=CCCO1"]
    return [embed_mol(smi) for smi in smis_flat]


@pytest.fixture()
def mol_awj_ideal():
    return MolFromMolFile("tests/conftest/mol_AWJ_ideal.sdf")


@pytest.fixture()
def mol_o2u_ideal():
    return MolFromMolFile("tests/conftest/mol_O2U_ideal.sdf")


@pytest.fixture()
def mol_pip_ideal():
    return MolFromMolFile("tests/conftest/mol_PIP_ideal.sdf")


@pytest.fixture()
def mol_pip_wrong():
    return MolFromMolFile("tests/conftest/mol_PIP_wrong.sdf")
