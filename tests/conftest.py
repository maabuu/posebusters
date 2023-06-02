from __future__ import annotations

import pytest
from rdkit.Chem.rdmolfiles import (
    MolFromMol2File,
    MolFromMolFile,
    MolFromPDBFile,
    SDMolSupplier,
)


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
def mol_calcium():
    return MolFromMolFile("tests/conftest/mol_calcium.sdf", sanitize=True)


@pytest.fixture
def mol_cholesterol():
    return MolFromMolFile("tests/conftest/mol_cholesterol.sdf", sanitize=True)


@pytest.fixture
def mol_true_1g9v():
    path = "tests/conftest/1G9V_RQ3/1G9V_RQ3_ligands.sdf"
    supplier = SDMolSupplier(str(path), sanitize=True, removeHs=True, strictParsing=True)
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
    supplier = SDMolSupplier(str(path), sanitize=True, removeHs=True, strictParsing=True)
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
