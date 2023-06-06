from __future__ import annotations

import pytest
from rdkit.Chem.rdmolfiles import MolFromMolFile, MolFromPDBFile


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
def mol_1a30_clash_2():
    return MolFromMolFile("tests/conftest/mol_1a30_clash_2.sdf")


@pytest.fixture
def mol_1a30_clash_3():
    return MolFromMolFile("tests/conftest/mol_1a30_clash_3.sdf")
