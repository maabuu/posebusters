import math

from rdkit.Chem.rdmolfiles import MolFromMolFile, MolFromPDBFile, MolFromSmiles

from posebusters import PoseBusters

mols_table = "tests/conftest/sample_bust_docks_table.csv"

mol_conf_2 = "tests/conftest/mol_1a30_clash_2.sdf"
mol_conf_3 = "tests/conftest/mol_1a30_clash_3.sdf"

mol_pred_1ia1 = "tests/conftest/1IA1_TQ3/1IA1_TQ3_ligand.sdf"
mol_true_1ia1 = "tests/conftest/1IA1_TQ3/1IA1_TQ3_ligands.sdf"
mol_cond_1ia1 = "tests/conftest/1IA1_TQ3/1IA1_TQ3_no_ligand_no_solvent.pdb"

mol_pred_1w1p = "tests/conftest/1W1P_GIO/1W1P_GIO_ligand.sdf"
mol_true_1w1p = "tests/conftest/1W1P_GIO/1W1P_GIO_ligands.sdf"
mol_cond_1w1p = "tests/conftest/1W1P_GIO/1W1P_GIO_protein.pdb"

mol_pred_7w2p = "tests/conftest/7W2P_8AI/7W2P_8AI_ligand.sdf"
mol_true_7w2p = "tests/conftest/7W2P_8AI/7W2P_8AI_ligands.sdf"
mol_cond_7w2p = "tests/conftest/7W2P_8AI/7W2P_8AI_no_ligand_no_solvent.pdb"

mol_single_h = "tests/conftest/mol_single_hydrogen.sdf"

mol_smaller = "tests/conftest/2HA2_SCK_2HA3_CHT/2HA2_SCK_2HA3_CHT_smaller_ligand.sdf"
mol_larger = "tests/conftest/2HA2_SCK_2HA3_CHT/2HA2_SCK_2HA3_CHT_larger_ligand.sdf"
mol_cond_smaller = "tests/conftest/2HA2_SCK_2HA3_CHT/2HA2_SCK_2HA3_CHT_smaller_receptor.pdb"

mol_true_5ze6 = "tests/conftest/5ze6/5ze6_true.mol2"
mol_pred_5ze6 = "tests/conftest/5ze6/5ze6_pred.sdf"
mol_cond_5ze6 = "tests/conftest/5ze6/5ze6_cond.pdb"

mol_true_sanity = "tests/conftest/sanity_error/true.sdf"
mol_pred_sanity = "tests/conftest/sanity_error/pred.sdf"
mol_cond_sanity = "tests/conftest/sanity_error/protein.pdb"

mol_true_sanity_2 = "tests/conftest/sanity_error_2/true.sdf"
mol_pred_sanity_2 = "tests/conftest/sanity_error_2/pred.sdf"
mol_cond_sanity_2 = "tests/conftest/sanity_error_2/protein.pdb"


def test_bust_redocks_1ia1() -> None:
    posebusters = PoseBusters("redock")
    df = posebusters.bust([mol_pred_1ia1], mol_true_1ia1, mol_cond_1ia1)
    assert df.all(axis=1).values[0]


def test_bust_redocks_1w1p() -> None:
    posebusters = PoseBusters("redock")
    df = posebusters.bust([mol_pred_1w1p], mol_true_1w1p, mol_cond_1w1p)
    assert df.all(axis=1).values[0]


def test_bust_redocks_5ze6() -> None:
    # check that mol2 files as true molecule can be loaded

    posebusters = PoseBusters("redock")
    df = posebusters.bust([mol_pred_5ze6], mol_true_5ze6, mol_cond_5ze6)
    assert df["mol_true_loaded"].all()


def test_bust_docks() -> None:
    posebusters = PoseBusters("dock")
    df = posebusters.bust([mol_pred_1ia1], mol_cond=mol_cond_1ia1)
    assert df.all(axis=1).values[0]


def test_bust_mols() -> None:
    posebusters = PoseBusters("mol")

    # pass one not in list
    df = posebusters.bust(mol_pred_1ia1)
    assert df.all(axis=1).values[0]

    # pass list
    df = posebusters.bust([mol_pred_1ia1])
    assert df.all(axis=1).values[0]


def test_bust_mol_rdkit() -> None:
    posebusters = PoseBusters(config="mol")
    mol = MolFromSmiles("C")

    df = posebusters.bust(mol)
    assert df.all(axis=1).values[0]

    df = posebusters.bust([mol])
    assert df.all(axis=1).values[0]


def test_bust_mols_hydrogen(threshold=8) -> None:
    posebusters = PoseBusters("mol")
    df = posebusters.bust([mol_single_h])
    assert df.sum(axis=1).values[0] >= threshold  # energy ratio test fails


def test_bust_mols_consistency(atol=1e-6) -> None:
    # check that running the same molecule twice gives the same result

    posebusters = PoseBusters("mol")
    result_2 = posebusters.bust([mol_conf_2])

    posebusters = PoseBusters("mol")
    result_3 = posebusters.bust([mol_conf_2])

    # check numbers are approximately equal
    for v1, v2 in zip(result_2.values, result_3.values):
        if v1[2] == v2[2] or math.isnan(v1[2]):
            continue
        assert abs(v1[2] - v2[2]) < atol, f"{v1[0], v1[1]}: {v1[2]} != {v2[2]}"


def test_bust_gen() -> None:
    posebusters = PoseBusters("gen")
    result = posebusters.bust(mol_pred=mol_larger, mol_true=mol_smaller, mol_cond=mol_cond_smaller, full_report=True)
    assert result["sucos"].iloc[0] > 0.2


def test_bust_loaded_mols() -> None:
    posebuster = PoseBusters("redock")

    mol_pred = MolFromMolFile(mol_pred_1ia1)
    mol_true = MolFromMolFile(mol_true_1ia1)
    mol_cond = MolFromPDBFile(mol_cond_1ia1)

    df = posebuster.bust([mol_pred], mol_true, mol_cond)
    assert df["mol_pred_loaded"].all()
    assert df["mol_true_loaded"].all()
    assert df["mol_cond_loaded"].all()


def test_check_energy_ratio_14gs_0(mol_pred_14gs_gen0):
    posebuster = PoseBusters("mol")
    df = posebuster.bust(mol_pred_14gs_gen0)

    assert df["mol_pred_loaded"].all()
    # assert that there is not a NA value given everyything else is correct
    assert not (df.eq(False) | df.isna()).all(axis=1).any()


def test_check_energy_ratio_1afs_87(mol_pred_1afs_gen87):
    posebuster = PoseBusters("mol")
    df = posebuster.bust(mol_pred_1afs_gen87)

    assert df["mol_pred_loaded"].all()
    # assert that there is not a NA value given everyything else is correct
    assert not (df.eq(False) | df.isna()).all(axis=1).any()


def test_check_energy_ratio_1afs_94(mol_pred_1afs_gen94):
    posebuster = PoseBusters("mol")
    df = posebuster.bust(mol_pred_1afs_gen94)

    assert df["mol_pred_loaded"].all()
    # assert that there is not a NA value given everyything else is correct
    assert not (df.eq(False) | df.isna()).all(axis=1).any()


def test_check_energy_ratio_1jn2_3(mol_pred_1jn2_gen3):
    posebuster = PoseBusters("mol")
    df = posebuster.bust(mol_pred_1jn2_gen3)

    assert df["mol_pred_loaded"].all()
    # assert that there is not a NA value given everyything else is correct
    assert not (df.eq(False) | df.isna()).all(axis=1).any()


def test_check_energy_ratio_1jn2_62(mol_pred_1jn2_gen62):
    posebuster = PoseBusters("mol")
    df = posebuster.bust(mol_pred_1jn2_gen62)

    assert df["mol_pred_loaded"].all()
    # assert that there is not a NA value given everyything else is correct
    assert not (df.eq(False) | df.isna()).all(axis=1).any()


def test_check_sanity():
    posebusters = PoseBusters("redock")
    df = posebusters.bust([mol_pred_sanity], mol_true_sanity, mol_cond_sanity, full_report=True)
    assert df["mol_true_loaded"].all()
    assert (df["rmsd"] < 3).all()


def test_check_sanity_2():
    posebusters = PoseBusters("redock")
    df = posebusters.bust([mol_pred_sanity_2], mol_true_sanity_2, mol_cond_sanity_2, full_report=True)
    assert df["mol_true_loaded"].all()
    assert (df["rmsd"] < 3).all()
