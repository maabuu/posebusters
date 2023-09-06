import math

from rdkit.Chem.rdmolfiles import MolFromSmiles

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

mol_single_h = "tests/conftest/single_hydrogen.sdf"


def test_bust_redocks_1ia1() -> None:
    posebusters = PoseBusters("redock")
    df = posebusters.bust([mol_pred_1ia1], mol_true_1ia1, mol_cond_1ia1)
    assert df.all(axis=1).values[0]


def test_bust_redocks_1w1p() -> None:
    posebusters = PoseBusters("redock")
    df = posebusters.bust([mol_pred_1w1p], mol_true_1w1p, mol_cond_1w1p)
    assert df.all(axis=1).values[0]


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


def test_bust_mols_hydrogen() -> None:
    posebusters = PoseBusters("mol")
    df = posebusters.bust([mol_single_h])
    assert df.sum(axis=1).values[0] >= 8  # energy ratio test fails


def test_bust_mols_consistency() -> None:
    # check that running the same molecule twice gives the same result

    posebusters = PoseBusters("mol")
    result_2 = posebusters.bust([mol_conf_2])

    posebusters = PoseBusters("mol")
    result_3 = posebusters.bust([mol_conf_2])

    # check numbers are approximately equal
    for v1, v2 in zip(result_2.values, result_3.values):
        if v1[2] == v2[2] or math.isnan(v1[2]):
            continue
        assert abs(v1[2] - v2[2]) < 1e-6, f"{v1[0], v1[1]}: {v1[2]} != {v2[2]}"
