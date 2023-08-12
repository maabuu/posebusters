from __future__ import annotations

from click.testing import CliRunner

from posebusters.cli import bust

runner = CliRunner()

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
    result = runner.invoke(
        bust, [mol_pred_1ia1, "-l", mol_true_1ia1, "-p", mol_cond_1ia1, "--full-report", "--outfmt", "csv"]
    )
    assert result.exit_code == 0


def test_bust_redocks_1w1p() -> None:
    result = runner.invoke(
        bust, [mol_pred_1w1p, "-l", mol_true_1w1p, "-p", mol_cond_1w1p, "--full-report", "--outfmt", "csv"]
    )
    assert result.exit_code == 0


def test_bust_docks() -> None:
    result = runner.invoke(bust, [mol_pred_1ia1, "-p", mol_cond_1ia1])
    assert result.exit_code == 0


def test_bust_mols() -> None:
    result = runner.invoke(bust, [mol_pred_1ia1])
    assert result.exit_code == 0


def test_bust_single_hydrogen() -> None:
    result = runner.invoke(bust, [mol_single_h])
    assert result.exit_code == 0


def test_consistency() -> None:
    result_2 = runner.invoke(bust, [mol_conf_2, "--outfmt", "long"]).output.splitlines()[2:]
    result_3 = runner.invoke(bust, [mol_conf_3, "--outfmt", "long"]).output.splitlines()[2:]

    # check numbers are approximately equal
    assert result_2 == result_3
