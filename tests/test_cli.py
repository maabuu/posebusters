from __future__ import annotations

from click.testing import CliRunner

from molbuster.cli import bust

runner = CliRunner()

mol_pred = "tests/conftest/1ia1/1ia1_ligand.sdf"
mol_true = "tests/conftest/1ia1/1ia1_ligands.sdf"
mol_cond = "tests/conftest/1ia1/1ia1_protein_one_lig_removed.pdb"
mols_table = "tests/conftest/sample_bust_docks_table.csv"

mol_conf_2 = "tests/conftest/mol_1a30_clash_2.sdf"
mol_conf_3 = "tests/conftest/mol_1a30_clash_3.sdf"


def test_bust_redocks() -> None:
    result = runner.invoke(bust, [mol_pred, "-l", mol_true, "-p", mol_cond])
    assert result.exit_code == 0


def test_bust_docks() -> None:
    result = runner.invoke(bust, [mol_pred, "-p", mol_cond])
    assert result.exit_code == 0


def test_bust_mols() -> None:
    result = runner.invoke(bust, [mol_pred])
    assert result.exit_code == 0


def test_bust_table() -> None:
    result = runner.invoke(bust, ["-t", mols_table])
    assert result.exit_code == 0


def test_consistency() -> None:
    result_2 = runner.invoke(bust, [mol_conf_2, "--outfmt", "long"]).output.splitlines()[2:]
    result_3 = runner.invoke(bust, [mol_conf_3, "--outfmt", "long"]).output.splitlines()[2:]

    # check numbers are approximately equal
    assert result_2 == result_3
