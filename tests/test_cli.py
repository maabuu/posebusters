from __future__ import annotations

from click.testing import CliRunner

from molbuster.cli import bust_dock, bust_mol, bust_redock, bust_table

runner = CliRunner()

mol_pred = "tests/conftest/1ia1/1ia1_ligand.sdf"
mol_true = "tests/conftest/1ia1/1ia1_ligands.sdf"
mol_cond = "tests/conftest/1ia1/1ia1_protein_one_lig_removed.pdb"
mols_table = "tests/conftest/sample_bust_docks_table.csv"


def test_bust_redocks() -> None:
    result = runner.invoke(bust_redock, [mol_pred, mol_true, mol_cond])
    assert result.exit_code == 0


def test_bust_docks() -> None:
    result = runner.invoke(bust_dock, [mol_pred, mol_cond])
    assert result.exit_code == 0


def test_bust_mols() -> None:
    result = runner.invoke(bust_mol, [mol_pred])
    assert result.exit_code == 0


def test_bust_table() -> None:
    result = runner.invoke(bust_table, [mols_table])
    assert result.exit_code == 0
