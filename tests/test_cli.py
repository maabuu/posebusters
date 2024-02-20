from __future__ import annotations

from pathlib import Path

from posebusters.cli import _parse_args, bust

mols_table = "tests/conftest/sample_bust_docks_table.csv"
mol_pred_1ia1 = "tests/conftest/1IA1_TQ3/1IA1_TQ3_ligand.sdf"
mol_true_1ia1 = "tests/conftest/1IA1_TQ3/1IA1_TQ3_ligands.sdf"
mol_cond_1ia1 = "tests/conftest/1IA1_TQ3/1IA1_TQ3_no_ligand_no_solvent.pdb"


def test_parse_args() -> None:
    args = _parse_args(
        [
            mol_pred_1ia1,
            "-l",
            mol_true_1ia1,
            "-p",
            mol_cond_1ia1,
            "--full-report",
            "--outfmt",
            "csv",
        ]
    )
    assert args.mol_pred[0] == Path(mol_pred_1ia1)
    assert args.mol_true == Path(mol_true_1ia1)
    assert args.mol_cond == Path(mol_cond_1ia1)
    assert args.full_report
    assert args.outfmt == "csv"

    args = _parse_args([mol_pred_1ia1, "--outfmt", "long"])
    assert args.mol_pred[0] == Path(mol_pred_1ia1)
    assert not args.full_report
    assert args.outfmt == "long"


def test_bust_mols() -> None:
    bust([Path(mol_pred_1ia1)], Path(mol_true_1ia1), Path(mol_cond_1ia1))


def test_bust_table() -> None:
    bust(table=Path(mols_table))


def test_output() -> None:
    output = Path(".test_output.csv")
    bust(table=Path(mols_table), outfmt="csv", output=output)
    assert output.exists()


def test_bust_none() -> None:
    # check that ValueError is raised when no molecules are provided
    try:
        bust()
    except ValueError:
        pass
    else:
        raise AssertionError("ValueError not raised.")
