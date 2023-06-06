"""Command line interface for MolBuster."""
from __future__ import annotations

from pathlib import Path

import click
import pandas as pd

from .molbuster import MolBuster
from .tools.formatting import _create_long_output, _create_short_output


def main():
    """Run MolBuster from the command line."""
    bust()


@click.command()
@click.argument("mol_pred", type=click.Path(exists=True, path_type=Path), required=False)
@click.option(
    "-l", "--mol_true", type=click.Path(exists=True, path_type=Path), required=False, default=None, help="True ligand."
)
@click.option(
    "-p",
    "--mol_cond",
    type=click.Path(exists=True, path_type=Path),
    required=False,
    default=None,
    help="Protein receptor.",
)
@click.option(
    "-t", "--table", type=click.Path(exists=True, path_type=Path), help="Run multiple inputs listed in a .csv file."
)
@click.option("-f", "--outfmt", type=click.Choice(["short", "long", "csv"]), default="short", help="Output format.")
@click.option(
    "-o", "--out", "output", type=click.File("w"), default="-", help="Output file. Prints to stdout by default."
)
@click.option("-c", "--config", type=click.File("r"), default=None, help="Configuration file.")
# @click.option("--debug", type=bool, default=False, is_flag=True, help="Enable debug output.")
@click.version_option()
def bust(table, outfmt, output, config, debug=False, **mol_args):
    """MolBuster: check generated 3D molecules with or without conditioning."""
    if debug:
        click.echo("Debug mode is on.")

    if table is None and mol_args.get("mol_pred") is None:
        # check that an input was provided
        click.echo("Provide either MOL_PRED or a table using the --table option.\n")
        click.echo(bust.get_help(click.Context(bust)))
        return None
    elif table is not None:
        # run on table
        file_paths = pd.read_csv(table, index_col=None)
        mode = _select_mode(file_paths.columns.tolist()) if config is None else config
        molbuster = MolBuster(mode)
        molbuster_results = molbuster.bust_table(file_paths)
    else:
        # run on file inputs
        mode = _select_mode([m for m, v in mol_args.items() if v is not None]) if config is None else config
        molbuster = MolBuster(mode)
        molbuster_results = molbuster.bust(**mol_args)

    for i, result in enumerate(molbuster_results):
        output.write(_print_results(result, outfmt, i))

    pass


def _print_results(df: pd.DataFrame, outfmt: str = "short", index: int = 0) -> str:
    if outfmt == "long":
        return _create_long_output(df)
    elif outfmt == "csv":
        header = True if index == 0 else False
        return df.to_csv(index=True, header=header)
    elif outfmt == "short":
        return _create_short_output(df)
    else:
        raise ValueError(f"Unknown output format {outfmt}")


def _select_mode(columns: list[str]) -> str:
    # decide on mode to run for provided input table
    if "mol_pred" in columns and "mol_true" in columns and "mol_cond" in columns:
        mode = "redock"
    elif "mol_pred" in columns and ("protein" in columns) or ("mol_cond" in columns):
        mode = "dock"
    elif any(column in columns for column in ("mol_pred", "ligand", "molecules", "molecule")):
        mode = "mol"
    else:
        raise NotImplementedError(f"No supported columns found in csv. Columns found are {columns}")
    return mode
