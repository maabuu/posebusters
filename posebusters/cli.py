"""Command line interface for PoseBusters."""
from __future__ import annotations

import logging
from pathlib import Path

import click
import pandas as pd

from .posebusters import PoseBusters
from .tools.formatting import create_long_output, create_short_output

logger = logging.getLogger(__name__)


def main():
    """Safe entry point for PoseBusters from the command line."""
    try:
        bust()
    except Exception as e:
        click.echo(e)


_path = click.Path(exists=True, path_type=Path)


@click.command(name="bust")
@click.argument("mol_pred", type=_path, required=False, default=None, nargs=-1)
@click.option("-l", "--mol_true", type=_path, required=False, default=None, help="True molecule, e.g. crystal ligand.")
@click.option("-p", "--mol_cond", type=_path, required=False, default=None, help="Conditioning molecule, e.g. protein.")
@click.option("-t", "--table", type=_path, help="Run multiple inputs listed in a .csv file.")
@click.option("-f", "--outfmt", type=click.Choice(["short", "long", "csv"]), default="short", help="Output format.")
@click.option(
    "-o", "--out", "output", type=click.File("w"), default="-", help="Output file. Prints to stdout by default."
)
@click.option("-c", "--config", type=click.File("r"), default=None, help="Configuration file.")
@click.option("--full-report", type=bool, default=False, is_flag=True, help="Print full report.")
@click.option("--no-header", type=bool, default=False, is_flag=True, help="Print without header.")
# @click.option("--print-header", type=bool, default=False, is_flag=True, help="Print header only.")
@click.option("--top-n", type=int, default=None, help="Run on top N results in MOL_PRED only.")
@click.option("--debug", type=bool, default=False, is_flag=True, help="Enable debug output.")
@click.version_option()
def bust(table, outfmt, output, config, debug, no_header, full_report, top_n, **mol_args):
    """PoseBusters: check generated 3D molecules with or without conditioning."""
    if debug:
        click.echo("Debug mode is on.")

    if outfmt == "short":
        if full_report:
            logger.warning("Full report is not available in short output format. Ignoring --full-report option.")
        full_report = False

    if table is None and mol_args.get("mol_pred") is None:
        # check that an input was provided
        click.echo("Provide either MOL_PRED or a table using the --table option.\n")
        click.echo(bust.get_help(click.Context(bust)))
        return None
    elif table is not None:
        # run on table
        file_paths = pd.read_csv(table, index_col=None)
        mode = _select_mode(file_paths.columns.tolist()) if config is None else config
        posebusters = PoseBusters(mode, top_n=top_n, debug=debug)
        posebusters_results = posebusters.bust_table(file_paths)
    else:
        # run on file inputs
        mode = _select_mode([m for m, v in mol_args.items() if v is not None]) if config is None else config
        posebusters = PoseBusters(mode, top_n=top_n, debug=debug)
        posebusters_results = posebusters.bust(**mol_args)

    for i, results_dict in enumerate(posebusters_results):
        results = _dataframe_from_output(results_dict, posebusters.config, full_report)
        output.write(_format_results(results, outfmt, no_header, i))


def _dataframe_from_output(results_dict, config, full_report: bool = False) -> pd.DataFrame:
    d = {id: {(module, output): value for module, output, value in results} for id, results in results_dict.items()}
    df = pd.DataFrame.from_dict(d, orient="index")

    test_columns = [(c["name"], n) for c in config["modules"] for n in c["chosen_binary_test_output"]]
    names_lookup = {(c["name"], k): v for c in config["modules"] for k, v in c["rename_outputs"].items()}
    suffix_lookup = {c["name"]: c["rename_suffix"] for c in config["modules"] if "rename_suffix" in c}

    available_columns = df.columns.tolist()
    missing_columns = [c for c in test_columns if c not in available_columns]
    extra_columns = [c for c in available_columns if c not in test_columns]
    columns = test_columns + extra_columns if full_report else test_columns

    df[missing_columns] = pd.NA
    df = df[columns]
    df.columns = [names_lookup.get(c, c[-1] + suffix_lookup.get(c[0], "")) for c in df.columns]

    return df


def _format_results(df: pd.DataFrame, outfmt: str = "short", no_header: bool = False, index: int = 0) -> str:
    if outfmt == "long":
        return create_long_output(df)
    elif outfmt == "csv":
        header = (not no_header) and (index == 0)
        return df.to_csv(index=True, header=header)
    elif outfmt == "short":
        return create_short_output(df)
    else:
        raise ValueError(f"Unknown output format {outfmt}")


def _select_mode(columns: list[str]) -> str:
    # decide on mode to run for provided input table
    if "mol_pred" in columns and "mol_true" in columns and "mol_cond" in columns:
        mode = "redock"
    elif "mol_pred" in columns and ("protein" in columns) or ("mol_cond" in columns):
        mode = "dock"
    elif any(column in columns for column in ("mol_pred", "molecule", "molecules", "molecule")):
        mode = "mol"
    else:
        raise NotImplementedError(f"No supported columns found in csv. Columns found are {columns}")
    return mode
