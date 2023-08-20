"""Command line interface for PoseBusters."""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd
from rdkit.Chem.rdchem import Mol

from . import __version__
from .posebusters import PoseBusters
from .tools.formatting import create_long_output, create_short_output

logger = logging.getLogger(__name__)


def main():
    """Safe entry point for PoseBusters from the command line."""
    parser = _parse_args(sys.argv[1:])
    try:
        bust(**vars(parser))
    except Exception as e:
        logger.error(e)


def bust(
    mols_pred: list[Mol | Path] = [],
    mol_true: Mol | Path | None = None,
    mol_cond: Mol | Path | None = None,
    table=None,
    outfmt="short",
    output=sys.stdout,
    config=None,
    debug=False,
    no_header=False,
    full_report=False,
    top_n=-1,
):
    """PoseBusters: Plausibility checks for generated molecule poses."""
    if table is None and len(mols_pred) == 0:
        raise ValueError("Provide either MOLS_PRED or TABLE.")
    elif table is not None:
        # run on table
        file_paths = pd.read_csv(table, index_col=None)
        mode = _select_mode(file_paths.columns.tolist()) if config is None else config
        posebusters = PoseBusters(mode, top_n=top_n, debug=debug)
        posebusters_results = posebusters.bust_table(file_paths)
    else:
        # run on single input
        d = {k for k, v in dict(mol_pred=mols_pred, mol_true=mol_true, mol_cond=mol_cond).items() if v}
        mode = _select_mode(d) if config is None else config
        posebusters = PoseBusters(mode, top_n=top_n, debug=debug)
        posebusters_results = posebusters.bust(mols_pred, mol_true, mol_cond)

    for i, results_dict in enumerate(posebusters_results):
        results = _dataframe_from_output(results_dict, posebusters.config, full_report)
        output.write(_format_results(results, outfmt, no_header, i))


def _parse_args(args):
    parser = argparse.ArgumentParser()

    # input
    parser.add_argument("mol_pred", default=[], type=argparse.FileType("r"), nargs="*", help="Predicted molecule(s).")
    parser.add_argument("-l", "--mol_true", type=argparse.FileType("r"), help="True molecule, e.g. crystal ligand.")
    parser.add_argument("-p", "--mol_cond", type=argparse.FileType("r"), help="Conditioning molecule, e.g. protein.")
    parser.add_argument("-t", "--table", type=argparse.FileType("r"), help="Run multiple inputs listed in a .csv file.")

    # output options
    parser.add_argument("-f", "--outfmt", choices=["short", "long", "csv"], default="short", help="Output format.")
    parser.add_argument("-o", "--output", type=argparse.FileType("w"), default=sys.stdout, help="Output file.")

    # config
    parser.add_argument("-c", "--config", type=argparse.FileType("r"), default=None, help="Configuration file.")
    parser.add_argument("--full-report", action="store_true", help="Print full report.")
    parser.add_argument("--no-header", action="store_true", help="Print without header.")
    parser.add_argument("--top-n", type=int, default=None, help="Run on top N results in MOL_PRED only.")

    # other
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")

    namespace = parser.parse_args(args)

    # check that either mol_pred or table was provided
    if namespace.table is None and len(namespace.mol_pred) == 0:
        parser.error("Provide either MOL_PRED or TABLE.\n")

    return namespace


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


def _select_mode(columns: Iterable[str]) -> str:
    # decide on mode to run for provided input table
    if "mol_pred" in columns and "mol_true" in columns and "mol_cond" in columns:
        mode = "redock"
    elif "mol_pred" in columns and ("protein" in columns) or ("mol_cond" in columns):
        mode = "dock"
    elif any(column in columns for column in ("mol_pred", "mols_pred", "molecule", "molecules", "molecule")):
        mode = "mol"
    else:
        raise NotImplementedError(f"No supported columns found in csv. Columns found are {columns}")
    return mode
