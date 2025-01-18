"""Command line interface for PoseBusters."""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Any, TextIO

import pandas as pd
from rdkit.Chem.rdchem import Mol
from yaml import safe_load

from . import __version__
from .posebusters import PoseBusters, _dataframe_from_output
from .tools.formatting import create_long_output, create_short_output

logger = logging.getLogger(__name__)


def main():
    """Safe entry point for PoseBusters from the command line."""
    parser = _parse_args(sys.argv[1:])
    try:
        bust(**vars(parser))
    except Exception as e:
        logger.error(e)


def bust(  # noqa: PLR0913
    mol_pred: list[Path | Mol] = [],
    mol_true: Path | Mol | None = None,
    mol_cond: Path | Mol | None = None,
    table: Path | None = None,
    outfmt: str = "short",
    output: Path | TextIO = sys.stdout,
    config: Path | None = None,
    no_header: bool = False,
    full_report: bool = False,
    top_n: int | None = None,
):
    """PoseBusters: Plausibility checks for generated molecule poses."""
    if table is None and len(mol_pred) == 0:
        raise ValueError("Provide either MOLS_PRED or TABLE.")

    if table is not None:
        # run on table
        file_paths = pd.read_csv(table, index_col=None)
        mode = _select_mode(config, file_paths.columns.tolist())
        posebusters = PoseBusters(mode, top_n=top_n)
        posebusters.file_paths = file_paths
        posebusters_results = posebusters._run()
    else:
        # run on single input
        d = {k for k, v in dict(mol_pred=mol_pred, mol_true=mol_true, mol_cond=mol_cond).items() if v}
        mode = _select_mode(config, d)
        posebusters = PoseBusters(mode, top_n=top_n)
        cols = ["mol_pred", "mol_true", "mol_cond"]
        posebusters.file_paths = pd.DataFrame([[mol_pred, mol_true, mol_cond] for mol_pred in mol_pred], columns=cols)
        posebusters_results = posebusters._run()

    if isinstance(output, Path):
        output = open(Path(output), "w", encoding="utf-8")

    for i, results_dict in enumerate(posebusters_results):
        results = _dataframe_from_output(results_dict, posebusters.config, full_report)
        output.write(_format_results(results, outfmt, no_header, i))


def _parse_args(args):
    desc = "PoseBusters: Plausibility checks for generated molecule poses."
    parser = argparse.ArgumentParser(description=desc, add_help=False)

    # Create two argument groups
    in_group = parser.add_argument_group(title="Input")
    out_group = parser.add_argument_group(title="Output")
    cfg_group = parser.add_argument_group(title="Configuration")
    inf_group = parser.add_argument_group(title="Information")

    # input
    help = "molecule(s) to check"
    in_group.add_argument("mol_pred", default=[], type=_path, nargs="*", help=help)
    in_group.add_argument("-l", dest="mol_true", type=_path, help="true molecule, e.g. crystal ligand")
    in_group.add_argument("-p", dest="mol_cond", type=_path, help="conditioning molecule, e.g. protein")
    help = "run multiple inputs listed in a .csv file"
    in_group.add_argument("-t", dest="table", type=_path, help=help)

    # output options
    out_group.add_argument("--outfmt", choices=["short", "long", "csv"], default="short", help="output format")
    out_group.add_argument("--output", type=Path, default=sys.stdout, help="output file (default: stdout)")
    # out_group.add_argument("--snake_case", action="store_false", help="use snake case for output columns")
    out_group.add_argument("--full-report", action="store_true", help="print details for each test")
    out_group.add_argument("--no-header", action="store_true", help="print output without header")

    # config
    cfg_group.add_argument("--config", type=_path, default=None, help="configuration file")
    cfg_group.add_argument(
        "--top-n", type=int, default=None, help="run on TOP_N results in MOL_PRED only (default: all)"
    )

    # other
    inf_group.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    inf_group.add_argument("-h", "--help", action="help", help="show this help message and exit")

    namespace = parser.parse_args(args)

    # check that either mol_pred or table was provided
    if namespace.table is None and len(namespace.mol_pred) == 0:
        parser.print_help()
        parser.exit(status=1, message="\nProvide either MOL_PRED or TABLE as input.\n")

    # full report only works with long and csv output
    if namespace.full_report and namespace.outfmt == "short":
        logger.warning("Option --full-report ignored. Please use --outfmt long or csv for --full-report.")
        namespace.full_report = False
    return namespace


def _format_results(df: pd.DataFrame, outfmt: str = "short", no_header: bool = False, index: int = 0) -> str:
    if outfmt == "long":
        return create_long_output(df)

    if outfmt == "csv":
        header = (not no_header) and (index == 0)
        df.index.names = ["file", "molecule"]
        df.columns = [c.lower().replace(" ", "_") for c in df.columns]
        return df.to_csv(index=True, header=header)

    if outfmt == "short":
        return create_short_output(df)

    raise ValueError(f"Unknown output format {outfmt}")


def _select_mode(config, columns: Iterable[str]) -> str | dict[str, Any]:
    # decide on mode to run

    # load config if provided
    if isinstance(config, Path):
        return dict(safe_load(open(config, encoding="utf-8")))

    # forward string if config provide
    if isinstance(config, str):
        return str(config)

    # select mode based on inputs
    if "mol_pred" in columns and "mol_true" in columns and "mol_cond" in columns:
        mode = "redock"
    elif "mol_pred" in columns and ("protein" in columns) or ("mol_cond" in columns):
        mode = "dock"
    elif any(column in columns for column in ("mol_pred", "mols_pred", "molecule", "molecules", "molecule")):
        mode = "mol"
    else:
        raise NotImplementedError(f"No supported columns found in csv. Columns found are {columns}")

    return mode


def _path(path_str: str):
    path = Path(path_str)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"File {path} not found!")
    return path
