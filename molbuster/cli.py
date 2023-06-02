"""Command line interface for MolBuster."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import click
import pandas as pd
from yaml import safe_load

from .molbuster import MolBuster


def main():
    """Run MolBuster from the command line."""
    bust()


def create_options_from_template_config() -> Callable:
    """Create click options from template configuration file.

    Returns:
        Decorator that adds click options to a function.
    """
    config_options = safe_load(open(Path(__file__).parent / "config" / "template.yml")).get("parameters")

    options = []
    for module_name, module_conf in config_options.items():
        for option_name, option_default in module_conf.items():
            if isinstance(option_default, dict):
                for sub_option_name, sub_option_default in option_default.items():
                    name = module_name + "_" + option_name + "_" + sub_option_name
                    options.append({"long": name, "required": False, "default": sub_option_default})
            # elif isinstance(option_default, list):
            #     name = module_name + "_" + option_name
            #     _type = list[type(option_default[0])]
            #     options.append({"long": name, "required": False,
            #          "default": option_default, "nargs": -1, "type": _type})
            # TODO: need to support lists ?
            else:
                name = module_name + "_" + option_name
                options.append({"long": name, "required": False, "default": option_default})

    def decorator(f):
        for option in reversed(options):
            args = ("--" + option["long"],)
            kwargs = dict(
                required=option["required"],
                default=option["default"],
                nargs=option.get("nargs", 1),
            )
            click.option(*args, **kwargs)(f)
        return f

    return decorator


@click.group(help="MolBuster: A tool for evaluating docking results.")
def bust():
    """Run MolBuster from the command line."""
    pass


@click.command("redock", help="Check docked ligands against crystal protein and ligand.")
@click.argument("mol_pred", type=click.Path(exists=True, path_type=Path))
@click.argument("mol_true", type=click.Path(exists=True, path_type=Path))
@click.argument("mol_cond", type=click.Path(exists=True, path_type=Path))
@click.option("-c", "--config", type=click.File("r"))
def bust_redock(mol_pred: Path, mol_true: Path, mol_cond: Path, config: dict[str, Any] = {}, __outfmt: str = "short"):
    """Run MolBuster on docked ligands with protein and native ligand."""
    molbuster = MolBuster("redock").bust(mol_pred, mol_true, mol_cond)
    _print_results(molbuster.results, __outfmt)


@click.command("dock", help="Check docked ligand with protein.")
@click.argument("mol_pred", nargs=1, type=click.Path(exists=True, path_type=Path))
@click.argument("protein", nargs=1, type=click.Path(exists=True, path_type=Path))
@click.argument("--outfmt", type=click.Choice(["short", "long", "csv"]), default="short")
def bust_dock(mol_pred: Path, protein: Path, config: dict[str, Any] = {}, __outfmt: str = "short"):
    """Run MolBuster on docked ligands and protein."""
    molbuster = MolBuster("dock").bust(mol_pred, None, protein)
    _print_results(molbuster.results, __outfmt)


@click.command("mol", help="Check molecules.")
@click.argument("molecule", type=click.Path(exists=True, path_type=Path))
@click.argument("--outfmt", type=click.Choice(["short", "long", "csv"]), default="short")
def bust_mol(molecule: list[Path], config: dict[str, Any] = {}, __outfmt: str = "short"):
    """Run MolBuster on docked ligands only."""
    molbuster = MolBuster("mol").bust(molecule, None, None)
    _print_results(molbuster.results, __outfmt)


@click.command("table", help="Run MolBuster on multiple inputs listed in a .csv file.")
@click.argument("table", type=click.Path(exists=True, path_type=Path))
@click.argument("--outfmt", type=click.Choice(["short", "long", "csv"]), default="short")
def bust_table(table: Path, config: dict[str, Any] = {}, __outfmt: str = "short"):
    """Run MolBuster on multiple inputs listed in a .csv file."""
    file_paths = pd.read_csv(table, index_col=None)
    mode = _select_mode(file_paths)
    molbuster = MolBuster(mode).bust_table(file_paths)
    _print_results(molbuster.results, __outfmt)


bust.add_command(bust_dock)
bust.add_command(bust_redock)
bust.add_command(bust_mol)
bust.add_command(bust_table)


def _print_results(df: pd.DataFrame, outfmt: str = "short"):
    if outfmt == "long":
        click.echo(df.to_string(index=True, header=True))
    elif outfmt == "csv":
        click.echo(df.to_csv(index=True, header=True))
    elif outfmt == "short":
        if df.shape[0] == 1:
            click.echo("MolBuster test summary:")
            click.echo(_create_single_output(df))
        else:
            click.echo("MolBuster test summary:")
            click.echo(_create_short_results(df))
    else:
        raise ValueError(f"Unknown output format {outfmt}")


def _create_single_output(df: pd.DataFrame) -> str:
    return df.T.to_string(index=True, header=False)


def _create_short_results(df: pd.DataFrame) -> str:
    results = df.copy()
    results.columns = results.columns.to_flat_index()
    columns = results.columns
    passes = results[columns].sum(axis=1)
    # results = results.replace(True, ".").replace(False, "x").replace(np.nan, " ")
    # results["Tests"] = results[columns].agg("".join, axis=1)
    results["Passed tests"] = passes.apply(lambda x: f"passes ({x} / {len(columns)})")
    results["Pass"] = results[columns].all(axis=1)
    # results.index.names = ("File", "Name")
    results.index.name = None
    results = results[["Passed tests"]]
    return results.to_string(index=True, header=False)


def _select_mode(file_paths: pd.DataFrame) -> str:
    # decide on mode to run for provided input table
    columns = set(file_paths.columns)
    if "mol_pred" in columns and "mol_true" in columns and "mol_cond" in columns:
        mode = "dock"
    elif "mol_pred" in columns and ("protein" in columns) or ("mol_cond" in columns):
        mode = "redock"
    elif any(column in columns for column in ("mol_pred", "ligand", "molecules", "molecule")):
        mode = "mol"
    else:
        raise NotImplementedError(f"No supported columns found in csv. Columns found are {columns}")
    return mode
