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
@click.version_option()
@click.option("--debug", default=False, is_flag=True, help="Show debug output.")
def bust(debug=False):
    """Run MolBuster from the command line."""
    if debug:
        click.echo(f"Debug mode is on.")
    # TODO: make debug show all output
    pass


@click.command("redock", help="Check docked ligand(s) against crystal protein and ligand.")
@click.argument("mol_pred", type=click.Path(exists=True, path_type=Path))
@click.argument("mol_true", type=click.Path(exists=True, path_type=Path))
@click.argument("mol_cond", type=click.Path(exists=True, path_type=Path))
@click.option("--outfmt", type=click.Choice(["short", "long", "csv"]), default="short")
@click.option("--output", type=click.File("w"), default="-")
@click.option("--config", type=click.File("r"))
def bust_redock(mol_pred, mol_true, mol_cond, outfmt, output, config):
    """Run MolBuster on docked ligands with protein and native ligand."""
    mode = "redock" if config is None else config
    molbuster = MolBuster(mode).bust(mol_pred, mol_true, mol_cond)
    _print_results(molbuster.results, outfmt)


@click.command("dock", help="Check docked ligand(s) with protein.")
@click.argument("mol_pred", nargs=1, type=click.Path(exists=True, path_type=Path))
@click.argument("mol_cond", nargs=1, type=click.Path(exists=True, path_type=Path))
@click.option("--outfmt", type=click.Choice(["short", "long", "csv"]), default="short")
@click.option("--output", type=click.File("w"), default="-")
@click.option("--config", type=click.File("r"))
def bust_dock(mol_pred, mol_cond, outfmt, output, config):
    """Run MolBuster on docked ligands and protein."""
    mode = "dock" if config is None else config
    molbuster = MolBuster(mode).bust(mol_pred, None, mol_cond)
    _print_results(molbuster.results, outfmt)


@click.command("mol", help="Check molecule(s).")
@click.argument("molecule", type=click.Path(exists=True, path_type=Path))
@click.option("--outfmt", type=click.Choice(["short", "long", "csv"]), default="short")
@click.option("--output", type=click.File("w"), default="-")
@click.option("--config", type=click.File("r"))
def bust_mol(molecule, outfmt, output, config):
    """Run MolBuster on docked ligands only."""
    mode = "mol" if config is None else config
    molbuster = MolBuster(mode).bust(molecule, None, None)
    _print_results(molbuster.results, outfmt)


@click.command("table", help="Check files listed in .csv file.")
@click.argument("table", type=click.Path(exists=True, path_type=Path))
@click.option("--outfmt", type=click.Choice(["short", "long", "csv"]), default="short")
@click.option("--output", type=click.File("w"), default="-")
@click.option("--config", type=click.File("r"))
def bust_table(table: Path, outfmt, output, config):
    """Run MolBuster on multiple inputs listed in a .csv file."""
    file_paths = pd.read_csv(table, index_col=None)
    mode = _select_mode(file_paths) if config is None else config
    molbuster = MolBuster(mode).bust_table(file_paths)
    _print_results(molbuster.results, outfmt)


bust.add_command(bust_dock)
bust.add_command(bust_redock)
bust.add_command(bust_mol)
bust.add_command(bust_table)


def _print_results(df: pd.DataFrame, outfmt: str = "short"):
    if outfmt == "long":
        _create_long_output(df)
    elif outfmt == "csv":
        click.echo(df.to_csv(index=True, header=True))
    elif outfmt == "short":
        _create_short_results(df)
    else:
        raise ValueError(f"Unknown output format {outfmt}")


def _create_long_output(df: pd.DataFrame):
    # return df.T.to_string(index=True, header=False)
    df = df.T
    cols = df.columns
    click.echo("MolBuster test summary:")
    for col in cols:
        click.echo("--> " + " ".join(col))
        click.echo(df[[col]].to_string(index=True, header=False))


def _create_short_results(df: pd.DataFrame):
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

    click.echo("MolBuster test summary:")
    click.echo(results.to_string(index=True, header=False))


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
