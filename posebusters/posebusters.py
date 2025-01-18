"""PoseBusters class for running all tests on a set of molecules."""

from __future__ import annotations

import inspect
import logging
from collections import defaultdict
from collections.abc import Generator, Iterable
from functools import partial
from pathlib import Path
from typing import Any, Callable

import pandas as pd
from rdkit.Chem.rdchem import Mol
from yaml import safe_load

from .modules.distance_geometry import check_geometry
from .modules.energy_ratio import check_energy_ratio
from .modules.flatness import check_flatness
from .modules.identity import check_identity
from .modules.intermolecular_distance import check_intermolecular_distance
from .modules.loading import check_loading
from .modules.rmsd import check_rmsd
from .modules.sanity import (
    check_all_atoms_connected,
    check_chemistry,
    check_chemistry_using_inchi,
    check_chemistry_using_rdkit,
)
from .modules.sucos import check_sucos
from .modules.volume_overlap import check_volume_overlap
from .tools.loading import safe_load_mol, safe_supply_mols

logger = logging.getLogger(__name__)


module_dict: dict[str, Callable] = {
    "loading": check_loading,
    "sanity": check_chemistry,  # for backwards compatibility
    "rdkit_sanity": check_chemistry_using_rdkit,
    "inchi_convertible": check_chemistry_using_inchi,
    "atoms_connected": check_all_atoms_connected,
    "identity": check_identity,
    "distance_geometry": check_geometry,
    "flatness": check_flatness,
    "energy_ratio": check_energy_ratio,
    "intermolecular_distance": check_intermolecular_distance,
    "volume_overlap": check_volume_overlap,
    "rmsd": check_rmsd,
    "sucos": check_sucos,
}
molecule_args = {"mol_cond", "mol_true", "mol_pred"}


class PoseBusters:
    """Class to run all tests on a set of molecules."""

    file_paths: pd.DataFrame
    module_name: list
    module_func: list
    module_args: list
    fname: list

    def __init__(self, config: str | dict[str, Any] = "redock", top_n: int | None = None):
        """Initialize PoseBusters object."""
        self.module_func: list  # dict[str, Callable]
        self.module_args: list  # dict[str, set[str]]

        if isinstance(config, str) and config in {"dock", "redock", "mol", "gen"}:
            logger.info("Using default configuration for mode %s.", config)
            with open(Path(__file__).parent / "config" / f"{config}.yml", encoding="utf-8") as config_file:
                self.config = safe_load(config_file)
        elif isinstance(config, dict):
            logger.info("Using configuration dictionary provided by user.")
            self.config = config
        else:
            logger.error("Configuration %s not valid. Provide 'dock', 'redock', 'mol', 'gen' or a dictionary.", config)
        assert len(set(self.config.get("tests", {}).keys()) - set(module_dict.keys())) == 0

        self.config["top_n"] = self.config.get("top_n", top_n)

        self.results: dict[tuple[str, str], list[tuple[str, str, Any]]] = defaultdict(list)

    def bust(
        self,
        mol_pred: Iterable[Mol | Path | str] | Mol | Path | str,
        mol_true: Mol | Path | str | None = None,
        mol_cond: Mol | Path | str | None = None,
        full_report: bool = False,
    ) -> pd.DataFrame:
        """Run tests on one or more molecules.

        Args:
            mol_pred: Generated molecule(s), e.g. de-novo generated molecule or docked ligand, with one or more poses.
            mol_true: True molecule, e.g. crystal ligand, with one or more poses.
            mol_cond: Conditioning molecule, e.g. protein.
            full_report: Whether to include all columns in the output or only the boolean ones specified in the config.

        Notes:
            - Molecules can be provided as rdkit molecule objects or file paths.

        Returns:
            Pandas dataframe with results.
        """
        mol_pred_list: Iterable[Mol | Path | str] = [mol_pred] if isinstance(mol_pred, (Mol, Path, str)) else mol_pred

        columns = ["mol_pred", "mol_true", "mol_cond"]
        self.file_paths = pd.DataFrame([[mol_pred, mol_true, mol_cond] for mol_pred in mol_pred_list], columns=columns)

        results_gen = self._run()

        df = pd.concat([_dataframe_from_output(d, self.config, full_report=full_report) for d in results_gen])
        df.index.names = ["file", "molecule"]
        df.columns = [c.lower().replace(" ", "_") for c in df.columns]

        return df

    def bust_table(self, mol_table: pd.DataFrame, full_report: bool = False) -> pd.DataFrame:
        """Run tests on molecules provided in pandas dataframe as paths or rdkit molecule objects.

        Args:
            mol_table: Pandas dataframe with columns "mol_pred", "mol_true", "mol_cond" containing paths to molecules.
            full_report: Whether to include all columns in the output or only the boolean ones specified in the config.

        Returns:
            Pandas dataframe with results.
        """
        self.file_paths = mol_table

        results_gen = self._run()

        df = pd.concat([_dataframe_from_output(d, self.config, full_report=full_report) for d in results_gen])
        df.index.names = ["file", "molecule"]
        df.columns = [c.lower().replace(" ", "_") for c in df.columns]

        return df

    def _run(self) -> Generator[dict, None, None]:
        """Run all tests on molecules provided in file paths.

        Yields:
            Generator of result dictionaries.
        """
        self._initialize_modules()

        for _, paths in self.file_paths.iterrows():
            mol_args = {}
            if "mol_cond" in paths and paths["mol_cond"] is not None:
                mol_cond_load_params = self.config.get("loading", {}).get("mol_cond", {})
                mol_args["mol_cond"] = safe_load_mol(path=paths["mol_cond"], **mol_cond_load_params)
            if "mol_true" in paths and paths["mol_true"] is not None:
                mol_true_load_params = self.config.get("loading", {}).get("mol_true", {})
                mol_args["mol_true"] = safe_load_mol(path=paths["mol_true"], **mol_true_load_params)

            mol_pred_load_params = self.config.get("loading", {}).get("mol_pred", {})
            for i, mol_pred in enumerate(safe_supply_mols(paths["mol_pred"], **mol_pred_load_params)):
                if self.config["top_n"] is not None and i >= self.config["top_n"]:
                    break

                mol_args["mol_pred"] = mol_pred

                results_key = (str(paths["mol_pred"]), self._get_name(mol_pred, i))

                for name, fname, func, args in zip(self.module_name, self.fname, self.module_func, self.module_args):
                    # pick needed arguments for module
                    args_needed = {k: v for k, v in mol_args.items() if k in args}
                    # loading takes all inputs
                    if fname == "loading":
                        args_needed = {k: args_needed.get(k, None) for k in args_needed}
                    # run module when all needed input molecules are valid Mol objects
                    if fname != "loading" and not all(args_needed.get(m, None) for m in args_needed):
                        module_output: dict[str, Any] = {"results": {}}
                    else:
                        module_output = func(**args_needed)

                    # save to object
                    self.results[results_key].extend([(name, k, v) for k, v in module_output["results"].items()])
                    # self.results[results_key]["details"].append(module_output["details"])

                # return results for this entry
                yield {results_key: self.results[results_key]}

    def _initialize_modules(self) -> None:
        self.module_name = []
        self.module_func = []
        self.module_args = []
        self.fname = []
        for module in self.config["modules"]:
            function = module_dict[module["function"]]
            parameters = module.get("parameters", {})
            module_args = set(inspect.signature(function).parameters).intersection({"mol_pred", "mol_true", "mol_cond"})

            self.module_name.append(module["name"])
            self.fname.append(module["function"])
            self.module_func.append(partial(function, **parameters))
            self.module_args.append(module_args)

    @staticmethod
    def _get_name(mol: Mol, i: int) -> str:
        if mol is None:
            return f"invalid_mol_at_pos_{i}"

        if not mol.HasProp("_Name") or mol.GetProp("_Name") == "":
            return f"mol_at_pos_{i}"

        return mol.GetProp("_Name")


def _dataframe_from_output(results_dict, config, full_report: bool = False) -> pd.DataFrame:
    d = {id: {(module, output): value for module, output, value in results} for id, results in results_dict.items()}
    df = pd.DataFrame.from_dict(d, orient="index")

    test_columns = [(c["name"], n) for c in config["modules"] for n in c.get("chosen_binary_test_output", [])]
    names_lookup = {(c["name"], k): v for c in config["modules"] for k, v in c.get("rename_outputs", {}).items()}
    suffix_lookup = {c["name"]: c["rename_suffix"] for c in config["modules"] if "rename_suffix" in c}

    available_columns = df.columns.tolist()
    missing_columns = [c for c in test_columns if c not in available_columns]
    extra_columns = [c for c in available_columns if c not in test_columns]
    columns = test_columns + extra_columns if full_report else test_columns

    df[missing_columns] = pd.NA
    df = df[columns]
    df.columns = [names_lookup.get(c, c[-1] + suffix_lookup.get(c[0], "")) for c in df.columns]

    return df
