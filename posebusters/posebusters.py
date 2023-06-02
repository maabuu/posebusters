"""PoseBusters class for running all tests on a set of molecules."""
from __future__ import annotations

import inspect
import logging
from collections import defaultdict
from functools import partial
from pathlib import Path
from typing import Any, Callable, Generator

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
from .modules.sanity import check_chemistry
from .modules.volume_overlap import check_volume_overlap
from .tools.loading import safe_load_mol, safe_supply_mols

logger = logging.getLogger(__name__)


module_dict: dict[str, Callable] = {
    "loading": check_loading,
    "sanity": check_chemistry,
    "identity": check_identity,
    "distance_geometry": check_geometry,
    "flatness": check_flatness,
    "energy_ratio": check_energy_ratio,
    "intermolecular_distance": check_intermolecular_distance,
    "volume_overlap": check_volume_overlap,
    "rmsd": check_rmsd,
}
molecule_args = {"mol_cond", "mol_true", "mol_pred"}


class PoseBusters:
    """Class to run all tests on a set of molecules."""

    def __init__(self, config: str | dict = "redock", debug: bool = False, top_n: int | None = None):
        """Initialize PoseBusters object."""
        self.module_func: list  # dict[str, Callable]
        self.module_args: list  # dict[str, set[str]]

        if isinstance(config, str) and config in {"dock", "redock", "mol"}:
            logger.info(f"Using default configuration for mode {config}.")
            self.config = safe_load(open(Path(__file__).parent / "config" / f"{config}.yml"))
        elif isinstance(config, dict):
            logger.info("Using configuration dictionary provided by user.")
            self.config = config
        else:
            logger.error(f"Configuration {config} not valid. Provide either 'dock', 'redock', 'mol' or a dictionary.")
        assert len(set(self.config.get("tests", {}).keys()) - set(module_dict.keys())) == 0

        self.config["debug"] = self.config.get("debug", debug)
        self.config["top_n"] = self.config.get("top_n", top_n)

        self.results: dict[tuple[str, str], list[tuple[str, str, Any]]] = defaultdict(list)

    def bust(
        self, mol_pred: tuple[Mol | Path], mol_true: Mol | Path | None, mol_cond: Mol | Path | None
    ) -> Generator[dict, None, None]:
        """Run all tests on one molecule.

        Args:
            mol_pred: Generated molecule, e.g. docked ligand, with one or more poses.
            mol_true: True molecule, e.g. crytal ligand, with one or more poses.
            mol_cond: Conditioning molecule, e.g. protein.

        Notes:
            - Molecules can be provided as rdkit molecule objects or file paths.

        Returns:
            PoseBusters object.
        """
        columns = ["mol_pred", "mol_true", "mol_cond"]
        self.file_paths = pd.DataFrame([[mol_pred, mol_true, mol_cond] for mol_pred in mol_pred], columns=columns)
        return self._run()

    def bust_table(self, mol_table: pd.DataFrame) -> Generator[dict, None, None]:
        """Run all tests on multiple molecules provided in pandas dataframe as paths or rdkit molecule objects.

        Args:
            mol_table: Pandas dataframe with columns "mol_pred", "mol_true", "mol_cond" containing paths to molecules.

        Returns:
            PoseBusters object.
        """
        self.file_paths = mol_table
        return self._run()

    def _run(self) -> Generator[dict, None, None]:
        """Run all tests on molecules provided in file paths.

        Yields:
            Generator of pandas dataframes with results.
        """
        self._initialize_modules()

        for i, paths in self.file_paths.iterrows():
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
                    args = {k: v for k, v in mol_args.items() if k in args}
                    # loading takes all inputs
                    if fname == "loading":
                        args = {k: args.get(k, None) for k in args}
                    # run module when all needed input molecules are valid Mol objects
                    if fname != "loading" and not all(args.get(m, None) for m in args):
                        module_output: dict[str, Any] = {"results": {}}
                    else:
                        module_output = func(**args)

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
        pass

    @staticmethod
    def _get_name(mol: Mol, i: int) -> str:
        if mol is None:
            return f"invalid_mol_at_pos_{i}"
        elif mol.GetProp("_Name") == "":
            return f"mol_at_pos_{i}"
        else:
            return mol.GetProp("_Name")
