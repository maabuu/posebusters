"""MolBuster class for running all tests on a set of molecules."""
from __future__ import annotations

import inspect
import logging
import operator
from collections import defaultdict
from functools import partial, reduce
from pathlib import Path
from typing import Callable

import numpy as np
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
from .tools.formatting import _create_long_output, _create_short_output
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


class MolBuster:
    """Class to run all tests on a set of molecules."""

    def __init__(self, config: str | dict = "dock"):
        """Initialize MolBuster object."""
        self.module_func: dict[str, Callable]
        self.module_args: dict[str, set[str]]

        if isinstance(config, str) and config in {"dock", "redock", "mol"}:
            logger.info(f"Using default configuration for mode {config}.")
            self.config = safe_load(open(Path(__file__).parent / "config" / f"{config}.yml"))
        elif isinstance(config, dict):
            logger.info("Using configuration dictionary provided by user.")
            self.config = config
        else:
            logger.error(f"Configuration {config} not valid. Provide either 'dock', 'redock', 'mol' or a dictionary.")
        assert len(set(self.config.get("tests", {}).keys()) - set(module_dict.keys())) == 0

    def bust(self, mol_pred: Mol | Path, mol_true: Mol | Path | None, mol_cond: Mol | Path | None):
        """Run all tests on one molecule.

        Args:
            mol_pred: _description_
            mol_true: _description_
            mol_cond: _description_

        Returns:
            MolBuster object.
        """
        self.file_paths = pd.DataFrame([[mol_pred, mol_true, mol_cond]], columns=["mol_pred", "mol_true", "mol_cond"])
        return self._run()

    def bust_table(self, mol_table: pd.DataFrame):
        """Run all tests on multiple molecules provided in pandas dataframe as paths or rdkit molecule objects.

        Args:
            mol_table: _description_

        Returns:
            MolBuster object.
        """
        self.file_paths = mol_table
        return self._run()

    def _run(self):
        """Run all tests on molecules provided in file paths.

        Yields:
            Generator of pandas dataframes with results.
        """
        self._initialize_modules()

        output = defaultdict(dict)
        for i, paths in self.file_paths.iterrows():
            mol_args = {}
            if "mol_cond" in paths and paths["mol_cond"] is not None:
                mol_cond_load_params = self.config.get("loading_options", {}).get("mol_cond", {})
                mol_args["mol_cond"] = safe_load_mol(path=paths["mol_cond"], **mol_cond_load_params)
            if "mol_true" in paths and paths["mol_true"] is not None:
                mol_true_load_params = self.config.get("loading_options", {}).get("mol_true", {})
                mol_args["mol_true"] = safe_load_mol(path=paths["mol_true"], **mol_true_load_params)

            mol_pred_load_params = self.config.get("loading_options", {}).get("mol_pred", {})
            for i, mol_pred in enumerate(safe_supply_mols(paths["mol_pred"], **mol_pred_load_params)):
                mol_args["mol_pred"] = mol_pred

                results_key = (str(paths["mol_pred"]), self._get_name(mol_pred, i))

                for module_name in self.config.get("tests", {}).keys():
                    # pick needed arguments for module
                    module_args = self.module_args[module_name]
                    args = {k: v for k, v in mol_args.items() if k in module_args}
                    # loading takes all inputs
                    if module_name == "loading":
                        args = {k: args.get(k, None) for k in module_args}
                    # check that no arguments are missing
                    if module_name != "loading" and not all(args.get(m, None) for m in module_args):
                        output[results_key][module_name] = {"results": {}}
                        continue
                    # run module
                    module_results = self.module_func[module_name](**args)

                    # add to all results
                    output[results_key][module_name] = module_results

                # return partial results to print progress
                yield self._get_partial_output(output, results_key)

    def _initialize_modules(self):
        self.module_func = {}
        self.module_args = {}
        for module_name in self.config.get("tests").keys():
            module_config = self.config.get("parameters", {}).get(module_name, {})
            self.module_func[module_name] = partial(module_dict[module_name], **module_config)
            module_args = set(inspect.signature(self.module_func[module_name]).parameters)
            self.module_args[module_name] = module_args.intersection({"mol_pred", "mol_true", "mol_cond"})

    def _apply_column_names_from_config(self, df):
        module_names = self.config["module_names"]
        test_names = self.config["test_names"]
        cols = df.columns
        cols = [(module_names.get(module, module), test_names.get(module, {}).get(test, test)) for module, test in cols]
        df.columns = pd.MultiIndex.from_tuples(cols)
        return df

    def _get_partial_output(self, output, results_key) -> pd.DataFrame:
        results = self._unfold_results_entry({results_key: output[results_key]}, "results")
        selected_columns = [(m, c) for m in self.config["tests"].keys() for c in self.config["tests"][m]]
        missing_columns = [c for c in selected_columns if c not in results]
        results[missing_columns] = np.nan
        results = self._apply_column_names_from_config(results[selected_columns])
        return results

    @staticmethod
    def _unfold_results_entry(d, field):
        # pick field
        ddd = {k: {(kk, kkk): vvv for kk, vv in v.items() for kkk, vvv in vv[field].items()} for k, v in d.items()}
        # make df
        return pd.DataFrame.from_dict(ddd, orient="index")

    @staticmethod
    def _get_name(mol: Mol, i: int) -> str:
        if mol is None:
            return f"invalid_mol_at_pos_{i}"
        elif mol.GetProp("_Name") == "":
            return f"mol_at_pos_{i}"
        else:
            return mol.GetProp("_Name")
