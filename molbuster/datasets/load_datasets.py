"""Load datasets from the datasets folder."""

from pathlib import Path

import pandas as pd

dataset_dir = Path(__file__).parent


def _prepend_path(file_path):
    return dataset_dir / file_path


def _load_dataset(cols):
    files = pd.read_csv(dataset_dir / "pdb.csv")
    files[cols] = files[cols].apply(_prepend_path)
    return files[["name"] + cols]


def four_redocks():
    """Load the four docks dataset."""
    return _load_dataset(["mol_pred", "mol_true", "mol_cond"])


def four_docks():
    """Load the four docks dataset."""
    return _load_dataset(["mol_pred", "mol_cond"])


def four_mols():
    """Load the four mols dataset."""
    return _load_dataset(["mol_pred"])
