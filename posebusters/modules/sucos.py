"""Module to calculate SuCOS score."""

from __future__ import annotations

import os

import numpy as np
from rdkit import RDConfig
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolChemicalFeatures import BuildFeatureFactory
from rdkit.Chem.rdmolops import AddHs, RemoveHs, SanitizeMol
from rdkit.Chem.rdShapeHelpers import ShapeProtrudeDist

FACTORY = BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))
PARAMETERS = {k: FeatMaps.FeatMapParams() for k in FACTORY.GetFeatureFamilies()}
KEEP = (
    "Donor",
    "Acceptor",
    "NegIonizable",
    "PosIonizable",
    "ZnBinder",
    "Aromatic",
    "Hydrophobe",
    "LumpedHydrophobe",
)


def get_feature_map_score(
    mol_small: Mol,
    mol_large: Mol,
    conf_id_small: int = -1,
    conf_id_large: int = -1,
) -> float:
    """Calculate the feature map score between two molecules.

    Good introduction:
        https://greglandrum.github.io/rdkit-blog/posts/2023-02-24-using-feature-maps.html
    """

    # raw features
    features_small = [f for f in FACTORY.GetFeaturesForMol(mol_small, confId=conf_id_small) if f.GetFamily() in KEEP]
    features_large = [f for f in FACTORY.GetFeaturesForMol(mol_large, confId=conf_id_large) if f.GetFamily() in KEEP]

    # feature map based on small molecule
    feature_map = FeatMaps.FeatMap(feats=features_small, weights=[1] * len(features_small), params=PARAMETERS)
    feature_map.scoreMode = FeatMaps.FeatMapScoreMode.Best

    # score feature in large molecule present in small molecule
    feature_map_score = feature_map.ScoreFeats(features_large) / min(feature_map.GetNumFeatures(), len(features_large))

    return feature_map_score


def get_sucos_score(
    mol_reference: Mol,
    mol_probe: Mol,
    conf_id_reference: int = -1,
    conf_id_probe: int = -1,
    heavy_only: bool | None = True,
) -> float:
    """Calculate the SuCOS score between a reference ligand and a list of probe ligands.

    Args:
        mol_reference: The smaller or reference molecule.
        mol_probe: The larger or query molecule.

    Returns:
        SuCOS score.

    Notes:
        SuCOS described in https://chemrxiv.org/engage/chemrxiv/article-details/60c741a99abda23230f8bed5
        Adapted from https://github.com/MarcMoesser/SuCOS/blob/master/calc_SuCOS_normalized.py
    """

    # explicit or implicit hydrogens should be same for both molecules
    mol_reference = handle_hydrogens(mol_reference, heavy_only=heavy_only)
    mol_probe = handle_hydrogens(mol_probe, heavy_only=heavy_only)

    feature_map_score = get_feature_map_score(
        mol_small=mol_reference,
        mol_large=mol_probe,
        conf_id_small=conf_id_reference,
        conf_id_large=conf_id_probe,
    )
    protrusion_distance = ShapeProtrudeDist(
        mol1=mol_reference,
        mol2=mol_probe,
        confId1=conf_id_reference,
        confId2=conf_id_probe,
        allowReordering=False,
        vdwScale=0.8,
        ignoreHs=True,
    )
    sucos_score = float(0.5 * feature_map_score + 0.5 * (1 - protrusion_distance))

    return sucos_score


def check_sucos(mol_pred: Mol, mol_true: Mol, sucos_threshold: float = 0.4) -> dict[str, dict[str, bool | float]]:
    """Calculate SuCOS and related metrics between predicted molecule and closest ground truth molecule.

    Args:
        mol_pred: Predicted molecule (docked ligand) with exactly one conformer.
        mol_true: Ground truth molecule (crystal ligand) with at least one conformer. If multiple conformers are
            present, the highest SuCOS will be reported.
        sucos_threshold: Threshold in angstrom for reporting whether SuCOS is within threshold. Defaults to 0.4.
        heavy_only: Whether to only consider heavy atoms for SuCOS calculation. Defaults to True.

    Returns:
        PoseBusters results dictionary.
    """
    assert isinstance(mol_true, Mol), "Ground truth molecule is missing."
    assert isinstance(mol_pred, Mol), "Predicted molecule is missing."
    num_conf = mol_true.GetNumConformers()
    assert num_conf > 0, "Ground truth molecule needs at least one conformer."
    assert mol_pred.GetNumConformers() == 1, "Predicted molecule should only have one conformer."

    try:
        SanitizeMol(mol_pred)
        SanitizeMol(mol_true)
    except Exception:
        return {"results": {"sucos": np.nan, "sucos_within_threshold": np.nan}}

    # iterate over all true molecules to find best sucos match
    sucos_scores = [get_sucos_score(mol_true, mol_pred, conf_id_reference=i) for i in range(num_conf)]
    best_sucos = max(sucos_scores)
    sucos_within_threshold = best_sucos >= sucos_threshold

    results = {"sucos": best_sucos, "sucos_within_threshold": sucos_within_threshold}
    return {"results": results}


def handle_hydrogens(mol: Mol, heavy_only: bool | None = True) -> Mol:
    """Remove, add or do not modify hydrogens in a molecule."""
    if heavy_only is None:
        return mol
    if heavy_only:
        return RemoveHs(mol)
    return AddHs(mol, addCoords=True)
