from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pytest
from rdkit.Chem import Conformer, MolFromPDBBlock, MolFromSmiles, SDWriter

from posebusters.cli import bust
from posebusters.diagnostics import create_diagnostic_output, diagnose_module
from posebusters.modules.distance_geometry import check_geometry
from posebusters.modules.flatness import check_flatness, nonflat
from posebusters.modules.intermolecular_distance import check_intermolecular_distance
from posebusters.posebusters import DiagnosticResultTuple, PoseBusters, ResultTuple

GEOMETRY_MODULE: dict[str, Any] = {
    "name": "Geometry",
    "function": "distance_geometry",
    "parameters": {
        "threshold_bad_bond_length": 0.25,
        "threshold_bad_angle": 0.25,
        "threshold_clash": 0.3,
    },
    "chosen_binary_test_output": [
        "bond_lengths_within_bounds",
        "bond_angles_within_bounds",
        "no_internal_clash",
    ],
    "rename_outputs": {
        "bond_lengths_within_bounds": "Bond lengths",
        "bond_angles_within_bounds": "Bond angles",
        "no_internal_clash": "Internal steric clash",
    },
}

NONFLAT_MODULE: dict[str, Any] = {
    "name": "Ring non-flatness",
    "function": "flatness",
    "parameters": {
        "flat_systems": nonflat,
        "check_nonflat": True,
        "threshold_flatness": 0.05,
    },
    "chosen_binary_test_output": ["flatness_passes"],
    "rename_outputs": {"flatness_passes": "Non-aromatic ring non-flatness"},
}

INTERMOLECULAR_MODULE: dict[str, Any] = {
    "name": "Distance to protein",
    "function": "intermolecular_distance",
    "parameters": {
        "radius_type": "vdw",
        "radius_scale": 1.0,
        "clash_cutoff": 0.75,
        "ignore_types": {"hydrogens", "organic_cofactors", "inorganic_cofactors", "waters"},
        "max_distance": 5.0,
    },
    "chosen_binary_test_output": ["not_too_far_away", "no_clashes"],
    "rename_outputs": {
        "not_too_far_away": "Protein-ligand maximum distance",
        "no_clashes": "Minimum distance to protein",
    },
}


def test_diagnose_geometry_uses_existing_failure_metrics(mol_pred_6yr2_t1c):
    out = check_geometry(
        mol_pred_6yr2_t1c,
        threshold_bad_bond_length=0.25,
        threshold_bad_angle=0.25,
        threshold_clash=0.3,
    )
    assert out["results"]["bond_lengths_within_bounds"] is False
    assert out["results"]["no_internal_clash"] is False

    diagnostics = diagnose_module(GEOMETRY_MODULE, out, {"mol_pred": mol_pred_6yr2_t1c})

    bond_diagnostics = [diagnostic for diagnostic in diagnostics if diagnostic.kind == "bond_length"]
    clash_diagnostics = [diagnostic for diagnostic in diagnostics if diagnostic.kind == "internal_clash"]

    assert len(bond_diagnostics) == 2
    assert tuple(atom.index for atom in bond_diagnostics[0].ligand_atoms) == (1, 3)
    bond_row = out["details"]["bonds"][
        out["details"]["bonds"]["atom_pair"].apply(lambda atoms: tuple(atoms) == (1, 3))
    ].iloc[0]
    assert bond_diagnostics[0].metrics["distance"] == pytest.approx(float(bond_row["distance"]))
    assert bond_diagnostics[0].metrics["upper_bound"] == pytest.approx(float(bond_row["upper_bound"]))
    assert bond_diagnostics[0].metrics["percent_error"] == pytest.approx(float(bond_row["percent_error"]))

    assert len(clash_diagnostics) == out["results"]["number_clashes"]
    assert tuple(atom.index for atom in clash_diagnostics[0].ligand_atoms) == (1, 5)
    clash_row = out["details"]["clash"][
        out["details"]["clash"]["atom_pair"].apply(lambda atoms: tuple(atoms) == (1, 5))
    ].iloc[0]
    assert clash_diagnostics[0].metrics["distance"] == pytest.approx(float(clash_row["distance"]))
    assert clash_diagnostics[0].metrics["lower_bound"] == pytest.approx(float(clash_row["lower_bound"]))

    summary_diagnostics = [diagnostic for diagnostic in diagnostics if diagnostic.kind == "summary"]

    assert not any(diagnostic.check == "Bond lengths" for diagnostic in summary_diagnostics)
    assert not any(diagnostic.check == "Internal steric clash" for diagnostic in summary_diagnostics)


def test_diagnose_angle_reports_atom_triple_and_dg_distance():
    molecule = MolFromSmiles("CCC")
    assert molecule is not None
    conformer = Conformer(molecule.GetNumAtoms())
    conformer.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conformer.SetAtomPosition(1, (1.54, 0.0, 0.0))
    conformer.SetAtomPosition(2, (0.77, 1.3337, 0.0))
    molecule.AddConformer(conformer)

    out = check_geometry(
        molecule,
        threshold_bad_bond_length=0.25,
        threshold_bad_angle=0.25,
        threshold_clash=0.3,
    )
    assert out["results"]["bond_angles_within_bounds"] is False

    diagnostics = diagnose_module(GEOMETRY_MODULE, out, {"mol_pred": molecule})
    angle_diagnostics = [diagnostic for diagnostic in diagnostics if diagnostic.kind == "bond_angle"]

    assert len(angle_diagnostics) == 1
    diagnostic = angle_diagnostics[0]
    assert tuple(atom.index for atom in diagnostic.ligand_atoms) == (0, 1, 2)
    angle_row = out["details"]["angles"][
        out["details"]["angles"]["angle"].apply(lambda atoms: tuple(atoms) == (0, 1, 2))
    ].iloc[0]
    assert diagnostic.metrics["distance_1_3"] == pytest.approx(float(angle_row["distance"]))
    assert diagnostic.metrics["lower_bound"] == pytest.approx(float(angle_row["lower_bound"]))
    assert diagnostic.metrics["absolute_percent_error"] == pytest.approx(
        float(angle_row["bound_absolute_percent_error"])
    )


def test_diagnose_flatness_reports_failed_atom_group(mol_pip_wrong):
    out = check_flatness(mol_pip_wrong, flat_systems=nonflat, check_nonflat=True, threshold_flatness=0.05)
    assert out["results"]["flatness_passes"] is False

    diagnostics = diagnose_module(NONFLAT_MODULE, out, {"mol_pred": mol_pip_wrong})

    # Several SMARTS patterns match this same ring; diagnostic output de-duplicates the atom group.
    assert len(diagnostics) == 1
    diagnostic = diagnostics[0]
    assert diagnostic.kind == "nonflatness"
    assert {atom.index for atom in diagnostic.ligand_atoms} == {0, 1, 2, 3, 4, 5}
    assert diagnostic.metrics["max_distance"] == pytest.approx(0.0027241122)
    assert diagnostic.metrics["threshold"] == 0.05


def test_diagnose_intermolecular_clash_reports_residue():
    ligand = MolFromSmiles("C")
    assert ligand is not None
    ligand_conf = Conformer(ligand.GetNumAtoms())
    ligand_conf.SetAtomPosition(0, (1.0, 0.0, 0.0))
    ligand.AddConformer(ligand_conf)

    protein = MolFromPDBBlock(
        "ATOM      1  CA  ALA A  10       0.000   0.000   0.000  1.00 20.00           C  \nEND\n",
        sanitize=False,
        removeHs=False,
        proximityBonding=False,
    )
    assert protein is not None

    parameters = INTERMOLECULAR_MODULE["parameters"]
    out = check_intermolecular_distance(ligand, protein, **parameters)
    assert out["results"]["no_clashes"] is False

    diagnostics = diagnose_module(INTERMOLECULAR_MODULE, out, {"mol_pred": ligand, "mol_cond": protein})
    clash_diagnostics = [diagnostic for diagnostic in diagnostics if diagnostic.kind == "intermolecular_clash"]

    assert len(clash_diagnostics) == 1
    diagnostic = clash_diagnostics[0]
    assert diagnostic.ligand_atoms[0].index == 0
    assert diagnostic.condition_atoms[0].index == 0
    assert diagnostic.condition_atoms[0].atom_name == "CA"
    assert diagnostic.condition_atoms[0].residue_name == "ALA"
    assert diagnostic.condition_atoms[0].residue_number == 10
    assert diagnostic.condition_atoms[0].chain_id == "A"
    assert diagnostic.metrics["distance"] == pytest.approx(1.0)
    assert diagnostic.metrics["relative_distance"] < 0.75


def test_create_diagnostic_output_contains_atom_ids_and_metrics(mol_pip_wrong):
    out = check_flatness(mol_pip_wrong, flat_systems=nonflat, check_nonflat=True, threshold_flatness=0.05)
    diagnostics = diagnose_module(NONFLAT_MODULE, out, {"mol_pred": mol_pip_wrong})

    text = create_diagnostic_output(("ligand.sdf", "PIP", 3), diagnostics)

    assert "Diagnostics for ligand.sdf | molecule: PIP | position: 3" in text
    assert "Non-aromatic ring non-flatness" in text
    assert "ligand atoms (RDKit indices)" in text
    assert "maximum distance from plane: 0.003 Å" in text
    assert "minimum required: 0.050 Å" in text


def test_cli_diagnose_multi_sdf(tmp_path, mol_pip_ideal, mol_pip_wrong):
    sdf_path = tmp_path / "poses.sdf"
    mol_pip_ideal.SetProp("_Name", "good_pose")
    mol_pip_wrong.SetProp("_Name", "bad_pose")
    with SDWriter(str(sdf_path)) as writer:
        writer.write(mol_pip_ideal)
        writer.write(mol_pip_wrong)

    output = StringIO()
    bust(
        [Path(sdf_path)],
        config=Path("posebusters/config/mol_fast.yml"),
        outfmt="diagnostic",
        output=output,
    )
    text = output.getvalue()

    assert "good_pose" in text
    assert "bad_pose" in text
    assert f"Diagnostics for {sdf_path} | molecule: bad_pose | position: 1" in text
    assert "Non-aromatic ring non-flatness" in text
    assert f"Diagnostics for {sdf_path} | molecule: good_pose | position: 0" not in text


def test_diagnose_mode_preserves_standard_results(mol_pip_wrong):
    buster = PoseBusters("mol_fast")
    buster.file_paths = pd.DataFrame([[mol_pip_wrong, None, None]], columns=["mol_pred", "mol_true", "mol_cond"])
    normal = list(buster._run())

    buster.file_paths = pd.DataFrame([[mol_pip_wrong, None, None]], columns=["mol_pred", "mol_true", "mol_cond"])
    diagnosed = list(buster._run(diagnose=True))

    assert len(normal) == len(diagnosed) == 1
    normal_key, normal_results = cast(ResultTuple, normal[0])
    diagnosed_key, diagnosed_results, diagnostics = cast(DiagnosticResultTuple, diagnosed[0])
    assert diagnosed_key == normal_key
    assert diagnosed_results == normal_results
    assert diagnostics


def test_diagnose_reports_failed_check_without_details():
    module = {
        "name": "Energy ratio",
        "function": "energy_ratio",
        "chosen_binary_test_output": ["energy_ratio_passes"],
        "rename_outputs": {
            "energy_ratio_passes": "Internal energy",
        },
    }

    module_output = {
        "results": {
            "energy_ratio_passes": False,
        }
    }

    diagnostics = diagnose_module(
        module,
        module_output,
        {},
    )

    assert len(diagnostics) == 1

    diagnostic = diagnostics[0]
    assert diagnostic.module == "Energy ratio"
    assert diagnostic.check == "Internal energy"
    assert diagnostic.kind == "summary"
    assert diagnostic.ligand_atoms == ()
    assert diagnostic.condition_atoms == ()
    assert diagnostic.metrics == {}

    text = create_diagnostic_output(
        ("ligand.sdf", "test_ligand", 0),
        diagnostics,
    )

    assert "Internal energy:" in text
    assert "failed (no additional diagnostic details)" in text


def test_diagnose_reports_summary_alongside_detailed_diagnostic(
    mol_pred_6yr2_t1c,
):
    out = check_geometry(
        mol_pred_6yr2_t1c,
        threshold_bad_bond_length=0.25,
        threshold_bad_angle=0.25,
        threshold_clash=0.3,
    )

    module = {
        **GEOMETRY_MODULE,
        "chosen_binary_test_output": [
            *GEOMETRY_MODULE["chosen_binary_test_output"],
            "other_check",
        ],
        "rename_outputs": {
            **GEOMETRY_MODULE["rename_outputs"],
            "other_check": "Other check",
        },
    }

    out["results"]["other_check"] = False

    diagnostics = diagnose_module(
        module,
        out,
        {"mol_pred": mol_pred_6yr2_t1c},
    )

    assert any(diagnostic.kind == "bond_length" for diagnostic in diagnostics)
    assert any(diagnostic.kind == "internal_clash" for diagnostic in diagnostics)

    summary = [diagnostic for diagnostic in diagnostics if diagnostic.kind == "summary"]

    assert len(summary) == 1
    assert summary[0].check == "Other check"


def test_diagnostic_output_format_includes_short_result_and_details(tmp_path, mol_pip_wrong):
    sdf_path = tmp_path / "bad_pose.sdf"
    mol_pip_wrong.SetProp("_Name", "bad_pose")
    with SDWriter(str(sdf_path)) as writer:
        writer.write(mol_pip_wrong)

    output = StringIO()
    bust(
        [Path(sdf_path)],
        config=Path("posebusters/config/mol_fast.yml"),
        outfmt="diagnostic",
        output=output,
    )
    text = output.getvalue()

    assert "bad_pose" in text
    assert "passes (" in text
    assert f"Diagnostics for {sdf_path} | molecule: bad_pose | position: 0" in text
    assert "Non-aromatic ring non-flatness" in text
