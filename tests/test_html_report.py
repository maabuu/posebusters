from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit.Chem import Conformer, MolFromPDBBlock, MolFromSmiles, SDWriter

from posebusters.cli import _parse_args, bust
from posebusters.diagnostics import AtomRef, Diagnostic, diagnose_module
from posebusters.html_report import (
    HtmlPoseResult,
    _atom_label_position,
    _measurement_label_position,
    _trim_line_to_circle_edges,
    create_html_report,
)
from posebusters.modules.flatness import check_flatness, nonflat
from posebusters.modules.intermolecular_distance import check_intermolecular_distance

NONFLAT_MODULE = {
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

INTERMOLECULAR_MODULE = {
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


def test_trim_line_to_circle_edges():
    point_1 = np.array([10.0, 20.0])
    point_2 = np.array([50.0, 20.0])

    line_1, line_2 = _trim_line_to_circle_edges(
        point_1,
        point_2,
        radius_1=5.0,
        radius_2=7.0,
        padding=1.0,
    )

    assert line_1 == pytest.approx(np.array([16.0, 20.0]))
    assert line_2 == pytest.approx(np.array([42.0, 20.0]))


def test_trim_line_to_circle_edges_keeps_zero_length_contact():
    point = np.array([10.0, 20.0])

    line_1, line_2 = _trim_line_to_circle_edges(
        point,
        point,
        radius_1=5.3,
        radius_2=5.3,
    )

    assert line_1 == pytest.approx(point)
    assert line_2 == pytest.approx(point)


def test_measurement_label_prefers_perpendicular_offset():
    point_1 = np.array([100.0, 100.0])
    point_2 = np.array([200.0, 100.0])
    occupied: list[np.ndarray] = []

    label_point = _measurement_label_position(
        point_1,
        point_2,
        avoid_points=None,
        occupied_labels=occupied,
        preferred_sign=1,
    )

    assert label_point[0] == pytest.approx(150.0)
    assert label_point[1] == pytest.approx(120.0)


def test_measurement_label_stays_close_to_midpoint_when_clear():
    point_1 = np.array([100.0, 100.0])
    point_2 = np.array([200.0, 100.0])
    occupied: list[np.ndarray] = []

    label_point = _measurement_label_position(
        point_1,
        point_2,
        avoid_points=None,
        occupied_labels=occupied,
        preferred_sign=1,
    )

    midpoint = (point_1 + point_2) / 2
    assert np.linalg.norm(label_point - midpoint) == pytest.approx(20.0)


def test_atom_label_is_placed_away_from_bonded_neighbours():
    point = np.array([100.0, 100.0])
    neighbour_points = np.array([[80.0, 100.0]])
    occupied: list[np.ndarray] = []

    label_point = _atom_label_position(
        point,
        neighbour_points,
        avoid_points=np.vstack([point, neighbour_points]),
        occupied_labels=occupied,
    )

    assert label_point[0] > point[0]


def test_html_report_contains_summary_and_flatness_visual(mol_pip_wrong):
    out = check_flatness(mol_pip_wrong, flat_systems=nonflat, check_nonflat=True, threshold_flatness=0.05)
    diagnostics = diagnose_module(NONFLAT_MODULE, out, {"mol_pred": mol_pip_wrong})
    results = pd.DataFrame([[False]], columns=pd.MultiIndex.from_tuples([("Ring non-flatness", "Flatness")]))

    html = create_html_report(
        [
            HtmlPoseResult(
                key=("ligand.sdf", "PIP", 0),
                results=results,
                diagnostics=tuple(diagnostics),
                mol_pred=mol_pip_wrong,
                mol_cond=None,
            )
        ]
    )

    assert "PoseBusters diagnostic report" in html
    assert "<strong>1</strong><span>Failing</span>" in html
    assert "Non-aromatic ring non-flatness" in html
    assert "3D geometry projected to 2D" in html
    assert "max plane distance: 0.003 A" in html
    assert "<svg" in html
    assert "Ligand overview" not in html


def test_html_report_all_passes_has_no_detail_sections(mol_pip_ideal):
    results = pd.DataFrame(
        [[True]],
        columns=pd.MultiIndex.from_tuples([("Ring non-flatness", "Flatness")]),
    )

    html = create_html_report(
        [
            HtmlPoseResult(
                key=("ligand.sdf", "passing_pose", 0),
                results=results,
                diagnostics=(),
                mol_pred=mol_pip_ideal,
                mol_cond=None,
            )
        ]
    )

    assert "<strong>1</strong><span>Passing</span>" in html
    assert "<strong>0</strong><span>Failing</span>" in html
    assert "All poses passed the selected PoseBusters checks." in html
    assert "<section id='pose-" not in html


def test_html_projection_colours_atoms_and_labels_only_highlights():
    ligand = MolFromSmiles("CON")
    assert ligand is not None
    conf = Conformer(ligand.GetNumAtoms())
    conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, (1.0, 0.0, 0.0))
    conf.SetAtomPosition(2, (2.0, 0.5, 0.0))
    ligand.AddConformer(conf)

    diagnostic = Diagnostic(
        module="Geometry",
        check="Internal steric clash",
        kind="internal_clash",
        ligand_atoms=(AtomRef(0, "C"), AtomRef(1, "O")),
        metrics={
            "distance": 1.0,
            "lower_bound": 2.0,
            "percent_error": -0.5,
            "threshold": 0.3,
        },
    )
    results = pd.DataFrame(
        [[False]],
        columns=pd.MultiIndex.from_tuples([("Geometry", "Internal steric clash")]),
    )

    html = create_html_report(
        [
            HtmlPoseResult(
                key=("ligand.sdf", "hetero", 0),
                results=results,
                diagnostics=(diagnostic,),
                mol_pred=ligand,
                mol_cond=None,
            )
        ]
    )

    assert "data-element='C'" in html
    assert "data-element='O'" in html
    assert "data-element='N'" in html
    assert "style='fill:#FF0000'" in html
    assert "style='fill:#0000FF'" in html
    assert ">C0<" in html
    assert ">O1<" in html
    assert ">N2<" not in html
    assert "Atom colours:" in html
    assert "involved in failure" in html
    assert "measure-bg" in html
    assert "rdkit-measurement-label" in html
    assert "class='measure-label'" in html
    assert "class='measurement'" in html
    assert "stroke-dasharray: 6 4" in html


def test_html_report_projects_protein_contact():
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

    out = check_intermolecular_distance(ligand, protein, **INTERMOLECULAR_MODULE["parameters"])
    diagnostics = diagnose_module(INTERMOLECULAR_MODULE, out, {"mol_pred": ligand, "mol_cond": protein})
    results = pd.DataFrame([[False]], columns=pd.MultiIndex.from_tuples([("Distance", "Minimum distance")]))

    html = create_html_report(
        [
            HtmlPoseResult(
                key=("ligand.sdf", "clash", 0),
                results=results,
                diagnostics=tuple(diagnostics),
                mol_pred=ligand,
                mol_cond=protein,
            )
        ]
    )

    assert "Minimum distance to protein" in html
    assert "ALA A:10 CA" in html
    assert "1.000 A" in html
    assert "3D geometry projected to 2D" in html
    assert "measurement" in html
    assert "measure-bg" in html
    assert "condition-callout" in html
    assert "label-leader" in html


def test_html_report_shows_one_overview_for_multiple_visual_failures(mol_pip_wrong):
    out = check_flatness(mol_pip_wrong, flat_systems=nonflat, check_nonflat=True, threshold_flatness=0.05)
    flatness = diagnose_module(NONFLAT_MODULE, out, {"mol_pred": mol_pip_wrong})[0]
    clash = Diagnostic(
        module="Geometry",
        check="Internal steric clash",
        kind="internal_clash",
        ligand_atoms=(AtomRef(0, "C"), AtomRef(1, "C")),
        metrics={
            "distance": 1.0,
            "lower_bound": 2.0,
            "percent_error": -0.5,
            "threshold": 0.3,
        },
    )
    results = pd.DataFrame(
        [[False, False]],
        columns=pd.MultiIndex.from_tuples(
            [
                ("Ring non-flatness", "Flatness"),
                ("Geometry", "Internal steric clash"),
            ]
        ),
    )

    html = create_html_report(
        [
            HtmlPoseResult(
                key=("ligand.sdf", "multi_failure", 0),
                results=results,
                diagnostics=(flatness, clash),
                mol_pred=mol_pip_wrong,
                mol_cond=None,
            )
        ]
    )

    assert html.count("Ligand overview (all failed checks)") == 1
    assert "Non-aromatic ring non-flatness" in html
    assert "Internal steric clash" in html


def test_cli_html_multi_sdf_is_single_report(tmp_path, mol_pip_ideal, mol_pip_wrong):
    sdf_path = tmp_path / "poses.sdf"
    report_path = tmp_path / "report.html"
    mol_pip_ideal.SetProp("_Name", "good_pose")
    mol_pip_wrong.SetProp("_Name", "bad_pose")
    with SDWriter(str(sdf_path)) as writer:
        writer.write(mol_pip_ideal)
        writer.write(mol_pip_wrong)

    bust(
        [Path(sdf_path)],
        config=Path("posebusters/config/mol_fast.yml"),
        outfmt="html",
        output=report_path,
    )
    html = report_path.read_text(encoding="utf-8")

    assert html.count("<!doctype html>") == 1
    assert "<strong>2</strong><span>Poses</span>" in html
    assert "<strong>1</strong><span>Passing</span>" in html
    assert "<strong>1</strong><span>Failing</span>" in html
    assert "good_pose" in html
    assert "bad_pose" in html
    assert html.count("<section id='pose-") == 1
    assert "Non-aromatic ring non-flatness" in html


def test_parse_args_html_requires_output():
    with pytest.raises(SystemExit):
        _parse_args(["tests/conftest/mol_PIP_wrong.sdf", "--outfmt", "html"])

    args = _parse_args(
        [
            "tests/conftest/mol_PIP_wrong.sdf",
            "--outfmt",
            "html",
            "--output",
            "report.html",
        ]
    )
    assert args.outfmt == "html"
    assert args.output == Path("report.html")
