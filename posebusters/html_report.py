"""Self-contained HTML reporting for PoseBusters diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from math import isfinite
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdchem import Atom, Mol

from .diagnostics import AtomRef, Diagnostic

_ATOM_PALETTE = rdMolDraw2D.MolDrawOptions().getAtomPalette()
_HIGHLIGHT_COLOUR = "#D95F02"
_MIN_VECTOR_NORM = 1e-6
_MIN_RING_ATOMS = 3
_MIN_DISTANCE_PAIR_ATOMS = 2
_MIN_ANGLE_ATOMS = 3


@dataclass(frozen=True)
class HtmlPoseResult:
    """Results and molecular context for one pose in an HTML report."""

    key: tuple[str, str, int]
    results: pd.DataFrame
    diagnostics: tuple[Diagnostic, ...]
    mol_pred: Mol | None
    mol_cond: Mol | None


def create_html_report(poses: list[HtmlPoseResult]) -> str:
    """Create a self-contained HTML diagnostic report for one or more poses."""
    summaries = [_summarise_pose(pose) for pose in poses]
    num_pass = sum(summary[0] for summary in summaries)
    num_fail = len(poses) - num_pass

    body = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'>",
        "<title>PoseBusters diagnostic report</title>",
        f"<style>{_CSS}</style></head><body>",
        "<h1>PoseBusters diagnostic report</h1>",
        "<p class='notice'><strong>Numbering:</strong> ligand and conditioning "
        "atom IDs are RDKit indices (0-based).</p>",
        "<div class='stats'>",
        _stat_card("Poses", len(poses)),
        _stat_card("Passing", num_pass),
        _stat_card("Failing", num_fail),
        "</div>",
        "<h2>Summary</h2>",
        _summary_table(poses, summaries),
    ]

    failing = [(pose, summary) for pose, summary in zip(poses, summaries) if not summary[0]]
    if failing:
        body.append("<h2>Failed poses</h2>")
        for index, (pose, summary) in enumerate(zip(poses, summaries)):
            if not summary[0]:
                body.append(_pose_section(index, pose, summary))
    else:
        body.append("<p class='pass-box'>All poses passed the selected PoseBusters checks.</p>")

    body.append("</body></html>")
    return "\n".join(body)


def _summarise_pose(pose: HtmlPoseResult) -> tuple[bool, int, int, tuple[str, ...]]:
    passes = int(pose.results.sum(axis=1).iloc[0])
    total = len(pose.results.columns)
    passed = bool(pose.results.all(axis=1).iloc[0])
    failed_checks = tuple(dict.fromkeys(diagnostic.check for diagnostic in pose.diagnostics))
    return passed, passes, total, failed_checks


def _summary_table(poses: list[HtmlPoseResult], summaries: list[tuple[bool, int, int, tuple[str, ...]]]) -> str:
    rows = []
    for index, (pose, summary) in enumerate(zip(poses, summaries)):
        passed, passes, total, failed_checks = summary
        file_name, molecule_name, position = pose.key
        label = molecule_name or f"pose {position}"
        href = "" if passed else f"<a href='#pose-{index}'>{escape(label)}</a>"
        molecule = escape(label) if passed else href
        status_class = "pass" if passed else "fail"
        status = "PASS" if passed else "FAIL"
        failed = ", ".join(escape(check) for check in failed_checks) or "-"
        rows.append(
            "<tr>"
            f"<td>{escape(file_name)}</td>"
            f"<td>{molecule}</td>"
            f"<td>{position}</td>"
            f"<td><span class='status {status_class}'>{status}</span></td>"
            f"<td>{passes} / {total}</td>"
            f"<td>{failed}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr><th>File</th><th>Molecule</th><th>Position</th><th>Status</th>"
        "<th>Passed</th><th>Failed checks</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
    )


def _pose_section(index: int, pose: HtmlPoseResult, summary: tuple[bool, int, int, tuple[str, ...]]) -> str:
    _, passes, total, _ = summary
    file_name, molecule_name, position = pose.key
    label = molecule_name or f"pose {position}"
    highlighted = sorted({atom.index for diagnostic in pose.diagnostics for atom in diagnostic.ligand_atoms})

    parts = [
        f"<section id='pose-{index}'>",
        f"<h3>{escape(label)}</h3>",
        f"<p class='meta'>{escape(file_name)} | position {position} | passes {passes} / {total}</p>",
    ]

    grouped_diagnostics = _group_diagnostics(pose.diagnostics)
    visual_groups = [
        diagnostics
        for _, diagnostics in grouped_diagnostics
        if any(diagnostic.ligand_atoms and diagnostic.kind != "summary" for diagnostic in diagnostics)
    ]
    if pose.mol_pred is not None and len(visual_groups) > 1:
        parts.extend(
            [
                "<div class='overview'>",
                "<div><h4>Ligand overview (all failed checks)</h4>",
                _ligand_svg(pose.mol_pred, highlighted),
                "</div></div>",
            ]
        )

    for check, diagnostics in grouped_diagnostics:
        parts.append(f"<div class='diagnostic'><h4>{escape(check)}</h4>")
        detailed = [diagnostic for diagnostic in diagnostics if diagnostic.kind != "summary"]
        if not detailed:
            parts.append("<p>Failed (no additional diagnostic details).</p></div>")
            continue

        visual = _diagnostic_visual(detailed, pose.mol_pred, pose.mol_cond)
        if visual:
            parts.append("<div class='grid'>")
            if pose.mol_pred is not None:
                parts.append("<div><h5>2D structure</h5>")
                parts.append(_ligand_svg(pose.mol_pred, _ligand_indices(detailed), detailed))
                parts.append("</div>")
            parts.append("<div><h5>3D geometry projected to 2D</h5>")
            parts.append(visual)
            parts.append("</div></div>")

        parts.append(_diagnostic_table(detailed))
        parts.append("</div>")

    parts.append("</section>")
    return "".join(parts)


def _group_diagnostics(diagnostics: tuple[Diagnostic, ...]) -> list[tuple[str, list[Diagnostic]]]:
    groups: dict[str, list[Diagnostic]] = {}
    for diagnostic in diagnostics:
        groups.setdefault(diagnostic.check, []).append(diagnostic)
    return list(groups.items())


def _ligand_indices(diagnostics: list[Diagnostic]) -> list[int]:
    return sorted({atom.index for diagnostic in diagnostics for atom in diagnostic.ligand_atoms})


def _ligand_svg(mol: Mol, highlighted: list[int], diagnostics: list[Diagnostic] | None = None) -> str:
    draw_mol = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(draw_mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(520, 340)
    draw_options = drawer.drawOptions()
    draw_options.addAtomIndices = True
    # Keep atom indices readable for larger ligands without fixing the overall
    # depiction scale. RDKit otherwise makes annotations quite small as the
    # molecule occupies more of the canvas.
    draw_options.annotationFontScale = 0.75
    draw_options.minFontSize = 8
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, draw_mol, highlightAtoms=highlighted)

    measurement_overlays: list[str] = []
    occupied_labels: list[np.ndarray] = []
    draw_points = np.array(
        [[drawer.GetDrawCoords(index).x, drawer.GetDrawCoords(index).y] for index in range(draw_mol.GetNumAtoms())],
        dtype=float,
    )
    if diagnostics:
        for measurement_index, diagnostic in enumerate(diagnostics):
            pair = _distance_pair(diagnostic)
            distance = _reported_distance(diagnostic)
            if pair is None or distance is None:
                continue
            atom_1, atom_2 = pair
            if atom_1 >= draw_mol.GetNumAtoms() or atom_2 >= draw_mol.GetNumAtoms():
                continue
            point_1 = drawer.GetDrawCoords(atom_1)
            point_2 = drawer.GetDrawCoords(atom_2)
            point_1_array = np.array([point_1.x, point_1.y])
            point_2_array = np.array([point_2.x, point_2.y])
            # Overlay measurements after RDKit has drawn the molecule so that
            # the line can use the same dashed styling as projected contacts.
            line_1, line_2 = _trim_line_to_circle_edges(
                point_1_array,
                point_2_array,
                radius_1=8.0,
                radius_2=8.0,
            )
            measurement_overlays.append(_svg_line(line_1, line_2, "measurement"))
            measurement_overlays.append(
                _svg_measurement_label(
                    point_1_array,
                    point_2_array,
                    f"{distance:.2f} A",
                    avoid_points=draw_points,
                    occupied_labels=occupied_labels,
                    preferred_sign=1 if measurement_index % 2 == 0 else -1,
                    width=520,
                    height=340,
                )
            )

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("<?xml version='1.0' encoding='iso-8859-1'?>", "")
    if measurement_overlays:
        svg = svg.replace("</svg>", "".join(measurement_overlays) + "</svg>")
    return f"<div class='svg-wrap'>{svg}</div>"


def _diagnostic_visual(  # noqa: PLR0911 -- dispatch is clearer as explicit early returns
    diagnostics: list[Diagnostic], mol_pred: Mol | None, mol_cond: Mol | None
) -> str:
    if mol_pred is None or mol_pred.GetNumConformers() == 0:
        return ""

    kind = diagnostics[0].kind
    if kind in {"flatness", "nonflatness"}:
        return _flatness_projection(mol_pred, diagnostics)
    if kind in {"bond_length", "bond_angle", "internal_clash"}:
        return _ligand_projection(mol_pred, diagnostics)
    if kind in {"intermolecular_clash", "intermolecular_max_distance"} and mol_cond is not None:
        if not any(diagnostic.condition_atoms for diagnostic in diagnostics):
            return ""
        if mol_cond.GetNumConformers() == 0:
            return ""
        return _contact_projection(mol_pred, mol_cond, diagnostics)
    return ""


def _ligand_projection(mol: Mol, diagnostics: list[Diagnostic]) -> str:
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())], dtype=float)
    xy = _project(coords)
    highlights = set(_ligand_indices(diagnostics))
    lines = []
    for bond in mol.GetBonds():
        lines.append(_svg_line(xy[bond.GetBeginAtomIdx()], xy[bond.GetEndAtomIdx()], "bond"))

    occupied_labels: list[np.ndarray] = []
    for index in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(index)
        label_point = None
        if index in highlights:
            neighbour_points = np.array([xy[neighbour.GetIdx()] for neighbour in atom.GetNeighbors()], dtype=float)
            label_point = _atom_label_position(xy[index], neighbour_points, xy, occupied_labels)
        lines.append(_svg_atom(xy[index], atom, index, index in highlights, label_point))

    for measurement_index, diagnostic in enumerate(diagnostics):
        pair = _distance_pair(diagnostic)
        distance = _reported_distance(diagnostic)
        if pair is not None and distance is not None:
            lines.append(
                _svg_measurement(
                    xy[pair[0]],
                    xy[pair[1]],
                    f"{distance:.2f} A",
                    avoid_points=xy,
                    occupied_labels=occupied_labels,
                    preferred_sign=1 if measurement_index % 2 == 0 else -1,
                )
            )
    return _svg_canvas(lines, _visible_atomic_numbers(mol))


def _flatness_projection(mol: Mol, diagnostics: list[Diagnostic]) -> str:
    atom_ids = _ligand_indices(diagnostics)
    if len(atom_ids) < _MIN_RING_ATOMS:
        return ""
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(index)) for index in atom_ids], dtype=float)
    centered = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    in_plane = vh[0]
    normal = vh[-1]
    projected = np.column_stack((centered @ in_plane, centered @ normal))
    xy = _scale_projection(projected)
    lookup = {atom_id: position for atom_id, position in zip(atom_ids, xy)}
    lines = []
    atom_set = set(atom_ids)
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in atom_set and end in atom_set:
            lines.append(_svg_line(lookup[begin], lookup[end], "bond"))
    occupied_labels: list[np.ndarray] = []
    avoid_points = np.array([lookup[atom_id] for atom_id in atom_ids], dtype=float)
    for atom_id in atom_ids:
        atom = mol.GetAtomWithIdx(atom_id)
        neighbour_points = np.array(
            [lookup[neighbour.GetIdx()] for neighbour in atom.GetNeighbors() if neighbour.GetIdx() in atom_set],
            dtype=float,
        )
        label_point = _atom_label_position(lookup[atom_id], neighbour_points, avoid_points, occupied_labels)
        lines.append(_svg_atom(lookup[atom_id], atom, atom_id, True, label_point))
    max_distance = diagnostics[0].metrics.get("max_distance")
    if isinstance(max_distance, (int, float)) and isfinite(float(max_distance)):
        lines.append(f"<text x='12' y='22' class='annotation'>max plane distance: {float(max_distance):.3f} A</text>")
    return _svg_canvas(lines, _visible_atomic_numbers(mol, atom_ids))


def _contact_projection(mol_pred: Mol, mol_cond: Mol, diagnostics: list[Diagnostic]) -> str:
    pred_conf = mol_pred.GetConformer()
    cond_conf = mol_cond.GetConformer()
    ligand_coords = np.array([list(pred_conf.GetAtomPosition(i)) for i in range(mol_pred.GetNumAtoms())], dtype=float)

    condition_ids = sorted({atom.index for diagnostic in diagnostics for atom in diagnostic.condition_atoms})
    residue_ids = _condition_context_atom_ids(mol_cond, condition_ids)
    condition_coords = np.array([list(cond_conf.GetAtomPosition(i)) for i in residue_ids], dtype=float)
    combined = ligand_coords if not residue_ids else np.vstack([ligand_coords, condition_coords])
    # Reserve a little more space around contact views for residue labels.
    xy = _project(combined, margin=60)
    ligand_xy = xy[: mol_pred.GetNumAtoms()]
    condition_xy = xy[mol_pred.GetNumAtoms() :]
    condition_lookup = {atom_id: position for atom_id, position in zip(residue_ids, condition_xy)}

    highlighted_ligand = set(_ligand_indices(diagnostics))
    highlighted_condition = set(condition_ids)
    lines = []
    for bond in mol_pred.GetBonds():
        lines.append(_svg_line(ligand_xy[bond.GetBeginAtomIdx()], ligand_xy[bond.GetEndAtomIdx()], "bond"))

    avoid_points = np.vstack([ligand_xy, condition_xy]) if len(condition_xy) else ligand_xy
    occupied_labels: list[np.ndarray] = []
    for index in range(mol_pred.GetNumAtoms()):
        atom = mol_pred.GetAtomWithIdx(index)
        label_point = None
        if index in highlighted_ligand:
            neighbour_points = np.array(
                [ligand_xy[neighbour.GetIdx()] for neighbour in atom.GetNeighbors()],
                dtype=float,
            )
            label_point = _atom_label_position(
                ligand_xy[index],
                neighbour_points,
                avoid_points,
                occupied_labels,
            )
        lines.append(_svg_atom(ligand_xy[index], atom, index, index in highlighted_ligand, label_point))

    for atom_id in residue_ids:
        atom = mol_cond.GetAtomWithIdx(atom_id)
        lines.append(_svg_condition_atom(condition_lookup[atom_id], atom, atom_id in highlighted_condition))

    label_targets: dict[int, list[np.ndarray]] = {}
    measurement_index = 0
    for diagnostic in diagnostics:
        if not diagnostic.ligand_atoms or not diagnostic.condition_atoms:
            continue
        ligand_id = diagnostic.ligand_atoms[0].index
        condition_id = diagnostic.condition_atoms[0].index
        if condition_id not in condition_lookup:
            continue
        distance = _reported_distance(diagnostic)
        label = f"{distance:.2f} A" if distance is not None else "contact"
        lines.append(
            _svg_measurement(
                ligand_xy[ligand_id],
                condition_lookup[condition_id],
                label,
                avoid_points=avoid_points,
                occupied_labels=occupied_labels,
                preferred_sign=1 if measurement_index % 2 == 0 else -1,
            )
        )
        label_targets.setdefault(condition_id, []).append(ligand_xy[ligand_id])
        measurement_index += 1

    # Draw residue labels last so that they remain legible above bonds and
    # distance markers. Place them away from the ligand contact where possible.
    for label_index, atom_id in enumerate(condition_ids):
        if atom_id not in condition_lookup:
            continue
        atom = mol_cond.GetAtomWithIdx(atom_id)
        label = _condition_label(atom_id, atom.GetSymbol(), atom.GetPDBResidueInfo())
        targets = label_targets.get(atom_id, [])
        target = np.mean(targets, axis=0) if targets else None
        lines.append(
            _svg_condition_label(
                condition_lookup[atom_id],
                label,
                target,
                label_index,
                avoid_points=avoid_points,
                occupied_labels=occupied_labels,
            )
        )
    atomic_numbers = _visible_atomic_numbers(mol_pred) | _visible_atomic_numbers(mol_cond, residue_ids)
    return _svg_canvas(lines, atomic_numbers)


def _condition_context_atom_ids(mol: Mol, atom_ids: list[int]) -> list[int]:
    selected = set(atom_ids)
    residue_keys: set[tuple[str, int, str, str]] = set()
    for atom_id in atom_ids:
        info = mol.GetAtomWithIdx(atom_id).GetPDBResidueInfo()
        if info is not None:
            residue_keys.add(
                (
                    info.GetChainId().strip(),
                    info.GetResidueNumber(),
                    info.GetInsertionCode().strip(),
                    info.GetResidueName().strip(),
                )
            )
    if not residue_keys:
        return sorted(selected)

    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        key = (
            info.GetChainId().strip(),
            info.GetResidueNumber(),
            info.GetInsertionCode().strip(),
            info.GetResidueName().strip(),
        )
        if key in residue_keys:
            selected.add(atom.GetIdx())
    return sorted(selected)


def _project(coords: np.ndarray, margin: float = 35) -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    if len(coords) == 1:
        projected = np.zeros((1, 2), dtype=float)
    else:
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        axes = vh[:2]
        if axes.shape[0] == 1:
            axes = np.vstack([axes, np.array([0.0, 1.0, 0.0])])
        projected = centered @ axes.T
    return _scale_projection(projected, margin=margin)


def _scale_projection(projected: np.ndarray, width: float = 520, height: float = 340, margin: float = 35) -> np.ndarray:
    if len(projected) == 0:
        return projected
    mins = projected.min(axis=0)
    maxs = projected.max(axis=0)
    spans = np.maximum(maxs - mins, 1e-6)
    scale = min((width - 2 * margin) / spans[0], (height - 2 * margin) / spans[1])
    scaled = (projected - (mins + maxs) / 2) * scale
    scaled[:, 0] += width / 2
    scaled[:, 1] = height / 2 - scaled[:, 1]
    return scaled


def _visible_atomic_numbers(mol: Mol, atom_ids: list[int] | None = None) -> set[int]:
    ids = range(mol.GetNumAtoms()) if atom_ids is None else atom_ids
    return {mol.GetAtomWithIdx(index).GetAtomicNum() for index in ids}


def _rdkit_atom_colour(atomic_number: int) -> str:
    colour = _ATOM_PALETTE.get(atomic_number, _ATOM_PALETTE[-1])
    red, green, blue = (round(255 * channel) for channel in colour[:3])
    return f"#{red:02X}{green:02X}{blue:02X}"


def _svg_canvas(contents: list[str], atomic_numbers: set[int] | None = None) -> str:
    legend = _atom_legend(atomic_numbers or set())
    return (
        "<div class='svg-wrap'><svg viewBox='0 0 520 340' xmlns='http://www.w3.org/2000/svg'>"
        + "".join(contents)
        + "</svg></div>"
        + legend
    )


def _atom_legend(atomic_numbers: set[int]) -> str:
    if not atomic_numbers:
        return ""
    periodic_table = Chem.GetPeriodicTable()
    entries = []
    for atomic_number in sorted(atomic_numbers):
        symbol = periodic_table.GetElementSymbol(atomic_number)
        colour = _rdkit_atom_colour(atomic_number)
        entries.append(
            "<span class='legend-item'>"
            f"<span class='legend-dot' style='background:{colour}'></span>{escape(symbol)}"
            "</span>"
        )
    entries.append(
        "<span class='legend-item'>"
        f"<span class='legend-dot problem-dot' style='background:white;border-color:{_HIGHLIGHT_COLOUR}'></span>"
        "involved in failure</span>"
    )
    return "<div class='atom-legend'><span>Atom colours:</span>" + "".join(entries) + "</div>"


def _trim_line_to_circle_edges(
    point_1: np.ndarray,
    point_2: np.ndarray,
    radius_1: float,
    radius_2: float,
    padding: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Trim a line so it meets circular atom markers at their edges."""
    vector = point_2 - point_1
    distance = float(np.linalg.norm(vector))
    trim_1 = radius_1 + padding
    trim_2 = radius_2 + padding
    if distance < _MIN_VECTOR_NORM or distance <= trim_1 + trim_2:
        return point_1, point_2

    direction = vector / distance
    return point_1 + direction * trim_1, point_2 - direction * trim_2


def _svg_line(point_1: np.ndarray, point_2: np.ndarray, css_class: str) -> str:
    return (
        f"<line x1='{point_1[0]:.1f}' y1='{point_1[1]:.1f}' x2='{point_2[0]:.1f}' y2='{point_2[1]:.1f}' "
        f"class='{css_class}'/>"
    )


def _svg_atom(
    point: np.ndarray,
    atom: Atom,
    index: int,
    highlighted: bool,
    label_point: np.ndarray | None = None,
) -> str:
    css_class = "atom highlight" if highlighted else "atom"
    radius = 5.3 if highlighted else 4.5
    colour = _rdkit_atom_colour(atom.GetAtomicNum())
    label = f"{atom.GetSymbol()}{index}" if highlighted else ""
    if label and label_point is None:
        label_point = point + np.array([7.0, -7.0])
    return (
        f"<circle cx='{point[0]:.1f}' cy='{point[1]:.1f}' r='{radius:.1f}' "
        f"class='{css_class}' style='fill:{colour}' data-element='{escape(atom.GetSymbol())}'/>"
        + (
            f"<text x='{label_point[0]:.1f}' y='{label_point[1]:.1f}' text-anchor='middle' "
            f"class='atom-label'>{escape(label)}</text>"
            if label and label_point is not None
            else ""
        )
    )


def _svg_condition_atom(point: np.ndarray, atom: Atom, highlighted: bool) -> str:
    css_class = "condition highlight" if highlighted else "condition"
    radius = 5.3 if highlighted else 4.5
    colour = _rdkit_atom_colour(atom.GetAtomicNum())
    return (
        f"<circle cx='{point[0]:.1f}' cy='{point[1]:.1f}' r='{radius:.1f}' "
        f"class='{css_class}' style='fill:{colour}' data-element='{escape(atom.GetSymbol())}'/>"
    )


def _unit_vector(vector: np.ndarray, fallback: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    return fallback if norm < _MIN_VECTOR_NORM else vector / norm


def _best_label_position(  # noqa: PLR0913 -- SVG layout inputs are intentionally explicit
    anchor: np.ndarray,
    preferred_direction: np.ndarray,
    avoid_points: np.ndarray | None,
    occupied_labels: list[np.ndarray],
    width: float = 520,
    height: float = 340,
    distances: tuple[float, ...] = (20.0, 30.0),
) -> np.ndarray:
    direction = _unit_vector(preferred_direction, np.array([0.0, -1.0]))
    perpendicular = np.array([-direction[1], direction[0]])
    directions = (
        direction,
        _unit_vector(direction + perpendicular, direction),
        _unit_vector(direction - perpendicular, direction),
        perpendicular,
        -perpendicular,
        _unit_vector(-direction + perpendicular, -direction),
        _unit_vector(-direction - perpendicular, -direction),
        -direction,
    )

    candidates = [
        np.array(
            [
                np.clip((anchor + candidate_direction * distance)[0], 36.0, width - 36.0),
                np.clip((anchor + candidate_direction * distance)[1], 18.0, height - 10.0),
            ]
        )
        for distance in distances
        for candidate_direction in directions
    ]

    def score(candidate: np.ndarray) -> float:
        penalty = 0.0
        if avoid_points is not None and len(avoid_points):
            atom_distances = np.linalg.norm(avoid_points - candidate, axis=1)
            penalty += float(np.maximum(24.0 - atom_distances, 0.0).sum())
        if occupied_labels:
            label_distances = np.array([np.linalg.norm(point - candidate) for point in occupied_labels])
            penalty += 3.0 * float(np.maximum(50.0 - label_distances, 0.0).sum())
        # Prefer a nearby candidate when several placements are otherwise
        # similarly clear. This keeps callouts tied visually to their atom.
        penalty += 0.12 * float(np.linalg.norm(candidate - anchor))
        return penalty

    label_point = min(candidates, key=score)
    occupied_labels.append(label_point)
    return label_point


def _atom_label_position(
    point: np.ndarray,
    neighbour_points: np.ndarray,
    avoid_points: np.ndarray | None,
    occupied_labels: list[np.ndarray],
    width: float = 520,
    height: float = 340,
) -> np.ndarray:
    """Place a highlighted atom label away from its local bonded geometry."""
    if len(neighbour_points):
        preferred_direction = point - neighbour_points.mean(axis=0)
    else:
        preferred_direction = point - np.array([width / 2.0, height / 2.0])
    return _best_label_position(
        point,
        preferred_direction,
        avoid_points,
        occupied_labels,
        width=width,
        height=height,
        distances=(14.0, 20.0, 28.0),
    )


def _measurement_label_position(  # noqa: PLR0913 -- SVG layout inputs are intentionally explicit
    point_1: np.ndarray,
    point_2: np.ndarray,
    avoid_points: np.ndarray | None,
    occupied_labels: list[np.ndarray],
    preferred_sign: int = 1,
    width: float = 520,
    height: float = 340,
) -> np.ndarray:
    """Place a distance label perpendicular to its measurement line."""
    midpoint = (point_1 + point_2) / 2
    line_vector = point_2 - point_1
    line_direction = _unit_vector(line_vector, np.array([1.0, 0.0]))
    perpendicular = np.array([-line_direction[1], line_direction[0]])

    signs = (1 if preferred_sign >= 0 else -1, -1 if preferred_sign >= 0 else 1)
    candidates: list[tuple[np.ndarray, float]] = []
    for side_index, sign in enumerate(signs):
        for offset_index, offset in enumerate((20.0, 28.0, 38.0)):
            for along in (0.0, 10.0, -10.0):
                raw = midpoint + perpendicular * sign * offset + line_direction * along
                point = np.array(
                    [
                        np.clip(raw[0], 36.0, width - 36.0),
                        np.clip(raw[1], 18.0, height - 10.0),
                    ]
                )
                # Prefer an exactly perpendicular placement on the requested
                # side, but allow small along-line shifts when the area is busy.
                preference_penalty = 3.0 * side_index + 0.55 * offset_index + 0.08 * abs(along)
                candidates.append((point, preference_penalty))

    def score(candidate: tuple[np.ndarray, float]) -> float:
        point, preference_penalty = candidate
        penalty = preference_penalty
        if avoid_points is not None and len(avoid_points):
            atom_distances = np.linalg.norm(avoid_points - point, axis=1)
            penalty += 1.4 * float(np.maximum(24.0 - atom_distances, 0.0).sum())
        if occupied_labels:
            label_distances = np.array([np.linalg.norm(label - point) for label in occupied_labels])
            penalty += 3.5 * float(np.maximum(50.0 - label_distances, 0.0).sum())
        # Among similarly clear positions, favour the label closest to the
        # measurement midpoint rather than letting it drift across the panel.
        penalty += 0.18 * float(np.linalg.norm(point - midpoint))
        return penalty

    label_point = min(candidates, key=score)[0]
    occupied_labels.append(label_point)
    return label_point


def _svg_condition_label(
    point: np.ndarray,
    label: str,
    target: np.ndarray | None,
    label_index: int,
    avoid_points: np.ndarray | None = None,
    occupied_labels: list[np.ndarray] | None = None,
) -> str:
    occupied = occupied_labels if occupied_labels is not None else []
    if target is None or float(np.linalg.norm(point - target)) < _MIN_VECTOR_NORM:
        fallback = ((1.0, -1.0), (1.0, 1.0), (-1.0, -1.0), (-1.0, 1.0))
        direction = np.array(fallback[label_index % len(fallback)])
    else:
        direction = point - target

    label_point = _best_label_position(
        point,
        direction,
        avoid_points,
        occupied,
        distances=(18.0, 28.0, 38.0),
    )
    text_width = max(48.0, min(150.0, 6.0 * len(label) + 12.0))
    rect_x = label_point[0] - text_width / 2
    rect_y = label_point[1] - 13.0

    return (
        "<g class='condition-callout'>" + f"<line x1='{point[0]:.1f}' y1='{point[1]:.1f}' x2='{label_point[0]:.1f}' "
        f"y2='{label_point[1] - 3:.1f}' class='label-leader'/>"
        + f"<rect x='{rect_x:.1f}' y='{rect_y:.1f}' width='{text_width:.1f}' height='19' "
        "rx='5' class='label-bg'/>" + f"<text x='{label_point[0]:.1f}' y='{label_point[1]:.1f}' text-anchor='middle' "
        f"class='condition-label'>{escape(label)}</text>" + "</g>"
    )


def _svg_measurement_label(  # noqa: PLR0913 -- SVG layout inputs are intentionally explicit
    point_1: np.ndarray,
    point_2: np.ndarray,
    label: str,
    avoid_points: np.ndarray | None = None,
    occupied_labels: list[np.ndarray] | None = None,
    preferred_sign: int = 1,
    width: float = 520,
    height: float = 340,
) -> str:
    """Render an offset label for a measurement already drawn by RDKit."""
    occupied = occupied_labels if occupied_labels is not None else []
    label_point = _measurement_label_position(
        point_1,
        point_2,
        avoid_points,
        occupied,
        preferred_sign=preferred_sign,
        width=width,
        height=height,
    )

    text_width = max(58.0, 6.5 * len(label) + 14.0)
    rect_x = label_point[0] - text_width / 2
    rect_y = label_point[1] - 13.0
    return (
        "<g class='rdkit-measurement-label'>"
        + f"<rect x='{rect_x:.1f}' y='{rect_y:.1f}' width='{text_width:.1f}' height='19' "
        "rx='5' class='measure-bg'/>" + f"<text x='{label_point[0]:.1f}' y='{label_point[1]:.1f}' text-anchor='middle' "
        f"class='measure-label'>{escape(label)}</text>" + "</g>"
    )


def _svg_measurement(
    point_1: np.ndarray,
    point_2: np.ndarray,
    label: str,
    avoid_points: np.ndarray | None = None,
    occupied_labels: list[np.ndarray] | None = None,
    preferred_sign: int = 1,
) -> str:
    occupied = occupied_labels if occupied_labels is not None else []
    label_point = _measurement_label_position(
        point_1,
        point_2,
        avoid_points,
        occupied,
        preferred_sign=preferred_sign,
    )

    text_width = max(58.0, 6.5 * len(label) + 14.0)
    rect_x = label_point[0] - text_width / 2
    rect_y = label_point[1] - 13.0
    line_1, line_2 = _trim_line_to_circle_edges(
        point_1,
        point_2,
        radius_1=5.3,
        radius_2=5.3,
    )
    return (
        _svg_line(line_1, line_2, "measurement")
        + f"<rect x='{rect_x:.1f}' y='{rect_y:.1f}' width='{text_width:.1f}' height='19' "
        "rx='5' class='measure-bg'/>" + f"<text x='{label_point[0]:.1f}' y='{label_point[1]:.1f}' text-anchor='middle' "
        f"class='measure-label'>{escape(label)}</text>"
    )


def _distance_pair(diagnostic: Diagnostic) -> tuple[int, int] | None:
    if (
        diagnostic.kind in {"bond_length", "internal_clash"}
        and len(diagnostic.ligand_atoms) >= _MIN_DISTANCE_PAIR_ATOMS
    ):
        return diagnostic.ligand_atoms[0].index, diagnostic.ligand_atoms[1].index
    if diagnostic.kind == "bond_angle" and len(diagnostic.ligand_atoms) >= _MIN_ANGLE_ATOMS:
        return diagnostic.ligand_atoms[0].index, diagnostic.ligand_atoms[2].index
    return None


def _reported_distance(diagnostic: Diagnostic) -> float | None:
    for key in ("distance", "distance_1_3"):
        value = diagnostic.metrics.get(key)
        if isinstance(value, (int, float)) and isfinite(float(value)):
            return float(value)
    return None


def _diagnostic_table(diagnostics: list[Diagnostic]) -> str:
    kind = diagnostics[0].kind
    headers, rows = _table_rows(kind, diagnostics)
    if not rows:
        return ""
    head = "".join(f"<th>{escape(header)}</th>" for header in headers)
    body = "".join("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>" for row in rows)
    return f"<table class='metrics'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _table_rows(  # noqa: PLR0911 -- each diagnostic kind has a distinct table schema
    kind: str, diagnostics: list[Diagnostic]
) -> tuple[list[str], list[list[str]]]:
    rows: list[list[str]] = []
    if kind == "bond_length":
        headers = ["Ligand atoms", "Distance", "DG bounds", "Deviation", "Allowed"]
        for diagnostic in diagnostics:
            metrics = diagnostic.metrics
            rows.append(
                [
                    _atoms_html(diagnostic.ligand_atoms),
                    _angstrom(metrics["distance"]),
                    f"{metrics['lower_bound']:.3f}-{metrics['upper_bound']:.3f} A",
                    f"{100 * metrics['percent_error']:.1f}%",
                    f"{100 * metrics['threshold']:.1f}%",
                ]
            )
        return headers, rows
    if kind == "bond_angle":
        headers = ["Ligand atoms", "1,3 distance", "DG bounds", "Deviation", "Allowed"]
        for diagnostic in diagnostics:
            metrics = diagnostic.metrics
            rows.append(
                [
                    _atoms_html(diagnostic.ligand_atoms),
                    _angstrom(metrics["distance_1_3"]),
                    f"{metrics['lower_bound']:.3f}-{metrics['upper_bound']:.3f} A",
                    f"{100 * metrics['absolute_percent_error']:.1f}%",
                    f"{100 * metrics['threshold']:.1f}%",
                ]
            )
        return headers, rows
    if kind == "internal_clash":
        headers = ["Ligand atoms", "Distance", "DG lower bound", "Deviation", "Allowed"]
        for diagnostic in diagnostics:
            metrics = diagnostic.metrics
            rows.append(
                [
                    _atoms_html(diagnostic.ligand_atoms),
                    _angstrom(metrics["distance"]),
                    _angstrom(metrics["lower_bound"]),
                    f"{100 * metrics['percent_error']:.1f}%",
                    f"{100 * metrics['threshold']:.1f}%",
                ]
            )
        return headers, rows
    if kind in {"flatness", "nonflatness"}:
        relation = "Minimum required" if kind == "nonflatness" else "Maximum allowed"
        headers = ["Ligand atoms", "Matched system", "Max plane distance", relation]
        for diagnostic in diagnostics:
            metrics = diagnostic.metrics
            rows.append(
                [
                    _atoms_html(diagnostic.ligand_atoms),
                    escape(str(metrics["system_type"])),
                    _angstrom(metrics["max_distance"]),
                    _angstrom(metrics["threshold"]),
                ]
            )
        return headers, rows
    if kind == "intermolecular_clash":
        headers = ["Ligand atom", "Conditioning atom", "Distance", "Scaled radii", "Relative distance", "Cutoff"]
        for diagnostic in diagnostics:
            metrics = diagnostic.metrics
            rows.append(
                [
                    _atoms_html(diagnostic.ligand_atoms),
                    _atoms_html(diagnostic.condition_atoms),
                    _angstrom(metrics["distance"]),
                    _angstrom(metrics["sum_radii_scaled"]),
                    f"{metrics['relative_distance']:.3f}",
                    f"{metrics['clash_cutoff']:.3f}",
                ]
            )
        return headers, rows
    if kind == "intermolecular_max_distance":
        headers = ["Ligand atom", "Conditioning atom", "Nearest distance", "Maximum allowed"]
        for diagnostic in diagnostics:
            metrics = diagnostic.metrics
            distance = _angstrom(metrics["distance"]) if "distance" in metrics else "not found within search distance"
            rows.append(
                [
                    _atoms_html(diagnostic.ligand_atoms),
                    _atoms_html(diagnostic.condition_atoms),
                    distance,
                    _angstrom(metrics["max_distance"]),
                ]
            )
        return headers, rows
    return [], []


def _atoms_html(atoms: tuple[AtomRef, ...]) -> str:
    return ", ".join(f"<code>{escape(_atom_label(atom))}</code>" for atom in atoms) or "-"


def _atom_label(atom: AtomRef) -> str:
    label = f"{atom.index} ({atom.element})"
    if atom.residue_name is None:
        return label
    residue = atom.residue_name
    if atom.chain_id:
        residue += f" {atom.chain_id}"
    if atom.residue_number is not None:
        residue += f":{atom.residue_number}"
    if atom.insertion_code:
        residue += atom.insertion_code
    if atom.atom_name:
        residue += f" {atom.atom_name}"
    return f"{label}; {residue}"


def _condition_label(atom_id: int, element: str, info: Any) -> str:
    if info is None:
        return f"{element}{atom_id}"
    residue = info.GetResidueName().strip()
    chain = info.GetChainId().strip()
    number = info.GetResidueNumber()
    atom_name = info.GetName().strip()
    prefix = f"{residue} {chain}:{number}" if chain else f"{residue}:{number}"
    return f"{prefix} {atom_name}".strip()


def _angstrom(value: Any) -> str:
    return f"{float(value):.3f} A"


def _stat_card(label: str, value: int) -> str:
    return f"<div class='stat'><strong>{value}</strong><span>{escape(label)}</span></div>"


_CSS = """
body { font-family: Arial, sans-serif; max-width: 1200px; margin: 30px auto; padding: 0 18px;
       line-height: 1.45; color: #222; }
h1, h2, h3, h4, h5 { margin-bottom: .45rem; }
table { border-collapse: collapse; width: 100%; margin: 14px 0 24px; font-size: 0.92rem; }
th, td { border: 1px solid #d5d5d5; padding: 7px 9px; text-align: left; vertical-align: top; }
th { background: #f3f3f3; }
code { background: #f5f5f5; padding: 2px 4px; border-radius: 3px; }
.notice { background: #fff3cd; padding: 11px 13px; border-left: 5px solid #e0a800; }
.pass-box { background: #eef8ef; padding: 12px; border-left: 5px solid #4f8f55; }
.stats { display: flex; gap: 12px; margin: 18px 0 26px; flex-wrap: wrap; }
.stat { border: 1px solid #ddd; border-radius: 8px; padding: 10px 18px; min-width: 100px; }
.stat strong { display: block; font-size: 1.4rem; }
.stat span { color: #555; }
.status { font-weight: bold; padding: 2px 7px; border-radius: 4px; }
.status.pass { background: #e8f5e9; }
.status.fail { background: #fdecea; }
section { border-top: 2px solid #ddd; padding-top: 12px; margin-top: 30px; }
.meta { color: #666; margin-top: 0; word-break: break-all; }
.diagnostic { margin: 22px 0 30px; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(380px, 1fr)); gap: 16px; align-items: start; }
.overview { max-width: 560px; }
.svg-wrap { width: 100%; border: 1px solid #ddd; border-radius: 6px; background: white; overflow: hidden; }
.svg-wrap svg { width: 100%; height: auto; display: block; }
.bond { stroke: #999; stroke-width: 2; }
.atom, .condition { stroke: white; stroke-width: 1; }
.highlight { stroke: #d95f02; stroke-width: 3; }
.measurement { stroke: #c43c39; stroke-width: 2.2; stroke-dasharray: 6 4; }
.measure-bg { fill: white; fill-opacity: 1; stroke: #ccc; stroke-width: 1; }
.atom-label, .condition-label { font-size: 11px; font-family: Arial, sans-serif; fill: #222; font-weight: bold; }
.atom-label { paint-order: stroke; stroke: white; stroke-width: 3px; stroke-linejoin: round; }
.condition-label { font-size: 10px; }
.label-bg { fill: white; fill-opacity: .94; stroke: #ddd; stroke-width: 1; }
.label-leader { stroke: #aaa; stroke-width: 1; }
.measure-label, .annotation { font-size: 11px; font-family: Arial, sans-serif; fill: #222; }
.measure-label { font-weight: bold; }
.atom-legend {
  display: flex; flex-wrap: wrap; gap: 8px 12px; align-items: center;
  margin: 7px 2px 0; color: #555; font-size: .82rem;
}
.legend-item { display: inline-flex; align-items: center; gap: 4px; }
.legend-dot {
  width: 10px; height: 10px; border-radius: 50%; border: 1px solid #bbb;
  display: inline-block; box-sizing: border-box;
}
.problem-dot { border-width: 2px; }
.metrics { margin-top: 14px; }
a { color: inherit; }
"""
