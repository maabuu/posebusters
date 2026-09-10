"""Diagnostic reporting for failed PoseBusters checks."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd
from rdkit.Chem.rdchem import Mol

from .modules.distance_geometry import col_bape, col_bpe, col_pe


@dataclass(frozen=True)
class AtomRef:
    """Reference to an atom involved in a diagnostic."""

    index: int
    element: str
    atom_name: str | None = None
    residue_name: str | None = None
    residue_number: int | None = None
    chain_id: str | None = None
    insertion_code: str | None = None


@dataclass(frozen=True)
class Diagnostic:
    """Structured information explaining one failed check."""

    module: str
    check: str
    kind: str
    ligand_atoms: tuple[AtomRef, ...] = ()
    condition_atoms: tuple[AtomRef, ...] = ()
    metrics: dict[str, Any] = field(default_factory=dict)


def diagnose_module(
    module: dict[str, Any],
    module_output: dict[str, Any],
    molecules: dict[str, Any],
) -> list[Diagnostic]:
    """Create diagnostics for failed outputs from one PoseBusters module."""
    results = module_output.get("results", {})
    chosen_outputs = module.get("chosen_binary_test_output", [])

    failed_outputs = [output for output in chosen_outputs if _is_failed(results.get(output))]
    if not failed_outputs:
        return []

    failed_output_set = set(failed_outputs)

    details = module_output.get("details")
    function_name = module["function"]
    diagnostics: list[Diagnostic] = []

    if details is not None:
        if function_name == "distance_geometry":
            diagnostics.extend(
                _diagnose_distance_geometry(
                    module,
                    details,
                    failed_output_set,
                    molecules["mol_pred"],
                )
            )
        elif function_name == "flatness":
            diagnostics.extend(
                _diagnose_flatness(
                    module,
                    details,
                    failed_output_set,
                    molecules["mol_pred"],
                )
            )
        elif function_name == "intermolecular_distance":
            diagnostics.extend(
                _diagnose_intermolecular_distance(
                    module,
                    details,
                    failed_output_set,
                    molecules["mol_pred"],
                    molecules["mol_cond"],
                )
            )

    diagnosed_checks = {diagnostic.check for diagnostic in diagnostics}

    for output in failed_outputs:
        check_name = _check_name(module, output)

        if check_name in diagnosed_checks:
            continue

        diagnostics.append(
            Diagnostic(
                module=module["name"],
                check=check_name,
                kind="summary",
            )
        )

    return diagnostics


def create_diagnostic_output(key: tuple[str, str, int], diagnostics: list[Diagnostic]) -> str:
    """Format diagnostics for one predicted pose as plain text."""
    if not diagnostics:
        return ""

    file_name, molecule_name, position = key
    header = f"Diagnostics for {file_name}"
    if molecule_name:
        header += f" | molecule: {molecule_name}"
    header += f" | position: {position}"

    lines = ["", header]
    for diagnostic in diagnostics:
        lines.extend(_format_diagnostic(diagnostic))
    return "\n".join(lines) + "\n"


def _diagnose_distance_geometry(
    module: dict[str, Any],
    details: dict[str, pd.DataFrame],
    failed_outputs: set[str],
    mol_pred: Mol,
) -> list[Diagnostic]:
    diagnostics: list[Diagnostic] = []
    parameters = module.get("parameters", {})

    if "bond_lengths_within_bounds" in failed_outputs:
        threshold = parameters.get("threshold_bad_bond_length", 0.2)
        bonds = details["bonds"]
        failures = bonds[(bonds[col_pe] < -threshold) | (bonds[col_pe] > threshold)]
        for _, row in failures.iterrows():
            atom_ids = tuple(int(i) for i in row["atom_pair"])
            diagnostics.append(
                Diagnostic(
                    module=module["name"],
                    check=_check_name(module, "bond_lengths_within_bounds"),
                    kind="bond_length",
                    ligand_atoms=_ligand_atom_refs(mol_pred, atom_ids),
                    metrics={
                        "distance": float(row["distance"]),
                        "lower_bound": float(row["lower_bound"]),
                        "upper_bound": float(row["upper_bound"]),
                        "percent_error": float(row[col_pe]),
                        "threshold": float(threshold),
                    },
                )
            )

    if "bond_angles_within_bounds" in failed_outputs:
        threshold = parameters.get("threshold_bad_angle", 0.2)
        angles = details["angles"]
        failures = angles[angles[col_bape] > threshold]
        for _, row in failures.iterrows():
            atom_ids = tuple(int(i) for i in row["angle"])
            diagnostics.append(
                Diagnostic(
                    module=module["name"],
                    check=_check_name(module, "bond_angles_within_bounds"),
                    kind="bond_angle",
                    ligand_atoms=_ligand_atom_refs(mol_pred, atom_ids),
                    metrics={
                        "distance_1_3": float(row["distance"]),
                        "lower_bound": float(row["lower_bound"]),
                        "upper_bound": float(row["upper_bound"]),
                        "absolute_percent_error": float(row[col_bape]),
                        "threshold": float(threshold),
                    },
                )
            )

    if "no_internal_clash" in failed_outputs:
        threshold = parameters.get("threshold_clash", 0.2)
        clashes = details["clash"]
        failures = clashes[clashes[col_bpe] < -threshold]
        for _, row in failures.iterrows():
            atom_ids = tuple(int(i) for i in row["atom_pair"])
            diagnostics.append(
                Diagnostic(
                    module=module["name"],
                    check=_check_name(module, "no_internal_clash"),
                    kind="internal_clash",
                    ligand_atoms=_ligand_atom_refs(mol_pred, atom_ids),
                    metrics={
                        "distance": float(row["distance"]),
                        "lower_bound": float(row["lower_bound"]),
                        "percent_error": float(row[col_bpe]),
                        "threshold": float(threshold),
                    },
                )
            )

    return diagnostics


def _diagnose_flatness(
    module: dict[str, Any],
    details: dict[str, Any],
    failed_outputs: set[str],
    mol_pred: Mol,
) -> list[Diagnostic]:
    if "flatness_passes" not in failed_outputs:
        return []

    parameters = module.get("parameters", {})
    threshold = float(parameters.get("threshold_flatness", 0.1))
    check_nonflat = bool(parameters.get("check_nonflat", False))
    diagnostics: list[Diagnostic] = []
    seen_groups: set[tuple[int, ...]] = set()

    for system_type, atom_group, max_distance, passes in zip(
        details["type"],
        details["planar_group"],
        details["max_distance"],
        details["flatness_passes"],
    ):
        if passes:
            continue
        atom_ids = tuple(int(i) for i in atom_group)
        canonical_group = tuple(sorted(atom_ids))
        if canonical_group in seen_groups:
            continue
        seen_groups.add(canonical_group)
        diagnostics.append(
            Diagnostic(
                module=module["name"],
                check=_check_name(module, "flatness_passes"),
                kind="nonflatness" if check_nonflat else "flatness",
                ligand_atoms=_ligand_atom_refs(mol_pred, atom_ids),
                metrics={
                    "system_type": str(system_type),
                    "max_distance": float(max_distance),
                    "threshold": threshold,
                },
            )
        )

    return diagnostics


def _diagnose_intermolecular_distance(
    module: dict[str, Any],
    details: pd.DataFrame,
    failed_outputs: set[str],
    mol_pred: Mol,
    mol_cond: Mol,
) -> list[Diagnostic]:
    diagnostics: list[Diagnostic] = []
    parameters = module.get("parameters", {})

    if "not_too_far_away" in failed_outputs:
        metrics = {
            "max_distance": float(parameters.get("max_distance", 5.0)),
            "search_distance": float(parameters.get("search_distance", 6.0)),
        }
        ligand_atoms: tuple[AtomRef, ...] = ()
        condition_atoms: tuple[AtomRef, ...] = ()
        if not details.empty:
            row = details.loc[details["distance"].idxmin()]
            ligand_idx = int(row["ligand_atom_id"])
            condition_idx = int(row["protein_atom_id"])
            ligand_atoms = _ligand_atom_refs(mol_pred, (ligand_idx,))
            condition_atoms = (_condition_atom_ref(mol_cond, condition_idx),)
            metrics["distance"] = float(row["distance"])

        diagnostics.append(
            Diagnostic(
                module=module["name"],
                check=_check_name(module, "not_too_far_away"),
                kind="intermolecular_max_distance",
                ligand_atoms=ligand_atoms,
                condition_atoms=condition_atoms,
                metrics=metrics,
            )
        )

    if "no_clashes" in failed_outputs:
        failures = details[details["clash"]]
        for _, row in failures.iterrows():
            ligand_idx = int(row["ligand_atom_id"])
            condition_idx = int(row["protein_atom_id"])
            diagnostics.append(
                Diagnostic(
                    module=module["name"],
                    check=_check_name(module, "no_clashes"),
                    kind="intermolecular_clash",
                    ligand_atoms=_ligand_atom_refs(mol_pred, (ligand_idx,)),
                    condition_atoms=(_condition_atom_ref(mol_cond, condition_idx),),
                    metrics={
                        "distance": float(row["distance"]),
                        "sum_radii": float(row["sum_radii"]),
                        "sum_radii_scaled": float(row["sum_radii_scaled"]),
                        "relative_distance": float(row["relative_distance"]),
                        "clash_cutoff": float(parameters.get("clash_cutoff", 0.75)),
                    },
                )
            )

    return diagnostics


def _check_name(module: dict[str, Any], output: str) -> str:
    return str(module.get("rename_outputs", {}).get(output, output))


def _is_failed(value: Any) -> bool:
    if value is None:
        return False
    try:
        if bool(pd.isna(value)):
            return False
    except (TypeError, ValueError):
        pass
    return not bool(value)


def _ligand_atom_refs(mol: Mol, atom_ids: tuple[int, ...]) -> tuple[AtomRef, ...]:
    return tuple(AtomRef(index=i, element=mol.GetAtomWithIdx(i).GetSymbol()) for i in atom_ids)


def _condition_atom_ref(mol: Mol, atom_id: int) -> AtomRef:
    atom = mol.GetAtomWithIdx(atom_id)
    info = atom.GetPDBResidueInfo()
    if info is None:
        return AtomRef(index=atom_id, element=atom.GetSymbol())

    return AtomRef(
        index=atom_id,
        element=atom.GetSymbol(),
        atom_name=info.GetName().strip() or None,
        residue_name=info.GetResidueName().strip() or None,
        residue_number=info.GetResidueNumber(),
        chain_id=info.GetChainId().strip() or None,
        insertion_code=info.GetInsertionCode().strip() or None,
    )


def _format_diagnostic(diagnostic: Diagnostic) -> list[str]:
    lines = [f"  {diagnostic.check}:"]
    ligand_atoms = _format_atoms(diagnostic.ligand_atoms)
    condition_atoms = _format_atoms(diagnostic.condition_atoms)
    metrics = diagnostic.metrics

    if ligand_atoms:
        lines.append(f"    ligand atoms (RDKit indices): {ligand_atoms}")
    if condition_atoms:
        lines.append(f"    conditioning atom: {condition_atoms}")

    if diagnostic.kind == "bond_length":
        lines.append(f"    distance: {metrics['distance']:.3f} Å")
        lines.append(f"    DG bounds: {metrics['lower_bound']:.3f}-{metrics['upper_bound']:.3f} Å")
        lines.append(f"    deviation: {100 * metrics['percent_error']:.1f}%")
        lines.append(f"    allowed deviation: {100 * metrics['threshold']:.1f}%")
    elif diagnostic.kind == "bond_angle":
        lines.append(f"    1,3 distance: {metrics['distance_1_3']:.3f} Å")
        lines.append(f"    DG bounds: {metrics['lower_bound']:.3f}-{metrics['upper_bound']:.3f} Å")
        lines.append(f"    deviation beyond bounds: {100 * metrics['absolute_percent_error']:.1f}%")
        lines.append(f"    allowed deviation: {100 * metrics['threshold']:.1f}%")
    elif diagnostic.kind == "internal_clash":
        lines.append(f"    distance: {metrics['distance']:.3f} Å")
        lines.append(f"    DG lower bound: {metrics['lower_bound']:.3f} Å")
        lines.append(f"    deviation below bound: {100 * metrics['percent_error']:.1f}%")
        lines.append(f"    allowed deviation: {100 * metrics['threshold']:.1f}%")
    elif diagnostic.kind in {"flatness", "nonflatness"}:
        lines.append(f"    matched system: {metrics['system_type']}")
        lines.append(f"    maximum distance from plane: {metrics['max_distance']:.3f} Å")
        relation = "minimum required" if diagnostic.kind == "nonflatness" else "maximum allowed"
        lines.append(f"    {relation}: {metrics['threshold']:.3f} Å")
    elif diagnostic.kind == "intermolecular_clash":
        lines.append(f"    distance: {metrics['distance']:.3f} Å")
        lines.append(f"    sum of radii: {metrics['sum_radii']:.3f} Å")
        lines.append(f"    scaled sum of radii: {metrics['sum_radii_scaled']:.3f} Å")
        lines.append(f"    relative distance: {metrics['relative_distance']:.3f}")
        lines.append(f"    minimum allowed relative distance: {metrics['clash_cutoff']:.3f}")
    elif diagnostic.kind == "intermolecular_max_distance":
        if "distance" in metrics:
            lines.append(f"    nearest retained distance: {metrics['distance']:.3f} Å")
        else:
            lines.append(f"    no conditioning atom found within search distance: {metrics['search_distance']:.3f} Å")
        lines.append(f"    maximum allowed distance: {metrics['max_distance']:.3f} Å")
    elif diagnostic.kind == "summary":
        lines.append("    failed (no additional diagnostic details)")

    return lines


def _format_atoms(atoms: tuple[AtomRef, ...]) -> str:
    return ", ".join(_format_atom(atom) for atom in atoms)


def _format_atom(atom: AtomRef) -> str:
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
