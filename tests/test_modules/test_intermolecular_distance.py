from __future__ import annotations

from posebusters.modules.intermolecular_distance import check_intermolecular_distance


def test_check_intermolecular_distance(protein_2bm2, ligand_2bm2):
    out = check_intermolecular_distance(ligand_2bm2, protein_2bm2)
    out["results"]["no_clashes"] is True


def test_check_intermolecular_distance_7ecr_sin(mol_lig_7ecr_sin, mol_cond_7ecr_sin):
    out = check_intermolecular_distance(
        mol_lig_7ecr_sin,
        mol_cond_7ecr_sin,
        radius_type="vdw",
        radius_scale=1.0,
        clash_cutoff=0.75,
        ignore_types={"organic_cofactors", "inorganic_cofactors"},
    )
    out["results"]["no_clashes"] is True
