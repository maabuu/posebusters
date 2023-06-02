from __future__ import annotations

from molbuster.modules.intermolecular_distance import check_intermolecular_distance


def test_check_intermolecular_distance(protein_2bm2, ligand_2bm2):
    out = check_intermolecular_distance(ligand_2bm2, protein_2bm2)
    out["results"]["no_clashes"] is True
