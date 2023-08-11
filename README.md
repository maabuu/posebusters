
# <img src="docs/source/_static/logo_banner.png" alt="drawing" height="75"/>

PoseBusters: Command line tool and Python package for checking the chemical
and physical sensibility of docked or generated molecules.


## Installation

```bash
# install with pip from PyPI
pip install posebusters
```
<!-- # install with conda from conda-forge
conda install posebusters -c conda-forge -->

## Usage

### Command line

```bash
# check re-docked ligand (new ligand into protein).
bust redock ligand_pred.sdf mol_true.sdf mol_cond.pdb

# check docked ligand (where crystal protein-ligand structure is known).
bust dock ligand_pred.sdf protein.pdb

# check generated molecule (only check ligand).
bust mol molecule_pred.sdf

# check multiple of the three above using a .csv input:
bust table file_table.csv
```

### Python API

```python
from dockbusters import DockBuster

# check re-docked ligand
DockBuster().bust_redock(ligand_pred_file, ligand_crystal_file, protein_crystal_file)

# check docked ligand
DockBuster().bust_dock(ligand_pred_file, protein_crystal_file)

# check molecule
DockBuster().bust_mol(ligand_pred_file, protein_crystal_file)
```

## Thanks

This program uses software written by other people. Notably:

- RDKit [https://github.com/rdkit/rdkit](https://github.com/rdkit/rdkit)
- Pandas [https://github.com/pandas-dev/pandas](https://github.com/pandas-dev/pandas)

