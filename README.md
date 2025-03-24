![posebusters_banner](https://github.com/maabuu/posebusters/blob/b5f4c2caed1499c2b73f5297a08e60ec7d263c61/docs/source/_static/logo_banner.png?raw=true "PoseBusters")

PoseBusters: Plausibility checks for generated molecule poses.

### Paper in [Chemical Science](https://pubs.rsc.org/en/content/articlelanding/2024/sc/d3sc04185a) and preprint on [arXiv](https://arxiv.org/abs/2308.05777)

## Installation

```bash
# install with pip from PyPI
pip install posebusters
```
<!-- # install with conda from conda-forge
conda install posebusters -c conda-forge -->

## Usage

<!-- ### Command line usage -->

```bash

# Check generated molecule pose.
bust molecule_pred.sdf
bust molecule_a.sdf molecule_b.sdf
bust molecule_*.sdf

# Check new ligand generated for a given protein.
bust ligand_pred.sdf -p mol_cond.pdb

# Check re-docked ligand (a pose that should recover the ligand in a given protein-ligand crystal complex).
bust ligand_pred.sdf -l mol_true.sdf -p protein.pdb

# Check any of the three by providing a csv with files to check together
bust -t file_table.csv
```

<!-- ### Python API

```python
from dockbusters import DockBuster

# check re-docked ligand
DockBuster().bust(ligand_pred_file, ligand_crystal_file, protein_crystal_file)

# check docked ligand
DockBuster().bust(ligand_pred_file, protein_crystal_file)

# check molecule
DockBuster().bust(ligand_pred_file, protein_crystal_file)
``` -->

## Documentation

Documentation is available at [https://posebusters.readthedocs.io](https://posebusters.readthedocs.io).

For more information about the tests and for a study using PoseBusters to compare docking methods, refer to our [paper](http://dx.doi.org/10.1039/D3SC04185A) or [preprint](https://arxiv.org/abs/2308.05777):

```
@article{buttenschoen2024posebusters,
  title = {{{PoseBusters}}: {{AI-based}} Docking Methods Fail to Generate Physically Valid Poses or Generalise to Novel Sequences},
  shorttitle = {{{PoseBusters}}},
  author = {Buttenschoen, Martin and Morris, Garrett M. and Deane, Charlotte M.},
  year = "2024",
  journal = "Chemical Science",
  volume = "15",
  issue = "9",
  pages = "3130-3139",
  publisher = "The Royal Society of Chemistry",
  doi = "10.1039/D3SC04185A",
  url = "http://dx.doi.org/10.1039/D3SC04185A",
}
```

The data used for the paper is available at [https://zenodo.org/record/8278563](https://zenodo.org/record/8278563).

## Feedback & Contact

We welcome all feedback. For code issues, please open an issue. For other inquiries contact us by email.

## Thanks

This program uses software written by other people. Notably:

- RDKit - [https://github.com/rdkit/rdkit](https://github.com/rdkit/rdkit)
- Pandas - [https://github.com/pandas-dev/pandas](https://github.com/pandas-dev/pandas)
