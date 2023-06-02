
====================================
Quick start
====================================

This quick start guide will show you how to install and use MolBuster to
check the quality of your docking result.

Install MolBuster
====================================

Installing ``molbuster`` will install the command-line tool and the Python library.
MolBuster is available via PyPI_, conda-forge_ and GitHub_.

.. code-block:: bash

   # install with pip from PyPI
   pip install molbuster

   # install with conda from conda-forge
   conda install molbuster -c conda-forge

   # install from source
   git clone https://github.com/maabuu/molbuster.git
   cd molbuster
   pip install flit
   flit install

.. todo::
   Create and publish pip and conda packages.

.. note::
    MolBuster requires Python 3.7 or higher and depends on RDKit_ and pandas_.


Check a docked molecule
====================================

MolBuster is called through the ``bust`` command-line tool and it provides several
sub-commands, for testing different sets of molecules, and each one accepts a different
set of arguments and options.

.. code-block:: bash

    # check re-docked ligand (new ligand into protein).
    bust redock mol_pred.sdf mol_true.sdf mol_cond.pdb

    # check docked ligand (where crystal protein-ligand structure is known).
    bust dock mol_pred.sdf protein.pdb

    # check generated molecule (only check ligand).
    bust mol molecule_pred.sdf

    # check multiple of the three above using a .csv input:
    bust table file_table.csv


MolBuster can also be used as a Python library by using the ``MolBuster`` class
which stores the test options, calls the tests, and combines the test results.

.. code-block:: python

    from molbuster import MolBuster

    # check re-docked ligand and use crystal ligand and crystal protein
    MolBuster().bust_redock(mol_pred, mol_true, mol_cond)

    # check docked ligand or generated molecule and use protein
    MolBuster().bust_dock(mol_pred, mol_cond)

    # check a generated molecule
    MolBuster().bust_dock(mol_pred)


.. _conda-forge: https://github.com/conda-forge/molbuster
.. _GitHub: https://github.com/maabuu/molbuster
.. _PyPI: https://pypi.org/project/molbuster

.. _pandas: https://pandas.pydata.org/
.. _RDKit: https://RDKit.org/
