
====================================
Quick start
====================================


Installation
====================================

MolBuster can be installed from PyPI.

.. code-block:: bash

   >>> pip install molbuster


Usage
====================================

MolBuster can be used as via the ``bust`` command or as a Python Library.

Use ``bust mol`` to check molecules without conditioning.

.. command-output:: bust mol generated_molecules.sdf --outfmt short
  :cwd: inputs
  :ellipsis: 12

.. command-output:: bust mol generated_molecules.sdf --outfmt long
  :cwd: inputs
  :ellipsis: 12

Use ``bust dock`` to check docked ligands or generated molecules conditioned on a protein.

.. command-output:: bust dock generated_ligands.sdf protein.pdb --outfmt short
  :cwd: inputs
  :ellipsis: 12

.. command-output:: bust dock generated_ligands.sdf protein.pdb --outfmt long
  :cwd: inputs
  :ellipsis: 12

Use ``bust redock`` to check a series of re-docked ligands against the crystal ligand and protein.

.. command-output:: bust redock redocked_ligand.sdf crystal_ligand.sdf protein.pdb --outfmt short
  :cwd: inputs
  :ellipsis: 12

.. command-output:: bust redock redocked_ligand.sdf crystal_ligand.sdf protein.pdb --outfmt long
  :cwd: inputs
  :ellipsis: 25

For more usage examples and bulk processing see the documentation for the `command line tool  <cli.rst>`_
and the `Python library <python_library.ipynb>`_.
