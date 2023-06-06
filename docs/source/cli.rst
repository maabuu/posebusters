.. highlight:: none

.. _ref_cli:

====================================
Command line tool
====================================

MolBuster provides the command ``bust`` for checking generated molecules
and optionally taking a conditioning protein or ligands into account.

You can run ``bust`` with no arguments to get more info.

.. command-output:: bust

Use ``bust`` to check a series of molecules within one ``.sdf`` file.

.. command-output:: bust generated_molecules.sdf --outfmt long
  :cwd: inputs
  :ellipsis: 21

Check a docked ligand or generated molecule conditioned on a protein.

.. command-output:: bust generated_ligands.sdf -p protein.pdb --outfmt long
  :cwd: inputs
  :ellipsis: 21

Check a series of re-docked ligands against the crystal ligand and protein.

.. command-output:: bust redocked_ligand.sdf -l crystal_ligand.sdf -p protein.pdb --outfmt long
  :cwd: inputs
  :ellipsis: 21

Use the `-t` option to bulk check multiple sets of files.

.. command-output:: bust -t molecule_table.csv --outfmt long
  :cwd: inputs
  :ellipsis: 21


Output format options
====================================

The short format is the default output format.

.. command-output:: bust generated_molecules.sdf --outfmt short
  :cwd: inputs

The long format lists each test result for each molecule/conformation.

.. command-output:: bust generated_molecules.sdf --outfmt long
  :cwd: inputs
  :ellipsis: 21

For copying and saving the output use the ``csv`` option.

.. command-output:: bust generated_molecules.sdf --outfmt csv
  :cwd: inputs


Saving the output
====================================

The ``--out`` option can be used to save the output to a file.

.. command-output:: bust generated_molecules.sdf --outfmt csv --out results.csv
  :cwd: inputs



Configuration settings
====================================

MolBuster will look for configuration parameters in a yaml file ``molbuster.yml``
in standard locations:

1. ``/etc/molbuster.cfg`` or ``c:\molbuster\molbuster.cfg`` (system-wide),
2. ``~/.config/molbuster.cfg`` (``$XDG_CONFIG_HOME``) and ``~/.molbuster.cfg`` (``$HOME``)
   for global (user-wide) settings, and
3. ``molbuster.cfg`` inside the working directory.
4. File location provided by the ``--config`` command line option.

Settings from these files are merged in the listed order of preference:
user-defined values have higher priority than system-wide defaults
and project-wide settings will override all others, when defined.
