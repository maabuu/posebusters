.. highlight:: none

.. _ref_cli:

====================================
Command line tool
====================================

MolBuster provides the command ``bust`` which has four subcommands ``mol``, ``dock``, ``redock``, and
``table`` for checking generated molecules with different deafult settings for different types
of conditioning and input.

You can run all commands with no arguments to get more info
about each command.


.. command-output:: bust


.. command-output:: bust redock
  :returncode: 2



Sub-commands
====================================

Use ``bust mol`` to check a series of molecules within one ``.sdf`` file.

.. command-output:: bust mol generated_molecules.sdf --outfmt long
  :cwd: inputs
  :ellipsis: 21

Use ``bust dock`` to check a docked ligand or generated molecule conditioned on a protein.

.. command-output:: bust dock generated_ligands.sdf protein.pdb --outfmt long
  :cwd: inputs
  :ellipsis: 21

Use ``bust redock`` to check a series of re-docked ligands against the crystal ligand and protein.

.. command-output:: bust redock redocked_ligand.sdf crystal_ligand.sdf protein.pdb --outfmt long
  :cwd: inputs
  :ellipsis: 21

Use ``bust table`` to bulk check multiple sets of files for any of the three modes above

.. command-output:: bust table molecule_table.csv --outfmt long
  :cwd: inputs
  :ellipsis: 21


Output format options
====================================


The long format lists each test result for each molecule/conformation.

.. command-output:: bust mol generated_molecules.sdf --outfmt long
  :cwd: inputs
  :ellipsis: 21

The short format is the default output format.

.. command-output:: bust mol generated_molecules.sdf --outfmt short
  :cwd: inputs

For copying and saving the output use the ``csv`` option.

.. command-output:: bust mol generated_molecules.sdf --outfmt csv
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
