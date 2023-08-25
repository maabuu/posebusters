.. Top level package

====================================
Python API
====================================

PoseBusters class
====================================

The PoseBusters class collects the molecules to test, runs the *modules*, and reports the test results.

.. autoclass:: posebusters.PoseBusters
    :members:
    :undoc-members:

Modules
====================================

A PoseBusters *module* is a function that takes one or more of ``mol_pred``, ``mol_true``, and ``mol_cond`` as input and returns one or more test results as a dictionary.

Inputs
    - Must take one of ``mol_pred``, ``mol_true``, and ``mol_cond`` provided as RDKit molecules.
    - Other inputs are parameters for which default values must be specified.
Outputs
    - The output must be a dictionary with at least the `results` entry.
    - The ``results`` entry contains a dictionary with keys corresponding to the test names and the test outcomes.
    - Other output entries to contain further results e.g. lengths and bound for all bonds in ligand.


Distance Geometry
------------------------------------
.. autofunction:: posebusters.modules.distance_geometry.check_geometry

Energy Ratio
------------------------------------
.. autofunction:: posebusters.modules.energy_ratio.check_energy_ratio

Flatness
------------------------------------
.. autofunction:: posebusters.modules.flatness.check_flatness

Identity
------------------------------------
.. autofunction:: posebusters.modules.identity.check_identity

Intermolecular Distance
------------------------------------
.. autofunction:: posebusters.modules.intermolecular_distance.check_intermolecular_distance

Loading
------------------------------------
.. autofunction:: posebusters.modules.loading.check_loading

RMSD
------------------------------------
.. autofunction:: posebusters.modules.rmsd.check_rmsd

Volume Overlap
------------------------------------
.. autofunction:: posebusters.modules.volume_overlap.check_volume_overlap
