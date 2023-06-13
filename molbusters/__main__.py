"""
Run molbusters module as a script.

Example: python -m molbusters mol_pred.sdf mol_true.sdf mol_cond.pdb
"""

from .cli import main

if __name__ == "__main__":
    main()
