"""
Run posebusters module as a script.

Example: python -m posebusters mol_pred.sdf -l mol_true.sdf -p mol_cond.pdb
"""

from .cli import main

if __name__ == "__main__":
    main()
