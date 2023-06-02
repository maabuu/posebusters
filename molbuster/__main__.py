"""
Run molbuster module as a script.

Example: python -m molbuster mol_pred.sdf mol_true.sdf mol_cond.pdb
"""

from .cli import main

if __name__ == "__main__":
    main()
