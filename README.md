# pymdtools
Some tools for the analysis and pre/post-processing of molecular dynamics data.

## Completing residues 

The script `complete_residues.py` is a friendly tool that completes a PDB 
structure with both missing residues and atoms. Example of usage is a folows:

    python complete_residues.py -i /path/to/input/PDB --pdbcode PDBcode

And it will return a pdb fill structure (see output for exact name).
