# pymdtools
Some tools for the analysis and pre/post-processing of molecular dynamics data.

## Requirements

It is recommended that you use an anaconda/conda environment, since it would allow to install
many of the needed libs and tols in an easier way (they are already included in some conda 
channels).

### Adding conda channels

Add `salilab` conda channel to install `Modeller` which is a required tool, as follows:

```bash
conda config --add channels salilab
```

### Building environment from file

The `environment.yml` file allows to build an anaconda environment that will be named `mdcli-env`
by default, with all the required packages. To build it just use:

```bash
conda env create -f environment.yml
```

## Completing residues 

The script `complete_residues.py` is a friendly tool that completes a PDB 
structure with both missing residues and atoms. Example of usage is a folows:

    python complete_residues.py -i /path/to/input/PDB --pdbcode PDBcode

And it will return a pdb fill structure (see output for exact name).
