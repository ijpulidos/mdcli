# mdcli
Some tools for the analysis and pre/post-processing of molecular dynamics data.

## Requirements

It is recommended that you use an anaconda/conda environment, since it would allow to install
many of the needed libs and tools in an easier way (they are already included in some conda 
channels).

### Adding conda channels

Add `salilab` conda channel to install `Modeller` which is a required tool and
`schrodinger` channel for `pymol` (as a python lib), as follows:

```bash
conda config --add channels salilab
conda config --add channels schrodinger
```

### Building environment from file

The `environment.yml` file allows to build an anaconda environment that will be
named `mdcli-env` by default, with all the required packages. To build it just
use:

```bash
conda env create -f environment.yml
```

## Features - capabilities
The current capabilities of the tools in this repository are the following.

### Completing residues 

The script `complete_residues.py` is a friendly tool that completes a PDB 
structure with both missing residues and atoms. Example of usage is as follows:

    python complete_residues.py -i /path/to/input/PDB --pdbcode PDBcode

And it will return a pdb fill structure (see output for exact name).

### Extract helical residues

Using the script `get_helices.py` you can extract the ids from residues
belonging to helical structure given a PDB file. It also outputs a PDB file with
the structure of only the helical residues. Example of its usage is:

    python get_helices.py --pdbfile ../../datafiles/conf.pdb -o /tmp/helices.pdb

For example using the `conf.pdb` file in the datafiles.