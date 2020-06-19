#!/bin/env python

"""
Script to help extracting the helical residues from a PDB structure or getting
the residues IDs from it.
"""

# Required modules
import argparse
import pymol


def get_helices(
    pdb_path: str,
    basename: str = "protein",
    outfile: str = "./helices.pdb"
) -> list:
    """
    Function that returns the helices residues indices from a PDB structure
    file. It uses PyMol to select by secondary structure.

    :param pdb_path: String with the path for the input (original) PDB file.
    :param basename: Selection string to be used for naming only (default
                     'protein').
    :param outfile: String to path of output file (default './helices.pdb')
    """
    # Loading pdb file
    # Selecting helices and saving to pdb file
    # TODO: Check if can be done without writing an extra pdb file
    # Delete everything first
    pymol.cmd.delete("all")
    pymol.cmd.load(pdb_path, basename)
    pymol.cmd.select(basename + "_helices", selection="ss H")
    pymol.cmd.save(
        filename=f"{outfile}",
        selection=basename + "_helices"
        )

    # Reading residues indices in helices
    # Get fourth (sixth?) column (indices column) from helices pdb
    helices_indices = [
        int(x.split()[4])
        for x in open(f"{outfile}").readlines()
        if "ATOM" in x.split()[0]
        ]
    return helices_indices


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get helical residues IDs.")
    parser.add_argument(
        "--pdbfile",
        "-pdb",
        type=str,
        help="Input file in PDB format (.pdb).",
        required=True
    )
    parser.add_argument(
        "--basename",
        type=str,
        help="Name for pymol selection.",
        required=False,
        default="protein"
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="Path for output PDB file.",
        required=False,
        default="./helices.pdb"
    )
    args = parser.parse_args()

    # Get helical residues indices
    res_ids = get_helices(
        args.pdbfile,
        basename=args.basename,
        outfile=args.output
    )
    print(f"Helical residues IDs for {args.pdbfile} are:")
    print(res_ids)
