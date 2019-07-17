#!/usr/bin/env python

"""
Script that builds an align file for MODELLER software to complete missing 
residues in a PDB structure file.

Requires BioPython>=1.74 and Modeller>=9.21.
"""

import argparse
import shutil

# Modules to read PDB and build align file
from Bio.PDB import PDBParser
from Bio.PDB.Selection import unfold_entities
import numpy as np
from Bio.SeqUtils import seq1
import textwrap

# Modules for completing missing residues
from modeller import *
from modeller.automodel import *    # Load the automodel class


def insert_gap(aaindex, residue_list):
    """
    Inserts a gap in the residue list in the position corresponding to the amino acid index.
    """
    tmplist = residue_list.copy()
    for index, res in enumerate(tmplist):
        if aaindex < res[0]:  # Given index less than current index in list
            residue_list.insert(index, (aaindex, "-"))
            return tmplist
    residue_list.append((aaindex, "-"))


def build_align_file(input_pdb, output_align_file, pdbcode):
    """
    Function that takes a PDB filepath, detects missing residues and builds a MODELLER align file with this information,
    to be later used for completing residues.
    :param input_pdb: PDB filepath
    :param output_align_file: Filepath for output align file.
    :param pdbcode: Code identifier for the PDB structure. Ex: 2y8d
    :return:
    """

    # Read structure and extract present and missing residues
    pdbparser = PDBParser()
    structure = pdbparser.get_structure(pdbcode, input_pdb)
    residues = unfold_entities(structure, "R")
    missing_residues = structure.header['missing_residues']

    residues_list = [(residue.id[1], seq1(residue.resname)) for residue in residues if residue.id[0]==" " ]
    for mis_res in missing_residues:
        insert_gap(mis_res['ssseq'], residues_list)

    cadena = "".join(np.array(residues_list)[:, 1])

    # Make the line width the correct/expected one for modeller align file
    textwrap.wrap(cadena, width=75, break_on_hyphens=False)

    full_seq = cadena
    for mis_res in missing_residues:
        full_seq = full_seq.replace("-", seq1(mis_res["res_name"]), 1)

    # For checking full_seq
    # print(full_seq)

    # Writing to file (test)
    # Remember sequences have to end with the * character
    with open(output_align_file, "w") as file:
        # Writing structure section
        file.write(">P1;" + structure.id + "\n")
        file.write("structureX:" + structure.id + 8*':.' + "\n")
        for line in textwrap.wrap(cadena+"*", width=75, break_on_hyphens=False):
            file.write("%s\n" % line)
        # Writing sequence section
        file.write(">P1;" + structure.id + "_fill\n")
        file.write("sequence:" + structure.id + 8*':.' + "\n")
        for line in textwrap.wrap(full_seq+"*", width=75, break_on_hyphens=False):
            file.write("%s\n" % line)


def complete_residues(code, align_file):
    """
    Function that completes residues based on an alignment file using MODELLER software.
    :param code: PDB code identifier of the structure with missing residues.
    :param align_file: Path to the align-formatted file with gaps as missing residues.
    :return:
    """
    # Get the sequence of the coded PDB file, and write to an alignment file
    e = environ()
    m = model(e, file=code)
    aln = alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code + '.seq')

    # Completing residues
    log.verbose()
    env = environ()

    # directories for input atom files (local dir)
    env.io.atom_files_directory = ['.']

    a = loopmodel(env, alnfile=align_file,
                  knowns=code, sequence=code+'_fill')

    a.starting_model = 1
    a.ending_model = 1

    a.loop.starting_model = 1
    a.loop.ending_model = 2
    a.loop.md_level = refine.fast

    a.make()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Builds alignment file for MODELLER.")
    parser.add_argument('--input', '-i', type=str, help='Input PDB file with missing residues.', required=True)
    parser.add_argument('--output', '-o', type=str, help='Output align (.ali) file for completing residues.',
                        required=True)
    parser.add_argument('--pdbcode', type=str, help='PDB code identifier for the structure.',
                        required=True)

    args = parser.parse_args()

    build_align_file(args.input, args.output, args.pdbcode)

    # the PDB file has to be in the same directory, copying and using the code as name.
    pdb_path = args.input
    code = args.pdbcode
    shutil.copy(pdb_path, "./"+code+".pdb")

    complete_residues(code, args.output)
