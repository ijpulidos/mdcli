#!/usr/bin/env python

"""
Script that completes a PDB structure with missing residues or atoms.

It uses MODELLER as backend, it automatically builds a MODELLER align file and
then performs the completion with a MODELLER model and its tools.

Requires BioPython>=1.74 and Modeller>=9.21.
"""

import os
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
from modeller.automodel import *  # Load the automodel class


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


def get_chain_missing_res(missing_residues, chainID):
    """
    Function that returns a list of missing residues from a given chain identifier/letter.
    """
    result = [residue for residue in missing_residues if residue["chain"] == chainID]
    return result


def build_align_file(input_pdb, pdbcode, output_align_file="protein.ali"):
    """
    Function that takes a PDB filepath, detects missing residues and builds a MODELLER align file with this information,
    to be later used for completing residues.
    :param input_pdb: PDB filepath
    :param pdbcode: Code identifier for the PDB structure. Ex: 2y8d
    :param output_align_file: Filepath for output align file.
    :return:
    """

    # Read structure and extract present and missing residues
    pdbparser = PDBParser()
    structure = pdbparser.get_structure(pdbcode, input_pdb)
    chains = unfold_entities(structure, "C")  # Get chains

    missing_residues = structure.header[
        "missing_residues"
    ]  # Get missing residues from whole structure

    # Remove alignment file if exists
    try:
        os.remove(output_align_file)
    except FileNotFoundError:
        pass

    # Where to store the sequences from structure separated by chains/index
    whole_gapped = []
    whole_full = []

    for chain in chains:
        chain_id = chain.get_id()
        residues = unfold_entities(chain, "R")  # Get residues of chain
        missing_res_chain = get_chain_missing_res(missing_residues, chain_id)

        # Residues with empty id[0] are the 'real' residues, others are solvent or different.
        residues_list = [
            (residue.id[1], seq1(residue.resname))
            for residue in residues
            if residue.id[0] == " "
        ]
        for mis_res in missing_res_chain:
            insert_gap(mis_res["ssseq"], residues_list)

        # Sequence with gaps
        gapped_seq = "".join(np.array(residues_list)[:, 1])
        # Make the line width the correct/expected one for modeller align file
        textwrap.wrap(gapped_seq, width=75, break_on_hyphens=False)

        # Full sequence without gaps by replacing gaps with the missing res
        full_seq = gapped_seq
        for mis_res in missing_residues:
            full_seq = full_seq.replace("-", seq1(mis_res["res_name"]), 1)

        whole_gapped.append(gapped_seq)
        whole_full.append(full_seq)

        # For checking full_seq
        # print(full_seq)

    # Building whole strings to write to file. "/" char separates chains.
    whole_gapped_str = "/".join(whole_gapped)
    whole_full_str = "/".join(whole_full)

    # Writing to file
    # Remember sequences have to end with the * character
    with open(output_align_file, "a+") as file:
        # Writing structure/gapped section
        file.write(">P1;" + structure.id + "\n")
        file.write("structureX:" + structure.id + ":FIRST:@ END:@" + 5 * ":." + "\n")
        for line in textwrap.wrap(
            whole_gapped_str + "*", width=75, break_on_hyphens=False
        ):
            file.write("%s\n" % line)
        # Writing full sequence section
        file.write(">P1;" + structure.id + "_fill\n")
        file.write("sequence:" + structure.id + ":FIRST:@ END:@" + 5 * ":." + "\n")
        for line in textwrap.wrap(
            whole_full_str + "*", width=75, break_on_hyphens=False
        ):
            file.write("%s\n" % line)


def complete_residues(pdbcode, align_file="protein.ali", loop_ref=False):
    """
    Function that completes residues based on an alignment file using MODELLER software.
    :param pdbcode: PDB code identifier of the structure with missing residues.
    :param align_file: Path to the align-formatted file with gaps as missing residues.
    :param loop_ref: (optional) Boolean for specifying loop refinement, doesn't always work.
    :return:
    """
    # Get the sequence of the coded PDB file, and write to an alignment file
    e = environ()
    m = model(e, file=pdbcode)
    aln = alignment(e)
    aln.append_model(m, align_codes=pdbcode)
    aln.write(file=pdbcode + ".seq")

    # Completing residues
    log.verbose()
    env = environ()

    # directories for input atom files (local dir)
    env.io.atom_files_directory = ["."]

    if loop_ref is True:
        # For loop refinement - Doesn't always work
        a = loopmodel(env, alnfile=align_file, knowns=pdbcode, sequence=code + "_fill")
        a.loop.starting_model = 1
        a.loop.ending_model = 2
        a.loop.md_level = refine.fast
    else:
        a = automodel(
            env, alnfile=align_file, knowns=pdbcode, sequence=pdbcode + "_fill"
        )

    a.starting_model = 1
    a.ending_model = 1

    a.make()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Builds alignment file for MODELLER.")
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        help="Input PDB file with missing residues.",
        required=True,
    )
    parser.add_argument(
        "--pdbcode",
        type=str,
        help="PDB code identifier for the structure.",
        required=True,
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="(Optional) Output align (.ali) file for completing residues.",
        required=False,
        default="protein.ali",
    )
    parser.add_argument(
        "--loop_ref",
        "-lr",
        dest="loop_ref",
        help="(Optional) Enables loop refinement. " "Does not always work.",
        required=False,
        action="store_true",
    )
    parser.set_defaults(feature=True)

    args = parser.parse_args()

    # the PDB file has to be in the same directory, copying and using the code as name.
    pdb_path = args.input
    code = args.pdbcode
    print("input: " + os.path.abspath(args.input))
    print("cwd: " + os.getcwd() + code + ".pdb")
    if os.path.abspath(args.input) == os.getcwd() + "/" + code + ".pdb":
        raise ValueError(
            "Input file comes from current working directory and cannot have "
            + code
            + ".pdb already as a name. Please change the name (or location) of input PDB file."
        )
    else:
        shutil.copy(pdb_path, "./" + code + ".pdb")

    build_align_file(args.input, args.pdbcode, output_align_file=args.output)
    complete_residues(code, align_file=args.output, loop_ref=args.loop_ref)
