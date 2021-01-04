#!/bin/env python

# Script that performas a CONAN analysis in a more friendly way, hopefully.

import argparse
from subprocess import Popen, PIPE, STDOUT


def build_conan_inp(
    trajectory,
    coordinate,
    begin,
    end,
    resolution,
    truncation,
    clusters,
    gnus_path,
    pdbfile,
    out_file="conan.inp",
    dt=None,
    trunc_inter=None,
    trunc_dr=None,
    domains="domains.txt",
):
    """
    Function that builds a CONAN .inp input file and saves it where specified.
    """

    # Compute nlevel with resolution and truncation
    nlevel = int(truncation / resolution + 1)
    # Compute trunc values if not given
    if trunc_inter:
        trunc_inter_ = trunc_inter
    else:
        trunc_inter_ = 0.5 * truncation
    if trunc_dr:
        trunc_dr_ = trunc_dr
    else:
        trunc_dr_ = truncation

    # write to file
    with open(out_file, "w") as inp_file:
        inp_file.write(f"TRAJ {trajectory}\n")
        inp_file.write(f"COORD {coordinate}\n")
        inp_file.write(f"NLEVEL {nlevel}\n")
        # take user's dt if given
        if dt:
            inp_file.write(f"DT {dt}\n")
        inp_file.write(f"TRUNC {truncation}\n")
        inp_file.write(f"TRUNC_INTER {trunc_inter_}\n")
        inp_file.write(f"TRUNC_DR {trunc_dr_}\n")
        # Check for begin and end times - else default values
        if begin:
            inp_file.write(f"BEGIN {begin}\n")
        if end:
            inp_file.write(f"END {end}\n")
        inp_file.write(f"DR_MODE init\n")
        inp_file.write(f"GNUS_PATH {gnus_path}\n")
        inp_file.write(f"RUN_MDMAT yes\n")
        # Don't need the matrices
        inp_file.write(f"MATRICES no\n")
        # Clean matrices to save up space
        inp_file.write(f"CLEAN_MATRICES yes\n")
        # Don't make a movie
        inp_file.write(f"MAKE_MOVIE no\n")
        # inp_file.write(f"DOMAINS {domains}\n")
        inp_file.write(f"COORD_PDB {pdbfile}\n")
        inp_file.write(f"K_TRAJ_CLUSTERS {clusters}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automatic CONAN analysis.")
    parser.add_argument(
        "--trajectory",
        "-t",
        type=str,
        help="Input trajectory file in .xtc format.",
        required=True,
    )
    parser.add_argument(
        "--coordinate",
        "-c",
        type=str,
        help="Input coordinate file (.tpr).",
        required=True,
    )
    parser.add_argument(
        "--gnus-path", type=str, help="Path to GNUPlot scripts.", required=True
    )
    parser.add_argument(
        "--pdbfile",
        "-pdb",
        type=str,
        help="Path to PDB file to obtain sequence.",
        required=True,
    )
    parser.add_argument(
        "--begin",
        "-b",
        type=int,
        help="Beginning time to do analysis, in ps.",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--end",
        "-e",
        type=int,
        help="Ending time to do analysis, in ps.",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--resolution",
        type=float,
        help="Spatial resolution in nm.",
        required=False,
        default=0.001,
    )
    parser.add_argument(
        "--truncation",
        "-trunc",
        type=float,
        help="Distance truncation for computation of contats, in nm.",
        required=False,
        default=1.0,
    )
    parser.add_argument(
        "--clusters",
        type=str,
        help="Cluster specification range. e.g. 1-5 for 1,2,3,4,5.",
        required=False,
        default="1-5",
    )
    parser.add_argument(
        "--inpfile",
        type=str,
        help="Path to output CONAN parameter .inp file to be created.",
        required=False,
        default="conan.inp",
    )
    parser.add_argument(
        "--conan",
        type=str,
        help="Path to CONAN executable. Default 'conan.py'.",
        required=False,
        default="conan.py",
    )
    parser.add_argument(
        "--logfile",
        type=str,
        help="Path to logfile where to store all output from runing CONAN. Prints to stdout if \
            not given.",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--dt", type=int, help="Time step in ps. Integer.", required=False, default=500
    )
    args = parser.parse_args()

    # build CONAN input file
    build_conan_inp(
        trajectory=args.trajectory,
        coordinate=args.coordinate,
        begin=args.begin,
        end=args.end,
        resolution=args.resolution,
        truncation=args.truncation,
        clusters=args.clusters,
        gnus_path=args.gnus_path,
        pdbfile=args.pdbfile,
        out_file=args.inpfile,
        dt=args.dt,
    )

    # Run CONAN subprocess
    conan_arg_list = ["python", args.conan, args.inpfile]
    if args.logfile:
        out = open(args.logfile, "w")
    else:
        out = PIPE
    process = Popen(conan_arg_list, stdin=PIPE, stdout=out, stderr=STDOUT)
    output = process.communicate(input=b"1\n")  # Wait for finish and get ouput
    print(output)
