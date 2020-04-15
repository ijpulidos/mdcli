#!/bin/bash

# This scripts performs the PCA computation and projection, using GROMACS tools, for the
# first 3 components for all of the domains, using only the helices residues.
# It assumes there's already a pdb with the helices called helices.pdb
# Expects command line argument for the directory where the helices.pdb file is located

# Exit on error
set -e

# Understanding command line arguments
# inspired by: https://stackoverflow.com/a/7069755
while test $# -gt 0; do
  case "$1" in
    -h|--help)
      echo "do_pca - attempt to compute PCA analysis for MD."
      echo "Check output PCA dir in output directory for results."
      echo " "
      echo "do_pca [arguments]"
      echo " "
      echo "arguments:"
      echo "-h, --help                 show brief help"
      echo "-f, --file=PATH            path to helices.pdb file"
      echo "-s, --structure=PATH       path to structure file (.pdb or .gro format)"
      echo "-t, --trajectory=PATH      path to trajectory file (.xtc format)"
      echo "-o, --output=PATH          path to output directory"
      exit 0
      ;;
    -f)
      shift
      if test $# -gt 0; then
        export HELICES=$1
      else
        echo "no helices.pdb file specified"
        exit 1
      fi
      shift
      ;;
    --file*)
      export HELICES=`echo $1 | sed -e 's/^[^=]*=//g'`
      shift
      ;;
    -s)
      shift
      if test $# -gt 0; then
        export STRUCTURE=$1
      else
        echo "no structure file specified"
        exit 1
      fi
      shift
      ;;
    --structure*)
      export STRUCTURE=`echo $1 | sed -e 's/^[^=]*=//g'`
      shift
      ;;
    -t)
      shift
      if test $# -gt 0; then
        export TRAJECTORY=$1
      else
        echo "no trajectory file specified"
        exit 1
      fi
      shift
      ;;
    --trajectory*)
      export TRAJECTORY=`echo $1 | sed -e 's/^[^=]*=//g'`
      shift
      ;;
    -o)
      shift
      if test $# -gt 0; then
        export OUTPUT=$1
      else
        echo "no output directory specified"
        exit 1
      fi
      shift
      ;;
    --output*)
      export OUTPUT=`echo $1 | sed -e 's/^[^=]*=//g'`
      shift
      ;;
    *)
      break
      ;;
  esac
done

# Basepath is current working directory
basepath=`pwd`

# Make output dir if doesn't exist
mkdir -p ${OUTPUT}
cd ${OUTPUT}
# get indices from pdb file passed as argument for script
# Note: Sometimes the f**king column number changes... 
cat ${HELICES} | awk '{print $5}' > helices.txt
# make index file with helices indices
rs=($(cat helices.txt)); { echo r ${rs[*]} ; echo q; } | gmx make_ndx -f ${STRUCTURE} -o helices.ndx
# make index file with backbone of the helices
printf "4 & 18\nq\n" | gmx make_ndx -f ${STRUCTURE} -n helices.ndx -o bb_helices.ndx 
# compute covariance matrix
printf "19\n19\nq\n" | gmx covar -s ${STRUCTURE} -f ${TRAJECTORY} -n bb_helices.ndx
# Analysis of first 3 eigenvalues
printf "19\n19\nq\n" | gmx anaeig -s ${STRUCTURE} -f ${TRAJECTORY} -n bb_helices.ndx -b 25000 -3d -first 1 -last 3 -extr -nframes 20
cd ${basepath}