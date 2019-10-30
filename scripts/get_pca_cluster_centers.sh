#!/bin/bash

# Script to get the closest structure to cluster centers from PCA.
# It uses R and k-means algorithm for clustering.

# Understanding command line arguments
# inspired by: https://stackoverflow.com/a/7069755
while test $# -gt 0; do
  case "$1" in
    -h|--help)
      echo "get_pca_cluster_centers - Perform clustering of PCA data and get centers of clusters."
      echo "Check output files in output directory for results."
      echo "It requires R language for performing the clustering."
      echo " "
      echo "get_pca_cluster_centers [arguments]"
      echo " "
      echo "arguments:"
      echo "-h, --help                 show brief help"
      echo "-f, --file=PATH            path to 3dproj.pdb file"
      echo "-k, --clusters=PATH        number of clusters"
      echo "-o, --output=PATH          path to output directory (must exist)"
      exit 0
      ;;
    -f)
      shift
      if test $# -gt 0; then
        export PROJ=$1
      else
        echo "no 3dproj.pdb file specified"
        exit 1
      fi
      shift
      ;;
    --file*)
      export PROJ=`echo $1 | sed -e 's/^[^=]*=//g'`
      shift
      ;;
    -k)
      shift
      if test $# -gt 0; then
        export CLUSTERS=$1
      else
        echo "number of clusters not specified"
        exit 1
      fi
      shift
      ;;
    --clusters*)
      export CLUSTERS=`echo $1 | sed -e 's/^[^=]*=//g'`
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


# 1. Get XYZ coordinate information from projected PCA coords PDB file and write to dat file
grep ATOM ${PROJ} | awk ' {printf("%f %f %f\n", $6,$7,$8) } ' > ${OUTPUT}/3dproj.dat

# 2. Clustering with K-means 
# 2.2. Kmeans cluster script to run in R - make script file
cat> ${OUTPUT}/kmeans_clustering.R <<EOF
# read data
read.table ("${OUTPUT}/3dproj.dat") -> data

# kmeans calculation (10 clusters)
fit <- kmeans(data,${CLUSTERS},10)

# geometrical centers of the clusters (to variable gc)
aggregate(data,by=list(fit\$cluster), FUN=mean) -> gc

# write geometrical centers to a file (only cols 3-4 and without row and column names)
write.table (gc[,c(2,3,4)], "${OUTPUT}/cluster_gcs.dat" , row.names= FALSE , col.names = FALSE)
EOF

# 2.3. Run the R script file - output in cluster_gcs.dat
R --no-save < ${OUTPUT}/kmeans_clustering.R >& out

# Remove copied 3dproj.dat file
rm -f 3dproj.dat

# 2.4. pdb with the geometrical centers of the cluster
awk ' 
# Read the pdb coordinates
{
X=$1; 
Y=$2; 
Z=$3; 
} 

# print the coordinates
{
printf("%-6s%5d %4s %3s %1s%4d%4s%8.3f%8.3f%8.3f%6.2f%6.2f\n", "ATOM ", NR, "GOC" , "GOC", " ", NR ," ",X ,
 Y , Z , 0 , 0) ;
}
' ${OUTPUT}/cluster_gcs.dat > ${OUTPUT}/clusters_gcs.pdb
