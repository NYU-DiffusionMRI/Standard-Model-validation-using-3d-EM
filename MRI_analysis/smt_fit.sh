#!/bin/bash
# Fit Spherical Mean Technique
# usage:
#      bash smt_fit.sh pathFolderIn pathFolderOut inputFilesBasename
# By Ricardo Coronado-Leija

# print date, command
now=$(date +%Y-%m-%d-%H-%M-%S)
echo $now
# print command
echo $*

# ================================================================================================== #
# Inputs

# Check number of inputs
if [ "$#" -ne 3 ]; then
echo "Illegal number of parameters (should be 3):"
exit 1
fi

# Set input/output variables
folderIn=${1}
folderOut=${2}
BasenameIn=${3}
# ================================================================================================== #
# Check input values and create output folders

# If input folder do not exists
if [ ! -d ${folderIn} ]; then
echo "Input folder " ${folderIn} " not found."
exit 2
fi

# Create output folder in case of needed
if [ ! -d ${folderOut} ]; then
echo "Creating Folder " ${folderOut}
mkdir ${folderOut} -p
fi
# ================================================================================================== #

SMT_NUM_THREADS=16

# SMT fit
fitmcmicro --bvals ${folderIn}/${BasenameIn}.bval    --bvecs ${folderIn}/${BasenameIn}.bvec  --maxdiff 2e-3 --b0 ${folderIn}/${BasenameIn}.nii.gz   ${folderOut}/smt.nii.gz --mask ${folderIn}/${BasenameIn}_dwi_mask.nii.gz


 
 






