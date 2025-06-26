#!/bin/bash
#Put this script into the namd folder
#VMD in HPC can't call the package catdcd in VMD, so download and specify the catdcd manually: https://www.ks.uiuc.edu/Development/MDTools/catdcd/files/catdcd-4.0b.tar.gz

# Set the number of equilibration phases
NUM_PHASES=6

# Build the list of input files
INPUT_FILES=""
for i in $(seq 1 $NUM_PHASES); do
    INPUT_FILES="$INPUT_FILES eq$i/final_eq$i.dcd"
done

# Run catdcd with all input files
/data01/genbiolab/mdanh/data/CatDCD/LINUXAMD64/bin/catdcd4.0/catdcd -o eq_all.dcd $INPUT_FILES