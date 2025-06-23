#!/bin/bash
#Put this script into the namd folder
#VMD in HPC can't call the package catdcd in VMD, so download and specify the catdcd manually: https://www.ks.uiuc.edu/Development/MDTools/catdcd/files/catdcd-4.0b.tar.gz

for i in {1..6}; do
    /data01/genbiolab/mdanh/data/CatDCD/LINUXAMD64/bin/catdcd4.0/catdcd -o eq_all.dcd eq$i/final_eq$i.dcd
done