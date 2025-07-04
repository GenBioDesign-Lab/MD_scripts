#!/bin/bash
#SBATCH --job-name=Packmol
#SBATCH --output=packmol.log
#SBATCH --error=packmol_error.log
#SBATCH --partition=cpus
##SBATCH --gres=gpu:a6000:1
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

ml load ambertools/2025.0
ml load vmd/2.0.0a4
ml load anaconda3/2024.10
conda activate mdpipe

cd /data01/genbiolab/mdanh/data/MD_scripts/working/test/scripts #where the namd config files are

./solvate_packmol_ionize.sh /data01/genbiolab/mdanh/data/MD_scripts/working/test/scripts/TIP3.pdb