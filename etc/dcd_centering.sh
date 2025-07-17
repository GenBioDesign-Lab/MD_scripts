#!/bin/bash
#SBATCH --job-name=dcdcentering
#SBATCH --output=dcdcentering_output.log
#SBATCH --error=dcdcentering_error.log
#SBATCH --partition=cpus
##SBATCH --gres=gpu:a6000:1
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

ml load anaconda3/2024.10
conda activate mdpipe

python /data01/genbiolab/mdanh/data/MD_scripts/md_analysis/dcdcentering.py \
       -structure /data01/genbiolab/mdanh/data/MD_scripts/working/test/cnt_mem_ionized.psf \
       -traj /data01/genbiolab/mdanh/data/MD_scripts/working/test/final_NPT.dcd \
       -sel "segid CNT1 or segid CNT2 or segid CNT3" \
       -center geometry \
       -o /data01/genbiolab/mdanh/data/MD_scripts/working/test/CNT_centered_wrapped.xtc