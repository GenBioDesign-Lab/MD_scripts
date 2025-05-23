#!/bin/bash
#SBATCH --job-name=24CPU2GPU_gpu01
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=gpus
#SBATCH --gres=gpu:8000:2
#SBATCH --nodelist=gpu01
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

ml load namd/3.0.1-GPU
cd /data01/genbiolab/mdanh/data/MD_pipeline/12CPU1GPU_gpu01 #where the namd config files are
#namd3 +p48 +devices 0 01_Minimization.namd > minimization.log
namd3 +p24 +setcpuaffinity +devices 0,1 02_Equilibration.namd > equilibration.log
#namd3 +p12 +idlepoll +setcpuaffinity +devices 0 03_Production_npt.namd > production_npt.log

