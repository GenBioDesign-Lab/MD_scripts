#!/bin/bash
#SBATCH --job-name=CNT_MEM_test
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=gpus
#SBATCH --gres=gpu:a6000:1
#SBATCH --nodelist=gpu03
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

ml load namd/3.0.1-GPU
cd /data01/genbiolab/mdanh/data/test_membrane/cnt_mem_ionized #where the namd config files are
namd3 +p24 +idlepoll +setcpuaffinity +devices 0 minimization/01_Minimization.namd > minimization/01_Minimization.log
#namd3 +p24 +idlepoll +setcpuaffinity +devices 0 equilibration/02_Equilibration.namd > equilibration/02_Equilibration.log
#namd3 +p24 +idlepoll +setcpuaffinity +devices 0 production/03_Production_npt.namd > production/03_Production.log