#!/bin/bash
#SBATCH --job-name=NAMD_CPU
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=cpus
##SBATCH --gres=gpu:a6000:1
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

ml load namd/3.0.1-CPU
cd /data01/genbiolab/mdanh/data/MD_scripts/test/2mlt/2mlt_ionized #where the namd config files are
namd3 +p48  minimization/01_Minimization.namd > minimization/01_Minimization.log
#namd3 +p48  equilibration/02_Equilibration.namd > equilibration/02_Equilibration.log
#namd3 +p48  production/03_Production_npt.namd > production/03_Production.log