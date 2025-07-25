#!/bin/bash
#SBATCH --job-name=DPPC_mem
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=gpus
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --nodelist=gpu04
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

ml load namd/3.0.1-GPU
cd /data01/genbiolab/mdanh/data/simulation/MEM_DPPC/new/DPPC_140/namd #where the namd config files are

# Loop through equilibration steps 6.1 to 6.6
for step in 6.1 6.2 6.3 6.4 6.5 6.6; do
    namd3 +p24 +idlepoll +setcpuaffinity +devices 0 step${step}_equilibration.inp > step${step}_equilibration.log
done

#namd3 +p24 +idlepoll +setcpuaffinity +devices 0 step7_production.inp > step7_production.log