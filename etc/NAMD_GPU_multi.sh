#!/bin/bash
#SBATCH --job-name=MEMB
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=gpus
#SBATCH --gres=gpu:a6000:1
#SBATCH --nodelist=gpu03
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

ml load namd/3.0.1-GPU
cd /data01/genbiolab/mdanh/data/MD_scripts/working/test/cnt6_50/namd #where the namd config files are

# Run minimization
namd3 +p24 +idlepoll +setcpuaffinity +devices 0 minimization/01_Minimization.namd > minimization/01_Minimization.log

# Array of equilibration folders
eq_folders=("eq1" "eq2" "eq3" "eq4" "eq5" "eq6")

# Base directory (current namd directory)
base_dir=$(pwd)

# Loop through each equilibration folder
for i in "${!eq_folders[@]}"; do
    eq_folder="${eq_folders[$i]}"
    eq_num=$((i + 1))

    # Change to the equilibration directory
    cd "${base_dir}/${eq_folder}"
    
    # Find the NAMD configuration file in this directory
    namd_config=$(ls *.namd 2>/dev/null | head -1)
    
    if [ -z "$namd_config" ]; then
        echo "ERROR: No NAMD configuration file found in ${eq_folder}"
        exit 1
    fi
        
    # Run NAMD simulation
    namd3 +p24 +idlepoll +setcpuaffinity +devices 0 "$namd_config" > "${namd_config%.namd}.log" 2>&1
done

# Return to base directory
cd "$base_dir"

# Run production
namd3 +p24 +idlepoll +setcpuaffinity +devices 0 production/03_Production_npt.namd > production/03_Production.log