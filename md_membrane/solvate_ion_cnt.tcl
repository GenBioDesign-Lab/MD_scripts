# Usage: vmd -e solvate_ion_cnt.tcl
# User-defined input files and output directory
set inp_psf "../clean/cnt_mem.psf"
set inp_pdb "../clean/cnt_mem.pdb"
set out_dir "tmp"
set ionized_dir "cnt_mem_ionized"

# Create temporary directory to store intermediate files
file mkdir $out_dir
file mkdir $ionized_dir

# Load required plugins
package require solvate

# 1. Load starting structure
mol new $inp_psf
mol addfile $inp_pdb

# 2. SOLVATE
#  -t 10   : add 10 Å padding in every direction
#  -x,-y,-z, +x,+y,+z : add padding in specific directions
#  -o ...  : prefix for output files
solvate $inp_psf $inp_pdb -t 10 -o $out_dir/${out_dir}_solv

# The call creates: <out_dir>/<out_dir>_solv.psf / <out_dir>/<out_dir>_solv.pdb
# Keep only the solvated molecule in memory
mol delete top

# ----‑‑‑ 3. REMOVE WATER OUTSIDE Z-RANGE ----------------------------------
# Remove water molecules outside z-range (-15 to 15)
mol new      $out_dir/${out_dir}_solv.psf
mol addfile  $out_dir/${out_dir}_solv.pdb 

# Select all atoms EXCEPT water outside z-range (atoms to keep)
set getwater [atomselect top "water and not (z>-15 and z < 15)"]

# Write structure with water outside z-range removed
$getwater writepdb  $out_dir/${out_dir}_water.pdb
$getwater writepsf $out_dir/${out_dir}_water.psf

# Clean up selections
$getwater delete
mol delete top

# ----‑‑‑ 4. Combine water and membrane ------------------------------------------

source ../combine.tcl

# ----‑‑‑ 4. REMOVE WATER NEAR CNT ------------------------------------------
# Remove water molecules within 5 Å of the carbon nanotube
mol new      $out_dir/${out_dir}_cnt_mem_water.psf
mol addfile  $out_dir/${out_dir}_cnt_mem_water.pdb

# Select all atoms EXCEPT water near CNT (atoms to keep)
set delcntwater [atomselect top "not (same residue as (water and within 10 of segname TUBE))"]

# Write structure with water removed
$delcntwater writepdb $ionized_dir/cnt_mem_fixed_wat.pdb
$delcntwater writepsf $ionized_dir/cnt_mem_fixed_wat.psf

# Clean up selections
$delcntwater delete
mol delete top

# ----‑‑‑ 5. AUTOIONIZE ------------------------------------------------------
#  -sc       : neutralize and set salt concentration (mol/L)
#  -cation 0.15        : add additional Na⁺/Cl⁻ to reach 0.15 M
#  -o ...            : prefix for ionized output
package require autoionize
autoionize -psf $ionized_dir/cnt_mem_fixed_wat.psf \
           -pdb $ionized_dir/cnt_mem_fixed_wat.pdb \
           -sc 0.15 \
           -cation SOD \
           -o $ionized_dir/cnt_mem_ionized

# Creates: <ionized_dir>/<ionized_dir>_ionized.psf / <ionized_dir>/<ionized_dir>_ionized.pdb

exit