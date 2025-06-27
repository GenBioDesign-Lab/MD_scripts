#!/usr/bin/env vmd
# DPPC membrane restraints generator
# Usage: vmd -dispdev text -e create_dppc_ref_files.tcl

# ====== CONFIGURATION - EDIT THESE PATHS ======
# Input files - MODIFY THESE PATHS FOR YOUR SYSTEM
set input_psf "/data01/genbiolab/mdanh/data/MD_scripts/test/DPPC_CNT/ionized/cnt_mem_ionized/cnt_mem_ionized.psf"
set input_pdb "/data01/genbiolab/mdanh/data/MD_scripts/test/DPPC_CNT/ionized/cnt_mem_ionized/cnt_mem_ionized.pdb"

# Output directory - MODIFY THIS PATH AS NEEDED
set output_dir "restraints"

# ====== DO NOT EDIT BELOW THIS LINE ======

# Check if input files exist
if {![file exists $input_psf]} {
    puts "Error: PSF file not found: $input_psf"
    puts "Please edit the input_psf variable in the script to point to your PSF file."
    exit 1
}

if {![file exists $input_pdb]} {
    puts "Error: PDB file not found: $input_pdb"
    puts "Please edit the input_pdb variable in the script to point to your PDB file."
    exit 1
}

# Create output directory if it doesn't exist
file mkdir $output_dir

# Output file names
set dppc_upper_ref "$output_dir/dppc_head_upper.ref"
set dppc_lower_ref "$output_dir/dppc_head_lower.ref"
set colvar_file "$output_dir/dppc_restraints.colvar"
set log_file "$output_dir/restraints_generation.log"

# Open log file
set log_fh [open $log_file w]

# Function to write to both terminal and log
proc log_message {message {terminal_only 0}} {
    global log_fh
    if {$terminal_only == 0} {
        puts $log_fh $message
        flush $log_fh
    }
    puts $message
}

# Function to write only to log file
proc log_only {message} {
    global log_fh
    puts $log_fh $message
    flush $log_fh
}

# Start logging
puts "DPPC Membrane Restraints Generator - Starting..."
puts "Log file: $log_file"

log_only "==================================================================================="
log_only "DPPC Membrane Restraints Generator"
log_only "==================================================================================="
log_only "Input PSF: $input_psf"
log_only "Input PDB: $input_pdb"
log_only "Output directory: $output_dir"
log_only "DPPC upper leaflet file: $dppc_upper_ref"
log_only "DPPC lower leaflet file: $dppc_lower_ref"
log_only "Colvar file: $colvar_file"
log_only "Log file: $log_file"
log_only "==================================================================================="

# ====== LOAD STRUCTURE ======
puts "Loading structure..."
log_only "Loading structure..."
mol new $input_psf
mol addfile $input_pdb
log_only "Structure loaded successfully."

# ====== FIND MEMBRANE CENTER ======
puts "Determining membrane center..."
log_only "Determining membrane center..."

# Try to find phosphate atoms to determine membrane center
set phosphate_atoms [atomselect top "resname DPPC and name P"]
if {[$phosphate_atoms num] > 0} {
    set z_coords [$phosphate_atoms get z]
    set membrane_z_center [expr {([tcl::mathfunc::max {*}$z_coords] + [tcl::mathfunc::min {*}$z_coords]) / 2.0}]
    log_only "Membrane center determined from phosphate atoms: Z = [format %.3f $membrane_z_center] Å"
} else {
    # Fallback: use all DPPC atoms
    set all_dppc [atomselect top "resname DPPC"]
    if {[$all_dppc num] > 0} {
        set z_coords [$all_dppc get z]
        set membrane_z_center [expr {([tcl::mathfunc::max {*}$z_coords] + [tcl::mathfunc::min {*}$z_coords]) / 2.0}]
        log_only "Membrane center determined from all DPPC atoms: Z = [format %.3f $membrane_z_center] Å"
        $all_dppc delete
    } else {
        puts "Error: No DPPC molecules found in the system!"
        log_only "Error: No DPPC molecules found in the system!"
        close $log_fh
        exit 1
    }
}
$phosphate_atoms delete

# ====== SELECT DPPC HEAD GROUP ATOMS ======
puts "Processing DPPC membrane..."
log_only "Selecting DPPC head group atoms..."

# Define head group atoms for DPPC (choline, phosphate, and glycerol backbone)
set headgroup_atoms "name N or name C13 or name C14 or name C15 or name C12 or name C11 or name P or name O13 or name O14 or name O12 or name O11 or name C1 or name HA or name HB or name C2 or name HS or name O21 or name C21 or name O22 or name C3 or name HX or name HY or name O31 or name C31 or name O32"

# Select all DPPC head group atoms
set all_heads [atomselect top "resname DPPC and ($headgroup_atoms)"]
log_only "Total DPPC head group atoms found: [$all_heads num]"

if {[$all_heads num] == 0} {
    puts "Error: No DPPC head group atoms found!"
    log_only "Error: No DPPC head group atoms found!"
    log_only "Make sure your structure contains DPPC lipids and uses standard CHARMM atom names."
    close $log_fh
    exit 1
}

# ====== SEPARATE UPPER AND LOWER LEAFLETS ======
log_only "Separating DPPC leaflets..."

# Select upper leaflet head groups (Z > membrane center)
set upper_heads [atomselect top "resname DPPC and ($headgroup_atoms) and z > $membrane_z_center"]
set upper_count [$upper_heads num]

# Select lower leaflet head groups (Z < membrane center)  
set lower_heads [atomselect top "resname DPPC and ($headgroup_atoms) and z < $membrane_z_center"]
set lower_count [$lower_heads num]

log_only "Upper leaflet head group atoms: $upper_count"
log_only "Lower leaflet head group atoms: $lower_count"

# ====== WRITE DPPC UPPER LEAFLET REFERENCE FILE ======
log_only "Writing DPPC upper leaflet reference file..."

# Set B-factors: 1.0 for upper leaflet atoms, 0.0 for lower leaflet atoms
$all_heads set beta 0.0
$upper_heads set beta 1.0

# Write the reference file for upper leaflet
$all_heads writepdb $dppc_upper_ref

log_only "DPPC upper leaflet reference file written: $dppc_upper_ref"
log_only "  - Upper leaflet atoms (B-factor = 1.0): $upper_count"
log_only "  - Lower leaflet atoms (B-factor = 0.0): $lower_count"

# ====== WRITE DPPC LOWER LEAFLET REFERENCE FILE ======
log_only "Writing DPPC lower leaflet reference file..."

# Set B-factors: 0.0 for upper leaflet atoms, 1.0 for lower leaflet atoms
$all_heads set beta 0.0
$lower_heads set beta 1.0

# Write the reference file for lower leaflet
$all_heads writepdb $dppc_lower_ref

log_only "DPPC lower leaflet reference file written: $dppc_lower_ref"
log_only "  - Upper leaflet atoms (B-factor = 0.0): $upper_count"
log_only "  - Lower leaflet atoms (B-factor = 1.0): $lower_count"

# ====== GENERATE COLVAR FILE ======
puts "Generating colvar configuration..."
log_only "Creating colvar configuration..."

set colvar_fh [open $colvar_file w]

puts $colvar_fh "# DPPC Membrane Restraints Colvar Configuration"
puts $colvar_fh "# Generated by create_dppc_ref_files.tcl"
puts $colvar_fh ""
puts $colvar_fh "Colvarstrajfrequency    100"
puts $colvar_fh "Colvarsrestartfrequency 100"
puts $colvar_fh ""

puts $colvar_fh "# ====== DPPC MEMBRANE RESTRAINTS ======"
puts $colvar_fh ""

# DPPC upper leaflet
puts $colvar_fh "# DPPC upper leaflet Z-position"
puts $colvar_fh "colvar {"
puts $colvar_fh "    name dppc_head_upper"
puts $colvar_fh "    distanceZ {"
puts $colvar_fh "        ref {"
puts $colvar_fh "            dummyAtom ( 0.000, 0.000, 0.000 )"
puts $colvar_fh "        }"
puts $colvar_fh "        main {"
puts $colvar_fh "            atomsFile      dppc_head_upper.ref"
puts $colvar_fh "            atomsCol       B"
puts $colvar_fh "            atomsColValue  1.0"
puts $colvar_fh "        }"
puts $colvar_fh "    }"
puts $colvar_fh "}"
puts $colvar_fh ""

# Calculate membrane positions with standard offset
set upper_center [expr {$membrane_z_center + 15.0}]
set lower_center [expr {$membrane_z_center - 15.0}]

puts $colvar_fh "harmonic {"
puts $colvar_fh "    colvars dppc_head_upper"
puts $colvar_fh "    centers [format %.3f $upper_center]"
puts $colvar_fh "    forceConstant \$membrane_fc"
puts $colvar_fh "}"
puts $colvar_fh ""

# DPPC lower leaflet
puts $colvar_fh "# DPPC lower leaflet Z-position"
puts $colvar_fh "colvar {"
puts $colvar_fh "    name dppc_head_lower"
puts $colvar_fh "    distanceZ {"
puts $colvar_fh "        ref {"
puts $colvar_fh "            dummyAtom ( 0.000, 0.000, 0.000 )"
puts $colvar_fh "        }"
puts $colvar_fh "        main {"
puts $colvar_fh "            atomsFile      dppc_head_lower.ref"
puts $colvar_fh "            atomsCol       B"
puts $colvar_fh "            atomsColValue  1.0"
puts $colvar_fh "        }"
puts $colvar_fh "    }"
puts $colvar_fh "}"
puts $colvar_fh ""
puts $colvar_fh "harmonic {"
puts $colvar_fh "    colvars dppc_head_lower"
puts $colvar_fh "    centers [format %.3f $lower_center]"
puts $colvar_fh "    forceConstant \$membrane_fc"
puts $colvar_fh "}"

close $colvar_fh
log_only "Colvar file written: $colvar_file"

# ====== WRITE SUMMARY TO LOG ======
log_only ""
log_only "==================================================================================="
log_only "SUMMARY:"
log_only "Membrane center: Z = [format %.3f $membrane_z_center] Å"
log_only "Upper leaflet center: Z = [format %.3f $upper_center] Å"
log_only "Lower leaflet center: Z = [format %.3f $lower_center] Å"
log_only "Total DPPC head group atoms: [$all_heads num]"
log_only "Upper leaflet atoms: $upper_count"
log_only "Lower leaflet atoms: $lower_count"
log_only ""
log_only "Generated files:"
log_only "  - DPPC upper leaflet: $dppc_upper_ref"
log_only "  - DPPC lower leaflet: $dppc_lower_ref"
log_only "  - Colvar config: $colvar_file"
log_only "  - Log file: $log_file"
log_only ""
log_only "Usage in NAMD:"
log_only "  colvars on"
log_only "  exec sed -e \"s/\\\$membrane_fc/10.0/g\" $colvar_file > restraints.colvars"
log_only "  colvarsConfig restraints.colvars"
log_only "==================================================================================="

# Clean up
$all_heads delete
$upper_heads delete
$lower_heads delete

# Final message
puts "DPPC membrane restraints generation completed successfully!"
puts "Check $log_file for detailed information."

# Close log file
close $log_fh

# Exit VMD
quit 