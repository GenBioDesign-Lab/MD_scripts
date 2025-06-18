#PSF builder for homo oligomers from pre-split PDB files
# Usage: vmd -dispdev text -e psf_homo.tcl

# ====== CONFIGURATION ======
# Input directory containing split PDB files
set input_dir "/data01/genbiolab/mdanh/data/MD_scripts/test/cnt_oh_mem70/clean/cleaned_mem"

# Topology directory path
set topdir "/data01/genbiolab/mdanh/data/CHARMMFF/toppar"

# Define list of topology files to use
set topology_files {
    "/top_all36_prot.rtf"
    "/top_all36_na.rtf"
    "/top_all36_carb.rtf"
    "/top_all36_lipid.rtf"
    "/top_all36_cgenff.rtf"
}

# Output files
set output_psf "rebuilt_mem.psf"
set output_pdb "rebuilt_mem.pdb"

# ====== READ INPUT FILES ======

# Check if input directory exists
if {![file exists $input_dir]} {
    puts "Error: Input directory $input_dir does not exist!"
    exit 1
}

# Get all PDB files with chain L pattern from input directory
set pdb_files [glob -nocomplain $input_dir/*_L_*.pdb]
if {[llength $pdb_files] == 0} {
    puts "Error: No PDB files matching pattern *_L_*.pdb found in $input_dir!"
    exit 1
}

puts "Found PDB files:"
foreach file $pdb_files {
    puts "  $file"
}

# ====== PSF GENERATION ======
puts "Building PSF structure..."
package require psfgen
resetpsf

# Load topology files
foreach topo_file $topology_files {
    topology "$topdir$topo_file"
}

# Process each PDB file as membrane component
set seg_num 1

foreach pdb_file $pdb_files {
    set filename [file tail $pdb_file]
    set basename [file rootname $filename]
    
    puts "Processing $filename"
    
    # Create membrane segment
    set segname "MEM$seg_num"
    puts "Building membrane segment $segname for $filename"
    segment $segname {
        pdb $pdb_file
    }
    coordpdb $pdb_file $segname
    incr seg_num
}

# Finalize structure
puts "Building missing atoms..."
guesscoord

# Write output files
puts "Writing output files..."
writepsf $output_psf
writepdb $output_pdb

puts "PSF generation completed successfully!"
puts "Output files: $output_psf, $output_pdb"

exit