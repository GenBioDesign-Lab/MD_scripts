# PSF Builder Script
#Usage: vmd -dispdev text -e psf_builder.tcl
# Must use vmd 2.0.4 because its supported new version of topology files


# Require packages
package require psfgen

# Set topology directory path
set input_path "/data01/genbiolab/mdanh/data/CHARMMFF/toppar"

topology "/data01/genbiolab/mdanh/data/CHARMMFF/cnt_ff/cnt_6_40_oh_CHARMM.rtf"

# Define list of topology files to use
set topology_files {
    "/top_all36_prot.rtf"
    "/top_all36_na.rtf"
    "/top_all36_carb.rtf"
    "/top_all36_cgenff.rtf"
}

# Load all topology files
foreach topo_file $topology_files {
    topology "$input_path$topo_file"
}
#open log file
psfgen_logfile "cnt_preparation.log"

#Build TUBE segment
# Segment name max 4 characters
segment TUBE {
    pdb centered_cnt.pdb
}

#Read protein coordinates from PDB file
coordpdb centered_cnt.pdb TUBE

guesscoord

# Write the PSF and PDB files
writepsf cnt_psf.psf
writepdb cnt_psf.pdb

psfgen_logfile close

resetpsf

#Open log file
psfgen_logfile "cnt_mem_preparation.log"
# Read the first structure (CNT)
readpsf cnt_psf.psf
coordpdb cnt_psf.pdb

# Read the second structure (membrane) and merge it
readpsf mem_8080_nowater.psf
coordpdb cleaned_mem.pdb

# Write the merged system
writepsf cnt_mem.psf
writepdb cnt_mem.pdb

psfgen_logfile close

exit
