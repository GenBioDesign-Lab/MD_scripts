# PSF Builder Script
#Usage: vmd -dispdev text -e psf_builder.tcl
# Must use vmd 2.0.4 because its supported new version of topology files


# Require packages
package require psfgen

# Set topology directory path
set input_path "/data01/genbiolab/mdanh/data/CHARMMFF/toppar"

# Define list of topology files to use
set topology_files {
    "/top_all36_prot.rtf"
    "/top_all36_na.rtf"
    "/top_all36_carb.rtf"
    "/top_all36_lipid.rtf"
    "/top_all36_cgenff.rtf"
}

# Load all topology files
foreach topo_file $topology_files {
    topology "$input_path$topo_file"
}
#open log file
psfgen_logfile "structure_preparation.log"

#change HIS to HSD
pdbalias residue HIS HSD

#change U to URA (for RNA)
#pdbalias residue U URA

#Build protein segment
# Segment name max 4 characters
segment CNT {
    pdb cnt_mem.pdb
}

#patch protein segment: patch [list] <patch residue name> <segid:resid>
#patch DISU BPTI:5 BPTI:55

pdbalias atom ILE CD1 CD

#Read protein coordinates from PDB file
coordpdb cnt_mem.pdb CNT

segment MEMB {
    pdb cnt_mem.pdb
}

coordpdb cnt_mem.pdb MEMB

#Build RNA segment
#segment RNA {
#    pdb 9bwz_aligned.pdb
#}

#coordpdb test/ligand.pdb RNA

# Build missing atoms and guess missing coordinates
guesscoord

# Write the PSF and PDB files
writepsf cnt_mem.psf
writepdb cnt_mem_psf.pdb

psfgen_logfile close
exit
