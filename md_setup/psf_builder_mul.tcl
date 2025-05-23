# PSF Builder Script
#Usage: vmd -dispdev text -e psf_builder.tcl
# Must use vmd 2.0.4 because its supported new version of topology files


# Require packages
package require psfgen

# Set topology directory path
set input_path "/data01/genbiolab/jiyun/data/toppar"

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
pdbalias residue U URA

# Fix atom naming for ILE
pdbalias atom ILE CD1 CD

#Build protein segment
# Segment name max 4 characters
segment PROA {
    pdb 9bwz_aligned.pdb
    first NTER
    last CTER
    chain A
}
coordpdb 9bwz_aligned.pdb PROA

segment PROB {
    pdb 9bwz_aligned.pdb
    first NTER
    last CTER
    chain B
}
coordpdb 9bwz_aligned.pdb PROB

segment PROC {
    pdb 9bwz_aligned.pdb
    first NTER
    last CTER
    chain C
}
coordpdb 9bwz_aligned.pdb PROC

segment PROD {
    pdb 9bwz_aligned.pdb
    first NTER
    last CTER
    chain D
}
coordpdb 9bwz_aligned.pdb PROD

#patch protein segment: patch [list] <patch residue name> <segid:resid>
#patch DISU BPTI:5 BPTI:55

#Build RNA segment
segment RNA {
    pdb 9bwz_aligned.pdb
    first 5TER
    last 3TER
    chain E
}

coordpdb 9bwz_aligned.pdb RNA

# Build missing atoms and guess missing coordinates
guesscoord

# Write the PSF and PDB files
writepsf 9bwz.psf
writepdb 9bwz_psf.pdb

psfgen_logfile close
exit
