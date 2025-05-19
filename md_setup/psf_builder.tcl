# PSF Builder Script

# Create working directory; remove old output files
#mkdir -p output
#rm -f output/output_protein.pdb

# Remove HETATM lines from input.pdb
#grep -v 'HETATM' 6pti.pdb > output/output_protein.pdb
#grep -v 'HOH' 6pti.pdb > output/output_water.pdb

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
#change HIS to HSD
pdbalias residue HIS HSD

#Build protein segment
# Segment name max 4 characters
segment 1BYI {
    pdb output/output_protein.pdb
}

#patch protein segment: patch [list] <patch residue name> <segid:resid>
#patch DISU BPTI:5 BPTI:55

#Read protein coordinates from PDB file
pdbalias atom ILE CD1 CD
coordpdb output/output_protein.pdb 1BYI

# Build missing atoms and guess missing coordinates
guesscoord

# Write the PSF and PDB files
writepsf 1BYI.psf
writepdb 1BYI.pdb

exit
