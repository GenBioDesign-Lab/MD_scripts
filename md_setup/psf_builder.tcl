# PSF Builder Script
#Usage: vmd -dispdev text -e psf_builder.tcl <.pdb file>

# Get input PDB file from command-line argument
set pdb_file [lindex $argv 0]
# Extract base name (remove path and extension)
set pdb_base [file rootname [file tail $pdb_file]]
# Use first 4 characters for segment name
set seg_name [string range $pdb_base 0 3]

# Create working directory; remove old output files
mkdir -p output
set out_protein "output/${pdb_base}_protein.pdb"
set out_hetatm "output/${pdb_base}_HETATM.pdb"
rm -f $out_protein $out_hetatm

# Remove HETATM lines from input PDB
exec grep -v HETATM $pdb_file > $out_protein
exec grep -v ATOM $pdb_file > $out_hetatm

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

# Build protein segment
segment $seg_name {
    pdb $out_protein
}

#patch protein segment: patch [list] <patch residue name> <segid:resid>
#patch DISU $seg_name:5 $seg_name:55

#Read protein coordinates from PDB file
pdbalias atom ILE CD1 CD
coordpdb $out_protein $seg_name

# Build missing atoms and guess missing coordinates
guesscoord

# Write the PSF and PDB files
writepsf ${seg_name}.psf
writepdb ${seg_name}.pdb

exit
