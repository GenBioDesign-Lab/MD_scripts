# Ultra-simple PSF builder for homo oligomers
# Usage: vmd -dispdev text -e psf_homo.tcl

# ====== CONFIGURATION ======
# Input PDB file
set input_pdb "9bwz_aligned.pdb"

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
set output_psf "9bwz.psf"
set output_pdb "9bwz_psf.pdb"

# Splitted chains directory
set chains_dir "splitted_chains"

# ====== CHAIN DETECTION AND SPLITTING ======

# Create directory for splitted chains if it doesn't exist
if {![file exists $chains_dir]} {
    file mkdir $chains_dir
    puts "Created directory: $chains_dir"
}

set mol [mol new $input_pdb]
set chains {}
foreach c {A B C D E F G H} {
    set sel [atomselect $mol "chain $c"]
    if {[$sel num] > 0} {
        lappend chains $c
        $sel writepdb $chains_dir/chain_$c.pdb
        puts "Found chain $c -> $chains_dir/chain_$c.pdb"
    }
    $sel delete
}
mol delete $mol
puts "Chains found: $chains"

# ====== PSF GENERATION ======
puts "Building PSF structure..."
package require psfgen
resetpsf

# Load topology files
foreach topo_file $topology_files {
    topology "$topdir$topo_file"
}

# Set residue and atom aliases
pdbalias residue HIS HSD
pdbalias residue U URA
pdbalias atom ILE CD1 CD

# Build segments
set protein_num 1
foreach c $chains {
    if {[string match {[A-D]} $c]} {
        # Protein chains (A-D)
        set segname "PRO$protein_num"
        puts "Building protein segment $segname for chain $c"
        segment $segname {
            pdb $chains_dir/chain_$c.pdb
            first NTER
            last CTER
        }
        coordpdb $chains_dir/chain_$c.pdb $segname
        incr protein_num
    } else {
        # RNA/DNA chains (E+)
        set segname "RNA$c"
        puts "Building RNA segment $segname for chain $c"
        segment $segname {
            pdb $chains_dir/chain_$c.pdb
            first 5TER
            last 3TER
        }
        coordpdb $chains_dir/chain_$c.pdb $segname
    }
}

# Finalize structure
puts "Building missing atoms..."
guesscoord

# Write output files
puts "Writing output files..."
writepsf $output_psf
writepdb $output_pdb

puts "Complete! Output files:"
puts "  - $output_psf"
puts "  - $output_pdb"

exit