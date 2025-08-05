package require psfgen
psfgen_logfile "{log_file}"

# Load topology files
{topology_commands}

# Define CNT segment
segment CNT {{
    pdb {cnt_pdb}
}}
coordpdb {cnt_pdb} CNT

# Define membrane segment  
segment MEMB {{
    pdb {membrane_pdb}
}}
coordpdb {membrane_pdb} MEMB

# Build structure
guesscoord
writepsf {output_psf}
writepdb {output_pdb}

psfgen_logfile close
quit 