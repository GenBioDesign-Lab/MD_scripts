#Usage: vmd -e combine.tcl
package require psfgen

# Reset PSFGen
resetpsf

# Read the cnt and membrane complex
readpsf $inp_psf
coordpdb $inp_pdb

# Merge the water molecules
readpsf $out_dir/${out_dir}_water.psf
coordpdb $out_dir/${out_dir}_water.pdb

# Write the merged system
writepsf $out_dir/${out_dir}_cnt_mem_water.psf
writepdb $out_dir/${out_dir}_cnt_mem_water.pdb
