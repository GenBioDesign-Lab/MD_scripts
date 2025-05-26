#Usage: vmd -e combine.tcl
package require psfgen

# Reset PSFGen
resetpsf

# Read the first structure (CNT)
readpsf 2mlt_mem.psf
coordpdb 2mlt_mem.pdb

# Read the second structure (membrane) and merge it
readpsf mem6060.psf
coordpdb mem6060.pdb

# Write the merged system
writepsf 2mlt_mem.psf
writepdb 2mlt_mem.pdb

exit