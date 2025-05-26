package require psfgen

# Reset PSFGen
resetpsf

# Read the first structure (CNT)
readpsf cnt666.psf
coordpdb cnt666_centered.pdb

# Read the second structure (membrane) and merge it
readpsf mem70.psf
coordpdb mem70_cleaned.pdb

# Write the merged system
writepsf cnt_mem.psf
writepdb cnt_mem.pdb

exit