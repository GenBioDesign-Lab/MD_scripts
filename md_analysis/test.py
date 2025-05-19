import MDAnalysis as mda
import warnings
# Load simulation results with a single line
u = mda.Universe('example_data/9BWZ_ionized.psf','example_data/final_eq.dcd')

# Select atoms
ag = u.select_atoms("segid 9BWZ and not (name H* or type OW)")

frames = [2,3,3,1]
selection  = u.trajectory[frames] 
positions = []
for ts in selection:
    print(ts)
    positions.append(ts.positions)

positions # each element of the list will be different. 

# Atom data made available as Numpy arrays
#ag.positions
#ag.velocities
#ag.forces

# Iterate through trajectories
#for ts in u.trajectory:
#    print(ag.center_of_mass())