import MDAnalysis as mda
import numpy as np
import argparse
import warnings

# Usage: python remove_lipid.py nanotube.pdb membrane.pdb -o cleaned_membrane.pdb -c complex.pdb -r 1.1 -b 2.0

warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.coordinates.PDB")

def remove_lipids_and_combine(nanotube_file, membrane_pdb, output_file, complex_file=None, radius_factor=1.1, buffer_distance=2.0):
    # Load structures
    mem = mda.Universe(membrane_pdb)
    tube = mda.Universe(nanotube_file)

    # Get nanotube dimensions
    tube_coords = tube.atoms.positions
    x_center, y_center = tube_coords[:,0].mean(), tube_coords[:,1].mean()
    z_min, z_max = tube_coords[:,2].min(), tube_coords[:,2].max()
    
    r_xy = np.sqrt((tube_coords[:,0] - x_center)**2 + (tube_coords[:,1] - y_center)**2)
    # Add buffer distance to the exclusion radius
    radius = r_xy.max() * radius_factor + buffer_distance
    
    print(f"Tube: center=({x_center:.2f}, {y_center:.2f}), z=({z_min:.2f}, {z_max:.2f})")
    print(f"CNT radius: {r_xy.max():.2f} Å, exclusion radius: {radius:.2f} Å (buffer: {buffer_distance:.1f} Å)")

    # Select atoms to process (try lipids and water first, then all atoms)
    try:
        atoms_to_check = mem.select_atoms("resname POPC or resname POPE or resname TIP3 or resname SOL or resname WAT")
        if len(atoms_to_check) == 0:
            atoms_to_check = mem.atoms
            print("Processing all membrane atoms")
        else:
            print(f"Processing lipids and water molecules: {len(atoms_to_check)} atoms")
    except:
        atoms_to_check = mem.atoms

    # Find atoms in nanotube region
    z_range_atoms = atoms_to_check.select_atoms(f"prop z > {z_min} and prop z < {z_max}")
    coords = z_range_atoms.positions
    dist_sq = (coords[:,0] - x_center)**2 + (coords[:,1] - y_center)**2
    inside_atoms = z_range_atoms[dist_sq < radius**2]
    
    # Remove any residue with ANY atom in the exclusion zone (zero tolerance)
    overlapping_residues = set()
    for atom in inside_atoms:
        overlapping_residues.add(atom.residue.resid)
    
    # Mark residues for removal
    bad_resids = [str(res_id) for res_id in overlapping_residues]
    
    if bad_resids:
        bad_residues = mem.select_atoms(f"resid {' '.join(bad_resids)}").residues
    else:
        bad_residues = mem.atoms.residues[0:0]
    
    print(f"{len(bad_residues)} residues removed")

    # Save cleaned membrane only if specified
    keep_atoms = mem.atoms - bad_residues.atoms
    if output_file:
        keep_atoms.write(output_file)
        print(f"Cleaned membrane: {len(keep_atoms)} atoms → {output_file}")
    
    # Create complex if requested
    if complex_file:
        combined = mda.Merge(tube.atoms, keep_atoms)
        combined.atoms.write(complex_file)
        print(f"Complex: {len(combined.atoms)} atoms → {complex_file}")

def main():
    parser = argparse.ArgumentParser(description='Remove overlapping residues and create complex structure with no overlap')
    parser.add_argument('nanotube', help='Nanotube PDB file')
    parser.add_argument('membrane_pdb', help='Membrane PDB file')
    parser.add_argument('-o', '--output', help='Cleaned membrane output (optional)')
    parser.add_argument('-c', '--complex', help='Complex structure output (optional)')
    parser.add_argument('-r', '--radius_factor', type=float, default=1.1, help='Radius scaling factor (default: 1.1)')
    parser.add_argument('-b', '--buffer_distance', type=float, default=2.0, help='Buffer distance in Å (default: 2.0)')
    
    args = parser.parse_args()
    remove_lipids_and_combine(args.nanotube, args.membrane_pdb, args.output, 
                             args.complex, args.radius_factor, args.buffer_distance)

if __name__ == "__main__":
    main()
