import MDAnalysis as mda
import numpy as np
import argparse
import sys
import warnings

#Usage: python remove_lipid.py nanotube.pdb membrane.psf membrane.pdb -o cleaned_membrane.pdb -r 1.0 -t 0.0

# Suppress MDAnalysis PDB format warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.coordinates.PDB")

def remove_lipids(nanotube_file, membrane_psf, membrane_pdb, output_file, radius_factor=1.0, overlap_threshold=0.0):
    mem = mda.Universe(membrane_psf, membrane_pdb)
    tube = mda.Universe(nanotube_file)

    tube_coords = tube.atoms.positions
    x_center = tube_coords[:,0].mean()
    y_center = tube_coords[:,1].mean()
    z_min = tube_coords[:,2].min()
    z_max = tube_coords[:,2].max()
    
    r_xy = np.sqrt((tube_coords[:,0] - x_center)**2 + (tube_coords[:,1] - y_center)**2)
    radius = r_xy.max() * radius_factor  # Apply scaling factor to radius
    
    print(f"Using tube dimensions: center=({x_center:.2f}, {y_center:.2f}), z=({z_min:.2f}, {z_max:.2f}), radius={radius:.2f}")

    lipid_atoms = mem.select_atoms("resname POPC or resname DOPE")
    
    # First select by z-range only
    z_range_atoms = lipid_atoms.select_atoms(f"prop z > {z_min} and prop z < {z_max}")
    
    # Then calculate distance from center for these atoms
    coords = z_range_atoms.positions
    dist_sq = (coords[:,0] - x_center)**2 + (coords[:,1] - y_center)**2
    
    # Get atoms within radius
    inside_atoms = z_range_atoms[dist_sq < radius**2]
    
    # Count how many atoms of each residue are inside
    residue_counts = {}
    for atom in inside_atoms:
        res_id = atom.residue.resid
        if res_id not in residue_counts:
            residue_counts[res_id] = {"inside": 0, "total": len(atom.residue.atoms)}
        residue_counts[res_id]["inside"] += 1
    
    # Only mark residues with more than threshold% atoms inside
    bad_resids = []
    for res_id, counts in residue_counts.items():
        if counts["inside"] / counts["total"] >= overlap_threshold:
            bad_resids.append(str(res_id))
    
    if bad_resids:
        bad_residues = mem.select_atoms(f"resid {' '.join(bad_resids)}").residues
    else:
        bad_residues = mem.atoms.residues[0:0]  # Empty residue group
    
    print(f"{len(bad_residues)} lipids marked for removal")

    keep_atoms = mem.atoms - bad_residues.atoms
    keep_atoms.write(output_file)
    print(f"Successfully wrote {len(keep_atoms)} atoms to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Remove lipids inside nanotube')
    parser.add_argument('nanotube', help='Nanotube PDB')
    parser.add_argument('membrane_psf', help='Membrane PSF')
    parser.add_argument('membrane_pdb', help='Membrane PDB')
    parser.add_argument('-o', '--output', default='cleaned_membrane.pdb')
    parser.add_argument('-r', '--radius_factor', type=float, default=1.0,
                        help='Factor to multiply nanotube radius (smaller = remove fewer lipids, default: 1.0)')
    parser.add_argument('-t', '--threshold', type=float, default=0.0,
                        help='Fraction of residue atoms that must be inside to remove (default: 0.0)')
    
    args = parser.parse_args()
    remove_lipids(args.nanotube, args.membrane_psf, args.membrane_pdb, args.output, 
                 args.radius_factor, args.threshold)

if __name__ == "__main__":
    main()
