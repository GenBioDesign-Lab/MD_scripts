#!/usr/bin/env python3
"""
Remove POPC/POPE lipid molecules within specified distance from protein.
Usage: python remove_lipid_near_protein.py <protein> <membrane> -o <output> [-d distance]
"""

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import argparse
import sys

def remove_lipids_near_protein(protein_file, membrane_file, output_file, cutoff_distance=2.5):
    """Remove lipid molecules within cutoff_distance of protein."""
    
    # Load structures
    print(f"Loading structures...")
    protein_u = mda.Universe(protein_file)
    membrane_u = mda.Universe(membrane_file)
    
    # Select POPC and POPE lipids
    lipids = membrane_u.select_atoms("resname POPC or resname POPE")
    if len(lipids) == 0:
        print("No POPC/POPE lipids found!")
        sys.exit(1)
    
    print(f"Found {len(lipids)} lipid atoms in {len(lipids.residues)} molecules")
    
    # Calculate distances
    print(f"Calculating distances (cutoff: {cutoff_distance} Å)...")
    dist_matrix = distance_array(protein_u.atoms.positions, lipids.positions)
    
    # Find close lipid atoms and their residues
    close_indices = np.where(dist_matrix <= cutoff_distance)[1]
    close_residues = set(lipids[close_indices].residues)
    
    if not close_residues:
        print("No lipids within cutoff distance")
        membrane_u.atoms.write(output_file)
        return
    
    # Remove close residues
    # Start with an empty AtomGroup from the membrane universe and
    # accumulate atoms from each residue using ``sum``.
    atoms_to_remove = sum(
        (res.atoms for res in close_residues),
        membrane_u.atoms[0:0]
    )
    atoms_to_keep = membrane_u.atoms - atoms_to_remove
    
    print(f"Removing {len(close_residues)} lipid molecules ({len(atoms_to_remove)} atoms)")
    print(f"Keeping {len(atoms_to_keep)}/{len(membrane_u.atoms)} atoms")
    
    # Write output
    atoms_to_keep.write(output_file)
    print(f"Wrote truncated membrane to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Remove POPC/POPE lipids near protein')
    parser.add_argument('protein', help='Protein structure file')
    parser.add_argument('membrane', help='Membrane structure file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-d', '--distance', type=float, default=3.0, help='Distance cutoff (Å)')
    
    args = parser.parse_args()
    remove_lipids_near_protein(args.protein, args.membrane, args.output, args.distance)

if __name__ == "__main__":
    main() 