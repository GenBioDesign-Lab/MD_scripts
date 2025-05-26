import MDAnalysis as mda
import numpy as np
import argparse
import warnings

# Usage: python center_in_membrane.py membrane.pdb protein_or_nanotube.pdb -o centered_structure.pdb

# Suppress MDAnalysis PDB format warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.coordinates.PDB")

def center_structure_in_membrane(membrane_file, structure_file, output_file, membrane_selection="name P"):
    """
    Center a protein or nanotube structure in the middle of a membrane.
    """
    
    # Load membrane and structure
    membrane = mda.Universe(membrane_file)
    structure = mda.Universe(structure_file)
    
    # Get membrane center
    membrane_atoms = membrane.select_atoms(membrane_selection)
    if len(membrane_atoms) == 0:
        print(f"Warning: No atoms found with selection '{membrane_selection}'. Using all atoms.")
        membrane_atoms = membrane.atoms
    
    mem_coords = membrane_atoms.positions
    mem_center = np.array([
        mem_coords[:, 0].mean(),
        mem_coords[:, 1].mean(),
        (mem_coords[:, 2].min() + mem_coords[:, 2].max()) / 2.0
    ])
    
    # Get current structure center
    struct_coords = structure.atoms.positions
    struct_center = struct_coords.mean(axis=0)
    
    # Calculate and apply translation
    translation = mem_center - struct_center
    structure.atoms.positions += translation
    
    # Print results
    print(f"Membrane center: ({mem_center[0]:.2f}, {mem_center[1]:.2f}, {mem_center[2]:.2f})")
    print(f"Structure moved by: ({translation[0]:.2f}, {translation[1]:.2f}, {translation[2]:.2f})")
    
    # Write output
    structure.atoms.write(output_file)
    print(f"Successfully wrote {len(structure.atoms)} atoms to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Center a protein or nanotube in a membrane')
    parser.add_argument('membrane', help='Membrane PDB file')
    parser.add_argument('structure', help='Protein/nanotube PDB file to center')
    parser.add_argument('-o', '--output', default='centered_structure.pdb',
                       help='Output PDB file (default: centered_structure.pdb)')
    parser.add_argument('-s', '--selection', default='name P',
                       help='Selection for membrane centering atoms (default: "name P")')
    
    args = parser.parse_args()
    center_structure_in_membrane(args.membrane, args.structure, args.output, args.selection)

if __name__ == "__main__":
    main() 