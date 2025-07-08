import MDAnalysis as mda
import numpy as np
import argparse
import warnings

# Usage: python center_in_membrane.py membrane.pdb protein_or_nanotube.pdb -o centered_structure.pdb

# Suppress MDAnalysis PDB format warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.coordinates.PDB")

# We need to rotate the CNT structure 90 degrees around the x-axis to make it
# parallel to the membrane (as the script creates the CNT lying in the x-y plane)

def rotate_structure(structure, x_angle, y_angle, z_angle):
    """
    Rotate a structure around x, y, and z axes by the specified angles in degrees.
    Rotations are applied in the order: X -> Y -> Z
    """
    if x_angle == 0 and y_angle == 0 and z_angle == 0:
        return
    
    # Convert to radians
    x_rad = np.radians(x_angle)
    y_rad = np.radians(y_angle)
    z_rad = np.radians(z_angle)
    
    # Rotation matrices for each axis
    # X-axis rotation
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(x_rad), -np.sin(x_rad)],
        [0, np.sin(x_rad), np.cos(x_rad)]
    ])
    
    # Y-axis rotation
    Ry = np.array([
        [np.cos(y_rad), 0, np.sin(y_rad)],
        [0, 1, 0],
        [-np.sin(y_rad), 0, np.cos(y_rad)]
    ])
    
    # Z-axis rotation
    Rz = np.array([
        [np.cos(z_rad), -np.sin(z_rad), 0],
        [np.sin(z_rad), np.cos(z_rad), 0],
        [0, 0, 1]
    ])
    
    # Combined rotation matrix (order: X -> Y -> Z)
    rotation_matrix = np.dot(Rz, np.dot(Ry, Rx))
    
    # Get center of mass for rotation
    center_of_mass = structure.atoms.center_of_mass()
    
    # Translate to origin, rotate, then translate back
    structure.atoms.positions -= center_of_mass
    structure.atoms.positions = np.dot(structure.atoms.positions, rotation_matrix.T)
    structure.atoms.positions += center_of_mass
    
    print(f"Structure rotated: X={x_angle}°, Y={y_angle}°, Z={z_angle}°")

def center_structure_in_membrane(membrane_file, structure_file, output_file, membrane_selection="name P", rotation_angles=None):
    """
    Center a protein or nanotube structure in the middle of a membrane.
    Optionally rotate the structure around x, y, z axes before centering.
    """
    
    # Load membrane and structure
    membrane = mda.Universe(membrane_file)
    structure = mda.Universe(structure_file)
    
    # Rotate structure if requested
    if rotation_angles:
        rotate_structure(structure, rotation_angles[0], rotation_angles[1], rotation_angles[2])
    
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
    parser.add_argument('-r', '--rotate', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       default=[90, 0, 0],
                       help='Rotate structure by X Y Z degrees around respective axes (default: 90 0 0, e.g., -r 0 0 90)')
    
    args = parser.parse_args()
    center_structure_in_membrane(args.membrane, args.structure, args.output, args.selection, args.rotate)

if __name__ == "__main__":
    main() 