import MDAnalysis as mda
import numpy as np
import argparse
import warnings

# Usage: python center_cnt.py cnt_structure.pdb -o centered_cnt.pdb -r 90 0 0

# Suppress MDAnalysis PDB format warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.coordinates.PDB")

def rotate_structure(structure, x_angle, y_angle, z_angle):
    """
    Rotate a structure around x, y, and z axes by the specified angles in degrees.
    Rotations are applied in the order: X -> Y -> Z
    """
    if x_angle == 0 and y_angle == 0 and z_angle == 0:
        print("No rotation applied (all angles are 0)")
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

def center_cnt_at_origin(cnt_file, output_file, rotation_angles=None):
    """
    Center a CNT structure at (0,0,0).
    Optionally rotate the structure around x, y, z axes before centering.
    """
    
    # Load CNT structure
    cnt = mda.Universe(cnt_file)
    
    print(f"Loaded CNT structure with {len(cnt.atoms)} atoms")
    
    # Get initial center of mass
    initial_center = cnt.atoms.center_of_mass()
    print(f"Initial CNT center: ({initial_center[0]:.3f}, {initial_center[1]:.3f}, {initial_center[2]:.3f})")
    
    # Rotate structure if requested (rotation happens around current center of mass)
    if rotation_angles and any(angle != 0 for angle in rotation_angles):
        rotate_structure(cnt, rotation_angles[0], rotation_angles[1], rotation_angles[2])
    
    # Get center of mass after rotation
    current_center = cnt.atoms.center_of_mass()
    
    # Calculate translation to move center to (0,0,0)
    translation = np.array([0.0, 0.0, 0.0]) - current_center
    
    # Apply translation
    cnt.atoms.positions += translation
    
    # Verify final center
    final_center = cnt.atoms.center_of_mass()
    
    # Print results
    print(f"Translation applied: ({translation[0]:.3f}, {translation[1]:.3f}, {translation[2]:.3f})")
    print(f"Final CNT center: ({final_center[0]:.3f}, {final_center[1]:.3f}, {final_center[2]:.3f})")
    
    # Write output
    cnt.atoms.write(output_file)
    print(f"Successfully wrote centered CNT with {len(cnt.atoms)} atoms to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Center a CNT structure at (0,0,0) with optional rotation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python center_cnt.py cnt.pdb                    # Default: 90° X rotation, then center at origin
  python center_cnt.py cnt.pdb -r 0 0 0          # No rotation, just center at origin
  python center_cnt.py cnt.pdb -r 90 45 0        # 90° X, 45° Y rotation, then center
  python center_cnt.py cnt.pdb -r 0 0 90 -o rotated_cnt.pdb  # 90° Z rotation, custom output
        """
    )
    
    parser.add_argument('cnt_file', help='CNT PDB file to center and rotate')
    parser.add_argument('-o', '--output', default='centered_cnt.pdb',
                       help='Output PDB file (default: centered_cnt.pdb)')
    parser.add_argument('-r', '--rotate', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       default=[90, 0, 0],
                       help='Rotate CNT by X Y Z degrees around respective axes before centering (default: 90 0 0)')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("CNT Centering and Rotation Script")
    print("=" * 60)
    print(f"Input file: {args.cnt_file}")
    print(f"Output file: {args.output}")
    print(f"Rotation: X={args.rotate[0]}°, Y={args.rotate[1]}°, Z={args.rotate[2]}°")
    print("-" * 60)
    
    center_cnt_at_origin(args.cnt_file, args.output, args.rotate)
    
    print("-" * 60)
    print("Process completed!")

if __name__ == "__main__":
    main() 