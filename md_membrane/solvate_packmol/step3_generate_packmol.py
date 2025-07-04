#!/usr/bin/env python3
"""
Step 3: Generate Packmol input file for placing water molecules in upper and lower leaflets
Usage: python3 step3_generate_packmol.py <waters_per_leaflet> <thickness> <water_pdb_file>
"""

import sys
import os

def parse_pdb_for_phosphorus_atoms(pdb_file):
    """Parse PDB file to find phosphorus atoms and their coordinates"""
    p_atoms = []
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_name = line[12:16].strip()
                    if atom_name == 'P':
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        p_atoms.append((x, y, z))
    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading PDB file {pdb_file}: {e}")
        return []
    
    return p_atoms

def get_system_dimensions_and_leaflet_bounds(solute_pdb):
    """Get system dimensions and determine leaflet boundaries from P atoms"""
    p_atoms = parse_pdb_for_phosphorus_atoms(solute_pdb)
    
    if not p_atoms:
        print("Warning: No phosphorus atoms found. Using default boundaries.")
        return {
            'x_min': -50, 'x_max': 50,
            'y_min': -50, 'y_max': 50,
            'upper_leaflet_start': 5,
            'lower_leaflet_end': -5
        }
    
    # Extract coordinates
    x_coords = [atom[0] for atom in p_atoms]
    y_coords = [atom[1] for atom in p_atoms]
    z_coords = [atom[2] for atom in p_atoms]
    
    # Find boundaries - limit to exact membrane boundaries (no padding)
    x_min, x_max = min(x_coords), max(x_coords)  # Exact membrane X boundaries
    y_min, y_max = min(y_coords), max(y_coords)  # Exact membrane Y boundaries
    
    # Use P atoms to determine leaflet boundaries
    upper_leaflet_start = max(z_coords)  # Maximum Z of P atoms
    lower_leaflet_end = min(z_coords)    # Minimum Z of P atoms
    
    print(f"Found {len(p_atoms)} phosphorus atoms")
    print(f"Membrane X range: {x_min:.2f} to {x_max:.2f}")
    print(f"Membrane Y range: {y_min:.2f} to {y_max:.2f}")
    print(f"P atom Z range: {lower_leaflet_end:.2f} to {upper_leaflet_start:.2f}")
    print(f"Upper leaflet will start at Z = {upper_leaflet_start:.2f}")
    print(f"Lower leaflet will end at Z = {lower_leaflet_end:.2f}")
    print("Water placement will be limited to exact membrane X,Y boundaries")
    
    return {
        'x_min': x_min, 'x_max': x_max,
        'y_min': y_min, 'y_max': y_max,
        'upper_leaflet_start': upper_leaflet_start,
        'lower_leaflet_end': lower_leaflet_end
    }

def generate_packmol_input(waters_per_leaflet, thickness, water_pdb_path):
    """Generate Packmol input file for water placement in leaflets"""
    
    # Get system dimensions and leaflet boundaries based on P atoms
    solute_pdb = "solvated/solute_only.pdb"
    dims = get_system_dimensions_and_leaflet_bounds(solute_pdb)
    
    # Calculate leaflet positions based on P atom positions
    # Upper leaflet: starts at max Z of P atoms, extends upward by thickness
    upper_z_min = dims['upper_leaflet_start']
    upper_z_max = upper_z_min + thickness
    
    # Lower leaflet: ends at min Z of P atoms, extends downward by thickness
    lower_z_max = dims['lower_leaflet_end']
    lower_z_min = lower_z_max - thickness
    
    packmol_input = f"""# Packmol input file for placing water in upper and lower leaflets
# Generated automatically by step3_generate_packmol.py
# Leaflet boundaries determined from phosphorus atom positions

tolerance 2.0
filetype pdb
output packmol_output/system_with_leaflets.pdb

# Add the solute (membrane/CNT system)
structure solvated/solute_only.pdb
  number 1
  fixed 0 0 0 0 0 0
end structure

# Upper leaflet water molecules (above maximum P atom Z position)
structure {water_pdb_path}
  number {waters_per_leaflet}
  inside box {dims['x_min']} {dims['y_min']} {upper_z_min} {dims['x_max']} {dims['y_max']} {upper_z_max}
end structure

# Lower leaflet water molecules (below minimum P atom Z position)
structure {water_pdb_path}
  number {waters_per_leaflet}
  inside box {dims['x_min']} {dims['y_min']} {lower_z_min} {dims['x_max']} {dims['y_max']} {lower_z_max}
end structure
"""

    # Write the Packmol input file
    with open('packmol_leaflets.inp', 'w') as f:
        f.write(packmol_input)
    
    print(f"Packmol input file generated: packmol_leaflets.inp")
    print(f"Using water template: {water_pdb_path}")
    print(f"Waters per leaflet: {waters_per_leaflet}")
    print(f"Leaflet thickness: {thickness} Å")
    print(f"Upper leaflet Z range: {upper_z_min:.2f} to {upper_z_max:.2f}")
    print(f"Lower leaflet Z range: {lower_z_min:.2f} to {lower_z_max:.2f}")



def main():
    if len(sys.argv) != 4:
        print("Usage: python3 step3_generate_packmol.py <waters_per_leaflet> <thickness> <water_pdb_file>")
        sys.exit(1)
    
    try:
        waters_per_leaflet = int(sys.argv[1])
        thickness = float(sys.argv[2])
        water_pdb_path = sys.argv[3]
    except (ValueError, IndexError):
        print("Error: waters_per_leaflet must be an integer, thickness must be a number")
        sys.exit(1)
    
    # Check if the water PDB file exists
    if not os.path.exists(water_pdb_path):
        print(f"Error: Water template file not found: {water_pdb_path}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs('packmol_output', exist_ok=True)
    
    # Generate Packmol input
    generate_packmol_input(waters_per_leaflet, thickness, water_pdb_path)

if __name__ == '__main__':
    main() 