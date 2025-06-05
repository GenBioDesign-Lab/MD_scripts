import MDAnalysis as mda
import numpy as np
import argparse
import warnings
import logging
import os
from datetime import datetime

# Usage: python remove_lipid.py nanotube.pdb membrane.pdb -o cleaned_membrane.pdb -c complex.pdb -r 1.1 -b 2.0

warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.coordinates.PDB")

def setup_logging(log_file):
    """
    Set up logging to write detailed information to a log file.
    
    Parameters:
    log_file (str): Path to the log file
    """
    # Create directory if it doesn't exist
    log_dir = os.path.dirname(log_file)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()  # Also print to console
        ]
    )
    
    # Log script start
    logging.info("="*60)
    logging.info("LIPID REMOVAL AND STRUCTURE COMBINATION LOG")
    logging.info("="*60)
    logging.info(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
def log_and_print(message):
    """
    Print message to console and log it.
    
    Parameters:
    message (str): Message to log and print
    """
    print(message)
    logging.info(message)

def remove_lipids_and_combine(nanotube_file, membrane_pdb, output_file, complex_file=None, radius_factor=1.0, buffer_distance=0.8):
    logging.info(f"Input parameters:")
    logging.info(f"  Nanotube file: {nanotube_file}")
    logging.info(f"  Membrane file: {membrane_pdb}")
    logging.info(f"  Output file: {output_file}")
    logging.info(f"  Complex file: {complex_file}")
    logging.info(f"  Radius factor: {radius_factor}")
    logging.info(f"  Buffer distance: {buffer_distance} Å")
    logging.info("-" * 40)
    
    # Load structures
    logging.info("Loading structures...")
    mem = mda.Universe(membrane_pdb)
    tube = mda.Universe(nanotube_file)
    
    logging.info(f"Membrane atoms: {len(mem.atoms)}")
    logging.info(f"Nanotube atoms: {len(tube.atoms)}")

    # Get nanotube dimensions
    tube_coords = tube.atoms.positions
    x_center, y_center = tube_coords[:,0].mean(), tube_coords[:,1].mean()
    z_min, z_max = tube_coords[:,2].min(), tube_coords[:,2].max()
    
    r_xy = np.sqrt((tube_coords[:,0] - x_center)**2 + (tube_coords[:,1] - y_center)**2)
    # Add buffer distance to the exclusion radius
    radius = r_xy.max() * radius_factor + buffer_distance
    
    message = f"Tube: center=({x_center:.2f}, {y_center:.2f}), z=({z_min:.2f}, {z_max:.2f})"
    log_and_print(message)
    message = f"CNT radius: {r_xy.max():.2f} Å, exclusion radius: {radius:.2f} Å (buffer: {buffer_distance:.1f} Å)"
    log_and_print(message)

    # Select atoms to process (try lipids and water first, then all atoms)
    try:
        atoms_to_check = mem.select_atoms("resname POPC or resname POPE or resname TIP3 or resname SOL or resname WAT")
        if len(atoms_to_check) == 0:
            atoms_to_check = mem.atoms
            message = "Processing all membrane atoms"
            log_and_print(message)
        else:
            message = f"Processing lipids and water molecules: {len(atoms_to_check)} atoms"
            log_and_print(message)
    except:
        atoms_to_check = mem.atoms
        logging.warning("Error selecting specific residues, processing all atoms")

    # Find atoms in nanotube region
    z_range_atoms = atoms_to_check.select_atoms(f"prop z > {z_min} and prop z < {z_max}")
    coords = z_range_atoms.positions
    dist_sq = (coords[:,0] - x_center)**2 + (coords[:,1] - y_center)**2
    inside_atoms = z_range_atoms[dist_sq < radius**2]
    
    logging.info(f"Atoms in z-range {z_min:.2f} to {z_max:.2f}: {len(z_range_atoms)}")
    logging.info(f"Atoms inside exclusion radius: {len(inside_atoms)}")
    
    # Remove any residue with ANY atom in the exclusion zone (zero tolerance)
    overlapping_residues = set()
    for atom in inside_atoms:
        overlapping_residues.add(atom.residue.resid)
    
    # Mark residues for removal
    bad_resids = [str(res_id) for res_id in overlapping_residues]
    
    if bad_resids:
        bad_residues = mem.select_atoms(f"resid {' '.join(bad_resids)}").residues
        logging.info(f"Residues marked for removal: {bad_resids}")
    else:
        bad_residues = mem.atoms.residues[0:0]
        logging.info("No residues found for removal")
    
    message = f"{len(bad_residues)} residues removed"
    log_and_print(message)

    # Save cleaned membrane only if specified
    keep_atoms = mem.atoms - bad_residues.atoms
    logging.info(f"Remaining atoms after cleanup: {len(keep_atoms)}")
    
    if output_file:
        keep_atoms.write(output_file)
        message = f"Cleaned membrane: {len(keep_atoms)} atoms → {output_file}"
        log_and_print(message)
        
        # Check if file exists and get size
        try:
            if os.path.exists(output_file):
                file_size = os.path.getsize(output_file)
                logging.info(f"Cleaned membrane file size: {file_size} bytes")
            else:
                logging.warning(f"Output file {output_file} was not found after writing")
        except Exception as e:
            logging.warning(f"Could not get file size for {output_file}: {e}")
    
    # Create complex if requested
    if complex_file:
        combined = mda.Merge(tube.atoms, keep_atoms)
        combined.atoms.write(complex_file)
        message = f"Complex: {len(combined.atoms)} atoms → {complex_file}"
        log_and_print(message)
        
        # Check if file exists and get size
        try:
            if os.path.exists(complex_file):
                file_size = os.path.getsize(complex_file)
                logging.info(f"Complex file size: {file_size} bytes")
            else:
                logging.warning(f"Complex file {complex_file} was not found after writing")
        except Exception as e:
            logging.warning(f"Could not get file size for {complex_file}: {e}")
    
    logging.info("-" * 40)
    logging.info(f"Process completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info("="*60)

def main():
    parser = argparse.ArgumentParser(description='Remove overlapping residues and create complex structure with no overlap')
    parser.add_argument('nanotube', help='Nanotube PDB file')
    parser.add_argument('membrane_pdb', help='Membrane PDB file')
    parser.add_argument('-o', '--output', help='Cleaned membrane output (optional)')
    parser.add_argument('-c', '--complex', help='Complex structure output (optional)')
    parser.add_argument('-r', '--radius_factor', type=float, default=1.0, help='Radius scaling factor (default: 1.0)')
    parser.add_argument('-b', '--buffer_distance', type=float, default=0.8, help='Buffer distance in Å (default: 0.8)')
    parser.add_argument('-l', '--log', help='Log file path (default: lipid_removal.log)', default='lipid_removal.log')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log)
    
    remove_lipids_and_combine(args.nanotube, args.membrane_pdb, args.output, 
                             args.complex, args.radius_factor, args.buffer_distance)

if __name__ == "__main__":
    main()
