import MDAnalysis as mda
import numpy as np
import argparse
import warnings
import logging
import os
from datetime import datetime

# Usage: python remove_lipid.py nanotube.pdb membrane.pdb -o -s cleaned_membrane.pdb -c complex.pdb

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

def remove_lipids_and_combine(nanotube_file, membrane_pdb, output_file, complex_file=None, radius_factor=1.0, buffer_distance=0.8, symmetric_removal=False):
    logging.info(f"Input parameters:")
    logging.info(f"  Nanotube file: {nanotube_file}")
    logging.info(f"  Membrane file: {membrane_pdb}")
    logging.info(f"  Output file: {output_file}")
    logging.info(f"  Complex file: {complex_file}")
    logging.info(f"  Radius factor: {radius_factor}")
    logging.info(f"  Buffer distance: {buffer_distance} Å")
    logging.info(f"  Symmetric removal: {symmetric_removal}")
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

    # Determine membrane center for symmetric removal
    membrane_z_center = None
    if symmetric_removal:
        # Try to find phosphate atoms (membrane headgroups) to determine membrane center
        try:
            phosphate_atoms = mem.select_atoms("name P")
            if len(phosphate_atoms) > 0:
                membrane_z_center = phosphate_atoms.positions[:, 2].mean()
                logging.info(f"Membrane center (from phosphate atoms): Z = {membrane_z_center:.2f} Å")
            else:
                # Fallback: use all membrane atoms
                membrane_z_center = mem.atoms.positions[:, 2].mean()
                logging.info(f"Membrane center (from all atoms): Z = {membrane_z_center:.2f} Å")
        except:
            membrane_z_center = mem.atoms.positions[:, 2].mean()
            logging.info(f"Membrane center (fallback): Z = {membrane_z_center:.2f} Å")

    # Select atoms to process (try lipids and water first, then all atoms)
    try:
        atoms_to_check = mem.select_atoms("resname POPC or resname POPE or resname DPPC or resname TIP3 or resname SOL or resname WAT")
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
    
    # Symmetric removal logic
    additional_removed_resids = []
    if symmetric_removal and membrane_z_center is not None and len(bad_residues) > 0:
        logging.info("Applying symmetric lipid removal...")
        
        # Separate lipids by leaflet
        try:
            lipid_residues = mem.select_atoms("resname POPC or resname POPE or resname DPPC").residues
            if len(lipid_residues) == 0:
                logging.warning("No POPC/POPE/DPPC lipids found for symmetric removal")
            else:
                # Classify removed lipids by leaflet
                upper_removed = []
                lower_removed = []
                
                for res in bad_residues:
                    # Check if this is a lipid residue
                    lipid_atoms = res.atoms.select_atoms("resname POPC or resname POPE or resname DPPC")
                    if len(lipid_atoms) > 0:
                        res_z_center = res.atoms.positions[:, 2].mean()
                        if res_z_center > membrane_z_center:
                            upper_removed.append(res)
                        else:
                            lower_removed.append(res)
                
                logging.info(f"Originally removed: {len(upper_removed)} upper leaflet, {len(lower_removed)} lower leaflet lipids")
                
                # Count non-lipid residues for reference
                non_lipid_removed = len(bad_residues) - len(upper_removed) - len(lower_removed)
                if non_lipid_removed > 0:
                    logging.info(f"Non-lipid residues removed (water, etc.): {non_lipid_removed}")
                
                # Balance the removal
                upper_count = len(upper_removed)
                lower_count = len(lower_removed)
                
                if upper_count > lower_count:
                    # Need to remove more from lower leaflet
                    needed = upper_count - lower_count
                    lower_leaflet_lipids = [res for res in lipid_residues 
                                          if res.atoms.positions[:, 2].mean() < membrane_z_center 
                                          and res.resid not in [r.resid for r in bad_residues]]
                    
                    if len(lower_leaflet_lipids) >= needed:
                        # Find closest lipids to the exclusion zone in XY plane
                        distances = []
                        for res in lower_leaflet_lipids:
                            res_center = res.atoms.positions.mean(axis=0)
                            xy_dist = np.sqrt((res_center[0] - x_center)**2 + (res_center[1] - y_center)**2)
                            distances.append((xy_dist, res))
                        
                        distances.sort(key=lambda x: x[0])
                        for i in range(needed):
                            additional_removed_resids.append(distances[i][1].resid)
                        
                        logging.info(f"Added {needed} lower leaflet lipids for symmetry")
                    else:
                        logging.warning(f"Not enough lower leaflet lipids to balance removal")
                        
                elif lower_count > upper_count:
                    # Need to remove more from upper leaflet
                    needed = lower_count - upper_count
                    upper_leaflet_lipids = [res for res in lipid_residues 
                                          if res.atoms.positions[:, 2].mean() > membrane_z_center 
                                          and res.resid not in [r.resid for r in bad_residues]]
                    
                    if len(upper_leaflet_lipids) >= needed:
                        # Find closest lipids to the exclusion zone in XY plane
                        distances = []
                        for res in upper_leaflet_lipids:
                            res_center = res.atoms.positions.mean(axis=0)
                            xy_dist = np.sqrt((res_center[0] - x_center)**2 + (res_center[1] - y_center)**2)
                            distances.append((xy_dist, res))
                        
                        distances.sort(key=lambda x: x[0])
                        for i in range(needed):
                            additional_removed_resids.append(distances[i][1].resid)
                        
                        logging.info(f"Added {needed} upper leaflet lipids for symmetry")
                    else:
                        logging.warning(f"Not enough upper leaflet lipids to balance removal")
                else:
                    logging.info("Lipid removal is already symmetric")
                    
        except Exception as e:
            logging.warning(f"Error during symmetric removal: {e}")
    
    # Combine all residues to remove
    all_bad_resids = bad_resids + [str(resid) for resid in additional_removed_resids]
    
    if all_bad_resids:
        all_bad_residues = mem.select_atoms(f"resid {' '.join(all_bad_resids)}").residues
        logging.info(f"Total residues marked for removal: {all_bad_resids}")
    else:
        all_bad_residues = mem.atoms.residues[0:0]
    
    # Detailed breakdown of removal
    original_removal = len(bad_residues)
    symmetry_removal = len(additional_removed_resids)
    total_removal = len(all_bad_residues)
    
    if symmetric_removal and symmetry_removal > 0:
        message = f"{total_removal} total residues removed: {original_removal} overlapping + {symmetry_removal} for leaflet symmetry"
        log_and_print(message)
        
        # Additional details about lipid symmetry
        if membrane_z_center is not None:
            try:
                final_lipids = mem.select_atoms("resname POPC or resname POPE or resname DPPC").residues
                final_lipids_remaining = [res for res in final_lipids if res.resid not in [int(r) for r in all_bad_resids]]
                
                upper_remaining = len([res for res in final_lipids_remaining 
                                     if res.atoms.positions[:, 2].mean() > membrane_z_center])
                lower_remaining = len([res for res in final_lipids_remaining 
                                     if res.atoms.positions[:, 2].mean() < membrane_z_center])
                
                logging.info(f"Final lipid count: {upper_remaining} upper leaflet, {lower_remaining} lower leaflet")
                if abs(upper_remaining - lower_remaining) <= 1:
                    logging.info("✓ Leaflet balance achieved (difference ≤ 1)")
                else:
                    logging.warning(f"⚠ Leaflet imbalance: difference = {abs(upper_remaining - lower_remaining)}")
            except Exception as e:
                logging.warning(f"Could not verify final leaflet balance: {e}")
    else:
        message = f"{total_removal} residues removed"
        log_and_print(message)

    # Save cleaned membrane only if specified
    keep_atoms = mem.atoms - all_bad_residues.atoms
    logging.info(f"Remaining atoms after cleanup: {len(keep_atoms)}")
    
    if output_file:
        try:
            # Get absolute path for better debugging
            abs_output_file = os.path.abspath(output_file)
            logging.info(f"Writing cleaned membrane to: {abs_output_file}")
            logging.info(f"Current working directory: {os.getcwd()}")
            
            # Ensure output directory exists
            output_dir = os.path.dirname(abs_output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)
                logging.info(f"Created output directory: {output_dir}")
            
            # Write the file
            keep_atoms.write(output_file)
            message = f"Cleaned membrane: {len(keep_atoms)} atoms → {output_file}"
            log_and_print(message)
            
            # Check if file exists and get size
            if os.path.exists(abs_output_file):
                file_size = os.path.getsize(abs_output_file)
                logging.info(f"Cleaned membrane file size: {file_size} bytes")
                logging.info(f"File successfully written to: {abs_output_file}")
            elif os.path.exists(output_file):
                file_size = os.path.getsize(output_file)
                logging.info(f"Cleaned membrane file size: {file_size} bytes")
                logging.info(f"File successfully written to: {output_file}")
            else:
                logging.error(f"Output file not found after writing!")
                logging.error(f"Attempted paths:")
                logging.error(f"  Relative: {output_file}")
                logging.error(f"  Absolute: {abs_output_file}")
                logging.error(f"  Files in current directory: {os.listdir('.')}")
                
        except Exception as e:
            logging.error(f"Error writing cleaned membrane file: {e}")
            logging.error(f"Exception type: {type(e).__name__}")
            import traceback
            logging.error(f"Traceback: {traceback.format_exc()}")
    
    # Create complex if requested
    if complex_file:
        try:
            # Get absolute path for better debugging
            abs_complex_file = os.path.abspath(complex_file)
            logging.info(f"Writing complex to: {abs_complex_file}")
            
            # Ensure output directory exists
            complex_dir = os.path.dirname(abs_complex_file)
            if complex_dir and not os.path.exists(complex_dir):
                os.makedirs(complex_dir)
                logging.info(f"Created complex directory: {complex_dir}")
            
            # Create and write complex
            combined = mda.Merge(tube.atoms, keep_atoms)
            combined.atoms.write(complex_file)
            message = f"Complex: {len(combined.atoms)} atoms → {complex_file}"
            log_and_print(message)
            
            # Check if file exists and get size
            if os.path.exists(abs_complex_file):
                file_size = os.path.getsize(abs_complex_file)
                logging.info(f"Complex file size: {file_size} bytes")
                logging.info(f"Complex file successfully written to: {abs_complex_file}")
            elif os.path.exists(complex_file):
                file_size = os.path.getsize(complex_file)
                logging.info(f"Complex file size: {file_size} bytes")
                logging.info(f"Complex file successfully written to: {complex_file}")
            else:
                logging.error(f"Complex file not found after writing!")
                logging.error(f"Attempted paths:")
                logging.error(f"  Relative: {complex_file}")
                logging.error(f"  Absolute: {abs_complex_file}")
                
        except Exception as e:
            logging.error(f"Error writing complex file: {e}")
            logging.error(f"Exception type: {type(e).__name__}")
            import traceback
            logging.error(f"Traceback: {traceback.format_exc()}")
    
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
    parser.add_argument('-s', '--symmetric', action='store_true', help='Enable symmetric lipid removal to maintain leaflet balance')
    parser.add_argument('-l', '--log', help='Log file path (default: lipid_removal.log)', default='lipid_removal.log')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log)
    
    remove_lipids_and_combine(args.nanotube, args.membrane_pdb, args.output, 
                             args.complex, args.radius_factor, args.buffer_distance, args.symmetric)

if __name__ == "__main__":
    main()
