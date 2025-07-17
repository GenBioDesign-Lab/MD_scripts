#!/usr/bin/env python3
"""
Example usage:
    python dcdcentering_update.py \
        -structure test/cnt_mem_ionized.psf \
        -traj test/final_NPT.dcd \
        -sel "segid CNT1 or segid CNT2 or segid CNT3" \
        -center geometry \
        -o CNT_centered_wrapped.xtc
        
Performance options:
    -stride 10          # Write every 10th frame
"""

# IMPORTS AND WARNING SUPPRESSION
import warnings
import os
import argparse

# Suppress all warnings comprehensively
warnings.simplefilter('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'

import MDAnalysis as mda
from MDAnalysis.transformations import unwrap, center_in_box, wrap
import MDAnalysis.coordinates.DCD

# Additional MDAnalysis specific warning suppression
warnings.filterwarnings('ignore', module='MDAnalysis')

# COMMAND LINE ARGUMENT PARSING

parser = argparse.ArgumentParser(description='Center and wrap trajectory')
parser.add_argument('-structure', required=True, 
                   help='Path to structure file (.psf)')
parser.add_argument('-traj', required=True, 
                   help='Path to trajectory file (.dcd)')
parser.add_argument('-sel', required=True, 
                   help='Selection string for atoms to center (e.g., "segid CNT1" or "resname PROT")')
parser.add_argument('-center', required=True, choices=['mass', 'geometry'], 
                   help='Center using mass or geometry')
parser.add_argument('-o', required=True, 
                   help='Output trajectory file name (.xtc)')
parser.add_argument('-stride', type=int, default=1,
                   help='Frame stride (default: 1, every frame)')

args = parser.parse_args()

# SYSTEM LOADING AND SETUP

print("Loading molecular system...")
print(f"  Structure: {args.structure}")
print(f"  Trajectory: {args.traj}")

# Load system with warnings suppressed
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    u = mda.Universe(args.structure, args.traj)
    u.transfer_to_memory(step=args.stride)

# Define atom selections
protein = u.select_atoms(args.sel)
all_atoms = u.atoms

# Validate selection
if len(protein) == 0:
    raise ValueError(f"Selection '{args.sel}' returned no atoms. Please check your selection string.")

print(f"  Selected {len(protein)} atoms for centering using: '{args.sel}'")
print(f"  Centering method: {args.center}")

# TRAJECTORY TRANSFORMATIONS

print("\nApplying transformations...")

# Transformation:
# 1. Unwrap molecules across periodic boundaries
# 2. Center the selected atoms in the box
# 3. Wrap all molecules back into the box
transformations = [
    unwrap(all_atoms),                                    # Fix PBC
    center_in_box(protein, wrap=False, center=args.center),  # Center selection
    wrap(all_atoms, compound='residues')                     # Wrap system
]

u.trajectory.add_transformations(*transformations)

# OUTPUT WRITING
total_frames = len(u.trajectory)
frames_to_write = total_frames  # Already strided in memory

print(f"\nWriting trajectory to: {args.o}")
print(f"  Total frames in memory: {total_frames}")
print(f"  Stride: {args.stride}")
print(f"  Frames to write: {frames_to_write}")

# Write XTC (trajectory already strided in memory)
with mda.Writer(args.o, all_atoms.n_atoms) as writer:
    frame_count = 0
    for i, ts in enumerate(u.trajectory):
        writer.write(all_atoms)
        frame_count += 1
        
        # Progress reporting every 100 frames or 10%
        if frame_count % 100 == 0 or frame_count % max(1, frames_to_write // 10) == 0:
            progress = (frame_count / frames_to_write) * 100
            print(f"  Progress: {frame_count}/{frames_to_write} frames ({progress:.1f}%)")

print("Trajectory processing completed")

