#!/usr/bin/env python3
"""
Usage: python split_pdb_chains.py input.pdb
"""

import os
import sys
from pathlib import Path
from collections import defaultdict


def split_pdb_chains(pdb_file):
    """Split PDB file into separate chain files, including renumbered segments."""
    # Check if file exists
    if not os.path.exists(pdb_file):
        print(f"Error: File '{pdb_file}' not found.")
        return False
    
    print(f"Processing: {pdb_file}")
    
    # Parse chains and detect renumbering
    chains = defaultdict(list)
    chain_counters = {}  # Track segment numbers for each chain
    last_residue_num = {}  # Track last residue number for each chain
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                # Extract chain ID and residue number
                chain_id = line[21] if len(line) > 21 and line[21] != ' ' else 'A'
                
                # Extract residue number (columns 23-26, but strip whitespace)
                try:
                    residue_num = int(line[22:26].strip())
                except (ValueError, IndexError):
                    residue_num = 1  # Default if can't parse
                
                # Initialize tracking for new chains
                if chain_id not in last_residue_num:
                    last_residue_num[chain_id] = residue_num
                    chain_counters[chain_id] = 1
                    segment_id = f"{chain_id}_{chain_counters[chain_id]}"
                else:
                    # Check if residue number has reset/renumbered
                    if residue_num < last_residue_num[chain_id]:
                        # Residue numbering has reset - start new segment
                        chain_counters[chain_id] += 1
                        print(f"  Detected renumbering in chain {chain_id}: {last_residue_num[chain_id]} -> {residue_num}")
                    
                    segment_id = f"{chain_id}_{chain_counters[chain_id]}"
                    last_residue_num[chain_id] = residue_num
                
                chains[segment_id].append(line)
    
    if not chains:
        print("No chains found in PDB file.")
        return False
    
    # Create output directory
    base_name = Path(pdb_file).stem
    output_dir = os.path.join(os.path.dirname(pdb_file), base_name)
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Found segments: {', '.join(sorted(chains.keys()))}")
    print(f"Output directory: {output_dir}")
    
    # Write separate files for each chain segment
    for segment_id in sorted(chains.keys()):
        output_file = os.path.join(output_dir, f"{base_name}_chain_{segment_id}.pdb")
        with open(output_file, 'w') as f:
            for line in chains[segment_id]:
                f.write(line)
            f.write("END\n")
        print(f"Created: {base_name}_chain_{segment_id}.pdb ({len(chains[segment_id])} atoms)")
    
    print(f"Successfully split {len(chains)} chain segments!")
    return True


def main():
    """Main function."""
    if len(sys.argv) != 2:
        print("Usage: python split_pdb_chains.py input.pdb")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    success = split_pdb_chains(pdb_file)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main() 