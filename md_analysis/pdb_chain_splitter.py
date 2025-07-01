#!/usr/bin/env python3
"""
PDB Chain Splitter
==================

This script automatically detects chains in a PDB file and splits them into separate files.
Each chain will be saved as a separate PDB file in the 'output/' directory.

Author: Generated for MD analysis workflow
Usage: python pdb_chain_splitter.py <pdb_file>
"""

import os
import sys
import argparse
from pathlib import Path
from collections import defaultdict

def create_output_directory(output_dir="output"):
    """Create output directory if it doesn't exist."""
    Path(output_dir).mkdir(exist_ok=True)
    return output_dir

def parse_pdb_file(pdb_file):
    """
    Parse PDB file and extract chains.
    Returns a dictionary with chain IDs as keys and list of lines as values.
    """
    chains = defaultdict(list)
    header_lines = []
    footer_lines = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            
            # Keep header information (HEADER, TITLE, COMPND, etc.)
            if line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 
                               'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK')):
                header_lines.append(line)
            
            # ATOM and HETATM records contain chain information
            elif line.startswith(('ATOM', 'HETATM')):
                if len(line) >= 22:  # Ensure line is long enough to contain chain ID
                    chain_id = line[21] if len(line) > 21 else ' '
                    # If chain ID is empty, use 'A' as default
                    if chain_id == ' ' or chain_id == '':
                        chain_id = 'A'
                    chains[chain_id].append(line)
            
            # Keep connectivity information
            elif line.startswith(('CONECT', 'MASTER', 'END')):
                footer_lines.append(line)
    
    return chains, header_lines, footer_lines

def write_chain_file(chain_id, chain_lines, header_lines, footer_lines, output_dir, base_filename):
    """Write individual chain to a separate PDB file."""
    output_filename = f"{base_filename}_chain_{chain_id}.pdb"
    output_path = os.path.join(output_dir, output_filename)
    
    with open(output_path, 'w') as f:
        # Write header information
        for line in header_lines:
            f.write(line + '\n')
        
        # Add a remark about the chain
        f.write(f"REMARK   1 CHAIN {chain_id} EXTRACTED FROM {base_filename}\n")
        
        # Write chain atoms
        for line in chain_lines:
            f.write(line + '\n')
        
        # Write footer information
        for line in footer_lines:
            f.write(line + '\n')
    
    return output_path

def get_chain_info(chains):
    """Get information about detected chains."""
    chain_info = {}
    for chain_id, lines in chains.items():
        atom_count = len([line for line in lines if line.startswith('ATOM')])
        hetatm_count = len([line for line in lines if line.startswith('HETATM')])
        chain_info[chain_id] = {
            'atoms': atom_count,
            'hetatms': hetatm_count,
            'total': len(lines)
        }
    return chain_info

def split_pdb_chains(pdb_file, output_dir="output"):
    """
    Main function to split PDB file into separate chain files.
    """
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    # Create output directory
    output_dir = create_output_directory(output_dir)
    
    # Parse PDB file
    print(f"Parsing PDB file: {pdb_file}")
    chains, header_lines, footer_lines = parse_pdb_file(pdb_file)
    
    if not chains:
        print("No chains detected in the PDB file.")
        return
    
    # Get chain information
    chain_info = get_chain_info(chains)
    
    # Display detected chains
    print(f"\nDetected {len(chains)} chain(s):")
    for chain_id, info in chain_info.items():
        print(f"  Chain {chain_id}: {info['atoms']} atoms, {info['hetatms']} HETATM records")
    
    # Split chains into separate files
    base_filename = Path(pdb_file).stem
    output_files = []
    
    print(f"\nSplitting chains into separate files in '{output_dir}/':")
    for chain_id, chain_lines in chains.items():
        output_path = write_chain_file(chain_id, chain_lines, header_lines, 
                                     footer_lines, output_dir, base_filename)
        output_files.append(output_path)
        print(f"  Created: {output_path}")
    
    print(f"\nSuccessfully split {len(chains)} chains from {pdb_file}")
    return output_files

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Split PDB file into separate chain files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pdb_chain_splitter.py protein.pdb
  python pdb_chain_splitter.py protein.pdb --output-dir my_chains
        """
    )
    
    parser.add_argument('pdb_file', help='Input PDB file to split')
    parser.add_argument('--output-dir', '-o', default='output',
                       help='Output directory for chain files (default: output)')
    
    args = parser.parse_args()
    
    try:
        split_pdb_chains(args.pdb_file, args.output_dir)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 