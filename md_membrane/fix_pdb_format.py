#!/usr/bin/env python3

import sys
import re

def fix_pdb_formatting(input_file, output_file):
    """
    Fix PDB file formatting issues by ensuring proper column alignment.
    """
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line_num, line in enumerate(infile, 1):
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                outfile.write(line)
                continue
                
            try:
                record_type = line[0:6].strip()
                atom_id = line[6:11].strip()
                
                remaining_line = line[11:].strip()
                fields = remaining_line.split()
                
                if len(fields) < 7:
                    raise ValueError("Not enough fields")
                
                atom_name = fields[0]
                res_name = fields[1] 
                chain_id = fields[2]
                res_id = fields[3]
                
                coord_start_idx = 4
                x = float(fields[coord_start_idx])
                y = float(fields[coord_start_idx + 1])
                z = float(fields[coord_start_idx + 2])
                
                occupancy = fields[coord_start_idx + 3] if len(fields) > coord_start_idx + 3 else "1.00"
                temp_factor = fields[coord_start_idx + 4] if len(fields) > coord_start_idx + 4 else "0.00"
                element = fields[coord_start_idx + 5] if len(fields) > coord_start_idx + 5 else ""
                
                if not element:
                    element_match = re.match(r'([A-Za-z]+)', atom_name)
                    element = element_match.group(1)[0].upper() if element_match else "C"
                
                if len(atom_name) <= 4:
                    formatted_atom_section = f"{atom_name:<4} "
                else:
                    formatted_atom_section = f"{atom_name:<5}"
                
                formatted_line = (
                    f"{record_type:<6}{atom_id:>5} {formatted_atom_section}"
                    f"{res_name:<3} {chain_id:<1}{res_id:>4}    "
                    f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6}{temp_factor:>6}"
                    f"           {element:>2}  \n"
                )
                
                outfile.write(formatted_line)
                
            except (ValueError, IndexError):
                try:
                    atom_id = line[6:11].strip() if len(line) > 11 else "1"
                    atom_section = line[12:18].strip() if len(line) > 18 else "C"
                    atom_name = atom_section.split()[0] if atom_section.split() else "C"
                    res_name = line[17:20].strip() if len(line) > 20 else "UNK"
                    chain_id = line[21:22].strip() if len(line) > 22 else "A"
                    res_id = line[22:26].strip() if len(line) > 26 else "1"
                    
                    coords = re.findall(r'-?\d+\.\d+', line)
                    if len(coords) >= 3:
                        x, y, z = float(coords[0]), float(coords[1]), float(coords[2])
                        
                        if len(atom_name) <= 4:
                            formatted_atom_section = f"{atom_name:<4} "
                        else:
                            formatted_atom_section = f"{atom_name:<5}"
                            
                        formatted_line = (
                            f"ATOM  {atom_id:>5} {formatted_atom_section}{res_name:<3} {chain_id:<1}{res_id:>4}    "
                            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C  \n"
                        )
                        outfile.write(formatted_line)
                    else:
                        outfile.write(line)
                except:
                    outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_pdb_format.py input.pdb output.pdb")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    fix_pdb_formatting(input_file, output_file)
    print("Done!") 