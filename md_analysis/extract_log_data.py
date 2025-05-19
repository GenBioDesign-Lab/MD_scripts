import re
import pandas as pd
import os
import sys

#Usage: python extract_log_data.py <log_file>

def extract_namd_data(log_file):
    # Get directory and base filename
    log_dir = os.path.dirname(log_file) or "."
    log_base = os.path.basename(log_file)
    base_name = os.path.splitext(log_base)[0]
    
    # Create output directory (output/ subfolder in same directory as log file)
    output_dir = os.path.join(log_dir, "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Set output path
    output_csv = os.path.join(output_dir, f"{base_name}.csv")
    
    print(f"Reading NAMD log file: {log_file}")
    
    # Read the log file
    with open(log_file, 'r') as f:
        log_content = f.readlines()
    
    # Find the ETITLE line to get column headers
    etitle_line = ''
    for line in log_content:
        if line.startswith('ETITLE:'):
            etitle_line = line
            break
    
    if not etitle_line:
        print("Error: ETITLE line not found in the log file")
        return False
    
    # Extract headers from ETITLE line
    headers = etitle_line.strip().split()[1:]
    
    # Extract data from ENERGY lines
    energy_data = []
    for line in log_content:
        if line.startswith('ENERGY:'):
            values = line.strip().split()[1:]
            # Convert string values to float
            values = [float(val) for val in values]
            energy_data.append(values)
    
    if not energy_data:
        print("Error: No ENERGY data found in the log file")
        return False
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(energy_data, columns=headers)
    
    # Save to CSV
    df.to_csv(output_csv, index=False)
    print(f"Success! Extracted {len(energy_data)} data points to {output_csv}")
    
    return True

def main():
    if len(sys.argv) != 2:
        print("Usage: python extract_log_data.py <log_file>")
        return
    
    log_file = sys.argv[1]
    extract_namd_data(log_file)

if __name__ == "__main__":
    main() 