import re
import pandas as pd
import argparse
import os
import sys

def extract_namd_data(log_file, output_csv=None):
    # If no output path specified, save in same directory as log file with .csv extension
    if output_csv is None:
        # Get directory and base filename
        log_dir = os.path.dirname(log_file) or "."
        log_base = os.path.basename(log_file)
        # Replace extension with .csv or add .csv if no extension
        base_name = os.path.splitext(log_base)[0]
        output_csv = os.path.join(log_dir, f"{base_name}.csv")
    
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
    
    # Create directory if it doesn't exist
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save to CSV
    df.to_csv(output_csv, index=False)
    print(f"Success! Extracted {len(energy_data)} data points to {output_csv}")
    
    # Suggest the next step to the user
    print("\nTo plot this data, run:")
    print(f"python plot_namd_data.py {output_csv}")
    
    return True

def main():
    # Simple usage
    if len(sys.argv) > 1:
        # Check if help is requested
        if sys.argv[1] in ['-h', '--help']:
            print("\nUsage:")
            print("  python extract_namd_data.py path/to/namd.log [-o output.csv]")
            print("\nOptions:")
            print("  -o, --output PATH    Specify output CSV file path (default: same name with .csv extension)")
            print("  -h, --help           Show this help message")
            return
        
        log_file = sys.argv[1]
        output_csv = None
        
        # Check for -o or --output option
        for i in range(2, len(sys.argv)):
            if sys.argv[i] in ['-o', '--output'] and i+1 < len(sys.argv):
                output_csv = sys.argv[i+1]
                break
        
        extract_namd_data(log_file, output_csv)
    else:
        print("\nUsage:")
        print("  python extract_namd_data.py path/to/namd.log [-o output.csv]")
        
        log_file = input("\nEnter the path to your NAMD log file: ")
        if log_file:
            output_choice = input("Save CSV to custom location? (y/n, default: n): ").lower()
            if output_choice.startswith('y'):
                output_csv = input("Enter the output CSV file path: ")
                extract_namd_data(log_file, output_csv)
            else:
                extract_namd_data(log_file)
        else:
            print("No file provided. Exiting.")

if __name__ == "__main__":
    main() 