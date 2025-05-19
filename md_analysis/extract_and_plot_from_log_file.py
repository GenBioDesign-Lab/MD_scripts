import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import re

#Usage: python extract_and_plot_from_log_file.py <log_file>

def select_column(available_columns):
    """Display column options and get user selection"""
    print("\n=== Available columns to plot ===")
    for i, col in enumerate(available_columns):
        print(f"{i}: {col}")
    
    print("\nType the number of the column you want to plot and press Enter")
    print("Type 'exit' to quit the program")
    
    while True:
        try:
            user_input = input("\nYour choice: ")
            if user_input.lower() == 'exit':
                return None
            
            column_index = int(user_input)
            if 0 <= column_index < len(available_columns):
                return column_index
            else:
                print(f"Please enter a number between 0 and {len(available_columns)-1}")
        except ValueError:
            print("Please enter a valid number")

def extract_and_plot_namd_data(log_file):
    try:
        print(f"\nReading NAMD log file: {log_file}")
        
        # Get the directory and base filename of the log file
        log_dir = os.path.dirname(log_file)
        if log_dir == '':
            log_dir = '.'
        log_basename = os.path.splitext(os.path.basename(log_file))[0]
        
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
            return
        
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
            return
        
        # Create DataFrame
        df = pd.DataFrame(energy_data, columns=headers)
        
        # Create output directory
        output_dir = os.path.join(log_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Save the extracted data to CSV
        csv_filename = os.path.join(output_dir, f"{log_basename}.csv")
        df.to_csv(csv_filename, index=False)
        print(f"\nExtracted data saved to: {csv_filename}")
        
        # Add Time column if TS exists
        if 'TS' in df.columns:
            df['Time (ns)'] = df['TS'] * 2e-6  # TS * 2fs / 1e6 = ns
            # Use Time as x-axis if available
            x_column = 'Time (ns)'
            # Get columns for plotting (excluding TS and Time since we use Time for x-axis)
            available_columns = [col for col in df.columns if col != 'TS' and col != 'Time (ns)']
        else:
            # If no TS column, use the first column as x-axis
            x_column = headers[0]
            available_columns = [col for col in df.columns if col != x_column]
        
        # Main plotting loop - continue until user exits
        while True:
            # Get column selection from user
            column_index = select_column(available_columns)
            
            # Exit if user chose to quit
            if column_index is None:
                print("Exiting...")
                return
            
            # Get the selected column name
            selected_column = available_columns[column_index]
            print(f"\nPlotting {selected_column} vs {x_column}")
            
            # Create the plot
            plt.figure(figsize=(10, 6))
            plt.plot(df[x_column], df[selected_column], linewidth=2)
            plt.xlabel(x_column, fontsize=12)
            plt.ylabel(selected_column, fontsize=12)
            plt.title(f'{selected_column} vs {x_column}', fontsize=14)
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.tight_layout()
            
            # Save the plot automatically to the same folder as the log file
            output_dir = os.path.join(log_dir, "output")
            os.makedirs(output_dir, exist_ok=True)
            plot_filename = os.path.join(output_dir, f"{log_basename}_{selected_column.replace(' ', '_')}.png")
            plt.savefig(plot_filename, dpi=300)
            print(f"\nPlot saved to: {plot_filename}")
            
            print("\nShowing plot... (close the plot window to return to menu)")
            # Show the plot (will return when user closes plot window)
            plt.show()
            
            print("\nReturning to column selection menu...")
    
    except FileNotFoundError:
        print(f"Error: Could not find the file '{log_file}'")
    except Exception as e:
        print(f"Error: {str(e)}")

def main():
    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
        print("\nNAMD Log Data Plotter")
        print("\nUsage:")
        print("  python extract_log_data.py path/to/logfile.log")
        print("\nOptions:")
        print("  -h, --help          Show this help message")
        return
    
    # Get log file from command line argument
    if len(sys.argv) > 1:
        log_file = sys.argv[1]
        extract_and_plot_namd_data(log_file)
    else:
        print("\nUsage:")
        print("  python extract_log_data.py path/to/logfile.log")
        
        log_file = input("\nEnter the path to your NAMD log file: ")
        if log_file:
            extract_and_plot_namd_data(log_file)
        else:
            print("No file provided. Exiting.")

if __name__ == "__main__":
    main() 