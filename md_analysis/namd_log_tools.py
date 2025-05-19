import re
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys

def extract_namd_data(log_file, output_csv=None):
    if output_csv is None:
        log_dir = os.path.dirname(log_file) or "."
        log_base = os.path.basename(log_file)
        base_name = os.path.splitext(log_base)[0]
        output_csv = os.path.join(log_dir, f"{base_name}.csv")
    print(f"Reading NAMD log file: {log_file}")
    with open(log_file, 'r') as f:
        log_content = f.readlines()
    etitle_line = ''
    for line in log_content:
        if line.startswith('ETITLE:'):
            etitle_line = line
            break
    if not etitle_line:
        print("Error: ETITLE line not found in the log file")
        return False
    headers = etitle_line.strip().split()[1:]
    energy_data = []
    for line in log_content:
        if line.startswith('ENERGY:'):
            values = line.strip().split()[1:]
            values = [float(val) for val in values]
            energy_data.append(values)
    if not energy_data:
        print("Error: No ENERGY data found in the log file")
        return False
    df = pd.DataFrame(energy_data, columns=headers)
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    df.to_csv(output_csv, index=False)
    print(f"Success! Extracted {len(energy_data)} data points to {output_csv}")
    print("\nTo plot this data, run:")
    print(f"python namd_log_tools.py plot {output_csv}")
    return True

def select_column(available_columns):
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

def plot_namd_data(csv_file, save_dir=None):
    try:
        print(f"\nReading data from: {csv_file}")
        df = pd.read_csv(csv_file)
        if 'TS' not in df.columns:
            print("Error: CSV file must contain a 'TS' column")
            return
        df['Time (ns)'] = df['TS'] * 2e-6
        available_columns = [col for col in df.columns if col != 'TS' and col != 'Time (ns)']
        while True:
            column_index = select_column(available_columns)
            if column_index is None:
                print("Exiting...")
                return
            selected_column = available_columns[column_index]
            print(f"\nPlotting {selected_column} vs Time (ns)")
            save_this_plot = False
            if save_dir is None:
                save_choice = input("\nDo you want to save this plot? (yes/no): ").lower()
                if save_choice.startswith('y'):
                    save_this_plot = True
                    save_path = input("Enter directory to save the plot (press Enter for current directory): ")
                    if not save_path:
                        save_path = "."
                else:
                    save_path = None
            else:
                save_path = save_dir
                save_this_plot = True
            plt.figure(figsize=(10, 6))
            plt.plot(df['Time (ns)'], df[selected_column], linewidth=2)
            plt.xlabel('Time (ns)', fontsize=12)
            plt.ylabel(selected_column, fontsize=12)
            plt.title(f'{selected_column} vs Time', fontsize=14)
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.tight_layout()
            if save_this_plot and save_path:
                os.makedirs(save_path, exist_ok=True)
                plot_filename = os.path.join(save_path, f"{selected_column.replace(' ', '_')}_vs_time.png")
                plt.savefig(plot_filename, dpi=300)
                print(f"\nPlot saved to: {plot_filename}")
            print("\nShowing plot... (close the plot window to return to menu)")
            plt.show()
            print("\nReturning to column selection menu...")
    except FileNotFoundError:
        print(f"Error: Could not find the file '{csv_file}'")
    except Exception as e:
        print(f"Error: {str(e)}")

def main():
    print("\nNAMD Log Tools")
    print("====================")
    print("Choose an option:")
    print("1: Extract NAMD log data to CSV")
    print("2: Plot data from CSV file")
    print("h: Help")
    print("q: Quit")
    while True:
        choice = input("\nEnter your choice (1/2/h/q): ").strip().lower()
        if choice == '1':
            log_file = input("Enter the path to your NAMD log file: ")
            if not log_file:
                print("No file provided. Returning to menu.")
                continue
            output_choice = input("Save CSV to custom location? (y/n, default: n): ").lower()
            if output_choice.startswith('y'):
                output_csv = input("Enter the output CSV file path: ")
                extract_namd_data(log_file, output_csv)
            else:
                extract_namd_data(log_file)
        elif choice == '2':
            csv_file = input("Enter the path to your CSV file: ")
            if not csv_file:
                print("No file provided. Returning to menu.")
                continue
            save_choice = input("Save plots to specific directory? (y/n, default: n): ").lower()
            if save_choice.startswith('y'):
                save_dir = input("Enter directory to save plots: ")
                plot_namd_data(csv_file, save_dir)
            else:
                plot_namd_data(csv_file)
        elif choice == 'h':
            print("\nNAMD Log Tools Help")
            print("1: Extract NAMD log data to CSV")
            print("2: Plot data from CSV file")
            print("You can also run this script with arguments:")
            print("  python namd_log_tools.py extract path/to/namd.log [-o output.csv]")
            print("  python namd_log_tools.py plot path/to/data.csv [-o save_directory]")
        elif choice == 'q':
            print("Exiting.")
            break
        else:
            print("Invalid choice. Please enter 1, 2, h, or q.")

def cli():
    # Command-line interface for direct calls
    if len(sys.argv) > 1:
        if sys.argv[1] in ['-h', '--help', 'help']:
            print("\nNAMD Log Tools Help")
            print("Usage:")
            print("  python namd_log_tools.py extract path/to/namd.log [-o output.csv]")
            print("  python namd_log_tools.py plot path/to/data.csv [-o save_directory]")
            print("Options:")
            print("  -o, --output PATH    Specify output file or directory")
            print("  -h, --help           Show this help message")
            return
        if sys.argv[1] == 'extract' and len(sys.argv) > 2:
            log_file = sys.argv[2]
            output_csv = None
            for i in range(3, len(sys.argv)):
                if sys.argv[i] in ['-o', '--output'] and i+1 < len(sys.argv):
                    output_csv = sys.argv[i+1]
                    break
            extract_namd_data(log_file, output_csv)
            return
        if sys.argv[1] == 'plot' and len(sys.argv) > 2:
            csv_file = sys.argv[2]
            save_dir = None
            for i in range(3, len(sys.argv)):
                if sys.argv[i] in ['-o', '--output'] and i+1 < len(sys.argv):
                    save_dir = sys.argv[i+1]
                    break
            plot_namd_data(csv_file, save_dir)
            return
    # If no valid CLI args, show menu
    main()

if __name__ == "__main__":
    cli() 