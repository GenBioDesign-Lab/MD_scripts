import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

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

def plot_namd_data(csv_file, save_dir=None):
    try:
        print(f"\nReading data from: {csv_file}")
        # Read the CSV file
        df = pd.read_csv(csv_file)
        
        # Check if 'TS' column exists
        if 'TS' not in df.columns:
            print("Error: CSV file must contain a 'TS' column")
            return
        
        # Convert TS to Time (ns)
        df['Time (ns)'] = df['TS'] * 2e-6  # TS * 2fs / 1e6 = ns
        
        # Get columns for plotting (excluding TS since we use it for x-axis)
        available_columns = [col for col in df.columns if col != 'TS' and col != 'Time (ns)']
        
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
            print(f"\nPlotting {selected_column} vs Time (ns)")
            
            # Check if user wants to save the plot
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
            
            # Create the plot
            plt.figure(figsize=(10, 6))
            plt.plot(df['Time (ns)'], df[selected_column], linewidth=2)
            plt.xlabel('Time (ns)', fontsize=12)
            plt.ylabel(selected_column, fontsize=12)
            plt.title(f'{selected_column} vs Time', fontsize=14)
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.tight_layout()
            
            # Save the plot if requested
            if save_this_plot and save_path:
                os.makedirs(save_path, exist_ok=True)
                plot_filename = os.path.join(save_path, f"{selected_column.replace(' ', '_')}_vs_time.png")
                plt.savefig(plot_filename, dpi=300)
                print(f"\nPlot saved to: {plot_filename}")
            
            print("\nShowing plot... (close the plot window to return to menu)")
            # Show the plot (will return when user closes plot window)
            plt.show()
            
            print("\nReturning to column selection menu...")
    
    except FileNotFoundError:
        print(f"Error: Could not find the file '{csv_file}'")
    except Exception as e:
        print(f"Error: {str(e)}")

def main():
    # Check for help flag first
    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
        print("\nNAMD Data Plotter")
        print("\nUsage:")
        print("  python plot_namd_data.py path/to/data.csv [-o save_directory]")
        print("\nOptions:")
        print("  -o, --output DIR    Directory to save plot images")
        print("  -h, --help          Show this help message")
        return
    
    # Simple usage with file name
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
        save_dir = None
        
        # Check for -o or --output option
        for i in range(2, len(sys.argv)):
            if sys.argv[i] in ['-o', '--output'] and i+1 < len(sys.argv):
                save_dir = sys.argv[i+1]
                break
        
        plot_namd_data(csv_file, save_dir)
    else:
        print("\nUsage:")
        print("  python plot_namd_data.py path/to/data.csv [-o save_directory]")
        
        csv_file = input("\nEnter the path to your CSV file: ")
        if csv_file:
            save_choice = input("Save plots to specific directory? (y/n, default: n, if no, just press enter): ").lower()
            if save_choice.startswith('y'):
                save_dir = input("Enter directory to save plots: ")
                plot_namd_data(csv_file, save_dir)
            else:
                plot_namd_data(csv_file)
        else:
            print("No file provided. Exiting.")

if __name__ == "__main__":
    main() 