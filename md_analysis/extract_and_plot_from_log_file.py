import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
#Usage: python extract_and_plot_from_log_file.py <log_file>
def select_column(available_columns):
    """Display column options and get user selection"""
    print("\nAvailable columns:")
    for i, col in enumerate(available_columns):
        print(f"{i}: {col}")
    
    while True:
        try:
            choice = input("\nEnter column number (or 'exit'): ")
            if choice.lower() == 'exit':
                return None
            idx = int(choice)
            if 0 <= idx < len(available_columns):
                return idx
            print(f"Enter number between 0 and {len(available_columns)-1}")
        except ValueError:
            print("Enter a valid number")

def extract_and_plot_namd_data(log_file):
    try:
        print(f"Reading: {log_file}")
        
        log_dir = os.path.dirname(log_file) or '.'
        log_basename = os.path.splitext(os.path.basename(log_file))[0]
        
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        # Find headers
        headers = []
        for line in lines:
            if line.startswith('ETITLE:'):
                headers = line.strip().split()[1:]
                break
        
        if not headers:
            print("Error: ETITLE line not found")
            return
        
        # Extract data
        data = []
        for line in lines:
            if line.startswith('ENERGY:'):
                values = [float(val) for val in line.strip().split()[1:]]
                data.append(values)
        
        if not data:
            print("Error: No ENERGY data found")
            return
        
        df = pd.DataFrame(data, columns=headers)
        
        # Save CSV
        output_dir = os.path.join(log_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        csv_file = os.path.join(output_dir, f"{log_basename}.csv")
        df.to_csv(csv_file, index=False)
        print(f"Data saved: {csv_file}")
        
        # Setup x-axis
        if 'TS' in df.columns:
            df['Time (ns)'] = df['TS'] * 2e-6
            x_col = 'Time (ns)'
            plot_cols = [col for col in df.columns if col not in ['TS', 'Time (ns)']]
        else:
            x_col = headers[0]
            plot_cols = [col for col in df.columns if col != x_col]
        
        # Plot loop
        while True:
            col_idx = select_column(plot_cols)
            if col_idx is None:
                break
            
            y_col = plot_cols[col_idx]
            print(f"Plotting {y_col} vs {x_col}")
            
            plt.figure(figsize=(10, 6))
            plt.plot(df[x_col], df[y_col])
            plt.xlabel(x_col)
            plt.ylabel(y_col)
            plt.title(f'{y_col} vs {x_col}')
            plt.grid(True, alpha=0.3)
            
            plot_file = os.path.join(output_dir, f"{log_basename}_{y_col.replace(' ', '_')}.png")
            plt.savefig(plot_file, dpi=300)
            print(f"Plot saved: {plot_file}")
            plt.show()
    
    except FileNotFoundError:
        print(f"File not found: {log_file}")
    except Exception as e:
        print(f"Error: {e}")

def main():
    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help']:
        print("Usage: python extract_and_plot_from_log_file.py <log_file>")
        return
    
    if len(sys.argv) > 1:
        extract_and_plot_namd_data(sys.argv[1])
    else:
        log_file = input("Enter NAMD log file path: ")
        if log_file:
            extract_and_plot_namd_data(log_file)

if __name__ == "__main__":
    main() 