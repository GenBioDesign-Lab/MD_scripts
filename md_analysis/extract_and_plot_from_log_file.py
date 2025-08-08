import argparse
import os
import sys
from typing import List

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract ENERGY table from NAMD log and save CSV and (optional) plots."
    )
    parser.add_argument("log_file", help="Path to NAMD log file")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory to write outputs. Default: '<log_dir>/output'",
    )
    parser.add_argument(
        "--x-column",
        default=None,
        help=(
            "X-axis column to use for plots. Default: 'time_ps' if TS present, otherwise first column"
        ),
    )
    parser.add_argument(
        "--plot-columns",
        nargs="*",
        default=None,
        help="Subset of columns to plot against the x-axis. Default: all except x",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not display plots interactively",
    )
    return parser.parse_args()


def read_namd_energy_table(log_path: str) -> pd.DataFrame:
    with open(log_path, "r") as file_handle:
        lines = file_handle.readlines()

    headers: List[str] = []
    for line in lines:
        if line.startswith("ETITLE:"):
            headers = line.strip().split()[1:]
            break
    if not headers:
        raise ValueError("ETITLE line not found in log file")

    data: List[List[float]] = []
    for line in lines:
        if line.startswith("ENERGY:"):
            values = [float(value) for value in line.strip().split()[1:]]
            data.append(values)
    if not data:
        raise ValueError("No ENERGY data found in log file")

    df = pd.DataFrame(data, columns=headers)

    # Add standardized time column in picoseconds if TS present.
    if "TS" in df.columns:
        df["time_ps"] = df["TS"] * 2e-3
    return df


def main() -> None:
    args = parse_args()

    log_file = args.log_file
    if not os.path.exists(log_file):
        print(f"Error: log file not found: {log_file}")
        sys.exit(1)

    log_dir = os.path.dirname(log_file) or "."
    output_dir = args.output_dir or os.path.join(log_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    log_basename = os.path.splitext(os.path.basename(log_file))[0]

    try:
        df = read_namd_energy_table(log_file)
    except Exception as error:
        print(f"Error: {error}")
        sys.exit(1)

    # Save CSV using snake_case columns where applicable
    csv_path = os.path.join(output_dir, f"{log_basename}.csv")
    df.to_csv(csv_path, index=False)
    print(f"Data saved: {csv_path}")

    # Determine x column
    if args.x_column is not None:
        x_col = args.x_column
    elif "time_ps" in df.columns:
        x_col = "time_ps"
    else:
        x_col = df.columns[0]

    # Determine columns to plot
    if args.plot_columns is not None and len(args.plot_columns) > 0:
        plot_columns = args.plot_columns
    else:
        plot_columns = [col for col in df.columns if col != x_col]

    for y_col in plot_columns:
        try:
            plt.figure(figsize=(10, 6))
            plt.plot(df[x_col], df[y_col])
            plt.xlabel(x_col)
            plt.ylabel(y_col)
            plt.title(f"{y_col} vs {x_col}")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()

            plot_file = os.path.join(
                output_dir, f"{log_basename}_{y_col.replace(' ', '_')}_plot.png"
            )
            plt.savefig(plot_file, dpi=300)
            print(f"Plot saved: {plot_file}")
            if not args.no_show:
                plt.show()
            plt.close()
        except Exception as plot_error:
            print(f"Warning: failed to plot column '{y_col}': {plot_error}")


if __name__ == "__main__":
    main()