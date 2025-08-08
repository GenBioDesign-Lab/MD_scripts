import argparse
import os
import sys
from typing import List

import MDAnalysis as mda
import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute radius of gyration from a trajectory and save CSV and plot."
    )
    parser.add_argument("topology", help="Path to topology file (e.g., PSF/PDB)")
    parser.add_argument("trajectory", help="Path to trajectory file (e.g., DCD/XTC)")
    parser.add_argument(
        "--selection",
        default="protein",
        help="Atom selection string (default: 'protein')",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory to write outputs. Default: '<trajectory_dir>/output'",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not display plots interactively",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.topology):
        print(f"Error: topology file not found: {args.topology}")
        sys.exit(1)
    if not os.path.exists(args.trajectory):
        print(f"Error: trajectory file not found: {args.trajectory}")
        sys.exit(1)

    traj_dir = os.path.dirname(args.trajectory)
    output_dir = args.output_dir or os.path.join(traj_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    traj_base = os.path.splitext(os.path.basename(args.trajectory))[0]

    u = mda.Universe(args.topology, args.trajectory)
    atom_group = u.select_atoms(args.selection)

    time_ps: List[float] = []
    rgyr_angstrom: List[float] = []

    for ts in u.trajectory:
        time_ps.append(ts.time)
        rgyr_angstrom.append(atom_group.radius_of_gyration())

    df = pd.DataFrame({
        "time_ps": time_ps,
        "rgyr_angstrom": rgyr_angstrom,
    })

    csv_path = os.path.join(output_dir, f"{traj_base}_rgyr.csv")
    df.to_csv(csv_path, index=False)

    plt.figure(figsize=(10, 6))
    plt.plot(df["time_ps"], df["rgyr_angstrom"])
    plt.xlabel("Time (ps)")
    plt.ylabel("Radius of Gyration (Å)")
    plt.title(f"Radius of Gyration - {traj_base}")
    plt.tight_layout()

    plot_path = os.path.join(output_dir, f"{traj_base}_rgyr_plot.png")
    plt.savefig(plot_path, dpi=300)
    if not args.no_show:
        plt.show()
    plt.close()


if __name__ == "__main__":
    main()