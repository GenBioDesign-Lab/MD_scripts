import argparse
import os
import sys

import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute per-residue RMSF from a trajectory and save CSV and plot."
    )
    parser.add_argument("topology", help="Path to topology file (e.g., PSF/PDB)")
    parser.add_argument("trajectory", help="Path to trajectory file (e.g., DCD/XTC)")
    parser.add_argument(
        "--selection",
        default="name CA",
        help="Atom selection for RMSF (default: 'name CA')",
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

    rmsf = RMSF(atom_group)
    rmsf.run()

    df = pd.DataFrame({
        "residue_index": atom_group.resids,
        "rmsf_angstrom": rmsf.results.rmsf,
    })

    csv_path = os.path.join(output_dir, f"{traj_base}_rmsf.csv")
    df.to_csv(csv_path, index=False)

    plt.figure(figsize=(10, 6))
    plt.plot(df["residue_index"], df["rmsf_angstrom"])
    plt.xlabel("Residue Index")
    plt.ylabel("RMSF (Å)")
    plt.title(f"Per-residue RMSF - {traj_base}")
    plt.tight_layout()

    plot_path = os.path.join(output_dir, f"{traj_base}_rmsf_plot.png")
    plt.savefig(plot_path, dpi=300)
    if not args.no_show:
        plt.show()
    plt.close()


if __name__ == "__main__":
    main()