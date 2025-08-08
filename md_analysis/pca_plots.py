import argparse
import os
import sys
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import pca


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run PCA on a trajectory selection and save CSVs and plots."
    )
    parser.add_argument("topology", help="Path to topology file (e.g., PSF/PDB)")
    parser.add_argument("trajectory", help="Path to trajectory file (e.g., DCD/XTC)")
    parser.add_argument(
        "--selection",
        default="resname CNT",
        help="Atom selection string for PCA (default: 'resname CNT')",
    )
    parser.add_argument(
        "--n-components",
        type=int,
        default=5,
        help="Number of principal components to transform and analyze (default: 5)",
    )
    parser.add_argument(
        "--no-align",
        action="store_true",
        help="Disable alignment before PCA (enabled by default)",
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


def maybe_show(no_show: bool) -> None:
    if not no_show:
        plt.show()
    plt.close()


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.topology):
        print(f"Error: topology file not found: {args.topology}")
        sys.exit(1)
    if not os.path.exists(args.trajectory):
        print(f"Error: trajectory file not found: {args.trajectory}")
        sys.exit(1)

    traj_dir = os.path.dirname(os.path.abspath(args.trajectory))
    output_dir = args.output_dir or os.path.join(traj_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(args.trajectory))[0]

    u = mda.Universe(args.topology, args.trajectory)

    align = not args.no_align
    pca_analysis = pca.PCA(
        u, select=args.selection, align=align, mean=None, n_components=None
    )
    pca_analysis.run()

    total_components = len(pca_analysis.results.p_components)
    n_components = min(args.n_components, total_components)

    atom_group = u.select_atoms(args.selection)
    transformed = pca_analysis.transform(atom_group, n_components=n_components)

    viz_components = min(3, n_components)
    df = pd.DataFrame(
        transformed[:, :viz_components],
        columns=[f"pc{i+1}" for i in range(viz_components)],
    )
    # Time in ps assuming trajectory.dt is in ps
    df["time_ps"] = df.index * u.trajectory.dt

    # Save transformed data
    df.to_csv(os.path.join(output_dir, f"{base_name}_pca_data.csv"), index=False)

    # Cumulative variance plot (first 10 PCs)
    plt.figure(figsize=(10, 6))
    plt.plot(pca_analysis.results.cumulated_variance[:10])
    plt.xlabel("Principal component")
    plt.ylabel("Cumulative variance")
    plt.title("Cumulative Variance by Principal Component")
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"{base_name}_pca_cumulative_variance_plot.png"),
        dpi=300,
        bbox_inches="tight",
    )
    maybe_show(args.no_show)

    # PairGrid of first 3 PCs colored by time
    g = sns.PairGrid(df[["time_ps"] + [f"pc{i+1}" for i in range(viz_components)]], hue="time_ps",
                     palette=sns.color_palette("viridis", n_colors=len(df)))
    g.map(plt.scatter, marker=".")
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"{base_name}_pca_pairgrid_plot.png"),
        dpi=300,
        bbox_inches="tight",
    )
    maybe_show(args.no_show)

    # PC1 vs PC2 scatter colored by time
    if viz_components >= 2:
        plt.figure(figsize=(10, 8))
        scatter = plt.scatter(df["pc1"], df["pc2"], c=df["time_ps"], cmap="viridis")
        cbar = plt.colorbar(scatter)
        cbar.set_label("Time (ps)")
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title(f"PCA of {base_name}")
        plt.tight_layout()
        plt.savefig(
            os.path.join(output_dir, f"{base_name}_pca_pc1_pc2_plot.png"),
            dpi=300,
            bbox_inches="tight",
        )
        maybe_show(args.no_show)

    # PC1 vs time
    plt.figure(figsize=(12, 6))
    plt.plot(df["time_ps"], df["pc1"])
    plt.xlabel("Time (ps)")
    plt.ylabel("PC1")
    plt.title("PC1 vs Time")
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"{base_name}_pca_pc1_vs_time_plot.png"),
        dpi=300,
        bbox_inches="tight",
    )
    maybe_show(args.no_show)

    # Cosine content per available component
    print("\n=== Measuring convergence with cosine content ===")
    print("Cosine content measures similarity to a cosine shape (0-1)")
    print("Values close to 1 may indicate poor sampling or random diffusion")
    print("Values below 0.7 do not indicate poor sampling")
    available_components = transformed.shape[1]
    cosine_contents: List[float] = []
    for i in range(available_components):
        cosine_contents.append(pca.cosine_content(transformed, i))

    cc_df = pd.DataFrame(
        {
            "pc": [f"pc{i+1}" for i in range(available_components)],
            "cosine_content": cosine_contents,
        }
    )
    cc_df.to_csv(
        os.path.join(output_dir, f"{base_name}_pca_cosine_content.csv"), index=False
    )

    plt.figure(figsize=(10, 6))
    plt.bar(cc_df["pc"], cc_df["cosine_content"])
    plt.axhline(y=0.7, color="r", linestyle="--", label="0.7 threshold")
    plt.xlabel("Principal Component")
    plt.ylabel("Cosine Content")
    plt.title("Cosine Content of Principal Components")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"{base_name}_pca_cosine_content_plot.png"),
        dpi=300,
        bbox_inches="tight",
    )
    maybe_show(args.no_show)

    # PC trajectories over time (first viz_components)
    melted = pd.melt(
        df,
        id_vars=["time_ps"],
        value_vars=[f"pc{i+1}" for i in range(viz_components)],
        var_name="pc",
        value_name="value",
    )

    g = sns.FacetGrid(melted, col="pc", height=4, aspect=1.2)
    g.map_dataframe(sns.lineplot, x="time_ps", y="value")
    g.set_axis_labels("Time (ps)", "Value")
    g.fig.suptitle("PC Trajectories Over Time", y=1.05)
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"{base_name}_pca_trajectories_plot.png"),
        dpi=300,
        bbox_inches="tight",
    )
    maybe_show(args.no_show)

    print(f"PCA analysis completed. Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
