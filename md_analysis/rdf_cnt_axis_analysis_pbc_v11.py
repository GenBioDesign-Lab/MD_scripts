#!/usr/bin/env python3
"""
CNT-Lipid Radial Distribution Function with Periodic Boundary Conditions
========================================================================

Calculates 2D RDF from CNT axis to lipid carbons using PBC.
Only applies periodic replication when analysis range exceeds box boundaries.

Version: 11
"""

import sys
import argparse
from math import ceil
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda


def principal_axis(coords):
    """Calculate principal axis using PCA."""
    centered = coords - coords.mean(0)
    _, vecs = np.linalg.eigh(np.cov(centered.T))
    axis = vecs[:, -1]
    return axis / np.linalg.norm(axis)


def radial_distance(points, center, axis):
    """Calculate radial distance from axis (perpendicular component)."""
    vec = points - center
    proj = (vec @ axis)[:, None] * axis
    perp = vec - proj
    return np.linalg.norm(perp, axis=1)


def compute_rdf(universe, cnt_sel="resname CNT*", lip_sel="resname *PC *PE and name C*",
                r_max=60.0, bins=300, stride=1, cnt_radius=5.0):
    # Select atoms
    cnt_atoms = universe.select_atoms(cnt_sel)
    lip_atoms = universe.select_atoms(lip_sel)
    
    if not cnt_atoms.n_atoms or not lip_atoms.n_atoms:
        raise ValueError("No atoms found - check selection strings")
    
    print(f"Selected {cnt_atoms.n_atoms} CNT atoms, {lip_atoms.n_atoms} lipid atoms")
    
    # Initialize histogram
    dr = r_max / bins
    bin_edges = np.linspace(0, r_max, bins + 1)
    bin_centers = bin_edges[:-1] + dr / 2
    hist_total = np.zeros(bins, dtype=int)
    
    # Process trajectory
    pbc_frames = 0
    for i, ts in enumerate(universe.trajectory[::stride]):
        Lx, Ly, Lz = ts.dimensions[:3]
        box_half = min(Lx, Ly) / 2
        
        # Smart PBC: only replicate when needed
        if r_max > box_half:
            # Calculate required replication
            nx = min(ceil((r_max - box_half) / Lx), 3)
            ny = min(ceil((r_max - box_half) / Ly), 3)
            
            # Replicate lipid positions in x,y plane
            replicated = []
            for dx in range(-nx, nx + 1):
                for dy in range(-ny, ny + 1):
                    shift = np.array([dx * Lx, dy * Ly, 0.0])
                    replicated.append(lip_atoms.positions + shift)
            lip_pos = np.vstack(replicated)
            pbc_frames += 1
        else:
            lip_pos = lip_atoms.positions
        
        # Calculate CNT axis and center
        cnt_center = cnt_atoms.positions.mean(0)
        axis = principal_axis(cnt_atoms.positions)
        
        # Compute radial distances and histogram
        distances = radial_distance(lip_pos, cnt_center, axis)
        distances = distances[distances <= r_max]
        
        hist, _ = np.histogram(distances, bins=bin_edges)
        hist_total += hist
    
    n_frames = len(range(0, len(universe.trajectory), stride))
    print(f"Processed {n_frames} frames ({pbc_frames} with PBC)")
    
    # Normalization
    box_area = universe.trajectory.ts.dimensions[0] * universe.trajectory.ts.dimensions[1]
    accessible_area = box_area - np.pi * cnt_radius**2
    density = lip_atoms.n_atoms / accessible_area
    
    # Calculate g(r)
    g_r = np.zeros_like(hist_total, dtype=float)
    for i, count in enumerate(hist_total):
        r = bin_centers[i]
        if r > 0:
            shell_area = 2 * np.pi * r * dr
            expected = density * shell_area * n_frames
            g_r[i] = count / expected if expected > 0 else 0.0
    
    print(f"RDF completed: max g(r) = {np.max(g_r):.2f}")
    return bin_centers, g_r, hist_total


def plot_results(r, g_r, raw_hist, output_prefix="rdf_results"):
    """Create publication-quality plots."""
    
    # Two-panel plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Raw histogram
    ax1.plot(r, raw_hist, 'b-', lw=2, color='steelblue')
    ax1.set_xlabel('Distance (Å)')
    ax1.set_ylabel('Counts')
    ax1.set_title('Raw Distance Histogram')
    ax1.grid(True, alpha=0.3)
    
    # Normalized RDF
    ax2.plot(r, g_r, 'r-', lw=2)
    ax2.axhline(1, color='k', ls='--', alpha=0.7, label='Bulk density')
    ax2.set_xlabel('Distance (Å)')
    ax2.set_ylabel('g(r)')
    ax2.set_title('Radial Distribution Function')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_combined.png", dpi=300, bbox_inches='tight')
    
    # Individual RDF plot
    plt.figure(figsize=(8, 6))
    plt.plot(r, g_r, 'r-', lw=2)
    plt.axhline(1, color='k', ls='--', alpha=0.7)
    plt.xlabel('Distance (Å)')
    plt.ylabel('g(r)')
    plt.title('CNT-Lipid Radial Distribution Function')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_rdf.png", dpi=300, bbox_inches='tight')
    
    print(f"Plots saved: {output_prefix}_combined.png, {output_prefix}_rdf.png")


def main():
    """Main analysis function."""
    parser = argparse.ArgumentParser(
        description="Calculate CNT-lipid RDF with smart PBC handling",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("topology", help="Topology file (PSF/GRO/PDB)")
    parser.add_argument("trajectory", help="Trajectory file (XTC/TRR/DCD)")
    parser.add_argument("-r", "--rmax", type=float, default=60.0,
                       help="Maximum distance (Å)")
    parser.add_argument("-b", "--bins", type=int, default=300,
                       help="Number of bins")
    parser.add_argument("-s", "--stride", type=int, default=1,
                       help="Frame stride")
    parser.add_argument("--cnt-sel", default="resname CNT*",
                       help="CNT selection string")
    parser.add_argument("--lip-sel", default="resname *PC *PE and name C*",
                       help="Lipid selection string")
    parser.add_argument("-o", "--output", default="rdf_results",
                       help="Output file prefix")
    
    args = parser.parse_args()
    
    # Load trajectory
    print(f"Loading: {args.topology}, {args.trajectory}")
    u = mda.Universe(args.topology, args.trajectory)
    print(f"System: {u.atoms.n_atoms} atoms, {len(u.trajectory)} frames")
    
    # Compute RDF
    r, g_r, raw_hist = compute_rdf(
        u, args.cnt_sel, args.lip_sel,
        r_max=args.rmax, bins=args.bins, stride=args.stride
    )
    
    # Save data
    data = np.column_stack([r, raw_hist, g_r])
    np.savetxt(f"{args.output}.dat", data,
               header="Distance(Å)  RawCounts  g(r)",
               fmt="%.6f")
    
    # Create plots
    plot_results(r, g_r, raw_hist, args.output)
    
    print(f"Analysis complete. Data saved to {args.output}.dat")


if __name__ == "__main__":
    main() 