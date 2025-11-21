"""
Test interpolation using actual AeoLiS configuration and run a single timestep comparison.
Loads configuration from external file and compares old vs new interpolation methods.
"""

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
import os

# Add aeolis to path
sys.path.insert(0, r'c:\Users\svries\Github\OE_aeolis-python')

import aeolis.inout
from aeolis.shear import WindShear

def interpolate_with_padding(x, y, z, xi, yi, z0_fill, x0, y0, istransect=False):
    """OLD method with padding"""
    # Rotation logic (simplified from shear.interpolate)
    dx = x[0,1] - x[0,0]
    dy = y[0,1] - y[0,0]
    angle = np.rad2deg(np.arctan(dy/dx))
    if dx <= 0 and dy<=0:
        angle += 180.
    
    # Simple rotation helper
    def rotate(x, y, alpha, origin):
        xr, yr = x - origin[0], y - origin[1]
        a = alpha / 180. * np.pi
        R = np.asmatrix([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        xy = np.concatenate((xr.reshape((-1,1)), yr.reshape((-1,1))), axis=1) * R
        return (np.asarray(xy[:,0].reshape(x.shape) + origin[0]),
                np.asarray(xy[:,1].reshape(y.shape) + origin[1]))
    
    x_rot, y_rot = rotate(x, y, angle, (x0, y0))
    xi_rot, yi_rot = rotate(xi, yi, angle, (x0, y0))
    
    if not np.all(sorted(y_rot[:,0]) == y_rot[:,0]) and not np.all(sorted(x_rot[0,:]) == x_rot[0,:]):
        x_rot, y_rot = rotate(x_rot, y_rot, 180, (x0, y0))
        xi_rot, yi_rot = rotate(xi_rot, yi_rot, 180, (x0, y0))
    
    xyi = np.concatenate((yi_rot.reshape((-1,1)), xi_rot.reshape((-1,1))), axis=1)
    
    if istransect:
        zi = np.interp(xi.flatten(), x.flatten(), z.flatten()).reshape(xi.shape)
    else:
        # OLD: with padding
        pad_w = max(x_rot.shape[0], x_rot.shape[1])
        x_pad = np.pad(x_rot, ((pad_w, pad_w), (pad_w, pad_w)), 'reflect', reflect_type='odd')
        y_pad = np.pad(y_rot, ((pad_w, pad_w), (pad_w, pad_w)), 'reflect', reflect_type='odd')
        z_pad = np.pad(z, ((pad_w, pad_w), (pad_w, pad_w)), 'edge')
        
        inter = scipy.interpolate.RegularGridInterpolator(
            (y_pad[:,0].copy(order='C'), x_pad[0,:].copy(order='C')), 
            z_pad, bounds_error=False, fill_value=z0_fill
        )
        zi = inter(xyi).reshape(xi.shape)
    
    return zi

def interpolate_no_padding(x, y, z, xi, yi, z0_fill, x0, y0, istransect=False):
    """NEW method without padding"""
    dx = x[0,1] - x[0,0]
    dy = y[0,1] - y[0,0]
    angle = np.rad2deg(np.arctan(dy/dx))
    if dx <= 0 and dy<=0:
        angle += 180.
    
    def rotate(x, y, alpha, origin):
        xr, yr = x - origin[0], y - origin[1]
        a = alpha / 180. * np.pi
        R = np.asmatrix([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        xy = np.concatenate((xr.reshape((-1,1)), yr.reshape((-1,1))), axis=1) * R
        return (np.asarray(xy[:,0].reshape(x.shape) + origin[0]),
                np.asarray(xy[:,1].reshape(y.shape) + origin[1]))
    
    x_rot, y_rot = rotate(x, y, angle, (x0, y0))
    xi_rot, yi_rot = rotate(xi, yi, angle, (x0, y0))
    
    if not np.all(sorted(y_rot[:,0]) == y_rot[:,0]) and not np.all(sorted(x_rot[0,:]) == x_rot[0,:]):
        x_rot, y_rot = rotate(x_rot, y_rot, 180, (x0, y0))
        xi_rot, yi_rot = rotate(xi_rot, yi_rot, 180, (x0, y0))
    
    xyi = np.concatenate((yi_rot.reshape((-1,1)), xi_rot.reshape((-1,1))), axis=1)
    
    if istransect:
        zi = np.interp(xi.flatten(), x.flatten(), z.flatten()).reshape(xi.shape)
    else:
        # NEW: no padding, fill_value=None for nearest extrapolation
        inter = scipy.interpolate.RegularGridInterpolator(
            (y_rot[:, 0].copy(order='C'), x_rot[0, :].copy(order='C')), 
            z.copy(order='C'), bounds_error=False, fill_value=None
        )
        zi = inter(xyi).reshape(xi.shape)
    
    return zi

# Load configuration
config_file = r"C:\Users\svries\Github\AeoLiS_Aquitaine\Carcans\Run_Batch\Batch_folder3\aeolis_U20_D225.txt"
print(f"Loading configuration from: {config_file}")

try:
    cfg = aeolis.inout.read_configfile(config_file)
    print("Configuration loaded successfully!")
    print(f"  Grid files:")
    print(f"    xgrid_file: {cfg.get('xgrid_file', 'N/A')}")
    print(f"    ygrid_file: {cfg.get('ygrid_file', 'N/A')}")
    print(f"    bed_file: {cfg.get('bed_file', 'N/A')}")
    print(f"  Shear parameters:")
    print(f"    dx: {cfg.get('dx', 'N/A')}")
    print(f"    dy: {cfg.get('dy', 'N/A')}")
    print(f"    L: {cfg.get('L', 'N/A')}")
    print(f"    l: {cfg.get('l', 'N/A')}")
    print(f"    z0: {cfg.get('z0', 'N/A')}")
    print(f"    buffer_width: {cfg.get('buffer_width', 100)}")
    
    # Load grids
    base_dir = os.path.dirname(config_file)
    x = np.loadtxt(os.path.join(base_dir, cfg['xgrid_file']))
    y = np.loadtxt(os.path.join(base_dir, cfg['ygrid_file']))
    z = np.loadtxt(os.path.join(base_dir, cfg['bed_file']))
    
    print(f"\nGrid loaded:")
    print(f"  Shape: {z.shape}")
    print(f"  X range: [{x.min():.2f}, {x.max():.2f}]")
    print(f"  Y range: [{y.min():.2f}, {y.max():.2f}]")
    print(f"  Z range: [{z.min():.2f}, {z.max():.2f}]")
    
    # Create WindShear instance to get computational grid
    ws = WindShear(
        x, y, z,
        dx=cfg.get('dx', 2.0),
        dy=cfg.get('dy', 2.0),
        L=cfg.get('L', 100.0),
        l=cfg.get('l', 10.0),
        z0=cfg.get('z0', 0.001),
        buffer_width=cfg.get('buffer_width', 100.0)
    )
    
    # Set computational grid for a wind direction
    udir = 225.0  # From filename
    ws.igrid = {'x': x, 'y': y, 'z': z}
    ws.set_computational_grid(udir)
    
    xi = ws.cgrid['xi']
    yi = ws.cgrid['yi']
    x0, y0 = ws.x0, ws.y0
    
    print(f"\nComputational grid:")
    print(f"  Shape: {xi.shape}")
    print(f"  Center: ({x0:.2f}, {y0:.2f})")
    
    # Interpolate with both methods
    print("\nInterpolating with OLD method (with padding)...")
    zi_old = interpolate_with_padding(x, y, z, xi, yi, 0.0, x0, y0, istransect=False)
    
    print("Interpolating with NEW method (no padding)...")
    zi_new = interpolate_no_padding(x, y, z, xi, yi, 0.0, x0, y0, istransect=False)
    
    # Compare
    diff = zi_new - zi_old
    abs_diff = np.abs(diff)
    rel_diff = np.abs(diff) / (np.abs(zi_old) + 1e-10) * 100
    
    print("\n" + "="*80)
    print("COMPARISON RESULTS")
    print("="*80)
    print(f"Max absolute difference: {abs_diff.max():.6e}")
    print(f"Mean absolute difference: {abs_diff.mean():.6e}")
    print(f"Max relative difference: {rel_diff.max():.2f}%")
    print(f"Mean relative difference: {rel_diff.mean():.2f}%")
    print(f"Points with >1% diff: {np.sum(rel_diff > 1.0)} / {rel_diff.size} ({100*np.sum(rel_diff > 1.0)/rel_diff.size:.2f}%)")
    print(f"Points with >5% diff: {np.sum(rel_diff > 5.0)} / {rel_diff.size} ({100*np.sum(rel_diff > 5.0)/rel_diff.size:.2f}%)")
    
    # Plot
    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    ax1 = fig.add_subplot(gs[0, 0])
    c1 = ax1.pcolormesh(x, y, z, cmap='terrain', shading='auto')
    ax1.set_title(f'Input Grid\n{z.shape}', fontsize=12)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    plt.colorbar(c1, ax=ax1, label='Elevation [m]')
    ax1.set_aspect('equal')
    
    ax2 = fig.add_subplot(gs[0, 1])
    c2 = ax2.pcolormesh(xi, yi, zi_old, cmap='terrain', shading='auto')
    ax2.set_title(f'OLD Method (with padding)\n{xi.shape}', fontsize=12)
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    plt.colorbar(c2, ax=ax2, label='Elevation [m]')
    ax2.set_aspect('equal')
    
    ax3 = fig.add_subplot(gs[0, 2])
    c3 = ax3.pcolormesh(xi, yi, zi_new, cmap='terrain', shading='auto')
    ax3.set_title(f'NEW Method (no padding)\n{xi.shape}', fontsize=12)
    ax3.set_xlabel('x [m]')
    ax3.set_ylabel('y [m]')
    plt.colorbar(c3, ax=ax3, label='Elevation [m]')
    ax3.set_aspect('equal')
    
    ax4 = fig.add_subplot(gs[1, 0])
    c4 = ax4.pcolormesh(xi, yi, diff, cmap='RdBu_r', shading='auto')
    ax4.set_title('Difference (NEW - OLD)', fontsize=12)
    ax4.set_xlabel('x [m]')
    ax4.set_ylabel('y [m]')
    plt.colorbar(c4, ax=ax4, label='Diff [m]')
    ax4.set_aspect('equal')
    
    ax5 = fig.add_subplot(gs[1, 1])
    c5 = ax5.pcolormesh(xi, yi, rel_diff, cmap='hot_r', shading='auto', vmin=0, vmax=min(5, rel_diff.max()))
    ax5.set_title('Relative Difference [%]', fontsize=12)
    ax5.set_xlabel('x [m]')
    ax5.set_ylabel('y [m]')
    plt.colorbar(c5, ax=ax5, label='% Diff')
    ax5.set_aspect('equal')
    
    ax6 = fig.add_subplot(gs[1, 2])
    stats_text = f"""STATISTICS

Max abs diff: {abs_diff.max():.6e} m
Mean abs diff: {abs_diff.mean():.6e} m

Max rel diff: {rel_diff.max():.2f}%
Mean rel diff: {rel_diff.mean():.4f}%

>1% diff: {100*np.sum(rel_diff > 1.0)/rel_diff.size:.2f}%
>5% diff: {100*np.sum(rel_diff > 5.0)/rel_diff.size:.2f}%

Identical: {100*np.sum(abs_diff < 1e-10)/rel_diff.size:.2f}%
"""
    ax6.text(0.1, 0.5, stats_text, ha='left', va='center', fontsize=11, family='monospace',
              bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    ax6.axis('off')
    
    plt.suptitle(f'Interpolation Comparison: {os.path.basename(config_file)}', fontsize=14, fontweight='bold')
    
    output_file = 'interpolation_comparison_carcans.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    plt.show()
    
    # Verdict
    print("\n" + "="*80)
    if abs_diff.max() < 1e-8:
        print("✓ IDENTICAL: Methods produce numerically identical results!")
    elif abs_diff.max() < 1e-4 and rel_diff.max() < 0.1:
        print("✓ EQUIVALENT: Methods produce effectively identical results")
    elif rel_diff.max() < 1.0:
        print("~ SIMILAR: Methods produce similar results with small differences")
    else:
        print("✗ DIFFERENT: Methods produce different results")
    print("="*80)
    
except FileNotFoundError as e:
    print(f"ERROR: File not found - {e}")
    print("\nPlease check that the configuration file exists and grid files are accessible.")
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
