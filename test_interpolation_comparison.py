"""
Test script to compare old padding-based interpolation vs new bounds-handling interpolation.
Uses barchan example data to verify identical output.
"""

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Read barchan grid data
def load_barchan_grid():
    """Load the barchan example grid files"""
    import os
    base_path = r"aeolis\examples\vanWesten2024\barchan_rotate"
    
    x = np.loadtxt(os.path.join(base_path, "x.grd"))
    y = np.loadtxt(os.path.join(base_path, "y.grd"))
    z = np.loadtxt(os.path.join(base_path, "z.grd"))
    
    return x, y, z


def interpolate_old_method(x, y, z, xi, yi, z0, x0, y0):
    """Original interpolation with padding (copied from old code)"""
    
    # First compute angle with horizontal
    dx = x[0,1] - x[0,0]
    dy = y[0,1] - y[0,0]

    angle = np.rad2deg(np.arctan(dy/dx))
    
    if dx <= 0 and dy<=0:
        angle += 180.
    
    # Rotate grids to align with horizontal
    xr = x - x0
    yr = y - y0
    a = angle / 180. * np.pi
    R = np.asmatrix([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    xy = np.concatenate((xr.reshape((-1,1)), yr.reshape((-1,1))), axis=1) * R
    x_rot = np.asarray(xy[:,0].reshape(x.shape) + x0)
    y_rot = np.asarray(xy[:,1].reshape(y.shape) + y0)
    
    xir = xi - x0
    yir = yi - y0
    xyi = np.concatenate((xir.reshape((-1,1)), yir.reshape((-1,1))), axis=1) * R
    xi_rot = np.asarray(xyi[:,0].reshape(xi.shape) + x0)
    yi_rot = np.asarray(xyi[:,1].reshape(yi.shape) + y0)
    
    # Rotate 180 deg if necessary
    if not np.all(sorted(y_rot[:,0]) == y_rot[:,0]) and not np.all(sorted(x_rot[0,:]) == x_rot[0,:]):
        xr2 = x_rot - x0
        yr2 = y_rot - y0
        a2 = 180 / 180. * np.pi
        R2 = np.asmatrix([[np.cos(a2), -np.sin(a2)], [np.sin(a2), np.cos(a2)]])
        xy2 = np.concatenate((xr2.reshape((-1,1)), yr2.reshape((-1,1))), axis=1) * R2
        x_rot = np.asarray(xy2[:,0].reshape(x.shape) + x0)
        y_rot = np.asarray(xy2[:,1].reshape(y.shape) + y0)
        
        xir2 = xi_rot - x0
        yir2 = yi_rot - y0
        xyi2 = np.concatenate((xir2.reshape((-1,1)), yir2.reshape((-1,1))), axis=1) * R2
        xi_rot = np.asarray(xyi2[:,0].reshape(xi.shape) + x0)
        yi_rot = np.asarray(xyi2[:,1].reshape(yi.shape) + y0)
    
    # Concatenate for interpolator
    xyi_coords = np.concatenate((yi_rot.reshape((-1,1)), xi_rot.reshape((-1,1))), axis=1)
    
    # OLD METHOD: With padding
    pad_w = np.maximum(np.shape(x_rot)[0], np.shape(x_rot)[1])
    x_pad = np.pad(x_rot, ((pad_w, pad_w), (pad_w, pad_w)), 'reflect', reflect_type='odd')
    y_pad = np.pad(y_rot, ((pad_w, pad_w), (pad_w, pad_w)), 'reflect', reflect_type='odd')
    z_pad = np.pad(z, ((pad_w, pad_w), (pad_w, pad_w)), 'edge')
    
    inter = scipy.interpolate.RegularGridInterpolator(
        (y_pad[:,0].copy(order='C'), x_pad[0,:].copy(order='C')), 
        z_pad, 
        bounds_error=False, 
        fill_value=z0
    )
    zi = inter(xyi_coords).reshape(xi.shape)
    
    return zi, x_rot, y_rot, xi_rot, yi_rot, x_pad, y_pad, z_pad


def interpolate_new_method(x, y, z, xi, yi, z0, x0, y0):
    """New interpolation without padding (current code)"""
    
    # First compute angle with horizontal
    dx = x[0,1] - x[0,0]
    dy = y[0,1] - y[0,0]

    angle = np.rad2deg(np.arctan(dy/dx))
    
    if dx <= 0 and dy<=0:
        angle += 180.
    
    # Rotate grids to align with horizontal
    xr = x - x0
    yr = y - y0
    a = angle / 180. * np.pi
    R = np.asmatrix([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    xy = np.concatenate((xr.reshape((-1,1)), yr.reshape((-1,1))), axis=1) * R
    x_rot = np.asarray(xy[:,0].reshape(x.shape) + x0)
    y_rot = np.asarray(xy[:,1].reshape(y.shape) + y0)
    
    xir = xi - x0
    yir = yi - y0
    xyi = np.concatenate((xir.reshape((-1,1)), yir.reshape((-1,1))), axis=1) * R
    xi_rot = np.asarray(xyi[:,0].reshape(xi.shape) + x0)
    yi_rot = np.asarray(xyi[:,1].reshape(yi.shape) + y0)
    
    # Rotate 180 deg if necessary
    if not np.all(sorted(y_rot[:,0]) == y_rot[:,0]) and not np.all(sorted(x_rot[0,:]) == x_rot[0,:]):
        xr2 = x_rot - x0
        yr2 = y_rot - y0
        a2 = 180 / 180. * np.pi
        R2 = np.asmatrix([[np.cos(a2), -np.sin(a2)], [np.sin(a2), np.cos(a2)]])
        xy2 = np.concatenate((xr2.reshape((-1,1)), yr2.reshape((-1,1))), axis=1) * R2
        x_rot = np.asarray(xy2[:,0].reshape(x.shape) + x0)
        y_rot = np.asarray(xy2[:,1].reshape(y.shape) + y0)
        
        xir2 = xi_rot - x0
        yir2 = yi_rot - y0
        xyi2 = np.concatenate((xir2.reshape((-1,1)), yir2.reshape((-1,1))), axis=1) * R2
        xi_rot = np.asarray(xyi2[:,0].reshape(xi.shape) + x0)
        yi_rot = np.asarray(xyi2[:,1].reshape(yi.shape) + y0)
    
    # Concatenate for interpolator
    xyi_coords = np.concatenate((yi_rot.reshape((-1,1)), xi_rot.reshape((-1,1))), axis=1)
    
    # NEW METHOD: Without padding, using fill_value=None for nearest extrapolation
    y_vec = y_rot[:, 0].copy(order='C')
    x_vec = x_rot[0, :].copy(order='C')
    z_c = z.copy(order='C')
    inter = scipy.interpolate.RegularGridInterpolator(
        (y_vec, x_vec), z_c, bounds_error=False, fill_value=None
    )
    zi = inter(xyi_coords).reshape(xi.shape)
    
    return zi, x_rot, y_rot, xi_rot, yi_rot


def create_computational_grid(x, y, buffer_width, dx, dy, udir=0):
    """Create a computational grid similar to WindShear.set_computational_grid"""
    x0, y0 = np.mean(x), np.mean(y)
    
    # Simple square grid covering input + buffer
    xmin, xmax = x.min() - buffer_width, x.max() + buffer_width
    ymin, ymax = y.min() - buffer_width, y.max() + buffer_width
    
    # Make it square
    width = max(xmax - xmin, ymax - ymin)
    
    xc = np.arange(x0 - width/2, x0 + width/2, dx)
    yc = np.arange(y0 - width/2, y0 + width/2, dy)
    xc, yc = np.meshgrid(xc, yc)
    
    return xc, yc, x0, y0


# Main comparison
print("Loading barchan grid data...")
x, y, z = load_barchan_grid()

print(f"Input grid shape: {z.shape}")
print(f"X range: [{x.min():.2f}, {x.max():.2f}]")
print(f"Y range: [{y.min():.2f}, {y.max():.2f}]")
print(f"Z range: [{z.min():.2f}, {z.max():.2f}]")

# Create computational grid
buffer_width = 50.0
dx_comp = 2.0
dy_comp = 2.0
xi, yi, x0, y0 = create_computational_grid(x, y, buffer_width, dx_comp, dy_comp)

print(f"\nComputational grid shape: {xi.shape}")
print(f"Grid center: ({x0:.2f}, {y0:.2f})")

# Interpolate with both methods
print("\nInterpolating with OLD method (with padding)...")
zi_old, x_rot, y_rot, xi_rot, yi_rot, x_pad, y_pad, z_pad = interpolate_old_method(
    x, y, z, xi, yi, z0=0.0, x0=x0, y0=y0
)

print(f"Padded grid shape: {z_pad.shape}")
print(f"Pad width: {(z_pad.shape[0] - z.shape[0]) // 2}")

print("\nInterpolating with NEW method (no padding)...")
zi_new, x_rot_new, y_rot_new, xi_rot_new, yi_rot_new = interpolate_new_method(
    x, y, z, xi, yi, z0=0.0, x0=x0, y0=y0
)

# Compare results
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

print(f"\nOld method - zi range: [{zi_old.min():.4f}, {zi_old.max():.4f}]")
print(f"New method - zi range: [{zi_new.min():.4f}, {zi_new.max():.4f}]")

# Create comprehensive plots
fig = plt.figure(figsize=(20, 12))
gs = GridSpec(3, 4, figure=fig, hspace=0.3, wspace=0.3)

# Row 1: Source grids
ax1 = fig.add_subplot(gs[0, 0])
c1 = ax1.pcolormesh(x, y, z, cmap='terrain', shading='auto')
ax1.set_title('Original Input Grid (z)')
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
plt.colorbar(c1, ax=ax1, label='Elevation [m]')
ax1.set_aspect('equal')

ax2 = fig.add_subplot(gs[0, 1])
c2 = ax2.pcolormesh(x_rot, y_rot, z, cmap='terrain', shading='auto')
ax2.plot(xi_rot[::5, ::5].flatten(), yi_rot[::5, ::5].flatten(), 'r.', ms=1, alpha=0.3)
ax2.set_title('Rotated Source Grid + Target Points')
ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
plt.colorbar(c2, ax=ax2, label='Elevation [m]')
ax2.set_aspect('equal')

ax3 = fig.add_subplot(gs[0, 2])
c3 = ax3.pcolormesh(x_pad, y_pad, z_pad, cmap='terrain', shading='auto')
ax3.plot(x_rot.flatten(), y_rot.flatten(), 'k.', ms=0.5, alpha=0.5, label='Original grid')
ax3.plot(xi_rot[::5, ::5].flatten(), yi_rot[::5, ::5].flatten(), 'r.', ms=1, alpha=0.3, label='Target points')
ax3.set_title('OLD: Padded Grid')
ax3.set_xlabel('x [m]')
ax3.set_ylabel('y [m]')
plt.colorbar(c3, ax=ax3, label='Elevation [m]')
ax3.legend(loc='upper right', fontsize=8)
ax3.set_aspect('equal')

ax4 = fig.add_subplot(gs[0, 3])
ax4.text(0.5, 0.5, f'Pad width: {(z_pad.shape[0] - z.shape[0]) // 2}\n\n'
                    f'Original: {z.shape}\n'
                    f'Padded: {z_pad.shape}\n\n'
                    f'Padding modes:\n'
                    f'  x, y: reflect (odd)\n'
                    f'  z: edge',
         ha='center', va='center', fontsize=11, family='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
ax4.axis('off')
ax4.set_title('Padding Info')

# Row 2: Interpolated results
ax5 = fig.add_subplot(gs[1, 0])
c5 = ax5.pcolormesh(xi, yi, zi_old, cmap='terrain', shading='auto')
ax5.set_title('OLD Method: Interpolated Result')
ax5.set_xlabel('x [m]')
ax5.set_ylabel('y [m]')
plt.colorbar(c5, ax=ax5, label='Elevation [m]')
ax5.set_aspect('equal')

ax6 = fig.add_subplot(gs[1, 1])
c6 = ax6.pcolormesh(xi, yi, zi_new, cmap='terrain', shading='auto')
ax6.set_title('NEW Method: Interpolated Result')
ax6.set_xlabel('x [m]')
ax6.set_ylabel('y [m]')
plt.colorbar(c6, ax=ax6, label='Elevation [m]')
ax6.set_aspect('equal')

ax7 = fig.add_subplot(gs[1, 2])
c7 = ax7.pcolormesh(xi, yi, diff, cmap='RdBu_r', shading='auto')
ax7.set_title('Difference (NEW - OLD)')
ax7.set_xlabel('x [m]')
ax7.set_ylabel('y [m]')
plt.colorbar(c7, ax=ax7, label='Diff [m]')
ax7.set_aspect('equal')

ax8 = fig.add_subplot(gs[1, 3])
c8 = ax8.pcolormesh(xi, yi, rel_diff, cmap='hot_r', shading='auto', vmin=0, vmax=5)
ax8.set_title('Relative Difference [%]')
ax8.set_xlabel('x [m]')
ax8.set_ylabel('y [m]')
plt.colorbar(c8, ax=ax8, label='% Diff')
ax8.set_aspect('equal')

# Row 3: Statistics and profiles
ax9 = fig.add_subplot(gs[2, 0:2])
mid_row = zi_old.shape[0] // 2
ax9.plot(xi[mid_row, :], zi_old[mid_row, :], 'b-', label='OLD method', linewidth=2)
ax9.plot(xi[mid_row, :], zi_new[mid_row, :], 'r--', label='NEW method', linewidth=2)
ax9.set_title(f'Cross-section at y={yi[mid_row, 0]:.1f} m')
ax9.set_xlabel('x [m]')
ax9.set_ylabel('Elevation [m]')
ax9.legend()
ax9.grid(True, alpha=0.3)

ax10 = fig.add_subplot(gs[2, 2])
ax10.hist(abs_diff.flatten(), bins=50, edgecolor='black', alpha=0.7)
ax10.set_xlabel('Absolute Difference [m]')
ax10.set_ylabel('Frequency')
ax10.set_title('Distribution of Abs. Differences')
ax10.set_yscale('log')
ax10.grid(True, alpha=0.3)

ax11 = fig.add_subplot(gs[2, 3])
stats_text = f"""STATISTICS

Max abs diff: {abs_diff.max():.6e} m
Mean abs diff: {abs_diff.mean():.6e} m
Median abs diff: {np.median(abs_diff):.6e} m

Max rel diff: {rel_diff.max():.2f}%
Mean rel diff: {rel_diff.mean():.4f}%

Points >1% diff: {100*np.sum(rel_diff > 1.0)/rel_diff.size:.2f}%
Points >5% diff: {100*np.sum(rel_diff > 5.0)/rel_diff.size:.2f}%

Identical points: {100*np.sum(abs_diff < 1e-10)/rel_diff.size:.2f}%
"""
ax11.text(0.1, 0.5, stats_text, ha='left', va='center', fontsize=10, family='monospace',
          bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
ax11.axis('off')

plt.suptitle('Interpolation Method Comparison: OLD (with padding) vs NEW (nearest extrapolation)', 
             fontsize=14, fontweight='bold')

output_file = 'interpolation_comparison.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\nPlot saved to: {output_file}")
plt.show()

# Final verdict
print("\n" + "="*80)
if abs_diff.max() < 1e-8:
    print("✓ IDENTICAL: Methods produce numerically identical results!")
elif abs_diff.max() < 1e-4 and rel_diff.max() < 0.1:
    print("✓ EQUIVALENT: Methods produce effectively identical results (differences negligible)")
elif rel_diff.max() < 1.0:
    print("~ SIMILAR: Methods produce similar results with small differences")
else:
    print("✗ DIFFERENT: Methods produce different results - review needed!")
print("="*80)
