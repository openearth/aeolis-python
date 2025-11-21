"""
Simple test to verify rotation and interpolation work correctly with 225° wind direction.
"""

import numpy as np
import sys
sys.path.insert(0, r'c:\Users\svries\Github\OE_aeolis-python')

from aeolis.shear import WindShear

# Create simple test grids
x = np.arange(0, 100, 10)
y = np.arange(0, 100, 10)
x, y = np.meshgrid(x, y)

# Simple elevation: a slope in x-direction
z = 0.1 * x + 5.0

print("="*80)
print("TESTING ROTATION WITH 225° WIND DIRECTION")
print("="*80)

# Create WindShear instance
ws = WindShear(x, y, z, dx=5.0, dy=5.0, L=100.0, l=10.0, z0=0.001, buffer_width=50.0)

# Set up for 225° wind direction
udir = 225.0
ws.igrid = {'x': x, 'y': y, 'z': z}
ws.set_computational_grid(udir)

print(f"\nOriginal grid:")
print(f"  Shape: {z.shape}")
print(f"  X range: [{x.min():.1f}, {x.max():.1f}]")
print(f"  Y range: [{y.min():.1f}, {y.max():.1f}]")
print(f"  Grid center (x0, y0): ({ws.x0:.1f}, {ws.y0:.1f})")

print(f"\nComputational grid (BEFORE rotation to wind direction):")
print(f"  Shape: {ws.cgrid['xi'].shape}")
print(f"  XI range: [{ws.cgrid['xi'].min():.1f}, {ws.cgrid['xi'].max():.1f}]")
print(f"  YI range: [{ws.cgrid['yi'].min():.1f}, {ws.cgrid['yi'].max():.1f}]")

# Now rotate computational grid to wind direction (as done in __call__)
u_angle = 270. - udir  # = 270 - 225 = 45°
print(f"\nWind direction: {udir}°")
print(f"Rotation angle (u_angle): {u_angle}° (= 270° - {udir}°)")

gc_x, gc_y = WindShear.rotate(ws.cgrid['xi'], ws.cgrid['yi'], -u_angle, origin=(ws.x0, ws.y0))

print(f"\nComputational grid (AFTER rotation by {-u_angle}°):")
print(f"  GC_X range: [{gc_x.min():.1f}, {gc_x.max():.1f}]")
print(f"  GC_Y range: [{gc_y.min():.1f}, {gc_y.max():.1f}]")

# Test interpolation
print(f"\n" + "="*80)
print("TESTING INTERPOLATION")
print("="*80)

zi = ws.interpolate(x, y, z, gc_x, gc_y, 0.0)

print(f"\nInterpolated grid:")
print(f"  Shape: {zi.shape}")
print(f"  Z range: [{zi.min():.2f}, {zi.max():.2f}]")
print(f"  Original Z range: [{z.min():.2f}, {z.max():.2f}]")

# Check a few sample points
print(f"\nSample interpolated values:")
ny, nx = zi.shape
for i, j in [(0, 0), (ny//2, nx//2), (ny-1, nx-1)]:
    print(f"  zi[{i},{j}] = {zi[i,j]:.2f} at (x={gc_x[i,j]:.1f}, y={gc_y[i,j]:.1f})")

print("\n" + "="*80)
print("✓ Rotation and interpolation completed successfully!")
print("  The computational grid was rotated by {:.1f}° for wind direction {}°".format(-u_angle, udir))
print("="*80)
