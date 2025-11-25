"""
Utility functions and constants for AeoLiS GUI.

This module contains:
- Constants for visualization parameters
- File path resolution utilities
- Time unit conversion utilities
- Data extraction utilities
- Hillshade computation
"""

import os
import numpy as np
import math


# ============================================================================
# Constants
# ============================================================================

# Hillshade parameters
HILLSHADE_AZIMUTH = 155.0
HILLSHADE_ALTITUDE = 5.0
HILLSHADE_AMBIENT = 0.35

# Time unit conversion thresholds (in seconds)
TIME_UNIT_THRESHOLDS = {
    'seconds': (0, 300),        # < 5 minutes
    'minutes': (300, 7200),     # 5 min to 2 hours
    'hours': (7200, 172800),    # 2 hours to 2 days
    'days': (172800, 7776000),  # 2 days to ~90 days
    'years': (7776000, float('inf'))  # >= 90 days
}

TIME_UNIT_DIVISORS = {
    'seconds': 1.0,
    'minutes': 60.0,
    'hours': 3600.0,
    'days': 86400.0,
    'years': 365.25 * 86400.0
}

# Visualization parameters
OCEAN_DEPTH_THRESHOLD = -0.5
OCEAN_DISTANCE_THRESHOLD = 200
SUBSAMPLE_RATE_DIVISOR = 25  # For quiver plot subsampling

# NetCDF coordinate and metadata variables to exclude from plotting
NC_COORD_VARS = {
    'x', 'y', 's', 'n', 'lat', 'lon', 'time', 'layers', 'fractions',
    'x_bounds', 'y_bounds', 'lat_bounds', 'lon_bounds', 'time_bounds',
    'crs', 'nv', 'nv2'
}

# Variable visualization configuration
VARIABLE_LABELS = {
    'zb': 'Elevation (m)',
    'zb+rhoveg': 'Vegetation-shaded Topography',
    'ustar': 'Shear Velocity (m/s)',
    'ustar quiver': 'Shear Velocity Vectors',
    'ustars': 'Shear Velocity S-component (m/s)',
    'ustarn': 'Shear Velocity N-component (m/s)',
    'zs': 'Surface Elevation (m)',
    'zsep': 'Separation Elevation (m)',
    'Ct': 'Sediment Concentration (kg/m²)',
    'Cu': 'Equilibrium Concentration (kg/m²)',
    'q': 'Sediment Flux (kg/m/s)',
    'qs': 'Sediment Flux S-component (kg/m/s)',
    'qn': 'Sediment Flux N-component (kg/m/s)',
    'pickup': 'Sediment Entrainment (kg/m²)',
    'uth': 'Threshold Shear Velocity (m/s)',
    'w': 'Fraction Weight (-)',
}

VARIABLE_TITLES = {
    'zb': 'Bed Elevation',
    'zb+rhoveg': 'Bed Elevation with Vegetation (Shaded)',
    'ustar': 'Shear Velocity',
    'ustar quiver': 'Shear Velocity Vector Field',
    'ustars': 'Shear Velocity (S-component)',
    'ustarn': 'Shear Velocity (N-component)',
    'zs': 'Surface Elevation',
    'zsep': 'Separation Elevation',
    'Ct': 'Sediment Concentration',
    'Cu': 'Equilibrium Concentration',
    'q': 'Sediment Flux',
    'qs': 'Sediment Flux (S-component)',
    'qn': 'Sediment Flux (N-component)',
    'pickup': 'Sediment Entrainment',
    'uth': 'Threshold Shear Velocity',
    'w': 'Fraction Weight',
}


# ============================================================================
# Utility Functions
# ============================================================================

def resolve_file_path(file_path, base_dir):
    """
    Resolve a file path relative to a base directory.
    
    Parameters
    ----------
    file_path : str
        The file path to resolve (can be relative or absolute)
    base_dir : str
        The base directory for relative paths
        
    Returns
    -------
    str
        Absolute path to the file, or None if file_path is empty
    """
    if not file_path:
        return None
    if os.path.isabs(file_path):
        return file_path
    return os.path.join(base_dir, file_path)


def make_relative_path(file_path, base_dir):
    """
    Make a file path relative to a base directory if possible.
    
    Parameters
    ----------
    file_path : str
        The absolute file path
    base_dir : str
        The base directory
        
    Returns
    -------
    str
        Relative path if possible and not too many levels up, otherwise absolute path
    """
    try:
        rel_path = os.path.relpath(file_path, base_dir)
        # Only use relative path if it doesn't go up too many levels
        parent_dir = os.pardir + os.sep + os.pardir + os.sep
        if not rel_path.startswith(parent_dir):
            return rel_path
    except (ValueError, TypeError):
        # Different drives on Windows or invalid path
        pass
    return file_path


def determine_time_unit(duration_seconds):
    """
    Determine appropriate time unit based on simulation duration.
    
    Parameters
    ----------
    duration_seconds : float
        Duration in seconds
        
    Returns
    -------
    tuple
        (time_unit_name, divisor) for converting seconds to the chosen unit
    """
    for unit_name, (lower, upper) in TIME_UNIT_THRESHOLDS.items():
        if lower <= duration_seconds < upper:
            return (unit_name, TIME_UNIT_DIVISORS[unit_name])
    # Default to years if duration is very large
    return ('years', TIME_UNIT_DIVISORS['years'])


def extract_time_slice(data, time_idx):
    """
    Extract a time slice from variable data, handling different dimensionalities.
    
    Parameters
    ----------
    data : ndarray
        Data array (3D or 4D with time dimension)
    time_idx : int
        Time index to extract
        
    Returns
    -------
    ndarray
        2D slice at the given time index
        
    Raises
    ------
    ValueError
        If data dimensionality is unexpected
    """
    if data.ndim == 4:
        # (time, n, s, fractions) - average across fractions
        return data[time_idx, :, :, :].mean(axis=2)
    elif data.ndim == 3:
        # (time, n, s)
        return data[time_idx, :, :]
    else:
        raise ValueError(f"Unexpected data dimensionality: {data.ndim}. Expected 3D or 4D array.")


def apply_hillshade(z2d, x1d, y1d, az_deg=HILLSHADE_AZIMUTH, alt_deg=HILLSHADE_ALTITUDE):
    """
    Compute a simple hillshade (0–1) for 2D elevation array.
    Uses safe gradient computation and normalization.
    Adapted from Anim2D_ShadeVeg.py
    
    Parameters
    ----------
    z2d : ndarray
        2D elevation data array
    x1d : ndarray
        1D x-coordinate array
    y1d : ndarray
        1D y-coordinate array
    az_deg : float, optional
        Azimuth angle in degrees (default: HILLSHADE_AZIMUTH)
    alt_deg : float, optional
        Altitude angle in degrees (default: HILLSHADE_ALTITUDE)
        
    Returns
    -------
    ndarray
        Hillshade values between 0 and 1
        
    Raises
    ------
    ValueError
        If z2d is not a 2D array
    """
    z = np.asarray(z2d, dtype=float)
    if z.ndim != 2:
        raise ValueError("apply_hillshade expects a 2D array")

    x1 = np.asarray(x1d).ravel()
    y1 = np.asarray(y1d).ravel()

    eps = 1e-8
    dx = np.mean(np.diff(x1)) if x1.size > 1 else 1.0
    dy = np.mean(np.diff(y1)) if y1.size > 1 else 1.0
    dx = 1.0 if abs(dx) < eps else dx
    dy = 1.0 if abs(dy) < eps else dy

    dz_dy, dz_dx = np.gradient(z, dy, dx)

    nx, ny, nz = -dz_dx, -dz_dy, np.ones_like(z)
    norm = np.sqrt(nx * nx + ny * ny + nz * nz)
    norm = np.where(norm < eps, eps, norm)
    nx, ny, nz = nx / norm, ny / norm, nz / norm

    az = math.radians(az_deg)
    alt = math.radians(alt_deg)
    lx = math.cos(alt) * math.cos(az)
    ly = math.cos(alt) * math.sin(az)
    lz = math.sin(alt)

    illum = np.clip(nx * lx + ny * ly + nz * lz, 0.0, 1.0)
    shaded = HILLSHADE_AMBIENT + (1.0 - HILLSHADE_AMBIENT) * illum  # ambient term
    return np.clip(shaded, 0.0, 1.0)
