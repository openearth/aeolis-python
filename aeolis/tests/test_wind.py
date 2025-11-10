"""
This file is part of AeoLiS test suite. It uses the Pytest framework.

Tests for wind module functionality, specifically for wind shear velocity components.
"""

import numpy as np
import pytest


class TestWindShearVelocityComponents:
    """Test wind shear velocity component calculations"""

    def test_ustars0_ustarn0_components_pure_s_direction(self):
        """Test that ustars0 and ustarn0 correctly represent directional components.
        
        For wind from 270 degrees (nautical, wind from west) with alfa=0:
        - Wind should be in pure s-direction (eastward)
        - ustars should be positive
        - ustarn should be 0
        - ustars0 should equal ustars (not the magnitude ustar)
        - ustarn0 should equal ustarn (which is 0, not the magnitude ustar)
        """
        from aeolis.wind import interpolate, initialize
        from aeolis.model import ModelState
        
        # Create minimal model state
        s = ModelState()
        
        # Grid setup - simple 2x2 grid
        nx, ny = 3, 3
        s['x'] = np.zeros((ny, nx))
        s['y'] = np.zeros((ny, nx))
        s['zb'] = np.zeros((ny, nx))
        s['uw'] = np.zeros((ny, nx))
        s['udir'] = np.zeros((ny, nx))
        s['uws'] = np.zeros((ny, nx))
        s['uwn'] = np.zeros((ny, nx))
        s['ustars'] = np.zeros((ny, nx))
        s['ustarn'] = np.zeros((ny, nx))
        s['ustar'] = np.zeros((ny, nx))
        s['ustars0'] = np.zeros((ny, nx))
        s['ustarn0'] = np.zeros((ny, nx))
        s['ustar0'] = np.zeros((ny, nx))
        s['tau'] = np.zeros((ny, nx))
        s['taus'] = np.zeros((ny, nx))
        s['taun'] = np.zeros((ny, nx))
        s['tau0'] = np.zeros((ny, nx))
        s['taus0'] = np.zeros((ny, nx))
        s['taun0'] = np.zeros((ny, nx))
        s['hveg'] = np.zeros((ny, nx))
        
        # Parameters
        p = {}
        p['alfa'] = 0.0  # Grid orientation (0 = aligned with north)
        p['kappa'] = 0.4  # von Karman constant
        p['z'] = 10.0  # Reference height (m)
        p['k'] = 0.001  # Roughness length (m)
        p['method_roughness'] = 'constant'
        p['rhoa'] = 1.25  # Air density (kg/m³)
        p['process_wind'] = False  # We'll set wind manually
        p['wind_file'] = None
        
        # Set constant wind from 270 degrees (from west, blowing east)
        # In nautical convention: 270° means wind FROM the west
        s['uw'][:, :] = 15.0  # 15 m/s
        s['udir'][:, :] = 270.0  # degrees nautical
        
        # Run wind interpolation (which calculates ustars, ustarn, etc.)
        s = interpolate(s, p, t=0.0)
        
        # For wind from 270° with alfa=0:
        # uws = -uw * sin((-alfa + udir) * pi/180) = -15 * sin(270°) = -15 * (-1) = 15
        # uwn = -uw * cos((-alfa + udir) * pi/180) = -15 * cos(270°) = -15 * 0 = 0
        
        # Therefore:
        # ustars should be positive (proportional to uws)
        # ustarn should be 0 (proportional to uwn)
        
        # Check that ustars is positive and ustarn is approximately 0
        assert np.all(s['ustars'] > 0), "ustars should be positive for wind from 270°"
        assert np.allclose(s['ustarn'], 0, atol=1e-10), "ustarn should be 0 for wind from 270°"
        
        # The key test: ustars0 should equal ustars, not ustar
        # and ustarn0 should equal ustarn (0), not ustar
        assert np.allclose(s['ustars0'], s['ustars'], rtol=1e-10), \
            "ustars0 should equal ustars, not ustar magnitude"
        assert np.allclose(s['ustarn0'], s['ustarn'], rtol=1e-10), \
            "ustarn0 should equal ustarn (0), not ustar magnitude"
        
        # Verify that ustar (magnitude) is not zero
        assert np.all(s['ustar'] > 0), "ustar magnitude should be positive"
        
        # Verify that ustarn0 is NOT equal to ustar (the bug we're fixing)
        assert not np.allclose(s['ustarn0'], s['ustar'], rtol=1e-10), \
            "ustarn0 should NOT equal ustar magnitude (this was the bug)"

    def test_ustars0_ustarn0_components_pure_n_direction(self):
        """Test ustars0 and ustarn0 for wind purely in n-direction.
        
        For wind from 0/360 degrees (nautical, wind from north) with alfa=0:
        - Wind should be in pure n-direction (southward)
        - ustars should be 0
        - ustarn should be positive
        """
        from aeolis.wind import interpolate
        from aeolis.model import ModelState
        
        # Create minimal model state
        s = ModelState()
        
        # Grid setup
        nx, ny = 3, 3
        s['x'] = np.zeros((ny, nx))
        s['y'] = np.zeros((ny, nx))
        s['zb'] = np.zeros((ny, nx))
        s['uw'] = np.zeros((ny, nx))
        s['udir'] = np.zeros((ny, nx))
        s['uws'] = np.zeros((ny, nx))
        s['uwn'] = np.zeros((ny, nx))
        s['ustars'] = np.zeros((ny, nx))
        s['ustarn'] = np.zeros((ny, nx))
        s['ustar'] = np.zeros((ny, nx))
        s['ustars0'] = np.zeros((ny, nx))
        s['ustarn0'] = np.zeros((ny, nx))
        s['ustar0'] = np.zeros((ny, nx))
        s['tau'] = np.zeros((ny, nx))
        s['taus'] = np.zeros((ny, nx))
        s['taun'] = np.zeros((ny, nx))
        s['tau0'] = np.zeros((ny, nx))
        s['taus0'] = np.zeros((ny, nx))
        s['taun0'] = np.zeros((ny, nx))
        s['hveg'] = np.zeros((ny, nx))
        
        # Parameters
        p = {}
        p['alfa'] = 0.0
        p['kappa'] = 0.4
        p['z'] = 10.0
        p['k'] = 0.001
        p['method_roughness'] = 'constant'
        p['rhoa'] = 1.25
        p['process_wind'] = False
        p['wind_file'] = None
        
        # Set wind from 0/360 degrees (from north, blowing south)
        s['uw'][:, :] = 15.0
        s['udir'][:, :] = 0.0  # or 360.0
        
        # Run wind interpolation
        s = interpolate(s, p, t=0.0)
        
        # For wind from 0° with alfa=0:
        # uws = -uw * sin(0) = 0
        # uwn = -uw * cos(0) = -15
        
        # Check components
        assert np.allclose(s['ustars'], 0, atol=1e-10), "ustars should be 0 for wind from 0°"
        assert np.all(s['ustarn'] < 0), "ustarn should be negative for wind from 0°"
        
        # Check that the 0 components are preserved
        assert np.allclose(s['ustars0'], s['ustars'], rtol=1e-10), \
            "ustars0 should equal ustars (0)"
        assert np.allclose(s['ustarn0'], s['ustarn'], rtol=1e-10), \
            "ustarn0 should equal ustarn"
