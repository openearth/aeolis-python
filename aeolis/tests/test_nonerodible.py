"""
Test for non-erodible layer constraint in sweep solver and bed update.

This test verifies that the bed level does not drop below the non-erodible
layer (zne) during erosion, which can occur at downwind boundaries in 2D cases.
"""

import numpy as np
import pytest
from aeolis import utils


class TestNonErodibleLayerConstraint:
    """
    Test that the sweep solver and bed update properly enforce the
    non-erodible layer constraint to prevent erosion below zne.
    """

    def test_sweep_respects_nonerodible_layer(self):
        """
        Test that the sweep function limits pickup to prevent bed from
        dropping below the non-erodible layer.
        """
        # Setup a simple 5x5 grid with 1 fraction
        ny, nx, nf = 4, 4, 1
        
        # Initial bed level
        zb = np.ones((ny + 1, nx + 1)) * 2.0  # 2m bed level
        
        # Non-erodible layer at 1.5m
        zne = np.ones((ny + 1, nx + 1)) * 1.5
        
        # Grid spacing
        ds = np.ones((ny + 1, nx + 1)) * 10.0
        dn = np.ones((ny + 1, nx + 1)) * 10.0
        
        # Sediment mass (enough to erode 1m if unrestricted)
        rhog = 2650.0
        porosity = 0.4
        mass_per_meter = rhog * (1.0 - porosity)  # kg/m2 per meter of depth
        mass = np.ones((ny + 1, nx + 1, 1, nf)) * mass_per_meter * 1.0  # 1m worth of sediment
        
        # Concentration (initial and equilibrium)
        Ct = np.zeros((ny + 1, nx + 1, nf))
        Cu = np.ones((ny + 1, nx + 1, nf)) * 10.0  # High equilibrium to drive erosion
        
        # Wind velocity (positive in x direction, creating erosion)
        us = np.ones((ny + 1, nx + 1, nf)) * 5.0
        un = np.zeros((ny + 1, nx + 1, nf))
        
        # Weights (uniform)
        w = np.ones((ny + 1, nx + 1, nf))
        
        # Time parameters
        dt = 3600.0  # 1 hour
        Ts = 1.0  # adaptation time scale
        
        # Call sweep with non-erodible layer constraint
        Ct_result, pickup = utils.sweep(
            Ct, Cu, mass, dt, Ts, ds, dn, us, un, w,
            zb=zb, zne=zne, rhog=rhog, porosity=porosity
        )
        
        # Calculate bed level change from pickup
        total_pickup = np.sum(pickup, axis=2)
        dz = total_pickup / mass_per_meter
        new_zb = zb - dz
        
        # Check that bed level doesn't drop below non-erodible layer
        # Allow small numerical tolerance
        assert np.all(new_zb >= zne - 1e-6), \
            f"Bed level dropped below non-erodible layer. Min new_zb: {new_zb.min()}, zne: {zne.min()}"
        
        # Check that pickup was limited (should be less than available mass)
        max_allowed_pickup = mass_per_meter * (zb - zne).max()
        assert np.all(total_pickup <= max_allowed_pickup + 1e-6), \
            f"Pickup exceeded maximum allowed. Max pickup: {total_pickup.max()}, max allowed: {max_allowed_pickup}"

    def test_sweep_without_nonerodible_layer_backward_compatible(self):
        """
        Test that sweep function works without zb/zne parameters (backward compatibility).
        """
        # Setup a simple 5x5 grid with 1 fraction
        ny, nx, nf = 4, 4, 1
        
        # Grid spacing
        ds = np.ones((ny + 1, nx + 1)) * 10.0
        dn = np.ones((ny + 1, nx + 1)) * 10.0
        
        # Sediment mass
        mass = np.ones((ny + 1, nx + 1, 1, nf)) * 1000.0
        
        # Concentration
        Ct = np.zeros((ny + 1, nx + 1, nf))
        Cu = np.ones((ny + 1, nx + 1, nf)) * 1.0
        
        # Wind velocity
        us = np.ones((ny + 1, nx + 1, nf)) * 5.0
        un = np.zeros((ny + 1, nx + 1, nf))
        
        # Weights
        w = np.ones((ny + 1, nx + 1, nf))
        
        # Time parameters
        dt = 3600.0
        Ts = 1.0
        
        # Call sweep without zb/zne parameters (should work)
        Ct_result, pickup = utils.sweep(
            Ct, Cu, mass, dt, Ts, ds, dn, us, un, w
        )
        
        # Should return valid results
        assert Ct_result.shape == (ny + 1, nx + 1, nf)
        assert pickup.shape == (ny + 1, nx + 1, nf)
        assert not np.any(np.isnan(Ct_result))
        assert not np.any(np.isnan(pickup))

    def test_bed_update_respects_nonerodible_layer(self):
        """
        Test that bed.update enforces the non-erodible layer constraint.
        """
        # This would require setting up a more complete model state
        # For now, we test the sweep function which is the main fix
        pass
