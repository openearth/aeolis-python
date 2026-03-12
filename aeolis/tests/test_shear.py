"""
Tests for aeolis.shear – specifically the spatial sweep alternative to the
FFT-based compute_shear method.

Run with::

    pytest aeolis/tests/test_shear.py
"""

import numpy as np
import pytest

from aeolis.shear import WindShear, _cauchy_hilbert_sweep


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_hill_grid(ny=10, nx=40, dx=2.0, dy=2.0, hill_height=3.0, hill_sigma=10.0):
    """Return a small 2-D grid with a centred Gaussian hill."""
    x1d = np.arange(nx) * dx - nx * dx / 2.0
    y1d = np.arange(ny) * dy - ny * dy / 2.0
    xx, yy = np.meshgrid(x1d, y1d)
    z = hill_height * np.exp(-(xx**2 + yy**2) / (2.0 * hill_sigma**2))
    return xx, yy, z


def _run_windshear(method, xx, yy, z, dx=2.0, dy=2.0,
                   L=100.0, l=10.0, z0=0.001, buffer_width=20.0,
                   u0=10.0, udir=270.0):
    """Instantiate WindShear with the given method and run it."""
    ws = WindShear(xx, yy, z, dx=dx, dy=dy, L=L, l=l, z0=z0,
                   buffer_width=buffer_width, method=method)
    taux0 = np.ones_like(z) * 0.1
    tauy0 = np.zeros_like(z)
    ws(x=xx, y=yy, z=z,
       taux=taux0, tauy=tauy0,
       u0=u0, udir=udir,
       process_separation=False, c=0.5, mu_b=30.0,
       taus0=0.1, taun0=0.0,
       sep_filter_iterations=0, zsep_y_filter=False)
    return ws.get_shear()


# ---------------------------------------------------------------------------
# Tests for _cauchy_hilbert_sweep
# ---------------------------------------------------------------------------

class TestCauchyHilbertSweep:
    """Unit tests for the module-level numba-compiled helper."""

    def test_flat_bed_gives_zero(self):
        """A flat bed has no slope, so the shear perturbation must be zero."""
        ny, nx = 5, 30
        field = np.zeros((ny, nx))
        result = _cauchy_hilbert_sweep(field, 1.0, ny, nx)
        np.testing.assert_allclose(result, 0.0, atol=1e-12)

    def test_output_shape(self):
        """Output shape must match input shape."""
        ny, nx = 7, 25
        field = np.random.rand(ny, nx)
        result = _cauchy_hilbert_sweep(field, 1.0, ny, nx)
        assert result.shape == (ny, nx)

    def test_antisymmetric_response_around_symmetric_hill(self):
        """For a symmetric Gaussian hill, the perturbation must be
        antisymmetric: positive on the windward side and negative on the
        leeward side of the crest."""
        nx = 60
        dx = 1.0
        x = np.arange(nx) * dx - nx * dx / 2.0
        z = 5.0 * np.exp(-x**2 / (2.0 * 8.0**2))
        dzdx = np.gradient(z, dx)
        field = dzdx.reshape(1, nx)

        result = _cauchy_hilbert_sweep(field, 1.0, 1, nx)
        row = result[0]

        # The hill is centred at x=0, which falls at index nx//2.
        # Indices 0..(nx//2 - 1) correspond to x < 0 (windward / upwind side).
        # Indices nx//2..nx-1 correspond to x >= 0 (leeward / downwind side).
        # The windward half should contain the peak perturbation (increased shear)
        # and the leeward half should contain the trough (reduced shear).
        windward = row[:nx // 2]   # x < 0, upwind of crest
        leeward  = row[nx // 2:]   # x >= 0, at/downwind of crest

        assert windward.max() > 0.0, "Expected positive perturbation on windward side"
        assert leeward.min() < 0.0, "Expected negative perturbation on leeward side"

    def test_amplitude_scales_with_alpha(self):
        """Output must scale linearly with alpha."""
        ny, nx = 4, 20
        field = np.random.rand(ny, nx) * 0.1
        r1 = _cauchy_hilbert_sweep(field, 1.0, ny, nx)
        r2 = _cauchy_hilbert_sweep(field, 3.0, ny, nx)
        np.testing.assert_allclose(r2, 3.0 * r1, rtol=1e-12)

    def test_row_independence(self):
        """Each y-row must be processed independently."""
        ny, nx = 4, 20
        # Make a field where every row is identical
        row = np.random.rand(nx) * 0.1
        field = np.tile(row, (ny, 1))
        result = _cauchy_hilbert_sweep(field, 1.0, ny, nx)
        for j in range(1, ny):
            np.testing.assert_allclose(result[j], result[0], rtol=1e-12)


# ---------------------------------------------------------------------------
# Tests for WindShear with method='spatial'
# ---------------------------------------------------------------------------

class TestWindShearSpatial:
    """Integration tests for WindShear.compute_shear_spatial."""

    @pytest.fixture
    def hill_grid(self):
        return _make_hill_grid()

    def test_init_accepts_spatial_method(self, hill_grid):
        """WindShear must accept method='spatial' without raising."""
        xx, yy, z = hill_grid
        ws = WindShear(xx, yy, z, dx=2.0, dy=2.0, L=100.0, l=10.0,
                       z0=0.001, buffer_width=20.0, method='spatial')
        assert ws.method == 'spatial'

    def test_init_rejects_invalid_method(self, hill_grid):
        """WindShear must raise ValueError for unknown method strings."""
        xx, yy, z = hill_grid
        with pytest.raises(ValueError, match="method must be"):
            WindShear(xx, yy, z, dx=2.0, dy=2.0, L=100.0, l=10.0,
                      z0=0.001, buffer_width=20.0, method='unknown')

    def test_zero_wind_gives_zero_perturbation(self, hill_grid):
        """With u0=0 the shear perturbation must be exactly zero."""
        xx, yy, z = hill_grid
        ws = WindShear(xx, yy, z, dx=2.0, dy=2.0, L=100.0, l=10.0,
                       z0=0.001, buffer_width=20.0, method='spatial')
        taux0 = np.ones_like(z) * 0.1
        tauy0 = np.zeros_like(z)
        ws(x=xx, y=yy, z=z,
           taux=taux0, tauy=tauy0,
           u0=0.0, udir=270.0,
           process_separation=False, c=0.5, mu_b=30.0,
           taus0=0.1, taun0=0.0,
           sep_filter_iterations=0, zsep_y_filter=False)
        taux, tauy = ws.get_shear()
        # With no wind the stresses are returned unchanged (no perturbation)
        np.testing.assert_allclose(taux, taux0, atol=1e-10)

    def test_spatial_produces_nonzero_perturbation(self, hill_grid):
        """The spatial method must produce a non-trivial shear perturbation
        for a Gaussian hill and non-zero wind."""
        xx, yy, z = hill_grid
        taux, tauy = _run_windshear('spatial', xx, yy, z)
        # There should be some cells above and below the flat-bed reference
        assert taux.max() > 0.1, "Expected positive taux somewhere"
        assert taux.min() >= 0.0, "taux must be clipped to non-negative"

    def test_spatial_and_fft_produce_correlated_results(self, hill_grid):
        """The spatial and FFT methods should produce qualitatively similar
        taux fields (correlation > 0.7)."""
        xx, yy, z = hill_grid
        taux_fft, _ = _run_windshear('fft', xx, yy, z)
        taux_sp, _  = _run_windshear('spatial', xx, yy, z)
        corr = np.corrcoef(taux_fft.flatten(), taux_sp.flatten())[0, 1]
        assert corr > 0.7, (
            f"Expected correlation > 0.7 between FFT and spatial methods, "
            f"got {corr:.4f}"
        )

    def test_spatial_output_shape_matches_input(self, hill_grid):
        """The output shear fields must have the same shape as the input grid."""
        xx, yy, z = hill_grid
        taux, tauy = _run_windshear('spatial', xx, yy, z)
        assert taux.shape == z.shape
        assert tauy.shape == z.shape

    def test_fft_method_still_works(self, hill_grid):
        """The original FFT method must continue to work after the refactor."""
        xx, yy, z = hill_grid
        taux, tauy = _run_windshear('fft', xx, yy, z)
        assert taux.max() > 0.1
        assert taux.shape == z.shape
