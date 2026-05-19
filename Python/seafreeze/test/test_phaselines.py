"""Unit tests for seafreeze.phaselines.

Run from this directory:
    python -m unittest test_phaselines.py
"""
import unittest as ut
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')  # headless

from seafreeze.phaselines import (
    PhaseLineResult, PhaseRange,
    phase_range, phase_lines, wpd,
)


class TestPhaseRange(ut.TestCase):
    def test_Ih(self):
        r = phase_range('Ih')
        self.assertIsInstance(r, PhaseRange)
        self.assertAlmostEqual(r.P[0], 0.0,   places=6)
        self.assertAlmostEqual(r.P[1], 400.0, places=6)
        self.assertAlmostEqual(r.T[0], 1.0,   places=6)
        self.assertAlmostEqual(r.T[1], 301.0, places=6)
        self.assertIsNone(r.m)

    def test_water1(self):
        r = phase_range('water1')
        self.assertIsNone(r.m)
        self.assertGreaterEqual(r.T[1], 500.0)

    def test_NaClaq_has_m(self):
        r = phase_range('NaClaq')
        self.assertIsNotNone(r.m)
        self.assertAlmostEqual(r.m[0], 0.0, places=6)
        self.assertGreater(r.m[1], 6.5)  # spline goes to ~7 mol/kg

    def test_unknown_raises(self):
        with self.assertRaises(ValueError):
            phase_range('NotARealPhase')


class TestPhaseLinesPureWater(ut.TestCase):
    """Smoke-and-shape tests for the 12 pure-water pairs."""

    def test_Ih_water1_basic(self):
        r = phase_lines('Ih', 'water1')
        self.assertIsInstance(r, PhaseLineResult)
        self.assertEqual((r.matA, r.matB), ('Ih', 'water1'))
        self.assertGreater(r.P.size, 100)
        self.assertEqual(r.P.shape, r.T.shape)
        self.assertEqual(r.P.shape, r.stable.shape)
        # Atmospheric melting near 273.15 K
        idx = np.argmin(r.P)
        self.assertAlmostEqual(r.T[idx], 273.15, delta=0.1)
        # Negative-slope curve: highest P should give lowest T
        self.assertLess(r.T[np.argmax(r.P)], r.T[np.argmin(r.P)])

    def test_VI_water1_high_pressure(self):
        r = phase_lines('VI', 'water1')
        # Curve should extend to high pressure
        self.assertGreater(r.P.max(), 2000.0)
        # Monotonic positive slope: VI-liquid melts hotter at higher P
        order = np.argsort(r.P)
        self.assertTrue(np.all(np.diff(r.T[order]) > -0.5))  # mostly monotonic

    def test_freezing_point_low_P(self):
        """Ih-water1 at near-zero pressure should be ~273.15 K."""
        r = phase_lines('Ih', 'water1')
        idx = np.argmin(r.P)
        self.assertLess(r.P[idx], 1.0)
        self.assertAlmostEqual(r.T[idx], 273.15, delta=0.05)

    def test_segment_stable(self):
        r = phase_lines('Ih', 'water1', segment='stable')
        self.assertTrue(r.stable.all())
        self.assertEqual(r.segment, 'stable')

    def test_segment_meta(self):
        r = phase_lines('Ih', 'water1', segment='meta')
        self.assertFalse(r.stable.any())
        self.assertEqual(r.segment, 'meta')

    def test_segment_all_partition(self):
        ra = phase_lines('Ih', 'water1', segment='all')
        rs = phase_lines('Ih', 'water1', segment='stable')
        rm = phase_lines('Ih', 'water1', segment='meta')
        self.assertEqual(ra.P.size, rs.P.size + rm.P.size)

    def test_II_water1_all_metastable(self):
        """II <-> water1 has no stable equilibrium curve in the standard
        phase diagram (II is surrounded by Ih/III/V/VI)."""
        r = phase_lines('II', 'water1', segment='all')
        self.assertEqual(r.stable.sum(), 0)
        # 'stable' segment should be empty
        rs = phase_lines('II', 'water1', segment='stable')
        self.assertEqual(rs.P.size, 0)

    def test_symmetric_lookup(self):
        r1 = phase_lines('Ih', 'water1')
        r2 = phase_lines('water1', 'Ih')
        # Should normalize to canonical (matA, matB) order
        self.assertEqual((r1.matA, r1.matB), (r2.matA, r2.matB))
        self.assertEqual(r1.P.size, r2.P.size)
        np.testing.assert_allclose(r1.P, r2.P)
        np.testing.assert_allclose(r1.T, r2.T)

    def test_invalid_segment(self):
        with self.assertRaises(ValueError):
            phase_lines('Ih', 'water1', segment='invalid')

    def test_unknown_pair(self):
        with self.assertRaises(ValueError):
            phase_lines('Ih', 'water_IAPWS95')

    def test_unknown_material(self):
        with self.assertRaises(ValueError):
            phase_lines('Ih', 'NotAPhase')

    def test_user_grid_override(self):
        # Custom T grid should be respected (within tolerance)
        T_user = np.arange(240.0, 273.0, 1.0)
        r = phase_lines('Ih', 'water1', T=T_user)
        # Resulting curve should stay within user T range
        self.assertGreaterEqual(r.T.min(), 239.0)
        self.assertLessEqual(r.T.max(), 274.0)


class TestPhaseLinesNaCl(ut.TestCase):
    def test_requires_m(self):
        with self.assertRaises(ValueError):
            phase_lines('Ih', 'NaClaq')

    def test_pure_pair_rejects_m(self):
        with self.assertRaises(ValueError):
            phase_lines('Ih', 'water1', m=0.5)

    def test_scalar_m_returns_single(self):
        r = phase_lines('Ih', 'NaClaq', m=0.5)
        self.assertIsInstance(r, PhaseLineResult)
        self.assertAlmostEqual(r.m, 0.5)
        self.assertGreater(r.P.size, 50)

    def test_vector_m_returns_list(self):
        rs = phase_lines('Ih', 'NaClaq', m=[0.5, 2.0])
        self.assertIsInstance(rs, list)
        self.assertEqual(len(rs), 2)
        self.assertAlmostEqual(rs[0].m, 0.5)
        self.assertAlmostEqual(rs[1].m, 2.0)

    def test_freezing_point_depression(self):
        """Higher molality lowers the Ih melting temperature at fixed P."""
        r0p5 = phase_lines('Ih', 'NaClaq', m=0.5)
        r2p0 = phase_lines('Ih', 'NaClaq', m=2.0)
        # Compare at P near 100 MPa
        def T_at(curve, P_target):
            order = np.argsort(curve.P)
            return np.interp(P_target, curve.P[order], curve.T[order])
        T_low_m  = T_at(r0p5, 100.0)
        T_high_m = T_at(r2p0, 100.0)
        self.assertGreater(T_low_m, T_high_m)
        # Should be a few K of depression
        self.assertGreater(T_low_m - T_high_m, 1.0)

    def test_m_out_of_range(self):
        with self.assertRaises(ValueError):
            phase_lines('Ih', 'NaClaq', m=100.0)  # > spline max ~7

    def test_NaClaq_marked_stable(self):
        r = phase_lines('Ih', 'NaClaq', m=1.0)
        # Per Matlab parity: ice-NaClaq curves are flagged 'all stable'
        self.assertTrue(r.stable.all())


class TestWPD(ut.TestCase):
    def test_runs(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fig = wpd()
        self.assertIsNotNone(fig)
        # Should have plotted at least 12 line artists (one per pair, more
        # when stable/meta segments are split)
        ax = fig.axes[0]
        self.assertGreaterEqual(len(ax.get_lines()), 6)

    def test_runs_with_NaCl_overlay(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fig = wpd(solute='NaCl', m=[0.5, 2.0])
        ax = fig.axes[0]
        self.assertGreater(len(ax.get_lines()), 6)

    def test_NaCl_requires_m(self):
        with self.assertRaises(ValueError):
            wpd(solute='NaCl')


if __name__ == '__main__':
    ut.main()
