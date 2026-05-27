"""
Tests for seafreeze.rho2P — pressure-from-density inversion.

Strategy: round-trip consistency.  Compute rho via getProp, invert to P with
rho2P, verify the residual is within tolerance.

Baptiste Journaux - 2026
"""

import warnings
import unittest
import numpy as np
import seafreeze as sf

TOL_P = 0.05    # MPa — acceptable round-trip error


def _scatter(P_list, T_list):
    """Build a scatter PTm array from lists."""
    n = len(P_list)
    PTm = np.empty(n, dtype=object)
    for i in range(n):
        PTm[i] = (float(P_list[i]), float(T_list[i]))
    return PTm


def _scatter_nacl(P_list, T_list, m_list):
    n = len(P_list)
    PTm = np.empty(n, dtype=object)
    for i in range(n):
        PTm[i] = (float(P_list[i]), float(T_list[i]), float(m_list[i]))
    return PTm


class TestRho2P_RoundTrip(unittest.TestCase):
    """Round-trip consistency tests."""

    # ------------------------------------------------------------------
    # Pure phases
    # ------------------------------------------------------------------

    def test_water1_scatter(self):
        P_ref = [0.1, 100., 500., 1000., 2000.]
        T_ref = [280., 300., 320.,  340.,  355.]
        out = sf.getProp(_scatter(P_ref, T_ref), 'water1', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_ref, 'water1')
        self.assertEqual(P_rec.shape, (5,))
        np.testing.assert_allclose(P_rec, P_ref, atol=TOL_P,
                                   err_msg='water1 round-trip failed')

    def test_iceVI_scatter(self):
        P_ref = [700., 900., 1100., 1500., 2000.]
        T_ref = [250., 255.,  260.,  270.,  280.]
        out = sf.getProp(_scatter(P_ref, T_ref), 'VI', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_ref, 'VI')
        np.testing.assert_allclose(P_rec, P_ref, atol=TOL_P,
                                   err_msg='ice VI round-trip failed')

    def test_iceIh_scatter(self):
        P_ref = [0.1, 50.,  100., 200., 350.]
        T_ref = [260., 265., 268., 270., 271.]
        out = sf.getProp(_scatter(P_ref, T_ref), 'Ih', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_ref, 'Ih')
        np.testing.assert_allclose(P_rec, P_ref, atol=TOL_P,
                                   err_msg='ice Ih round-trip failed')

    def test_water2_high_T(self):
        """water2 is monotone at high T (T=1000 K), wider tolerance allowed."""
        P_ref = [200., 400., 600., 800., 1000.]
        T_ref = [1000.] * 5
        out = sf.getProp(_scatter(P_ref, T_ref), 'water2', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_ref, 'water2')
        np.testing.assert_allclose(P_rec, P_ref, atol=1.0,
                                   err_msg='water2 round-trip failed')

    # ------------------------------------------------------------------
    # NaClaq
    # ------------------------------------------------------------------

    def test_NaClaq_scatter_m1(self):
        P_ref = [0.1, 200., 500., 1500., 5000.]
        T_ref = [280., 300., 320.,  350.,  400.]
        m_ref = [1.0] * 5
        out = sf.getProp(_scatter_nacl(P_ref, T_ref, m_ref),
                         'NaClaq', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_ref, 'NaClaq', m=1.0)
        np.testing.assert_allclose(P_rec, P_ref, atol=TOL_P,
                                   err_msg='NaClaq m=1 round-trip failed')

    def test_NaClaq_5GPa_2024_scatter_m2(self):
        P_ref = [100., 500., 1000., 2000., 3000.]
        T_ref = [280., 300.,  320.,  350.,  400.]
        m_ref = [2.0] * 5
        out = sf.getProp(_scatter_nacl(P_ref, T_ref, m_ref),
                         'NaClaq_5GPa_2024', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_ref, 'NaClaq_5GPa_2024', m=2.0)
        np.testing.assert_allclose(P_rec, P_ref, atol=TOL_P,
                                   err_msg='NaClaq_5GPa_2024 m=2 round-trip failed')


class TestRho2P_OutputShape(unittest.TestCase):
    """Output shape and broadcast tests."""

    def test_scalar_in_scalar_out(self):
        out = sf.getProp(_scatter([200.], [300.]), 'water1', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(float(out.rho[0]), 300.0, 'water1')
        self.assertEqual(P_rec.shape, ())  # 0-d (scalar input → scalar output)

    def test_1d_array_preserved(self):
        P_ref = [100., 300., 600., 900., 1200.]
        T_ref = [280., 290., 300., 310., 320.]
        out = sf.getProp(_scatter(P_ref, T_ref), 'water1', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_ref, 'water1')
        self.assertEqual(P_rec.shape, (5,))

    def test_scalar_T_broadcast(self):
        P_ref = [100., 300., 600., 900., 1200.]
        T_fix = 300.
        out = sf.getProp(_scatter(P_ref, [T_fix] * 5), 'water1', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_fix, 'water1')  # scalar T
        self.assertEqual(P_rec.shape, (5,))
        np.testing.assert_allclose(P_rec, P_ref, atol=TOL_P)

    def test_2d_input_shape_preserved(self):
        P_ref = np.array([[100., 500.], [1000., 1500.]])  # 2×2
        T_ref = np.array([[280., 300.], [320.,  340.]])
        out = sf.getProp(_scatter(P_ref.ravel(), T_ref.ravel()),
                         'water1', sf.seafreeze.defpath, 'rho')
        rho_2d = out.rho.reshape(2, 2)
        P_rec = sf.rho2P(rho_2d, T_ref, 'water1')
        self.assertEqual(P_rec.shape, (2, 2))
        np.testing.assert_allclose(P_rec, P_ref, atol=TOL_P)


class TestRho2P_NaN(unittest.TestCase):
    """NaN propagation for out-of-range inputs."""

    def test_rho_below_min_is_nan(self):
        # Hot water at low P has density ~958 kg/m³; target far below → NaN
        out = sf.getProp(_scatter([0.1], [373.]), 'water1', sf.seafreeze.defpath, 'rho')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            P_rec = sf.rho2P(out.rho[0] - 200., 373., 'water1')
        self.assertTrue(np.isnan(P_rec).all())

    def test_rho_above_max_is_nan(self):
        # 1600 kg/m³ is well above any water1 density (max ~1400 at high P) → NaN
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            P_rec = sf.rho2P(1600., 300., 'water1')
        self.assertTrue(np.isnan(P_rec).all())

    def test_nan_rho_target_is_nan(self):
        P_rec = sf.rho2P(np.nan, 300., 'water1')
        self.assertTrue(np.isnan(P_rec).all())


class TestRho2P_Physics(unittest.TestCase):
    """Physical sanity checks."""

    def test_higher_rho_higher_P_isothermal(self):
        """At fixed T, density increases with P → inversion must be monotone."""
        T_fix = 300.
        P_grid = [100., 500., 1000., 1500.]
        out = sf.getProp(_scatter(P_grid, [T_fix]*4), 'water1', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_fix, 'water1')
        self.assertTrue(np.all(np.diff(P_rec) > 0))

    def test_ambient_water_near_01_MPa(self):
        """Density from getProp at (0.1 MPa, 298 K) should round-trip to ≈ 0.1 MPa."""
        import numpy as _np
        PT = _np.empty(1, dtype=object); PT[0] = (0.1, 298.0)
        rho_ref = sf.getProp(PT, 'water1', sf.seafreeze.defpath, 'rho').rho[0]
        P_rec = sf.rho2P(rho_ref, 298.0, 'water1')
        self.assertAlmostEqual(float(_np.asarray(P_rec).flat[0]), 0.1, delta=0.5)


class TestRho2P_Options(unittest.TestCase):
    """Keyword options: P0, tol, max_iter."""

    def test_custom_P0(self):
        P_ref = [800., 900., 1000.]
        T_ref = [255., 258.,  262.]
        out = sf.getProp(_scatter(P_ref, T_ref), 'VI', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho, T_ref, 'VI', P0=1000.)
        np.testing.assert_allclose(P_rec, P_ref, atol=TOL_P)

    def test_tight_tol(self):
        out = sf.getProp(_scatter([500.], [300.]), 'water1', sf.seafreeze.defpath, 'rho')
        P_rec = sf.rho2P(out.rho[0], 300., 'water1', tol=1e-4)
        self.assertAlmostEqual(float(np.asarray(P_rec).flat[0]), 500., delta=0.001)


class TestRho2P_Errors(unittest.TestCase):
    """Input validation errors."""

    def test_unknown_phase_raises(self):
        with self.assertRaises(ValueError):
            sf.rho2P(1000., 300., 'IceX')

    def test_NaClaq_without_m_raises(self):
        with self.assertRaises(ValueError):
            sf.rho2P(1000., 300., 'NaClaq')

    def test_T_size_mismatch_raises(self):
        with self.assertRaises(ValueError):
            sf.rho2P([1000., 1010.], [280., 290., 300.], 'water1')


if __name__ == '__main__':
    unittest.main()
