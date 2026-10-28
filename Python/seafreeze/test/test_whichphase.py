import warnings, unittest as ut
import numpy as np
import seafreeze.seafreeze as sf


class TestWhichPhase(ut.TestCase):
    def setup(self):
        warnings.simplefilter('ignore', category=ImportWarning)
    def tearDown(self):
        pass
    def test_scatter_singlePt(self):
        PT = np.empty((1,), object)
        PT[0] = (2000, 334)
        self.assertEqual(6, sf.whichphase(PT, 'water1')[0])
    def test_scatter_multiPt(self):
        PT = np.empty((3,), object)
        PT[0] = (100, 200)
        PT[1] = (400, 250)
        PT[2] = (1000, 300)
        out = sf.whichphase(PT)
        self.assertTrue(np.array_equal([1, 5, 6], out))
    def test_scatter_all_extrapolations(self):
        PT = np.empty((3,), object)
        PT[0] = (1000001, 250)  # bad P, ok T
        PT[1] = (400, 5000)    # ok P, bad T
        PT[2] = (2000, 2000)   # bad P and T
        np.testing.assert_array_equal([np.nan, np.nan, np.nan], sf.whichphase(PT))
    def test_grid(self):
        P = np.arange(0, 1001, 200)
        T = np.arange(200, 351, 50)
        exp = np.array([[1, 1, 0, 0], [2, 1, 0, 0], [2, 5, 0, 0], [2, 5, 0, 0], [6, 6, 0, 0], [6, 6, 6, 0]], dtype=object)
        act = sf.whichphase(np.array([P, T], dtype=object), 'water1')
        np.testing.assert_array_equal(exp, act)
    def test_grid_all_extrapolations(self):
        P = np.arange(0, 1001, 200)
        T = np.arange(200, 351, 50)
        exp = np.array([[1, 1, 0, 0], [2, 1, 0, 0], [2, 5, 0, 0], [2, 5, 0, 0], [6, 6, 0, 0], [6, 6, 6, 0]], dtype=object)
        act = sf.whichphase(np.array([P, T], dtype=object), 'water1')
        np.testing.assert_array_equal(exp, act)
    def test_phases(self):
        PT = np.empty((7,), object)
        PT[0] = (500, 300)   # liquid water
        PT[1] = (1, 200)     # ice Ih
        PT[2] = (300, 200)   # ice II
        PT[3] = (300, 250)   # ice III
        PT[4] = (500, 250)   # ice V
        PT[5] = (800, 250)   # ice VI
        PT[6] = (2000, 334)  # high-P — ice VI stable at this condition
        act = sf.whichphase(PT, 'water1')
        np.testing.assert_array_equal([0, 1, 2, 3, 5, 6, 6], act)

    def test_NaCl_grid(self):
        P = np.arange(0, 1001, 200)
        T = np.arange(200, 351, 50)
        M = np.arange(1, 2, 1)
        out = sf.whichphase(np.array([P, T, M], dtype=object), 'NaClaq')
        self.assertEqual(1, out[0][0])

    def test_NaCl_singlePt(self):
        PTm = np.empty((1,), dtype=object)
        PTm[0] = (100, 250, 3)
        out = sf.whichphase(PTm, 'NaClaq')
        self.assertEqual(1, out)

    def test_NaCl_singlePt2(self):
        PTm = np.empty((1,), dtype=object)
        PTm[0] = (900, 285, 3)
        out = sf.whichphase(PTm, 'NaClaq')
        self.assertEqual(0, out)


class TestWhichPhaseSolutes(ut.TestCase):
    """Solute validation and agreement between liquid EOS (added 2026-10)."""

    def _pts(self):
        PT = np.empty((4,), object)
        PT[0] = (100.0, 280.0)   # liquid
        PT[1] = (600.0, 250.0)   # ice V
        PT[2] = (1.0, 200.0)     # ice Ih
        PT[3] = (1500.0, 300.0)  # ice VI
        return PT

    def test_all_pure_water_eos_agree(self):
        ref = sf.whichphase(self._pts(), 'water1')
        for sol in ('water2', 'water_IAPWS95'):
            np.testing.assert_array_equal(ref, sf.whichphase(self._pts(), sol), err_msg=sol)

    def test_pure_water_eos_agree_on_grid(self):
        P = np.arange(0, 1001, 200.0); T = np.arange(200, 351, 50.0)
        G = np.array([P, T], dtype=object)
        ref = sf.whichphase(G, 'water1')
        np.testing.assert_array_equal(ref, sf.whichphase(G, 'water2'))

    def test_every_nacl_variant_runs(self):
        PTm = np.empty((2,), object)
        PTm[0] = (100.0, 250.0, 3.0)   # ice Ih
        PTm[1] = (900.0, 285.0, 3.0)   # liquid
        for sol in ('NaClaq', 'NaClaq_LP', 'NaClaq_5GPa_2024'):
            np.testing.assert_array_equal([1, 0], sf.whichphase(PTm, sol), err_msg=sol)

    def test_ice_as_solute_rejected(self):
        with self.assertRaises(ValueError):
            sf.whichphase(self._pts(), 'Ih')

    def test_unknown_solute_rejected(self):
        with self.assertRaises(ValueError):
            sf.whichphase(self._pts(), 'brine')

    def test_agrees_with_phase_lines_at_Ih_melting(self):
        """Just below the Ih–water1 melting curve the answer is Ih, just above it is liquid."""
        from seafreeze import phase_lines
        line = phase_lines('Ih', 'water1', segment='stable')
        for P_target in (0.1, 100.0, 180.0):
            i = int(np.argmin(np.abs(line.P - P_target)))
            P_eq, T_eq = float(line.P[i]), float(line.T[i])
            PT = np.empty((2,), object)
            PT[0] = (P_eq, T_eq - 1.0)
            PT[1] = (P_eq, T_eq + 1.0)
            np.testing.assert_array_equal([1, 0], sf.whichphase(PT, 'water1'),
                                          err_msg=f'at P={P_eq:.1f} MPa, T_eq={T_eq:.2f} K')

    def test_salt_depresses_freezing(self):
        """A point that is ice Ih in pure water is liquid in 3 mol/kg brine."""
        PT = np.empty((1,), object); PT[0] = (0.1, 268.0)
        self.assertEqual(1, sf.whichphase(PT, 'water1')[0])
        PTm = np.empty((1,), object); PTm[0] = (0.1, 268.0, 3.0)
        self.assertEqual(0, sf.whichphase(PTm, 'NaClaq')[0])


class TestPhasenum2phase(ut.TestCase):
    def test_round_trip_ices(self):
        for name, desc in sf.phases.items():
            if desc.shear_mod_parms and name != 'VII_X_French':
                self.assertEqual(name, sf.phasenum2phase(desc.phase_num))

    def test_zero_is_liquid(self):
        self.assertEqual('water1', sf.phasenum2phase(0))
        self.assertEqual('NaClaq', sf.phasenum2phase(0, 'NaClaq'))

    def test_nan_defaults_to_liquid(self):
        self.assertEqual('water1', sf.phasenum2phase(np.nan))
