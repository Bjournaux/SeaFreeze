import warnings, unittest as ut
import numpy as np
import seafreeze as sf


class TestGetPhaseThermodynamics(ut.TestCase):
    def setup(self):
        warnings.simplefilter('ignore', category=ImportWarning)
    def tearDown(self):
        pass
    #########################################
    ## _is_scatter
    #########################################
    # def test_is_scatter_singlept_simple(self):
    #     PT = (1,2)
    #     self.assertTrue(sf._is_scatter(PT))
    # def test_is_scatter_singlept_list1(self):
    #     PT = np.array([(1,2)])
    #     self.assertTrue(sf._is_scatter(PT))
    # def test_is_scatter_single_grid(self):
    #     PT = np.array([7, 8])
    #     self.assertFalse(sf._is_scatter(PT))
    def test_is_scatter_singlept_preallocated(self):
        PT = np.empty((1,), np.object)
        PT[0] = (1,2)
        self.assertTrue(sf._is_scatter(PT))
    def test_is_scatter_scatter(self):
        PT = np.empty((3,), np.object)
        PT[0] = (1, 2)
        PT[1] = (3, 4)
        PT[2] = (5, 6)
        self.assertTrue(sf._is_scatter(PT))
    def test_is_scatter_grid(self):
        P = np.arange(0.1, 1000.2, 10)
        T = np.arange(240, 501, 2)
        PT = np.array([P, T])
        self.assertFalse(sf._is_scatter(PT))
    #########################################
    ## _get_T
    #########################################
    def test_get_T_singlept(self):
        PT = np.empty((1,), np.object)
        PT[0] = (1, 2)
        self.assertEqual(2, sf._get_T(PT, True))
    def test_get_T_multipt(self):
        PT = np.empty((3,), np.object)
        PT[0] = (1, 2)
        PT[1] = (3, 4)
        PT[2] = (5, 6)
        self.assertTrue(np.array_equal(np.array([2,4,6]), sf._get_T(PT, True)))
    def test_get_T_grid(self):
        P = np.arange(0.1, 1000.2, 10)
        T = np.arange(240, 501, 2)
        PT = np.array([P, T])
        self.assertTrue((np.array_equal(T, sf._get_T(PT, False))))
    #########################################
    ## _get_shear_mod_GPa
    #########################################
    def test_get_shear_mod_GPa_singlept(self):
        PT = np.empty((1,), np.object)
        PT[0] = (900, 255)
        rho = 1.356072490993616e+03
        sm = sf._get_shear_mod_GPa(sf.phases['VI'].shear_mod_parms, rho, sf._get_T(PT, True))
        self.assertAlmostEqual(7.303268592388283, sm[0], places=10)
    def test_get_shear_mod_GPa_multipt(self):
        PT = np.empty((3,), np.object)
        PT[0] = (900, 255)
        PT[1] = (910, 265)
        PT[2] = (920, 275)
        rho = np.array([1356.072490993616, 1354.106807053618, 1352.057926458423])
        sm = sf._get_shear_mod_GPa(sf.phases['VI'].shear_mod_parms, rho, sf._get_T(PT, True))
        expected = np.array([7.303268592388283, 7.128869123438319, 6.953013713022407])
        self.assertTrue(np.allclose(sm, expected))
    def test_get_shear_mod_GPa_grid(self):
        P = np.arange(900, 921, 10)
        T = np.arange(255, 276, 10)
        PT = np.array([P, T])
        rho = np.array([[1356.072490993616, 1353.307249697806, 1350.440903858578],
                        [1356.862314715232, 1354.106807053618, 1351.250712225762],
                        [1357.649726085489, 1354.903864412315, 1352.057926458423]])
        sm = sf._get_shear_mod_GPa(sf.phases['VI'].shear_mod_parms, rho, sf._get_T(PT, False))
        expected = np.array([[7.303268592388283, 7.114876869711612, 6.924715817525114],
                             [7.317090507516553, 7.128869123438319, 6.938887463950844],
                             [7.330870206496061, 7.142817627215505, 6.953013713022400]])
        self.assertTrue(np.allclose(sm, expected))
    #########################################
    ## __get_Vp
    #########################################
    def test_get_Vp_singlept(self):
        smg = 7.303268592388283
        rho = 1356.072490993616
        Ks = 18323.49756691741
        Vp = sf._get_Vp(smg, rho, Ks)
        self.assertAlmostEqual(4548.954381485812, Vp)
    def test_get_Vp_multipt(self):
        smg = np.array([7.303268592388283, 7.128869123438319, 6.953013713022407])
        rho = np.array([1356.072490993616, 1354.106807053618, 1352.057926458423])
        Ks = np.array([18323.49756691741, 18221.54159041028, 18115.56318084243])
        Vp = sf._get_Vp(smg, rho, Ks)
        expected = np.array([4548.954381485812, 4525.042206621253, 4500.581390585054])
        self.assertTrue(np.allclose(Vp, expected))
    def test_get_Vp_grid(self):
        smg = np.array([[7.303268592388283, 7.114876869711612, 6.924715817525114],
                        [7.317090507516553, 7.128869123438319, 6.938887463950844],
                        [7.330870206496061, 7.142817627215505, 6.953013713022400]])
        rho = np.array([[1356.072490993616, 1353.307249697806, 1350.440903858578],
                        [1356.862314715232, 1354.106807053618, 1351.250712225762],
                        [1357.649726085489, 1354.903864412315, 1352.057926458423]])
        Ks = np.array([[18323.49756691741, 18159.55544116658, 17990.81337195430],
                       [18385.04086077863, 18221.54159041028, 18053.25505449579],
                       [18446.46522308364, 18283.40144175383, 18115.56318084243]])
        Vp = sf._get_Vp(smg, rho, Ks)
        expected = np.array([[4548.954381485813, 4519.791517059500, 4489.896438587049],
                             [4554.105836046733, 4525.042206621253, 4495.251116541151],
                             [4559.235383167323, 4530.269763950370, 4500.581390585055]])
        self.assertTrue(np.allclose(Vp, expected))
    #########################################
    ## __get_Vs
    #########################################
    def test_get_Vs_singlept(self):
        smg = 7.303268592388283
        rho = 1356.072490993616
        Vs = sf._get_Vs(smg, rho)
        self.assertAlmostEqual(2320.690281146717, Vs)
    def test_get_Vs_multipt(self):
        smg = np.array([7.303268592388283, 7.128869123438319, 6.953013713022407])
        rho = np.array([1356.072490993616, 1354.106807053618, 1352.057926458423])
        Vs = sf._get_Vs(smg, rho)
        expected = np.array([2320.690281146717, 2294.477800749785, 2267.717197934159])
        self.assertTrue(np.allclose(Vs, expected))
    def test_get_Vs_grid(self):
        smg = np.array([[7.303268592388283, 7.114876869711612, 6.924715817525114],
                        [7.317090507516553, 7.128869123438319, 6.938887463950844],
                        [7.330870206496061, 7.142817627215505, 6.953013713022400]])
        rho = np.array([[1356.072490993616, 1353.307249697806, 1350.440903858578],
                        [1356.862314715232, 1354.106807053618, 1351.250712225762],
                        [1357.649726085489, 1354.903864412315, 1352.057926458423]])
        Vs = sf._get_Vs(smg, rho)
        expected = np.array([[2320.690281146717, 2292.901984107111, 2264.452345731802],
                             [2322.209103249143, 2294.477800749785, 2266.088955400825],
                             [2323.720540877948, 2296.045764570984, 2267.717197934159]])
        self.assertTrue(np.allclose(Vs, expected))
    # #########################################
    # ## get_phase_thermodynamics
    # #########################################
    # def test_get_phase_thermodynamics_singlept(self):
    #     PT = np.empty((1,), np.object)
    #     PT[0] = (900, 255)
    #     out = sf.get_phase_thermodynamics('VI', PT, '../../SeaFreeze_Gibbs.mat')
    #     self.assertAlmostEqual(1.356072490993616e+03, out.rho[0], places=7)
    #     self.assertAlmostEqual(7.303268592388283e+03, out.shear[0], places=7)
    #     self.assertAlmostEqual(4.548954381485812e+03, out.Vp[0], places=7)
    #     self.assertAlmostEqual(2.320690281146717e+03, out.Vs[0], places=7)
