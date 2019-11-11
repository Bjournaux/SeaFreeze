import warnings, unittest as ut
import numpy as np
import seafreeze.seafreeze as sf


class TestWhichPhase(ut.TestCase):
    def setup(self):
        warnings.simplefilter('ignore', category=ImportWarning)
    def tearDown(self):
        pass
    def test_scatter_singlePt(self):
        PT = np.empty((1,), np.object)
        PT[0] = (2000, 334) #(950, 255)
        self.assertEqual(6, sf.whichphase(PT, '../../../Matlab/SeaFreeze_Gibbs.mat')[0])
    def test_scatter_multiPt(self):
        PT = np.empty((3,), np.object)
        PT[0] = (100, 200)
        PT[1] = (400, 250)
        PT[2] = (1000, 300)
        self.assertTrue(np.array_equal([1, 5, 6], sf.whichphase(PT, '../../../Matlab/SeaFreeze_Gibbs.mat')))
    def test_scatter_all_extrapolations(self):
        PT = np.empty((3,), np.object)
        PT[0] = (5000, 250)  # bad P, ok T
        PT[1] = (400, 5000)  # ok P, bad T
        PT[2] = (2000, 2000) # bad P and T
        np.testing.assert_array_equal([np.nan, np.nan, np.nan], sf.whichphase(PT, '../../../Matlab/SeaFreeze_Gibbs.mat'))
    def test_grid(self):
        P = np.arange(0, 1001, 200)
        T = np.arange(200, 351, 50)
        exp = np.array([[1, 1, 0, 0], [2, 1, 0, 0], [2, 5, 0, 0], [2, 5, 0, 0], [6, 6, 0, 0], [6, 6, 6, 0]])
        act = sf.whichphase(np.array([P, T]), '../../../Matlab/SeaFreeze_Gibbs.mat')
        np.testing.assert_array_equal(exp, act)
    def test_grid_all_extrapolations(self):
        P = np.arange(0, 1001, 200)
        T = np.arange(200, 351, 50)
        exp = np.array([[1, 1, 0, 0], [2, 1, 0, 0], [2, 5, 0, 0], [2, 5, 0, 0], [6, 6, 0, 0], [6, 6, 6, 0]])
        act = sf.whichphase(np.array([P, T]), '../../../Matlab/SeaFreeze_Gibbs.mat')
        np.testing.assert_array_equal(exp, act)
    def test_phases(self):
        PT = np.empty((7,), np.object)
        PT[0] = (500, 300)  # liquid water
        PT[1] = (1, 200)    # ice Ih
        PT[2] = (300, 200)  # ice II
        PT[3] = (300, 250)  # ice III
        PT[4] = (500, 250)  # ice V
        PT[5] = (800, 250)  # ice VI
        PT[6] = (2000, 334) # extreme pressures at low temps should map to ice VI
        act = sf.whichphase(PT, '../../../Matlab/SeaFreeze_Gibbs.mat')
        np.testing.assert_array_equal([0, 1, 2, 3, 5, 6, 6], act)  # will throw AssertionError if not equal

