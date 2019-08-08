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
        PT[0] = (950, 255)
        self.assertEqual(6, sf.whichphase(PT, '../../../Matlab/SeaFreeze_Gibbs.mat')[0])
    def test_scatter_multiPt(self):
        PT = np.empty((3,), np.object)
        PT[0] = (100, 200)
        PT[1] = (400, 250)
        PT[2] = (1000, 300)
        self.assertTrue(np.array_equal([1, 5, 6], sf.whichphase(PT, '../../../Matlab/SeaFreeze_Gibbs.mat')))
    def test_grid(self):
        P = np.arange(0, 1001, 200)
        T = np.arange(200, 351, 50)
        exp = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [3, 5, 0, 0], [5, 5, 0, 0], [6, 6, 0, 0], [6, 6, 6, 0]])
        self.assertTrue(np.array_equal(exp, sf.whichphase(np.array([P, T]), '../../../Matlab/SeaFreeze_Gibbs.mat')))
