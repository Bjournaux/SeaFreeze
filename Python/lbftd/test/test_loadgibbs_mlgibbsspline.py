import warnings, unittest as ut
import numpy as np
import scipy.io as sio
from lbftd import loadGibbs as lg


class TestLoadGibbsMLSpline(ut.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
    def tearDown(self):
        pass
    def test_loadgibbs_3d_ScalarGo(self):
        gsp = lg.loadGibbsSpline('gsp_1scalarGo.mat')
        # just do a spot check of a field in the spline since it calls load.getSplineDict which is tested elsewhere
        self.assertTrue(np.array_equal(gsp['sp']['number'], np.array([29, 20, 14])))
        self.assertTrue(np.array_equal(gsp['MW'], np.array([0.01801528, 0.05844])))
        self.assertEqual(gsp['nu'], 2)
    def test_loadgibbs_2d(self):
        gsp = lg.loadGibbsSpline('gsp_puresubstance.mat')
        self.assertTrue(np.array_equal(gsp['sp']['number'], np.array([50, 40])))
        self.assertTrue(np.array_equal(gsp['MW'], np.array([0.01801528])))
        self.assertFalse('nu' in gsp)


if __name__ == '__main__':
    ut.main()
