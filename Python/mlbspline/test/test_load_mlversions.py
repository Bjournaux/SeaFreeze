import warnings, unittest as ut
import numpy as np
from mlbspline import load


class TestLoadSplineWithMatlabVersions(ut.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
    def tearDown(self):
        pass
    def check1dspline(self, sp):
        self.assertEqual(sp['form'], 'B-', 'form not properly loaded')
        self.assertEqual(sp['dim'], 1, 'dim not properly loaded')
        self.assertTrue(np.array_equal(sp['number'], np.array([131])), 'number not properly loaded')
        self.assertTrue(np.array_equal(sp['order'], np.array([4])), 'order not properly loaded')
        self.assertEqual(sp['knots'].ndim, 1, 'knots loaded with wrong nesting')
        self.assertEqual(sp['knots'][0].shape, (135,), 'knots not all loaded')
        # spot check a few values in knots
        self.assertEqual(sp['knots'][0][0], 240, 'first knot incorrect value')
        self.assertEqual(sp['knots'][0][49], 334, 'fiftieth knot incorrect value')
        self.assertEqual(sp['knots'][0][134], 500, 'last knot incorrect value')
        self.assertEqual(sp['coefs'].ndim, 1, 'coefs not properly nested')
        self.assertEqual(sp['coefs'].shape, (131,), 'coefs not all loaded')
        # spot check a few coefs - precision is number of decimal points  shown minus power of 10
        self.assertAlmostEqual(sp['coefs'][0], 7.934961866962184e+03, 12, 'first coef value not equal')
        self.assertAlmostEqual(sp['coefs'][9], 5.008173328406900e+03, 12, 'tenth coef value not equal')
        self.assertAlmostEqual(sp['coefs'][130], -1.833688771741746e+04, 11, 'last coef value not equal')
    def test1dsplinev7(self):
        sp = load.loadSpline('spline1d_v7.mat')
        self.check1dspline(sp)
    def test1dsplinev73(self):
        sp = load.loadSpline('spline1d_v73.mat')
        self.check1dspline(sp)
    def check2dspline(self, sp):
        self.assertEqual(sp['form'], 'B-', 'form not properly loaded')
        self.assertEqual(sp['dim'], 1, 'dim not properly loaded')
        self.assertTrue(np.array_equal(sp['number'], np.array([210, 203])), 'number not properly loaded')
        self.assertTrue(np.array_equal(sp['order'], np.array([6, 6])), 'order not properly loaded')
        self.assertEqual(sp['knots'].size, 2, 'knots not properly loaded')
        self.assertEqual(sp['knots'][0].shape, (216, ), 'knots first dim size incorrect')
        self.assertEqual(sp['knots'][1].shape, (209, ), 'knots second dim size incorrect')
        # spot check a few values in knots
        self.assertAlmostEqual(sp['knots'][0][0], 0., 15, 'first dim, first knot incorrect value')
        self.assertAlmostEqual(sp['knots'][0][49], 1.850000000000007e+03, 12, 'first dim, fiftieth knot incorrect value')
        self.assertAlmostEqual(sp['knots'][0][215], 10000., 15, 'first dim, last knot incorrect value')
        self.assertEqual(sp['knots'][1][0], 240., 'second dim, first knot incorrect value')
        self.assertEqual(sp['knots'][1][59], 520., 'second dim, sixtieth knot incorrect value')
        self.assertEqual(sp['knots'][1][208], 1250., 'second dim, last knot incorrect value')
        self.assertEqual(sp['coefs'].ndim, 2, 'coefs not properly loaded')
        self.assertEqual(sp['coefs'].shape, (210, 203), 'coefs dim sizes incorrect')
        # spot check a few coefs - precision is number of decimal points shown minus power of 10
        self.assertAlmostEqual(sp['coefs'][0, 0], 1.178305489097509e+03, 12, 'coef value a not equal')
        self.assertAlmostEqual(sp['coefs'][9, 9], 1.594732018089757e+03, 12, 'coef value b not equal')
        self.assertAlmostEqual(sp['coefs'][76, 77], 3.825051421428625e+03, 12, 'coef value c not equal')
        self.assertAlmostEqual(sp['coefs'][209, 202], 5.382233818477050e+03, 12, 'coef value d not equal')
    def test2dsplinev7(self):
        sp = load.loadSpline('spline2d_v7.mat')
        self.check2dspline(sp)
    def test2dsplinev73(self):
        sp = load.loadSpline('spline2d_v73.mat')
        self.check2dspline(sp)
    def check3dspline(self, sp):
        self.assertEqual(sp['form'], 'B-', 'form not properly loaded')
        self.assertEqual(sp['dim'], 1, 'dim not properly loaded')
        self.assertTrue(np.array_equal(sp['number'], np.array([29, 30, 14])), 'number not properly loaded')
        self.assertTrue(np.array_equal(sp['order'], np.array([6, 6, 4])), 'order not properly loaded')
        self.assertEqual(sp['knots'].size, 3, 'knots not properly loaded')
        self.assertEqual(sp['knots'][0].shape, (35, ), 'knots first dim size incorrect')
        self.assertEqual(sp['knots'][1].shape, (36, ), 'knots second dim size incorrect')
        self.assertEqual(sp['knots'][2].shape, (18, ), 'knots third dim size incorrect')
        # spot check a few values in knots
        self.assertAlmostEqual(sp['knots'][0][0], 0.1, 15, 'first dim, first knot incorrect value')
        self.assertAlmostEqual(sp['knots'][0][7], 40.1, 14, 'first dim, eighth knot incorrect value')
        self.assertAlmostEqual(sp['knots'][0][34], 8.000999999999998e+03, 12, 'first dim, last knot incorrect value')
        self.assertAlmostEqual(sp['knots'][1][0], 239., 15,  'second dim, first knot incorrect value')
        self.assertAlmostEqual(sp['knots'][1][18], 3.745172413793104e+02, 12, 'second dim, nineteenth knot incorrect value')
        self.assertAlmostEqual(sp['knots'][1][34], 501., 15, 'second dim, last knot incorrect value')
        self.assertEqual(sp['knots'][2][0], 0., 'third dim, first knot incorrect value')
        self.assertEqual(sp['knots'][2][9], 1.2, 'third dim, tenth knot incorrect value')
        self.assertEqual(sp['knots'][2][17], 8., 'third dim, last knot incorrect value')
        self.assertEqual(sp['coefs'].ndim, 3, 'coefs not properly loaded')
        self.assertEqual(sp['coefs'].shape, (29, 30, 14), 'coefs dim sizes incorrect')
        # spot check a few coefs - precision is number of decimal points shown minus power of 10
        self.assertAlmostEqual(sp['coefs'][0, 0, 0], -9.559719527277608e+03, 12, 'coef value a not equal')
        self.assertAlmostEqual(sp['coefs'][9, 9, 9], 7.851744906010371e+04, 11, 'coef value b not equal')
        self.assertAlmostEqual(sp['coefs'][11, 12, 13], 1.144636179511138e+05, 10, 'coef value c not equal')
        self.assertAlmostEqual(sp['coefs'][28, 29, 13], 0.007521531955442, 15, 'coef value d not equal')
    def test3dsplinev7(self):
        sp = load.loadSpline('spline3d_v7.mat')
        self.check3dspline(sp)
    def test3dsplinev73(self):
        sp = load.loadSpline('spline3d_v73.mat')
        self.check3dspline(sp)


if __name__ == '__main__':
    ut.main()

