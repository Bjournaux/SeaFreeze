import warnings, unittest as ut
import numpy as np
import scipy.io as sio
from mlbspline import load, eval


class TestEval3DSpline(ut.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.spline = load.loadSpline('spline3d_v7.mat')
        self.mlout = sio.loadmat('spline3d_out.mat')['sp3d_out']
    def tearDown(self):
        pass
    def test3dsplineeval_gridWithExtrapolations(self):
        x = np.empty(3, object)
        x[0] = np.arange(.1, 8001, 10)
        x[1] = np.arange(239, 502, 5)
        x[2] = np.arange(0, 8.1, .5)
        out = eval.evalMultivarSpline(self.spline, x)
        self.assertEqual(out.shape, self.mlout.shape, 'shapes not equal')
        # unfortunately, the three rounds of estimations that go on with a 3d spline result in more error
        if not (np.allclose(out, self.mlout, rtol=0, atol=2e-9)  # check both abs and rel differences
                and np.allclose(out, self.mlout, rtol=1e-10, atol=0)):
            absDiffs = np.absolute(out - self.mlout)
            relDiffs = absDiffs / np.absolute(self.mlout)
            self.fail('Output has absolute differences as large as ' + str(np.max(absDiffs)) +
                      ' and relative differences as large as ' + str(np.max(relDiffs)) + '.\n')
    def test3dsplineeval_singlepoint_tuple(self):
        out = eval.evalMultivarSpline(self.spline, (.1, 239, 0))
        self.assertEqual(1, out.size, 'Output has too many values')
        mlout = self.mlout[0][0][0]
        if not (np.allclose(out, mlout, rtol=0, atol=2e-9)  # check both abs and rel differences
                and np.allclose(out, mlout, rtol=1e-10, atol=0)):
            absDiffs = np.absolute(out - mlout)
            relDiffs = absDiffs / np.absolute(mlout)
            self.fail('Output has absolute differences as large as ' + str(np.max(absDiffs)) +
                      ' and relative differences as large as ' + str(np.max(relDiffs)) + '.\n')
    def test3dsplineeval_singlepoint_arrofarr(self):
        x = np.empty(3, object)
        for i in np.arange(0, 3):
            x[i] = np.empty((1,), float)
        x[0][0] = 0.1
        x[1][0] = 244
        x[2][0] = 1
        out = eval.evalMultivarSpline(self.spline, x)
        self.assertEqual(1, out.size, 'Output has too many values')
        mlout = self.mlout[0][1][2]
        if not (np.allclose(out, mlout, rtol=0, atol=2e-9)  # check both abs and rel differences
                and np.allclose(out, mlout, rtol=1e-10, atol=0)):
            absDiffs = np.absolute(out - mlout)
            relDiffs = absDiffs / np.absolute(mlout)
            self.fail('Output has absolute differences as large as ' + str(np.max(absDiffs)) +
                      ' and relative differences as large as ' + str(np.max(relDiffs)) + '.\n')

    def test3dsplineeval_gridDisallowExtrapolations(self):
        x = np.empty(3, object)
        x[0] = np.arange(.1, 8001, 10)
        x[1] = np.arange(239, 502, 5)
        x[2] = np.arange(0, 8.1, .5)
        # try to eval a couple of invalid points at both ends of all dimensions
        x[0][0] = .001; # knots for this dim go from .1 - 8001
        x[0][-1] = 9000
        x[1][0] = 200; # knots for this dim go from 239 - 501
        x[1][-1] = 1400
        x[2][0] = -.5; # knots for this dim go from 0 - 8
        x[2][-1] = 10
        out = eval.evalMultivarSpline(self.spline, x, allowExtrapolations=False)
        self.assertEqual(out.shape, self.mlout.shape, 'shapes not equal')
        # confirm nan for invalid points
        badslc = [0, -1]
        self.assertTrue(np.all(np.isnan(out[badslc, :, :])), 'some extrapolated values remain for first dim')
        self.assertTrue(np.all(np.isnan(out[:, badslc, :])), 'some extrapolated values remain for second dim')
        self.assertTrue(np.all(np.isnan(out[:, :, badslc])), 'some extrapolated values remain for third dim')
        # confirm values for valid points
        goodslc = tuple([slice(1, -1)] * 3)
        if not (np.allclose(out[goodslc], self.mlout[goodslc], rtol=0, atol=2e-9)  # check both abs and rel differences
                and np.allclose(out[goodslc], self.mlout[goodslc], rtol=1e-10, atol=0)):
            absDiffs = np.absolute(out[goodslc] - self.mlout[goodslc])
            relDiffs = absDiffs / np.absolute(self.mlout[goodslc])
            self.fail('Output has absolute differences as large as ' + str(np.max(absDiffs)) + \
                      ' and relative differences as large as ' + str(np.max(relDiffs)) + '.\n')


if __name__ == '__main__':
    ut.main()
