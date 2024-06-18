import warnings, unittest as ut
import numpy as np
import scipy.io as sio
from mlbspline import load, eval


class TestEval1DSpline(ut.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.spline = load.loadSpline('spline1d_v7.mat')
        # TODO: shouldn't have to squeeze
        self.mlout = sio.loadmat('spline1d_out.mat')['sp1d_out'].squeeze()
    def tearDown(self):
        pass
    def test1dsplineeval_grid(self):
        x = np.empty(1, object)              # sadly, must do this so numpy won't nest or unnest against your will
        x[0] = np.arange(240, 501, 20)          # in Matlab, a = 240:20:500;
        out = eval.evalMultivarSpline(self.spline, x)
        self.assertEqual(out.shape, self.mlout.shape, 'shapes not equal')
        if not (np.allclose(out, self.mlout, rtol=0, atol=1e-10)  # check both abs and rel differences
                and np.allclose(out, self.mlout, rtol=1e-10, atol=0)):
            absDiffs = np.absolute(out - self.mlout)
            relDiffs = absDiffs / np.absolute(self.mlout)
            self.fail('Output has absolute differences as large as ' + str(np.max(absDiffs)) + \
                      ' and relative differences as large as ' + str(np.max(relDiffs)) + '.\n')
    def test1dsplineeval_singlepoint(self):
        out = eval.evalMultivarSpline(self.spline, (240,))
        self.assertEqual(1, out.size, 'Output has too many values')
        mlout = self.mlout[0]   # pick out just the upper left corner of the output
        self.assertTrue(np.allclose(out, mlout, rtol=0, atol=1e-11), 'output not within absolute tolerances')
        self.assertTrue(np.allclose(out, mlout, rtol=1e-12, atol=0), 'output not within relative tolerances')

    def test1dspline_isExtrapolation_GridAllValid(self):
        knots = np.arange(240, 501, 20)
        x = np.arange(250, 501, 10)
        extrapIdx = eval._isExtrapolation(knots, x)
        self.assertEqual(extrapIdx.size, 0, 'valid values considered extrapolations')

    def test1dspline_isExtrapolation_GridAllExtrap(self):
        knots = np.arange(240, 500, 20)
        x = np.arange(30, 200, 30)
        self.assertTrue(np.all(eval._isExtrapolation(knots, x) == np.arange(x.size)),
                        'extrapolative x values considered valid')

    def test1dspline_isExtrapolation_GridSomeInvalidOutsideKnotRange(self):
        knots = np.arange(240, 500, 20) # doesn't include 500
        x = np.arange(150, 601, 50)
        expected = np.array([0, 1, 7, 8, 9])
        actual = eval._isExtrapolation(knots, x)
        self.assertTrue(np.logical_and(expected.size == actual.size, np.all(expected == actual)),
                         'extrapolations incorrectly identified in grid')

    def test1dspline_isExtrapolation_ScatterSomeInvalidValues(self):
        knots = np.arange(240, 501, 20)
        x = np.array([425, 36, 250, 64, 600, 499])
        expected = np.array([1, 3, 4])
        actual = eval._isExtrapolation(knots, x)
        self.assertTrue(np.logical_and(expected.size == actual.size, np.all(expected == actual)),
                         'extrapolations incorrectly identified in scatter')

    def test1dsplineeval_noextrap_grid(self):
        x = np.empty(1, object)
        x[0] = np.arange(100, 601, 20)
        out = eval.evalMultivarSpline(self.spline, x, allowExtrapolations=False)
        self.assertEqual(out.shape, x[0].shape, 'incorrect length')
        invalidXIndices = eval._isExtrapolation(self.spline['knots'][0], x[0])
        validXMask = np.ones(len(x[0]), bool)
        validXMask[invalidXIndices] = 0
        self.assertTrue(np.logical_or(np.allclose(out[validXMask], self.mlout, rtol=0, atol=1e-11),
                                       np.allclose(out[validXMask], self.mlout, rtol=1e-12, atol=0)),
                        'output values within range of spline are not correct')
        self.assertTrue(np.all(np.isnan(out[invalidXIndices])), 'invalid value of x does not return nan')


if __name__ == '__main__':
    ut.main()




