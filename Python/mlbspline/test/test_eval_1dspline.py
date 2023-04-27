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
        x = np.empty(1, np.object)              # sadly, must do this so numpy won't nest or unnest against your will
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


if __name__ == '__main__':
    ut.main()




