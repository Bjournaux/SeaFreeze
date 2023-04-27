import warnings, unittest as ut
import numpy as np
import scipy.io as sio
from mlbspline import load, eval


class TestEval2DSpline(ut.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.spline = load.loadSpline('spline2d_v7.mat')
        self.mlout = sio.loadmat('spline2d_out.mat')['sp2d_out']
    def tearDown(self):
        pass
    def test2dsplineeval(self):
        x = np.empty(2, np.object)
        x[0] = np.logspace(-1, 3.7, 50)
        x[1] = np.linspace(250, 800, 20)
        out = eval.evalMultivarSpline(self.spline, x)
        self.assertEqual(out.shape, self.mlout.shape, 'shapes not equal')
        if not (np.allclose(out, self.mlout, rtol=0, atol=1e-10)  # check both abs and rel differences
                and np.allclose(out, self.mlout, rtol=1e-10, atol=0)):
            absDiffs = np.absolute(out - self.mlout)
            relDiffs = absDiffs / np.absolute(self.mlout)
            self.fail('Output has absolute differences as large as ' + str(np.max(absDiffs)) + \
                      ' and relative differences as large as ' + str(np.max(relDiffs)) + '.\n')
    def test2dsplineeval_singlepoint(self):
        out = eval.evalMultivarSpline(self.spline, (0.1, 250))
        self.assertEqual(1, out.size, 'Output has too many values')
        mlout = self.mlout[0][0]
        if not (np.allclose(out, mlout, rtol=0, atol=1e-10)  # check both abs and rel differences
                and np.allclose(out, mlout, rtol=1e-10, atol=0)):
            absDiffs = np.absolute(out - mlout)
            relDiffs = absDiffs / np.absolute(mlout)
            self.fail('Output has absolute differences as large as ' + str(np.max(absDiffs)) + \
                      ' and relative differences as large as ' + str(np.max(relDiffs)) + '.\n')

if __name__ == '__main__':
    ut.main()
