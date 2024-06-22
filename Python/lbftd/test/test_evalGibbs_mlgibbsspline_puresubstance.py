from random import randint
import warnings, unittest as ut
import numpy as np
import scipy.io as sio
from mlbspline import load
from lbftd import loadGibbs as lg, evalGibbs as eg


class TestEvalGibbsPureSubstance(ut.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.spline = lg.loadGibbsSpline('gsp_puresubstance.mat')['sp']
        self.mlout = load._stripNestingToFields(sio.loadmat('gsp2d_out.mat')['gsp2d_out'])
        # knot ranges for this spline are P [0-3001], T[0-401]
        self.P = np.arange(0, 3001, 200).astype(float)
        self.T = np.arange(0, 401, 50).astype(float)
    def tearDown(sThermodyelf):
        pass
    def test_evalgibbs_puresubstance_grid_allmeasures(self):
        out = eg.evalSolutionGibbsGrid(self.spline, np.array([self.P, self.T], dtype=object))
        valErrs = ''
        # check all values and output just one error for all of them
        for tdv in vars(out).keys():
            outfield = getattr(out, tdv)
            if tdv not in self.mlout.dtype.fields:
                warnings.warn('Matlab output does not include tdv ' + tdv)
            else:
                mloutfield = self.mlout[tdv]
                self.assertEqual(outfield.shape, mloutfield.shape, tdv+' output not the same shape as MatLab output')
                if not (np.allclose(outfield, mloutfield, rtol=relTolerance, atol=0)  # check both abs and rel differences
                        and np.allclose(outfield, mloutfield, rtol=0, atol=absTolerance)):
                    absDiffs = np.absolute(outfield - mloutfield)
                    relDiffs = absDiffs / np.absolute(mloutfield)
                    valErrs = valErrs + 'Output for '+tdv+' has absolute differences as large as '+str(np.max(absDiffs)) +\
                              ' and relative differences as large as '+str(np.max(relDiffs))+'.\n'
        if valErrs:
            self.fail(valErrs)

    def test_evalGibbs_puresubstance_grid_confirmExtrapolationBehavior(self):
        P = np.array([self.P[1], 3002], float)
        T = np.array([self.T[1], 402], float)
        extrapOut = eg.evalSolutionGibbsGrid(self.spline, np.array([P, T], dtype=object), allowExtrapolations=True)
        noExtrapsOut = eg.evalSolutionGibbsGrid(self.spline, np.array([P, T], dtype=object), allowExtrapolations=False)
        # check a TDV that directly requires allowExtrapolations
        self.assertTrue(np.all(np.logical_not(np.isnan(extrapOut.G))),
                        'extrapolations are missing although they were allowed')
        self.assertTrue(np.all(np.logical_not(np.isnan(noExtrapsOut.G[0, 0]))),
                        'valid values are missing when extrapolations not allowed')
        self.assertTrue(np.all(np.isnan(noExtrapsOut.G[1, :])), 'extrapolations are present in the first dim')
        self.assertTrue(np.all(np.isnan(noExtrapsOut.G[:, 1])), 'extrapolations are present in the second dim')
        # check a TDV that depends on derivatives only
        self.assertTrue(np.all(np.logical_not(np.isnan(extrapOut.S))),
                        'extrapolations are missing although they were allowed')
        self.assertTrue(np.all(np.logical_not(np.isnan(noExtrapsOut.S[0, 0]))),
                        'valid values are missing when extrapolations not allowed')
        self.assertTrue(np.all(np.isnan(noExtrapsOut.S[1, :])), 'extrapolations are present in the first dim')
        self.assertTrue(np.all(np.isnan(noExtrapsOut.S[:, 1])), 'extrapolations are present in the second dim')

    def test_evalgibbs_puresubstance_scatter_singlepoint_allmeasures(self):
        pidx = 0; tidx = 0
        PTM = np.empty((1,), object)
        PTM[0] = (self.P[pidx], self.T[tidx])
        out = eg.evalSolutionGibbsScatter(self.spline, PTM)
        valErrs = ''
        # check all values and output just one error for all of them
        for tdv in vars(out).keys():
            outfield = getattr(out, tdv)
            self.assertEqual(1, outfield.size, 'Output for '+tdv+' has too many values')
            if tdv not in self.mlout.dtype.fields:
                warnings.warn('Matlab output does not include tdv ' + tdv)
            else:
                mloutfield = self.mlout[tdv][pidx][tidx]
                if not (np.allclose(outfield, mloutfield, rtol=relTolerance,
                                    atol=0)  # check both abs and rel differences
                        and np.allclose(outfield, mloutfield, rtol=0, atol=absTolerance)):
                    absDiffs = np.absolute(outfield - mloutfield)
                    relDiffs = absDiffs / np.absolute(mloutfield)
                    valErrs = valErrs + 'Output for ' + tdv + ' has absolute differences as large as ' + str(
                        np.max(absDiffs)) + \
                              ' and relative differences as large as ' + str(np.max(relDiffs)) + '.\n'
        if valErrs:
            self.fail(valErrs)
    def test_evalgibbs_puresubstance_scatter_multipoint_allmeasures(self):
        numpts = 4
        ptindices = np.empty((numpts,), object)
        PTM = np.empty((numpts,), object)
        for i in np.arange(0, numpts):        # get random data points (pick P/T separately)
            pidx = randint(0, len(self.P)-1)
            tidx = randint(0, len(self.T)-1)
            ptindices[i] = (pidx, tidx)
            PTM[i] = (self.P[pidx], self.T[tidx])
        out = eg.evalSolutionGibbsScatter(self.spline, PTM)
        valErrs = ''
        # check all values and output just one error for all of them
        for tdv in vars(out).keys():
            outfield = getattr(out, tdv)
            self.assertEqual(numpts, outfield.size, 'Output for ' + tdv + ' has too many values')
            if tdv not in self.mlout.dtype.fields:
                warnings.warn('Matlab output does not include tdv ' + tdv)
            else:
                # get randomly scattered outputs
                mlout = np.empty((numpts,), float)
                for i in np.arange(0, numpts):
                    mlout[i] = self.mlout[tdv][ptindices[i]]
                if not (np.allclose(outfield, mlout, rtol=relTolerance,
                                    atol=0)  # check both abs and rel differences
                        and np.allclose(outfield, mlout, rtol=0, atol=absTolerance)):
                    absDiffs = np.absolute(outfield - mlout)
                    relDiffs = absDiffs / np.absolute(mlout)
                    valErrs = valErrs + 'Output for ' + tdv + ' has absolute differences as large as ' + str(
                        np.max(absDiffs)) + \
                              ' and relative differences as large as ' + str(np.max(relDiffs)) + '.\n'
        if valErrs:
            self.fail(valErrs)



relTolerance = 5e-9
absTolerance = 1e-6

if __name__ == '__main__':
    ut.main()