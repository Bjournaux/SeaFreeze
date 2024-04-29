from random import randint
import warnings, unittest as ut
import numpy as np
import scipy.io as sio
from mlbspline import load
from lbftd import loadGibbs as lg, evalGibbs as eg, statevars as sv


class TestEvalGibbsSingleSolute(ut.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.spline = lg.loadGibbsSpline('gsp_singlesolute.mat', 'sp_NaCl')
        self.spline = self.spline['sp']
        self.mlout = load._stripNestingToFields(sio.loadmat('gsp3d_out.mat')['out'])
        self.P = np.arange(0.1, 8000, 200).astype(float)
        self.T = np.arange(239, 501, 50).astype(float)
        self.M = np.arange(0.5, 8, 2).astype(float)
        self.ssM = 1e-3  # standard state solution concentration (1 mol/Kg = 0.001 mol/cc)
        self.stP = 0.1 # STP pressure in MPa
    def tearDown(sThermodyelf):
        pass
    def test_evalgibbs_singlesolute_grid_allmeasures(self):
        out = eg.evalSolutionGibbsGrid(self.spline, np.array([self.P, self.T, self.M], dtype=object))
        valErrs = ''
        # check all values and output just one error for all of them
        for tdv in vars(out).keys():
            outfield = getattr(out, tdv)
            if tdv not in self.mlout.dtype.fields:
                warnings.warn('Matlab output does not include tdv ' + tdv)
            else:
                mloutfield = self.mlout[tdv]
                self.assertEqual(outfield.shape, mloutfield.shape, tdv+' output not the same shape as MatLab output')
                if not np.allclose(outfield, mloutfield, rtol=relTolerance, atol=0):
                    absDiffs = np.absolute(outfield - mloutfield)
                    relDiffs = absDiffs / np.absolute(mloutfield)
                    valErrs = valErrs + 'Output for '+tdv+' has relative differences as large as '+str(np.max(relDiffs))+'.\n'
        if valErrs:
            self.fail(valErrs)
    def test_evalgibbs_singlesolute_grid_Cpa_no0M(self):
        PTM = np.array([self.P, self.T, self.M[1:]], dtype=object)
        out = eg.evalSolutionGibbsGrid(self.spline, PTM)

        valErrs = ''
        for tdv in vars(out).keys():
            outfield = getattr(out, tdv)
            if tdv in self.mlout.dtype.fields:
                mloutfield = self.mlout[tdv][:, :, 1:]   # get everything but the 0 concentration in the original output
                self.assertEqual(outfield.shape, mloutfield.shape, tdv + ' output not the same shape as MatLab output')
                if not np.allclose(outfield, mloutfield, rtol=relTolerance, atol=0):
                    absDiffs = np.absolute(outfield - mloutfield)
                    relDiffs = absDiffs / np.absolute(mloutfield)
                    valErrs = valErrs + 'Output for ' + tdv + ' has relative differences as large as ' + str(np.max(relDiffs)) + '.\n'
        if valErrs:
            self.fail(valErrs)
    def test_evalgibbs_singlesolute_scatter_singlepoint_allmeasures(self):
        pidx = 0; tidx = 0; midx = 1;
        PTM = np.empty((1,), object)
        PTM[0] = (self.P[pidx], self.T[tidx], self.M[midx])
        out = eg.evalSolutionGibbsScatter(self.spline, PTM)
        valErrs = ''
        # check all values and output just one error for all of them
        for tdv in vars(out).keys():
            outfield = getattr(out, tdv)
            self.assertEqual(1, outfield.size, 'Output for '+tdv+' has too many values')
            if tdv not in self.mlout.dtype.fields:
                warnings.warn('Matlab output does not include tdv ' + tdv)
            else:
                mloutfield = self.mlout[tdv][pidx][tidx][midx]
                self.assertEqual(1, outfield.size, tdv + ' output not the same shape as MatLab output')
                if not np.allclose(outfield, mloutfield, rtol=relTolerance, atol=0):
                    absDiffs = np.absolute(outfield - mloutfield)
                    relDiffs = absDiffs / np.absolute(mloutfield)
                    valErrs = valErrs + 'Output for ' + tdv + ' has relative differences as large as ' + str(np.max(relDiffs)) + '.\n'
        if valErrs:
            self.fail(valErrs)
    def test_evalgibbs_singlesolute_scatter_singlepoint_Va_no0M(self):
        pidx = 1; tidx = 2; midx = 3;
        PTM = np.empty((1,), object)
        PTM[0] = (self.P[pidx], self.T[tidx], self.M[midx])
        out = eg.evalSolutionGibbsScatter(self.spline, PTM)
        valErrs = ''
        for tdv in vars(out).keys():
            outfield = getattr(out, tdv)
            self.assertEqual(1, outfield.size, 'Output for ' + tdv + ' has too many values')
            if tdv in self.mlout.dtype.fields:
                mloutfield = self.mlout[tdv][pidx][tidx][midx]
                if not np.allclose(outfield, mloutfield, rtol=relTolerance, atol=0):
                    absDiffs = np.absolute(outfield - mloutfield)
                    relDiffs = absDiffs / np.absolute(mloutfield)
                    valErrs = valErrs + 'Output for ' + tdv + ' has relative differences as large as ' + str(np.max(relDiffs)) + '.\n'
        if valErrs:
            self.fail(valErrs)
    def test_evalgibbs_singlesolute_scatter_multipoint_allmeasures(self):
        numpts = 4
        ptindices = np.empty((numpts,), object)
        PTM = np.empty((numpts,), object)
        for i in np.arange(0, numpts):        # get random data points (pick P/T separately)
            pidx = randint(0, len(self.P)-1)
            tidx = randint(0, len(self.T)-1)
            midx = randint(0, len(self.M)-1)
            ptindices[i] = (pidx, tidx, midx)
            PTM[i] = (self.P[pidx], self.T[tidx], self.M[midx])
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
                if not np.allclose(outfield, mlout, rtol=relTolerance, atol=0):
                    absDiffs = np.absolute(outfield - mlout)
                    relDiffs = absDiffs / np.absolute(mlout)
                    valErrs = valErrs + 'Output for ' + tdv + ' has relative differences as large as ' + str(np.max(relDiffs)) + '.\n'
        if valErrs:
            self.fail(valErrs)

    def test_G0_standard_state(self):
        out = sv._getG0ss(self.T, self.spline['Go'])
        mlGoss = sio.loadmat('Gss1.mat')['gss1']
        self.assertEqual(out.shape, mlGoss.shape,' output not the same shape as MatLab output')
        if not np.allclose(out, mlGoss, rtol=relTolerance, atol=0):
            absDiffs = np.absolute(out - mlGoss)
            relDiffs = absDiffs / np.absolute(mlGoss)
            self.fail('Output for G0ss has relative differences as large as ' + str(np.max(relDiffs)))

    def test_dGss_standard_state(self):
        out = sv._getdGss(self.P, self.T, self.spline, self.spline['MW'][1])
        mldGss = sio.loadmat('dGss.mat')['dgss']
        self.assertEqual(out.shape, mldGss.shape,' output not the same shape as MatLab output')
        if not np.allclose(out, mldGss, rtol=relTolerance, atol=0):
            absDiffs = np.absolute(out - mldGss)
            relDiffs = absDiffs / np.absolute(mldGss)
            self.fail('Output for dGss has relative differences as large as ' + str(np.max(relDiffs)))

    def test_G2(self):
        ss_PTM = np.array([self.P, self.T, np.array([self.ssM])], dtype=object)
        out = sv._getG2(ss_PTM, self.spline)
        mlG2 = sio.loadmat('g2.mat')['G2']
        self.assertEqual(out.shape, mlG2.shape,' output not the same shape as MatLab output')
        if not np.allclose(out, mlG2, rtol=relTolerance, atol=0):
            absDiffs = np.absolute(out - mlG2)
            relDiffs = absDiffs / np.absolute(mlG2)
            self.fail('Output for G2 has relative differences as large as ' + str(np.max(relDiffs)))

    def test_dGm2(self):
        ss_PTM = np.array([self.P, self.T, np.array([self.ssM])], dtype=object)
        out = sv._getdGm2(ss_PTM, self.spline)
        mldgm2 = sio.loadmat('dgdm2.mat')['dGdm2']
        self.assertEqual(out.shape, mldgm2.shape,' output not the same shape as MatLab output')
        if not np.allclose(out, mldgm2, rtol=relTolerance, atol=0):
            absDiffs = np.absolute(out - mldgm2)
            relDiffs = absDiffs / np.absolute(mldgm2)
            self.fail('Output for mldgm2 has relative differences as large as ' + str(np.max(relDiffs)))

    def test_G1b(self):
        PTM_ss_1_bar = np.array([np.array([self.stP]), self.T, np.array([self.ssM])], dtype=object)
        out = sv._getG1b(PTM_ss_1_bar, self.spline)
        mlg1b = sio.loadmat('g1b.mat')['G1b']
        self.assertEqual(out.shape, mlg1b.shape,' output not the same shape as MatLab output')
        if not np.allclose(out, mlg1b, rtol=relTolerance, atol=0):
            absDiffs = np.absolute(out - mlg1b)
            relDiffs = absDiffs / np.absolute(mlg1b)
            self.fail('Output for mlg1b has relative differences as large as ' + str(np.max(relDiffs)))

    def test_dGdm1(self):
        PTM_ss_1_bar = np.array([np.array([self.stP]), self.T, np.array([self.ssM])], dtype=object)
        out = sv._getdGdm1(PTM_ss_1_bar, self.spline)
        mldGdm1 = sio.loadmat('dgdm1.mat')['dGdm1']
        self.assertEqual(out.shape, mldGdm1.shape,' output not the same shape as MatLab output')
        if not np.allclose(out, mldGdm1, rtol=relTolerance, atol=0):
            absDiffs = np.absolute(out - mldGdm1)
            relDiffs = absDiffs / np.absolute(mldGdm1)
            self.fail('Output for mldGdm1 has relative differences as large as ' + str(np.max(relDiffs)))


relTolerance = 5e-6

if __name__ == '__main__':
    ut.main()

