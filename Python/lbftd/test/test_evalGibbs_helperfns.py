import warnings, unittest as ut
from collections.__init__ import namedtuple

import numpy as np

from lbftd.statevars import _getTDVSpec
from lbftd import loadGibbs as lg, evalGibbs as eg


class TestEvalGibbsHelperFns(ut.TestCase):
    '''
    Tests functions that are called by evalSolutionGibbs*
    '''
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
        self.puresubstancespline = lg.loadGibbsSpline('gsp_puresubstance.mat')
    def tearDown(self):
        pass
    def evalFoo(self):
        ''' a fake eval fn used for testing _buildEvalArgs'''
        pass
    #########################################
    ## 2d spline _checkInputs tests
    #########################################
    # ignore req0X checks for 2d spline
    # no tdv yet defined requires X but not f
    # no tdv defined that requires solvent molecular weight for a pure substance
    # no tdv yet defined that requires MWu but not f or X
    def test_checkinputs_grid_2d_reqFOnly(self):
        P = np.arange(0, 3001, 100)
        T = np.arange(0, 401, 50)
        with self.assertRaisesRegex(ValueError, 'You cannot calculate {\'mus\'} with a spline ' +
                         'that does not include concentration. Remove those statevars and all their dependencies, ' +
                         'or supply a spline that includes concentration.'):
            eg.evalSolutionGibbsGrid(self.puresubstancespline['sp'], np.array([P, T], dtype=object), 'mus')


    def test_checkinputs_grid_2d_reqFandX(self):
        P = np.arange(0, 3001, 100)
        T = np.arange(0, 401, 50)
        with self.assertRaisesRegex(ValueError, 'You cannot calculate {\'muw\'} with a spline ' +
                         'that does not include concentration. Remove those statevars and all their dependencies, ' +
                         'or supply a spline that includes concentration.'):
            eg.evalSolutionGibbsGrid(self.puresubstancespline['sp'], np.array([P, T], dtype=object), 'muw')

    def test_checkinputs_grid_2d_extrapPOnly(self):
        P = np.arange(0, 5001, 100)
        T = np.arange(0, 401, 50)
        with self.assertRaisesRegex(ValueError, 'Dimensions {\'P\'} ' +
                'contain values that fall outside the knot sequence for the given spline.'):
            eg.evalSolutionGibbsGrid(self.puresubstancespline['sp'], np.array([P, T], dtype=object), 'G',
                                     failOnExtrapolate=True)

    def test_checkinputs_grid_2d_extrapTOnly(self):
        P = np.arange(0, 3001, 100)
        T = np.arange(0, 601, 50)
        with self.assertRaisesRegex(ValueError, 'Dimensions {\'T\'} ' +
                'contain values that fall outside the knot sequence for the given spline.'):
            eg.evalSolutionGibbsGrid(self.puresubstancespline['sp'], np.array([P, T], dtype=object), 'G',
                                     failOnExtrapolate=True)

    def test_checkinputs_grid_2d_extrapAllDims(self):
        P = np.arange(0, 5001, 100)
        T = np.arange(0, 601, 50)
        # no order imposed on list of failed dimensions so accommodate either order
        with self.assertRaisesRegex(ValueError, 'Dimensions {\'[PT]\', \'[PT]\'} ' +
                'contain values that fall outside the knot sequence for the given spline.'):
            eg.evalSolutionGibbsGrid(self.puresubstancespline['sp'], np.array([P, T], dtype=object), 'G',
                                     failOnExtrapolate=True)


    def test_checkinputs_scatter_2d_reqFOnly(self):
        PTM = np.empty((1,), object)
        PTM[0] = (15, 40)
        with self.assertRaisesRegex(ValueError, 'You cannot calculate {\'mus\'} with a spline ' +
                         'that does not include concentration. Remove those statevars and all their dependencies, ' +
                         'or supply a spline that includes concentration.'):
            eg.evalSolutionGibbsScatter(self.puresubstancespline['sp'], PTM, 'mus')

    def test_checkinputs_scatter_2d_reqFandX(self):
        PTM = np.empty((1,), object)
        PTM[0] = (15, 40)
        with self.assertRaisesRegex(ValueError, 'You cannot calculate {\'muw\'} with a spline ' +
                         'that does not include concentration. Remove those statevars and all their dependencies, ' +
                         'or supply a spline that includes concentration.'):
            eg.evalSolutionGibbsScatter(self.puresubstancespline['sp'], PTM, 'muw')

    def test_checkinputs_scatter_2d_extrapPOnly(self):
        PTM = np.empty((1,), object)
        PTM[0] = (5000, 400)
        with self.assertRaisesRegex(ValueError, 'Dimensions {\'P\'} ' +
                'contain values that fall outside the knot sequence for the given spline.'):
            eg.evalSolutionGibbsScatter(self.puresubstancespline['sp'], PTM, failOnExtrapolate=True)

    def test_checkinputs_grid_2d_extrapTOnly(self):
        PTM = np.empty((1,), object)
        PTM[0] = (3000, 600)
        with self.assertRaisesRegex(ValueError, 'Dimensions {\'T\'} ' +
                'contain values that fall outside the knot sequence for the given spline.'):
            eg.evalSolutionGibbsScatter(self.puresubstancespline['sp'], PTM, failOnExtrapolate=True)

    def test_checkinputs_scatter_2d_extrapAllDims(self):
        PTM = np.empty((1,), object)
        PTM[0] = (5000, 600)
        # no order imposed on list of failed dimensions so accommodate either order
        with self.assertRaisesRegex(ValueError, 'Dimensions {\'[PT]\', \'[PT]\'} ' +
                'contain values that fall outside the knot sequence for the given spline.'):
            eg.evalSolutionGibbsScatter(self.puresubstancespline['sp'], PTM, failOnExtrapolate=True)
    #########################################
    ## 2d spline _buildEvalArgs tests
    #########################################
    def test_buildEvalArgs_derivsCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqDerivs=['d1P'], parmderivs='fooderivs')
        fooargs = eg._buildEvalArgs(tdv, derivs='foo', PTM=None, gPTM=None, tdvout=None,
                                    gibbsSp=None, f=None, allowExtrapolations=None)
        self.assertEqual(1, len(fooargs))
        self.assertFalse('derivs' in fooargs)
        self.assertTrue(fooargs['fooderivs'] == 'foo')
    def test_buildEvalArgs_gridCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqGrid=True, parmgrid='foogrid')
        fooargs = eg._buildEvalArgs(tdv, derivs=None, PTM=None, gPTM='foo', tdvout=None,
                                    gibbsSp=None, f=None, allowExtrapolations=None)
        self.assertEqual(1, len(fooargs))
        self.assertFalse('gPTM' in fooargs)
        self.assertTrue(fooargs['foogrid'] == 'foo')
    def test_buildEvalArgs_MWvCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqMWv=True, parmMWv='foomwv')
        fooargs = eg._buildEvalArgs(tdv, derivs=None, PTM=None, gPTM=None, tdvout=None,
                                    gibbsSp={'MW': list(['foo', 'foo2'])}, f=None, allowExtrapolations=None)
        self.assertEqual(1, len(fooargs))
        self.assertFalse('MWv' in fooargs)
        self.assertTrue(fooargs['foomwv'] == 'foo')
    def test_buildEvalArgs_MWuCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqMWu=True, parmMWu='foomwu')
        fooargs = eg._buildEvalArgs(tdv, derivs=None, PTM=None, gPTM=None, tdvout=None,
                                    gibbsSp={'MW': list(['foo1', 'foo'])}, f=None, allowExtrapolations=None)
        self.assertEqual(1, len(fooargs))
        self.assertFalse('MWu' in fooargs)
        self.assertTrue(fooargs['foomwu'] == 'foo')
    def test_buildEvalArgs_nuCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqNu=True, parmNu='foonu')
        fooargs = eg._buildEvalArgs(tdv, derivs=None, PTM=None, gPTM=None, tdvout=None,
                                    gibbsSp={'nu': 'foo'}, f=None, allowExtrapolations=None)
        self.assertEqual(1, len(fooargs))
        self.assertFalse('nu' in fooargs)
        self.assertTrue(fooargs['foonu'] == 'foo')
    def test_buildEvalArgs_TDVCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqTDV=['a'], parmtdv='footdv')
        fooargs = eg._buildEvalArgs(tdv, derivs=None, PTM=None, gPTM=None, tdvout='foo',
                                    gibbsSp=None, f=None, allowExtrapolations=None)
        self.assertEqual(1, len(fooargs))
        self.assertFalse('tdv' in fooargs)
        self.assertTrue(fooargs['footdv'] == 'foo')
    def test_buildEvalArgs_splineCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqSpline=True, parmspline='foospline')
        fooargs = eg._buildEvalArgs(tdv, derivs=None, PTM=None, gPTM=None, tdvout=None,
                                    gibbsSp='foo', f=None, allowExtrapolations=None)
        self.assertEqual(1, len(fooargs))
        self.assertFalse('gibbsSp' in fooargs)
        self.assertTrue(fooargs['foospline'] == 'foo')
    def test_buildEvalArgs_PTMCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqPTM=True, parmptm='fooptm')
        fooargs = eg._buildEvalArgs(tdv, derivs=None, PTM='foo', gPTM=None, tdvout=None,
                                    gibbsSp=None, f=None, allowExtrapolations=None)
        self.assertEqual(1, len(fooargs))
        self.assertFalse('PTM' in fooargs)
        self.assertTrue(fooargs['fooptm'] == 'foo')
    def test_buildEvalArgs_FCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqF=True, parmf='foof')
        fooargs = eg._buildEvalArgs(tdv, derivs=None, PTM=None, gPTM=None, tdvout=None,
                                    gibbsSp=None, f='foo', allowExtrapolations=None)
        self.assertEqual(1, len(fooargs))
        self.assertFalse('f' in fooargs)
        self.assertTrue(fooargs['foof'] == 'foo')
    def test_buildEvalArgs_AllowExtrapCustomParm(self):
        tdv = _getTDVSpec('foo', self.evalFoo, reqAllowExtrap=True, parmAllowExtrap='fooae')
        fooargs = eg._buildEvalArgs(tdv, derivs=None, PTM=None, gPTM=None, tdvout=None,
                                    gibbsSp=None, f=None, allowExtrapolations='foo')
        self.assertEqual(1, len(fooargs))
        self.assertFalse('allowExtrapolations' in fooargs)
        self.assertTrue(fooargs['fooae'] == 'foo')

    def test_isPointExtrapolation_2DTrue(self):
        knots = np.array([np.arange(2), np.arange(3)], object)
        self.assertTrue(eg._isPointExtrapolation(2, knots, np.array([np.array([3.1]), np.array([1.2])], object)))
        self.assertTrue(eg._isPointExtrapolation(2, knots, np.array([np.array([.5]), np.array([2.2])], object)))

    def test_isPointExtrapolation_2DFalse(self):
        knots = np.array([np.arange(2), np.arange(3)], object)
        self.assertFalse(eg._isPointExtrapolation(2, knots, np.array([np.array([0.1]), np.array([1.2])], object)))

    def test_isPointExtrapolation_3DTrue(self):
        knots = np.array([np.arange(3), np.arange(2), np.arange(4)], object)
        self.assertTrue(eg._isPointExtrapolation(3, knots, np.array([np.array([3.1]), np.array([1.2]), np.array([2.4])], object)))
        self.assertTrue(eg._isPointExtrapolation(3, knots, np.array([np.array([.5]), np.array([2.2]), np.array([2.4])], object)))
        self.assertTrue(eg._isPointExtrapolation(3, knots, np.array([np.array([.5]), np.array([1.2]), np.array([4.6])], object)))

    def test_isPointExtrapolation_3DFalse(self):
        knots = np.array([np.arange(3), np.arange(2), np.arange(4)], object)
        self.assertFalse(
            eg._isPointExtrapolation(3, knots, np.array([np.array([1.1]), np.array([0.2]), np.array([2.4])], object)))


#TODO: write additional test for if the length of MW is greater or less than 2 (only if reqMWu), as well as a test for if MW is a number if reqMWv


if __name__ == '__main__':
    ut.main()