import warnings, unittest as ut
import numpy as np
from mlbspline import load, eval


class TestEvalUtilFns(ut.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
    def tearDown(self):
        pass

    def testIsExtrapolation_gridAllValid(self):
        knots = np.arange(240, 501, 20)
        x = np.arange(250, 501, 10)
        extrapIdx = eval._isExtrapolation(knots, x)
        self.assertEqual(extrapIdx.size, 0, 'valid values considered extrapolations')

    def testIsExtrapolation_gridAllExtrap(self):
        knots = np.arange(240, 500, 20)
        x = np.arange(30, 200, 30)
        self.assertTrue(np.all(eval._isExtrapolation(knots, x) == np.arange(x.size)),
                        'extrapolative x values considered valid')

    def testIsExtrapolation_gridSomeInvalidOutsideKnotRange(self):
        knots = np.arange(240, 500, 20)  # doesn't include 500
        x = np.arange(150, 601, 50)
        expected = np.array([0, 1, 7, 8, 9])
        actual = eval._isExtrapolation(knots, x)
        self.assertTrue(np.logical_and(expected.size == actual.size, np.all(expected == actual)),
                        'extrapolations incorrectly identified in grid')

    def testIsExtrapolation_scatterSomeInvalidValues(self):
        knots = np.arange(240, 501, 20)
        x = np.array([425, 36, 250, 64, 600, 499])
        expected = np.array([1, 3, 4])
        actual = eval._isExtrapolation(knots, x)
        self.assertTrue(np.logical_and(expected.size == actual.size, np.all(expected == actual)),
                        'extrapolations incorrectly identified in scatter')

    def testSetExtrapolationsToNan_1dgrid(self):
        # _setExtrapolationsToNan(spd, x, y)
        spd = {'number': np.array([3], dtype=int),
               'knots': np.array([np.array([2, 3, 5], float)], dtype=object)}
        x = np.empty(1, object)
        x[0] = np.array([1, 2.4, 5, 5.5])
        y = np.random.random(x[0].size)
        out = eval._setExtrapolationsToNan(spd, x, y)
        # ensure out is nan for invalid values of x
        self.assertTrue(np.all(np.isnan(out[[0, 3]])), 'invalid values have not been set to nan')
        # ensure y is untouched
        self.assertTrue(np.all(np.logical_not(np.isnan(y[[0, 3]]))), 'y values were changed')
        # ensure all other values match y
        self.assertTrue(np.all(y[[1, 2]] == out[[1, 2]]), 'some out values do not match y input')

    def testSetExtrapolationsToNan_3dgrid(self):
        # _setExtrapolationsToNan(spd, x, y)
        spd = {'number': np.array([3, 4, 2], dtype=int),
               'knots': np.array([np.array([2.1, 3, 5]), np.array([8, 13.2, 21, 34]), np. array([5.5, 7])], object)}
        x = np.array([np.array([1, 2.4, 5, 5.5]), np.array([7.1, 8., 35]), np.array([1.2, 5.7, 7])], object)
        badslc = [[0, 3], [0, 2], [0]]
        goodslc = tuple([[1, 2], [1], [1, 2]])
        y = np.random.random((x[0].size, x[1].size, x[2].size))
        out = eval._setExtrapolationsToNan(spd, x, y)
        # ensure out is nan for invalid values of x
        self.assertTrue(np.all(np.isnan(out[badslc[0], :, :])), 'invalid values in the first dim have not been set to nan')
        self.assertTrue(np.all(np.isnan(out[:, badslc[1], :])), 'invalid values in the second dim have not been set to nan')
        self.assertTrue(np.all(np.isnan(out[:, :, badslc[2]])), 'invalid values in the third dim have not been set to nan')
        # ensure y is untouched
        self.assertTrue(np.all(np.logical_not(np.isnan(y[badslc[0], :, :]))), 'y values were changed in the first dim')
        self.assertTrue(np.all(np.logical_not(np.isnan(y[:, badslc[1], :]))), 'y values were changed in the second dim')
        self.assertTrue(np.all(np.logical_not(np.isnan(y[:, :, badslc[2]]))), 'y values were changed in the third dim')
        # ensure all other values match y
        self.assertTrue(np.all(y[goodslc] == out[goodslc]), 'some out values do not match y input')

if __name__ == '__main__':
    ut.main()

