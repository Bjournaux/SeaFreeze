import warnings, unittest as ut
from collections import Counter

import lbftd.statevars


class TestStatevarsExpandTDVSpecs(ut.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)
    def tearDown(self):
        pass
    #########################################
    ## 2d spline expandTDVSpec tests
    #########################################
    def test_evalgibbs_2d_expandEmptySpec(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec((), 2)]
        self.assertEqual(Counter(tdvs),
                         Counter(['G', 'rho', 'vel', 'Cp', 'Cv', 'alpha', 'S', 'U', 'A', 'H', 'Kt', 'Kp', 'Ks', 'V']))
    def test_evalgibbs_2d_errOnUnknown(self):
        with self.assertRaisesRegex(ValueError, 'One or more unsupported statevars have been requested: {\'foo\'}'):
            lbftd.statevars.expandTDVSpec(('foo',), 2)
    def test_evalgibbs_2d_expandG(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('G', 2)]
        self.assertEqual(Counter(tdvs), Counter(['G']))
    def test_evalgibbs_2d_expandrho(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('rho', 2)]
        self.assertEqual(Counter(tdvs), Counter(['rho']))
    def test_evalgibbs_2d_expandvel(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('vel', 2)]
        self.assertEqual(Counter(tdvs), Counter(['vel']))
    def test_evalgibbs_2d_expandCp(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('Cp', 2)]
        self.assertEqual(Counter(tdvs), Counter(['Cp']))
    def test_evalgibbs_2d_expandCv(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('Cv', 2)]
        self.assertEqual(Counter(tdvs), Counter(['Cv', 'Cp']))
    def test_evalgibbs_2d_expandalpha(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('alpha', 2)]
        self.assertEqual(Counter(tdvs), Counter(['alpha', 'rho']))
    def test_evalgibbs_2d_expandU(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('U', 2)]
        self.assertEqual(Counter(tdvs), Counter(['U', 'G', 'rho', 'S']))
    def test_evalgibbs_2d_expandA(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('A', 2)]
        self.assertEqual(Counter(tdvs), Counter(['A', 'U', 'S', 'G', 'rho']))
    def test_evalgibbs_2d_expandH(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('H', 2)]
        self.assertEqual(Counter(tdvs), Counter(['H', 'S', 'G']))
    def test_evalgibbs_2d_expandS(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('S', 2)]
        self.assertEqual(Counter(tdvs), Counter(['S']))
    def test_evalgibbs_2d_expandKt(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('Kt', 2)]
        self.assertEqual(Counter(tdvs), Counter(['Kt']))
    def test_evalgibbs_2d_expandKp(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('Kp', 2)]
        self.assertEqual(Counter(tdvs), Counter(['Kp']))
    def test_evalgibbs_2d_expandKs(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('Ks', 2)]
        self.assertEqual(Counter(tdvs), Counter(['Ks', 'rho', 'vel']))
    def test_evalgibbs_2d_expandV(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('V', 2)]
        self.assertEqual(Counter(tdvs), Counter(['V', 'rho']))
    def test_evalgibbs_2d_expandmultiple(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec(('V', 'Cp', 'U'), 2)]
        self.assertEqual(Counter(tdvs), Counter(['V', 'rho', 'Cp', 'U', 'G', 'S']))
    def test_statevars_2d_list(self):
        orderedTdvs = ['G', 'U', 'A', 'H', 'S', 'rho', 'V', 'Cp', 'Cv', 'Kt', 'Kp', 'Ks', 'vel', 'alpha']
        self.assertListEqual(orderedTdvs, [v.name for v in lbftd.statevars.tdvsPTOnly])
    #########################################
    ## 3d spline expandTDVSpec tests
    #########################################
    def test_evalgibbs_3d_expandEmptySpec(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec((), 3)]
        self.assertEqual(Counter(tdvs),
                         Counter(['G', 'rho', 'vel', 'Cp', 'Cv', 'alpha', 'S', 'U', 'A', 'H', 'Kt', 'Kp', 'Ks', 'V', 'mus', 'muw', 'Vm', 'Cpm', 'Cpa', 'Va']))
    def test_evalgibbs_3d_expandmus(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('mus', 3)]
        self.assertEqual(Counter(tdvs), Counter(['mus', 'G']))
    def test_evalgibbs_3d_expandmuw(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('muw', 3)]
        self.assertEqual(Counter(tdvs), Counter(['muw', 'G']))
    def test_evalgibbs_3d_expandVm(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('Vm', 3)]
        self.assertEqual(Counter(tdvs), Counter(['Vm']))
    def test_evalgibbs_3d_expandCpm(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('Cpm', 3)]
        self.assertEqual(Counter(tdvs), Counter(['Cpm', 'Cp']))
    def test_evalgibbs_3d_expandCpa(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('Cpa', 3)]
        self.assertEqual(Counter(tdvs), Counter(['Cpa', 'Cp']))
    def test_evalgibbs_3d_expandVa(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec('Va', 3)]
        self.assertEqual(Counter(tdvs), Counter(['Va', 'V', 'rho']))
    def test_evalgibbs_3d_expandmultiple(self):
        tdvs = [t.name for t in lbftd.statevars.expandTDVSpec(('muw', 'Va', 'Vm'), 3)]
        self.assertEqual(Counter(tdvs), Counter(['muw', 'G', 'Va', 'V', 'rho', 'Vm']))
    def test_statevars_3d_list(self):
        orderedTdvs = ['G', 'U', 'A', 'H', 'S', 'rho', 'V', 'Cp', 'Cv', 'Kt', 'Kp', 'Ks', 'vel', 'alpha', 'mus', 'muw', 'Vm', 'Cpm', 'Cpa', 'Va']
        self.assertListEqual(orderedTdvs, [tdv.name for tdv in lbftd.statevars.statevars])



if __name__ == '__main__':
    ut.main()