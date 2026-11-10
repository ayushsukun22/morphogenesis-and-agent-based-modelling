import unittest
import numpy as np
from branching_traction_budget import allocate,persistent


class AllocationTests(unittest.TestCase):
    def test_force_budget_and_anchors(self):
        d=np.array([.8,.1,.4,.9]); fixed=np.array([True,False,False,False])
        for b in [0,.25,1]:
            f,fallback=allocate(d,fixed,7.5,b)
            self.assertFalse(fallback); self.assertAlmostEqual(float(f[:,0].sum()),7.5)
            np.testing.assert_array_equal(f[fixed],0); np.testing.assert_array_equal(f[:,1:],0)
            self.assertTrue((f[:,0]>=0).all())

    def test_uniform_control_and_activity_weighting(self):
        d=np.array([.1,.5,.9]); fixed=np.zeros(3,bool)
        uniform,_=allocate(d,fixed,7.5,1)
        np.testing.assert_allclose(uniform[:,0],2.5)
        selective,_=allocate(d,fixed,7.5,0)
        self.assertGreater(selective[-1,0],uniform[-1,0])
        self.assertLess(selective[0,0],uniform[0,0])

    def test_total_scaling_and_zero_activity_fallback(self):
        d=np.array([.2,.8]); fixed=np.zeros(2,bool)
        a,_=allocate(d,fixed,2.5,.25); b,_=allocate(d,fixed,7.5,.25)
        np.testing.assert_allclose(b,3*a)
        f,fallback=allocate(np.zeros(2),fixed,2.5,0)
        self.assertTrue(fallback); np.testing.assert_allclose(f[:,0],1.25)

    def test_persistence_requires_three_matching_samples(self):
        self.assertEqual(persistent([[{1,2,3}],[{2,3,4}],[{3,4,5}]]),1)
        self.assertEqual(persistent([[{1,2,3}],[],[{1,2,3}]]),0)


if __name__=='__main__': unittest.main()
