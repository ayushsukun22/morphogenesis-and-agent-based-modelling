import unittest
import numpy as np
from branching_junction_sprouting import mechanical_interval,parent_attachment,observe,prepattern
from branching_sprouting_assay import initial,contact_graph,regulate


class CouplingTests(unittest.TestCase):
    def test_fixed_anchors_and_passive_interval(self):
        x,fixed,_=initial(42); rest=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2)
        y,b,stats=mechanical_interval(x,fixed,contact_graph(x),rest,np.zeros_like(x),32)
        np.testing.assert_array_equal(y[fixed],x[fixed])
        self.assertEqual(stats['increases'],0); self.assertEqual(stats['broken'],0)
        self.assertEqual(stats['min_material_attachment'],1)

    def test_junction_break_is_symmetric(self):
        x=np.array([[0.,0,0],[1.59,0,0]])
        fixed=np.array([True,False]); bonds=np.array([[False,True],[True,False]])
        rest=np.array([[0.,1.],[1.,0.]])
        active=np.array([[0.,0,0],[100.,0,0]])
        y,b,stats=mechanical_interval(x,fixed,bonds,rest,active,32)
        self.assertEqual(stats['broken'],1); self.assertFalse(b.any())
        self.assertEqual(float(parent_attachment(b,fixed)[1]),0)

    def test_prepattern_replays_reference_not_new_circuit(self):
        x,_,d=initial(42); n=np.zeros(36)
        a,b=prepattern(x,d.copy(),n.copy(),True)
        for _ in range(400): d,n=regulate(x,d,n,True,.05)
        np.testing.assert_array_equal(a,d); np.testing.assert_array_equal(b,n)

    def test_attachment_and_high_activity_are_separate(self):
        x,fixed,d=initial(42); n=np.zeros(36); graph=contact_graph(x)
        r,c=observe(x,fixed,d,n,graph)
        self.assertEqual(r['material_attachment'],1)
        self.assertEqual(r['mobile_high'],0); self.assertEqual(c,[])


if __name__=='__main__': unittest.main()
