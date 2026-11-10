import unittest
import numpy as np
from branching_prerequisites import signaling_step,mechanical_force_energy,attachment
from branching_sprouting_assay import initial,contact_graph
from branching_adhesion_assay import forces_energy


class PrerequisiteTests(unittest.TestCase):
    def test_reference_force_matches(self):
        x,_,_=initial(42)
        f,u=mechanical_force_energy(x,4,2.5)
        g,v=forces_energy(x,np.zeros(36),1,1)
        # Decimal .49 versus the original .7**2 differs at floating roundoff.
        # This diagnostic promises the same force law, not bitwise trajectory replay.
        np.testing.assert_allclose(f,g,rtol=1e-13,atol=2e-14)
        self.assertAlmostEqual(u,v,places=11)

    def test_short_range_removes_noncontact_attraction(self):
        x=np.array([[0.,0,0],[1.8,0,0]])
        f,_=mechanical_force_energy(x,40,1.2)
        np.testing.assert_array_equal(f,0)

    def test_material_springs_rest_force_and_gradient(self):
        x=np.array([[0.,0,0],[1.05,0,0]])
        rest=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2); bonds=contact_graph(x)
        f,_=mechanical_force_energy(x,40,1.2,bonds,rest)
        np.testing.assert_array_equal(f,0)
        x[1,0]+=.1
        f,_=mechanical_force_energy(x,40,1.2,bonds,rest)
        self.assertAlmostEqual(f[0,0],.4)
        plus=x.copy(); minus=x.copy(); plus[1,0]+=1e-6; minus[1,0]-=1e-6
        gradient=(mechanical_force_energy(plus,40,1.2,bonds,rest)[1]-mechanical_force_energy(minus,40,1.2,bonds,rest)[1])/2e-6
        self.assertAlmostEqual(f[1,0],-gradient,places=7)

    def test_inhibition_and_neighbor_locality(self):
        graph=np.array([[False,True,False],[True,False,False],[False,False,False]])
        d=np.array([.9,.1,.5]); n=np.array([0.,1.,0.])
        new_d,new_n=signaling_step(d,n,graph,.7,.05)
        self.assertLess(new_d[1],d[1]); self.assertEqual(new_n[2],0)
        a,_=signaling_step(d,n,graph,.7,.05,inhibition=False)
        self.assertGreater(a[1],d[1])

    def test_attachment_respects_gap(self):
        x=np.array([[0.,0,0],[1.,0,0],[4.,0,0]])
        self.assertEqual(attachment(x,np.array([True,False,False]),1.2),.5)


if __name__=='__main__': unittest.main()
