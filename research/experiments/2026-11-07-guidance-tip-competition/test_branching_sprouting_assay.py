import unittest
import numpy as np
from branching_sprouting_assay import initial,cue,regulate,traction,readout
from branching_adhesion_assay import forces_energy


class SproutingTests(unittest.TestCase):
    def test_cue_gradient(self):
        x,_,_=initial(42); v,g=cue(x)
        plus=x.copy(); minus=x.copy(); plus[:,0]+=1e-5; minus[:,0]-=1e-5
        np.testing.assert_allclose(g[:,0],(cue(plus)[0]-cue(minus)[0])/2e-5,atol=1e-10)
        self.assertTrue((g[:,0]>0).all())

    def test_notch_suppresses_own_tip_activity(self):
        x,_,d=initial(42)
        low,_=regulate(x,d,np.zeros(36),True,.05)
        high,_=regulate(x,d,np.ones(36),True,.05)
        self.assertTrue((high<low).all())
        a,_=regulate(x,d,np.zeros(36),False,.05)
        b,_=regulate(x,d,np.ones(36),False,.05)
        np.testing.assert_array_equal(a,b)

    def test_guidance_changes_direction_not_strength(self):
        x,_,d=initial(42)
        a=traction(x,d,True); b=traction(x,d,False)
        np.testing.assert_allclose(np.linalg.norm(a,axis=1),np.linalg.norm(b,axis=1))
        self.assertTrue((a[:,0]>0).all()); np.testing.assert_array_equal(a[:,1:],0)

    def test_anchors_and_active_work(self):
        x,fixed,d=initial(42); old=x.copy(); active=traction(x,d,True); active[fixed]=0
        f,u=forces_energy(x,np.zeros(36),1,1); w=u-np.sum(active*x)
        total=f+active; total[fixed]=0; x+=.0001*total/2
        _,after=forces_energy(x,np.zeros(36),1,1)
        self.assertLess(after-np.sum(active*x),w)
        np.testing.assert_array_equal(x[fixed],old[fixed])

    def test_detached_cells_not_counted_as_sprout(self):
        x=np.array([[0.,0,0],[4.,0,0],[4.5,0,0],[5.,0,0]])
        r,c=readout(x,np.array([True,False,False,False]),np.ones(4))
        self.assertEqual(r['detached_mobile_cells'],3); self.assertEqual(c,[])


if __name__=='__main__': unittest.main()
