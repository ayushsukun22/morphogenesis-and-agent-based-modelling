import unittest
import numpy as np
from branching_adhesion_assay import forces_energy,initial
from branching_baseline import load_simulation
from branching_passive_refinement import pair_step


class AdhesionAssayTests(unittest.TestCase):
    def test_equal_matches_original_passive_forces(self):
        sim=load_simulation()(); prototype=next(iter(sim.cells.values()))
        cls=type(prototype)
        x,_=initial(42); x=x[:8].copy()
        cells=[cls(cell_id=i,position=p.copy(),radius=.5,adhesion_strength=.7) for i,p in enumerate(x)]
        expected,_=forces_energy(x,np.arange(8)%2,1,1)
        pair_step(cells,sim.mechanics,.001)
        np.testing.assert_allclose(np.array([c.position for c in cells]),x+.001*expected/2,atol=1e-14,rtol=1e-14)

    def test_cross_factor_only_changes_unlike_pair(self):
        x=np.array([[0.,0,0],[1.5,0,0]])
        full,_=forces_energy(x,np.array([0,1]),1,1)
        weak,_=forces_energy(x,np.array([0,1]),1,.2)
        same,_=forces_energy(x,np.array([0,0]),1,.2)
        np.testing.assert_allclose(weak,.2*full)
        np.testing.assert_array_equal(same,full)
        np.testing.assert_allclose(full.sum(axis=0),0,atol=1e-15)

    def test_force_is_negative_potential_gradient(self):
        x=np.array([[0.,0,0],[.8,.1,0],[1.8,.3,.2]])
        labels=np.array([0,1,0]); force,_=forces_energy(x,labels,1,.2)
        gradient=np.zeros_like(x)
        for i in range(3):
            for j in range(3):
                plus=x.copy(); minus=x.copy(); plus[i,j]+=1e-6; minus[i,j]-=1e-6
                gradient[i,j]=(forces_energy(plus,labels,1,.2)[1]-forces_energy(minus,labels,1,.2)[1])/2e-6
        np.testing.assert_allclose(force,-gradient,atol=1e-8)

    def test_initial_is_balanced_reproducible(self):
        a,l=initial(42); b,m=initial(42)
        np.testing.assert_array_equal(a,b); np.testing.assert_array_equal(l,m)
        self.assertEqual(int(l.sum()),32)


if __name__=='__main__':
    unittest.main()
