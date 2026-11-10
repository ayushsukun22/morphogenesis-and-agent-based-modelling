"""Check complete-state checkpoint fidelity and deterministic noise control."""
import json
import unittest
import copy

import numpy as np

from branching_baseline import load_simulation
from branching_shared_state import (canonical_state, decode, encode, registry, MidpointRandom,
                                    scaled_memory, install_time_scaled_polarity,
                                    install_frozen_regulation, regulatory_arrays)
from branching_baseline import install_growth_arrest, install_protrusion_ablation
from types import SimpleNamespace
from branching_passive_report import pair_potential
from branching_passive_refinement import pair_step
from branching_refinement_report import current_forces


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.sim = load_simulation()()
        self.sim.update(.05)
        self.classes = registry(self.sim)
        self.payload = canonical_state(self.sim, self.classes)
        self.clone = decode(json.loads(self.payload), self.classes)

    def test_complete_roundtrip_and_integer_cell_keys(self):
        self.assertEqual(canonical_state(self.clone, self.classes), self.payload)
        self.assertTrue(all(isinstance(k, int) for k in self.clone.cells))
        self.assertEqual(set(vars(self.clone)), set(vars(self.sim)))

    def test_future_stochastic_update_matches_exactly(self):
        self.sim.update(.05)
        self.clone.update(.05)
        self.assertEqual(canonical_state(self.clone, self.classes), canonical_state(self.sim, self.classes))

    def test_clone_mutations_do_not_affect_original(self):
        cell_id = next(iter(self.clone.cells))
        before = self.sim.cells[cell_id].position.copy()
        self.clone.cells[cell_id].position[0] += 2
        self.clone.cells[cell_id].grn.snail = 1.5
        self.clone.oxygen_field.grid[:] = 0
        self.clone.morphogen_field.grid[:] = 8
        np.testing.assert_array_equal(self.sim.cells[cell_id].position, before)
        self.assertNotEqual(self.sim.cells[cell_id].grn.snail, 1.5)
        self.assertTrue(np.any(self.sim.oxygen_field.grid > 0))
        self.assertFalse(np.all(self.sim.morphogen_field.grid == 8))

    def test_midpoint_cancels_only_existing_grn_noise_term(self):
        noise = MidpointRandom()
        self.assertEqual(.05 * (noise.random() - .5), 0)
        self.assertEqual(noise.calls, 1)
        with self.assertRaises(AttributeError):
            noise.normal()

    def test_unplanned_checkpoint_types_rejected(self):
        with self.assertRaises(TypeError):
            encode(lambda: None, self.classes)

    def test_scalar_memory_retention_composes(self):
        self.assertEqual(scaled_memory(.2, .05, .05), .2)
        self.assertAlmostEqual(scaled_memory(.2, .025, .05)**2, .2)
        self.assertAlmostEqual(scaled_memory(.2, .0125, .05)**4, .2)
        self.assertEqual(scaled_memory(0., .025, .05), 0.)
        self.assertEqual(scaled_memory(1., .025, .05), 1.)
        with self.assertRaises(ValueError):
            scaled_memory(.2, 0., .05)

    def test_reference_step_reproduces_original_update(self):
        install_time_scaled_polarity(self.clone, .05)
        for _ in range(10):
            self.sim.update(.05)
            self.clone.update(.05)
        for key in self.sim.cells:
            a, b = self.sim.cells[key], self.clone.cells[key]
            for field in ["position", "polarity", "velocity"]:
                np.testing.assert_array_equal(getattr(a, field), getattr(b, field))
            self.assertEqual(vars(a.grn), vars(b.grn))
        np.testing.assert_array_equal(self.sim.oxygen_field.grid,self.clone.oxygen_field.grid)
        self.assertEqual(self.clone.polarity_scaled_steps,10)

    def test_scaled_memory_restored_after_failure(self):
        def fail(cells, dt):
            self.assertAlmostEqual(sim.mechanics.polarity_memory, np.sqrt(.2))
            raise RuntimeError("expected test exception")
        sim = SimpleNamespace(mechanics=SimpleNamespace(polarity_memory=.2,update=fail))
        install_time_scaled_polarity(sim,.05)
        with self.assertRaises(RuntimeError):
            sim.mechanics.update([], .025)
        self.assertEqual(sim.mechanics.polarity_memory,.2)

    def test_frozen_regulation_retains_nodes_while_motion_and_fields_continue(self):
        sim = self.sim
        cells = list(sim.cells.values())
        nodes, outputs = regulatory_arrays(cells)
        positions = np.array([c.position.copy() for c in cells])
        oxygen = sim.oxygen_field.grid.copy()
        install_growth_arrest(sim,0.)
        install_frozen_regulation(sim)
        sim.rng = MidpointRandom()
        for _ in range(5):
            sim.update(.05)
        new_nodes,new_outputs = regulatory_arrays(cells)
        np.testing.assert_array_equal(nodes,new_nodes)
        np.testing.assert_array_equal(outputs,new_outputs)
        self.assertFalse(np.array_equal(positions,[c.position for c in cells]))
        self.assertFalse(np.array_equal(oxygen,sim.oxygen_field.grid))
        self.assertEqual(sim.rng.calls,0)
        self.assertEqual(sim.frozen_regulation_calls,5*len(cells))

    def test_frozen_regulation_preserves_cell_specific_strengths(self):
        cells = list(self.sim.cells.values())
        cells[0].grn.adhesion=.2
        cells[1].grn.adhesion=.8
        for cell in cells:
            cell.sync_phenotype()
        install_growth_arrest(self.sim,0.)
        install_frozen_regulation(self.sim)
        self.sim.update(.05)
        self.assertEqual([c.adhesion_strength for c in cells],[.2,.8])

    def test_frozen_nodes_ignore_updated_pressure_input(self):
        cell = next(iter(self.sim.cells.values()))
        yap = cell.grn.yap
        cell.pressure=7.
        install_growth_arrest(self.sim,0.)
        install_frozen_regulation(self.sim)
        self.sim.update(.05)
        self.assertEqual(cell.grn.pressure,7.)
        self.assertEqual(cell.grn.yap,yap)

    def test_two_cell_passive_repulsion_matches_expected_motion(self):
        sim=self.sim
        cells=list(sim.cells.values())
        for i,cell in enumerate(cells):
            cell.position=np.array([-.4 if i==0 else .4,0.,0.])
            cell.radius=.5
            cell.grn.adhesion=0.
            cell.sync_phenotype()
        install_growth_arrest(sim,0.)
        install_frozen_regulation(sim)
        install_protrusion_ablation(sim,0.)
        sim.update(.05)
        # overlap=.2; repulsion=4*.2=.8; dx=.05*.8/2=.02.
        np.testing.assert_allclose([c.position for c in cells],[[-.42,0,0],[.42,0,0]],rtol=0,atol=1e-14)
        np.testing.assert_allclose(np.mean([c.position for c in cells],axis=0),0,atol=1e-14)
        self.assertTrue(all(c.protrusion_strength==c.grn.protrusion for c in cells))

    def test_two_cell_passive_adhesion_matches_expected_motion(self):
        sim=self.sim
        cells=list(sim.cells.values())
        for i,cell in enumerate(cells):
            cell.position=np.array([-.75 if i==0 else .75,0.,0.])
            cell.radius=.5
            cell.grn.adhesion=.5 if i==0 else .8
            cell.sync_phenotype()
        install_growth_arrest(sim,0.)
        install_frozen_regulation(sim)
        install_protrusion_ablation(sim,0.)
        sim.update(.05)
        # No overlap. Adhesion=.5*.8=.4; dx=.05*.4/2=.01.
        np.testing.assert_allclose([c.position for c in cells],[[-.74,0,0],[.74,0,0]],rtol=0,atol=1e-14)
        np.testing.assert_allclose(np.mean([c.position for c in cells],axis=0),0,atol=1e-14)

    def test_polarity_memory_does_not_change_passive_translation(self):
        a,b=self.sim,self.clone
        for sim,memory in [(a,0.),(b,1.)]:
            sim.mechanics.polarity_memory=memory
            for cell in sim.cells.values():
                cell.polarity=np.array([1.,0.,0.])
            install_growth_arrest(sim,0.)
            install_frozen_regulation(sim)
            install_protrusion_ablation(sim,0.)
        for _ in range(10):
            a.update(.05)
            b.update(.05)
        for key in a.cells:
            np.testing.assert_array_equal(a.cells[key].position,b.cells[key].position)
            self.assertEqual(a.cells[key].pressure,b.cells[key].pressure)
        self.assertTrue(any(not np.array_equal(a.cells[k].polarity,b.cells[k].polarity) for k in a.cells))

    def test_pair_potential_gradient_matches_existing_mechanical_force(self):
        positions=np.array([[0.,0,0],[.8,.2,0],[1.7,-.1,.3]])
        radii=np.array([.5,.6,.55])
        adhesion=np.array([.3,.5,.8])
        cells=[SimpleNamespace(cell_id=i,position=p.copy(),radius=r,adhesion_strength=a,
                               protrusion_strength=0.,polarity=np.array([1.,0,0]),
                               velocity=np.zeros(3),pressure=0.,
                               grn=SimpleNamespace(snail=.1,polarity_bias=.5))
               for i,(p,r,a) in enumerate(zip(positions,radii,adhesion))]
        mechanics=self.sim.mechanics
        gradient=np.zeros_like(positions)
        epsilon=1e-6
        for i in range(3):
            for j in range(3):
                plus,minus=positions.copy(),positions.copy()
                plus[i,j]+=epsilon
                minus[i,j]-=epsilon
                args=(radii,adhesion,mechanics.interaction_radius,mechanics.repulsion_strength)
                gradient[i,j]=(pair_potential(plus,*args)[0]-pair_potential(minus,*args)[0])/(2*epsilon)
        mechanics.update(cells,.001)
        force=np.array([c.velocity for c in cells])*mechanics.damping
        np.testing.assert_allclose(force,-gradient,rtol=0,atol=1e-7)

    def test_pair_potential_continuous_at_cutoff(self):
        # Also cover contact distance greater than interaction cutoff.
        for radii in [np.array([.5,.5]),np.array([2.,2.])]:
            args=(radii,np.array([.5,.8]),2.5,4.)
            at,near=pair_potential(np.array([[0.,0,0],[2.5,0,0]]),*args)
            inside,_=pair_potential(np.array([[0.,0,0],[2.5-1e-8,0,0]]),*args)
            outside,_=pair_potential(np.array([[0.,0,0],[2.5+1e-8,0,0]]),*args)
            self.assertEqual(at,0.)
            self.assertEqual(outside,0.)
            self.assertLess(abs(inside),1e-6)
            self.assertEqual(near,0)

    def test_reduced_pair_step_exactly_matches_original_across_updates(self):
        rng=np.random.default_rng(123)
        cell_type=type(next(iter(self.sim.cells.values())))
        template=[cell_type(cell_id=i,position=rng.normal(0,1,3),radius=float(rng.uniform(.3,.8))) for i in range(8)]
        for c in template:
            c.adhesion_strength=float(rng.uniform(.2,.9))
            c.protrusion_strength=0.
            c.grn.polarity_bias=float(rng.uniform(0,1))
        for dt in [.05,.0125,.003125]:
            a,b=copy.deepcopy(template),copy.deepcopy(template)
            for _ in range(20):
                self.sim.mechanics.update(a,dt)
                pair_step(b,self.sim.mechanics,dt)
                for x,y in zip(a,b):
                    np.testing.assert_array_equal(x.position,y.position)
                    np.testing.assert_array_equal(x.velocity,y.velocity)
                    self.assertEqual(x.pressure,y.pressure)

    def test_reduced_pair_step_leaves_nonmechanical_state_unchanged(self):
        sim=self.sim
        cells=list(sim.cells.values())
        before_nodes,before_outputs=regulatory_arrays(cells)
        oxygen=sim.oxygen_field.grid.copy()
        morphogen=sim.morphogen_field.grid.copy()
        age=[c.age for c in cells]
        center=np.mean([c.position for c in cells],axis=0)
        for _ in range(20):
            pair_step(cells,sim.mechanics,.0125)
        nodes,outputs=regulatory_arrays(cells)
        np.testing.assert_array_equal(nodes,before_nodes)
        np.testing.assert_array_equal(outputs,before_outputs)
        np.testing.assert_array_equal(sim.oxygen_field.grid,oxygen)
        np.testing.assert_array_equal(sim.morphogen_field.grid,morphogen)
        self.assertEqual([c.age for c in cells],age)
        np.testing.assert_allclose(np.mean([c.position for c in cells],axis=0),center,rtol=0,atol=1e-14)

    def test_common_time_force_probe_matches_source_without_advancing_positions(self):
        cells=copy.deepcopy(list(self.sim.cells.values()))
        for cell in cells:
            cell.protrusion_strength=0.
        positions=np.array([c.position for c in cells])
        frame=dict(ids=np.array([c.cell_id for c in cells]),positions=positions.copy(),
                   radii=np.array([c.radius for c in cells]),
                   mechanical_strengths=np.array([[c.adhesion_strength,0.] for c in cells]))
        settings=dict(force_order_ids=[c.cell_id for c in cells])
        force=current_forces(frame,settings,self.sim.mechanics)
        self.sim.mechanics.update(cells,0.)
        np.testing.assert_array_equal(force,np.array([c.velocity for c in cells])*self.sim.mechanics.damping)
        np.testing.assert_array_equal(frame['positions'],positions)


if __name__ == "__main__":
    unittest.main()
