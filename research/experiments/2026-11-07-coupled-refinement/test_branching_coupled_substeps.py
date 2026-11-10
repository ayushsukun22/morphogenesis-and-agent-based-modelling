"""Verify coupled splitting, fixed outer timing and feedback buffers."""
import copy
import json
import unittest

import numpy as np

from branching_baseline import load_simulation,install_growth_arrest
from branching_shared_state import canonical_state,decode,registry,MidpointRandom
from branching_coupled_substeps import (split_update,orient_and_signal,install_substeps,
                                        instrument_outer_updates)


class CoupledSubstepTests(unittest.TestCase):
    def setUp(self):
        self.sim=load_simulation()()
        self.classes=registry(self.sim)
        self.clone=decode(json.loads(canonical_state(self.sim,self.classes)),self.classes)

    def test_one_substep_full_outer_update_matches_original(self):
        for sim in [self.sim,self.clone]:
            install_growth_arrest(sim,0.)
            sim.rng=MidpointRandom()
        install_substeps(self.clone,1)
        for _ in range(10):
            self.sim.update(.05)
            self.clone.update(.05)
            for key in self.sim.cells:
                a,b=self.sim.cells[key],self.clone.cells[key]
                for name in ['position','velocity','polarity']:
                    np.testing.assert_array_equal(getattr(a,name),getattr(b,name))
                self.assertEqual(a.pressure,b.pressure)
                self.assertEqual(vars(a.grn),vars(b.grn))
            np.testing.assert_array_equal(self.sim.oxygen_field.grid,self.clone.oxygen_field.grid)
            np.testing.assert_array_equal(self.sim.morphogen_field.grid,self.clone.morphogen_field.grid)

    def test_orientation_and_signal_held_across_inner_steps(self):
        cells=list(self.sim.cells.values())
        expected=copy.deepcopy(cells)
        mask=orient_and_signal(expected,self.sim.mechanics)
        _,_,actual_mask=split_update(cells,self.sim.mechanics,.05,8)
        self.assertEqual(mask,actual_mask)
        for a,b in zip(cells,expected):
            np.testing.assert_array_equal(a.polarity,b.polarity)
            self.assertEqual(a.grn.neighbor_signal,b.grn.neighbor_signal)

    def test_pressure_buffer_preserves_first_force_not_last(self):
        cells=list(self.sim.cells.values())
        reference=copy.deepcopy(cells)
        self.sim.mechanics.update(reference,.05)
        expected=[c.pressure for c in reference]
        first,last,_=split_update(cells,self.sim.mechanics,.05,8)
        np.testing.assert_array_equal(first,expected)
        np.testing.assert_array_equal([c.pressure for c in cells],first)
        self.assertTrue(any(a!=b for a,b in zip(first,last)))
        np.testing.assert_allclose([np.linalg.norm(c.velocity*self.sim.mechanics.damping) for c in cells],last,rtol=0,atol=1e-14)

    def test_substeps_do_not_add_grn_field_or_orientation_updates(self):
        sim=self.sim
        old_pressure=np.array([c.pressure for c in sim.cells.values()])
        old_cycle=np.array([c.grn.cell_cycle for c in sim.cells.values()])
        radii=np.array([c.radius for c in sim.cells.values()])
        install_growth_arrest(sim,0.)
        sim.rng=MidpointRandom()
        stats=install_substeps(sim,8)
        counts=instrument_outer_updates(sim)
        sim.update(.05)
        self.assertEqual(counts,dict(grn_calls=2,oxygen_steps=1,morphogen_steps=1))
        self.assertEqual(stats['inner_calls'],8)
        self.assertEqual(stats['orientation_cell_evaluations'],2)
        self.assertEqual(sim.rng.calls,2)
        np.testing.assert_array_equal([c.grn.pressure for c in sim.cells.values()],old_pressure)
        np.testing.assert_array_equal([c.radius for c in sim.cells.values()],radii)
        self.assertFalse(np.array_equal([c.grn.cell_cycle for c in sim.cells.values()],old_cycle))

    def test_symmetric_neighbor_and_isolated_guards_match_source(self):
        cell_type=type(next(iter(self.sim.cells.values())))
        cells=[cell_type(cell_id=i,position=np.array(p,dtype=float)) for i,p in
               enumerate([[0,0,0],[-1,0,0],[1,0,0],[20,0,0]])]
        reference=copy.deepcopy(cells)
        self.sim.mechanics.update(reference,.05)
        _,_,mask=split_update(cells,self.sim.mechanics,.05,1)
        self.assertFalse(mask[0])
        self.assertFalse(mask[3])
        for a,b in zip(cells,reference):
            np.testing.assert_array_equal(a.position,b.position)
            np.testing.assert_array_equal(a.polarity,b.polarity)
            self.assertEqual(a.pressure,b.pressure)

    def test_invalid_substep_count_rejected(self):
        for value in [0,-1,1.5]:
            with self.assertRaises(ValueError):
                install_substeps(self.sim,value)


if __name__=='__main__':
    unittest.main()
