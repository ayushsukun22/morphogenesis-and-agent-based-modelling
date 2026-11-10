"""Check complete-state checkpoint fidelity and deterministic noise control."""
import json
import unittest

import numpy as np

from branching_baseline import load_simulation
from branching_shared_state import (canonical_state, decode, encode, registry, MidpointRandom,
                                    scaled_memory, install_time_scaled_polarity)
from types import SimpleNamespace


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


if __name__ == "__main__":
    unittest.main()
