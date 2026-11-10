"""Verify experiment switches without requiring Omniverse."""
import unittest
from types import SimpleNamespace

import numpy as np

from branching_baseline import install_growth_arrest, install_protrusion_ablation, load_simulation


class InterventionTests(unittest.TestCase):
    def test_growth_arrest_blocks_ready_cells_but_advances_age_and_grn(self):
        sim = load_simulation()()
        for cell in sim.cells.values():
            cell.radius = 0.95  # already ready: growth-rate zero alone is insufficient
        old_ids = set(sim.cells)
        old_radii = [cell.radius for cell in sim.cells.values()]
        old_cycle = [cell.grn.cell_cycle for cell in sim.cells.values()]
        rate = sim.growth_rate
        install_growth_arrest(sim, 0.0)
        sim.update(0.05)
        self.assertEqual(set(sim.cells), old_ids)
        self.assertEqual([c.radius for c in sim.cells.values()], old_radii)
        self.assertTrue(all(c.age == 0.05 for c in sim.cells.values()))
        self.assertNotEqual([c.grn.cell_cycle for c in sim.cells.values()], old_cycle)
        self.assertEqual(sim.growth_rate, rate)
        self.assertTrue(all("ready_for_division" not in c.__dict__ for c in sim.cells.values()))
        self.assertEqual(sim.growth_arrest_steps, 1)

    def test_growth_arrest_starts_on_interval_boundary(self):
        sim = load_simulation()()
        install_growth_arrest(sim, 0.05)
        radii = np.array([c.radius for c in sim.cells.values()])
        sim.update(0.05)
        grown = np.array([c.radius for c in sim.cells.values()])
        self.assertTrue(np.all(grown > radii))
        self.assertEqual(sim.growth_arrest_steps, 0)
        sim.update(0.05)
        np.testing.assert_array_equal([c.radius for c in sim.cells.values()], grown)
        self.assertEqual(sim.growth_arrest_steps, 1)

    def test_growth_arrest_restores_parameters_after_exception(self):
        cell = SimpleNamespace()
        cell.ready_for_division = lambda threshold: True
        original = cell.ready_for_division

        def fail(dt):
            self.assertEqual(sim.growth_rate, 0)
            self.assertFalse(cell.ready_for_division(0))
            raise RuntimeError("expected test exception")

        sim = SimpleNamespace(time=10., growth_rate=0.18, cells={1: cell}, update=fail)
        install_growth_arrest(sim, 10.)
        with self.assertRaises(RuntimeError):
            sim.update(0.05)
        self.assertIs(cell.ready_for_division, original)
        self.assertEqual(sim.growth_rate, 0.18)

    def test_combined_switches_apply_together(self):
        sim = load_simulation()()
        install_protrusion_ablation(sim, 0.)
        install_growth_arrest(sim, 0.)
        radii = [c.radius for c in sim.cells.values()]
        sim.update(0.05)
        self.assertEqual(sim.ablation_steps, 1)
        self.assertEqual(sim.growth_arrest_steps, 1)
        self.assertEqual([c.radius for c in sim.cells.values()], radii)
        self.assertTrue(all(c.protrusion_strength == c.grn.protrusion for c in sim.cells.values()))


if __name__ == "__main__":
    unittest.main()
