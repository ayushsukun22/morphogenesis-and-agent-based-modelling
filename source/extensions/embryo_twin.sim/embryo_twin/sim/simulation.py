"""
simulation.py

First prototype of the embryonic digital twin.

This is an engineering prototype, NOT a validated
embryological model.
"""

import numpy as np
from .mechanics import TissueMechanics

from .cell import Cell
from .chemical_field import ChemicalField


class EmbryoSimulation:
    MAX_CELLS = 1500

    def __init__(
        self,
        bounds_min=(-10.0, -10.0, -10.0),
        bounds_max=(10.0, 10.0, 10.0),

    ):

        self.mechanics = TissueMechanics(
        interaction_radius=2.5,
        repulsion_strength=4.0,
        adhesion_strength=1.8,
        protrusion_strength=0.12,
        polarity_memory=0.2,
        damping=2.0
                    )
        self.time = 0.0

        self.next_cell_id = 2

        self.cells = {}
        self.rng = np.random.default_rng(42)

        self.bounds_min = bounds_min
        self.bounds_max = bounds_max

        # -------------------------
        # Chemical fields
        # -------------------------

        self.oxygen_field = ChemicalField(
            bounds_min=bounds_min,
            bounds_max=bounds_max,

            diffusion_rate=0.8,
            decay_rate=0.0,      # only depleted by consumption
            baseline=1.0,
        )

        self.morphogen_field = ChemicalField(
            bounds_min=bounds_min,
            bounds_max=bounds_max,

            diffusion_rate=0.3,
            decay_rate=0.05,
            baseline=0.0,
        )

        # -------------------------
        # Simulation parameters
        # -------------------------

        self.growth_rate = 0.18
        self.division_radius = 0.9

        self.oxygen_consumption_rate = 0.05
        self.morphogen_production_rate = 0.02

        self.oxygen_division_threshold = 0.2   # below this, division halts
        self.fate_commit_threshold = 0.6       # min time in a fate before it can flip

        # -------------------------
        # Seed cell
        # -------------------------

        self._spawn_seed_cell()

    # ======================================================
    # Setup / reset
    # ======================================================

    def _spawn_seed_cell(self, count=2):
        self.cells = {}
        for i in range(count):
            position = self.rng.normal(
                loc=0.0,
                scale=0.65,
                size=3, )

            cell = Cell(
                cell_id=i + 1,
                position=position,
                radius=0.5,
                generation=0,
                lineage_id=str(i + 1),
            )
            cell.sync_phenotype()
            self.cells[cell.cell_id] = cell
            self.next_cell_id = count + 1

    def reset(self):
        self.time = 0.0
        self.oxygen_field.grid[:] = 1.0
        self.morphogen_field.grid[:] = 0.0
        self.rng = np.random.default_rng(42)
        self._spawn_seed_cell()

    # ======================================================
    # Queries
    # ======================================================

    def alive_cells(self):
        return [cell for cell in self.cells.values() if cell.alive]


    def log_metrics(self, tag=""):
        cells = self.alive_cells()
        pos = np.array([c.position for c in cells])
        snail = np.array([c.grn.snail for c in cells])
        # newborn daughters have not felt a mechanics step yet,
        # so leave them out of the pressure statistics
        settled = [c for c in cells if c.age > 0.3] or cells
        press = np.array([c.pressure for c in settled])

        centre = pos.mean(axis=0)
        dist = np.linalg.norm(pos - centre, axis=1)

        eig = np.sort(np.linalg.eigvalsh(np.cov((pos - centre).T)))
        roundness = float(np.sqrt(eig[-1] / max(eig[0], 1e-9)))

        tips = snail > 0.8
        tip_share = float(tips.mean())
        tip_reach = (
            float(dist[tips].mean() / dist.mean()) if tips.any() else float("nan")
        )
        p10, p50, p90 = np.percentile(press, [10, 50, 90])

        print(
            f"[METRICS{tag}] t={self.time:.1f}h cells={len(cells)} "
            f"tip_share={tip_share:.2f} "
            f"pressure p10/p50/p90={p10:.1f}/{p50:.1f}/{p90:.1f} "
            f"roundness={roundness:.2f} tip_reach={tip_reach:.2f}",
            flush=True,
        )


    # ======================================================
    # Main update loop
    # ======================================================

    def update(self, dt):

        self.time += dt

        new_cells = []
        pending_divisions = 0

        current_cells = list(self.alive_cells())

        for cell in current_cells:

            # -------------------------
            # Environment
            # -------------------------

            cell.oxygen = self.oxygen_field.sample(
                cell.position
            )

            cell.morphogen_a = self.morphogen_field.sample(
                cell.position
            )

            # -------------------------
            # GRN update (reads last step's pressure, which
            # TissueMechanics.update() wrote at the end of the
            # previous frame)
            # -------------------------

            cell.grn.oxygen = cell.oxygen
            cell.grn.morphogen = cell.morphogen_a
            cell.grn.pressure = cell.pressure
            cell.grn.update(dt, self.rng)
            cell.sync_phenotype()

            # -------------------------
            # Growth
            # -------------------------

            if cell.oxygen > self.oxygen_division_threshold:
                cell.grow(
                    dt,
                    self.growth_rate
                )

            # -------------------------
            # Consumption / production
            # -------------------------

            self.oxygen_field.consume(
                cell.position,
                amount=(
                    self.oxygen_consumption_rate
                    * cell.radius
                    * dt
                ),
            )

            self.morphogen_field.deposit(
                cell.position,
                amount=(
                    self.morphogen_production_rate
                    * dt
                ),
            )

            # -------------------------
            # Fate / GRN
            # -------------------------

            self.update_fate(
                cell,
                dt
            )

            # -------------------------
            # Division
            # -------------------------

            if (
                cell.ready_for_division(
                    self.division_radius
                )
                and
                cell.oxygen >
                self.oxygen_division_threshold
                and len(self.cells) + pending_divisions < self.MAX_CELLS
            ):

                daughter_a, daughter_b = (
                    cell.divide(
                        self.next_cell_id,
                        self.rng,
                    )
                )

                self.next_cell_id += 2
                pending_divisions += 1

                new_cells.append(
                    daughter_a
                )

                new_cells.append(
                    daughter_b
                )

        # ==================================================
        # Mechanics (also writes cell.pressure for next step's GRN)
        # ==================================================

        self.mechanics.update(
            current_cells,
            dt,
        )

        # ==================================================
        # Add newly created daughters
        # ==================================================

        for cell in new_cells:
            self.cells[cell.cell_id] = cell

        # ==================================================
        # Chemical fields
        # ==================================================

        self.oxygen_field.step(dt)

        self.morphogen_field.step(dt)

        # ==================================================
        # Cleanup
        # ==================================================

        self._prune_dead_cells()
        if int(self.time/5) != int((self.time-dt)/5):
            self.log_metrics()

    # ======================================================
    # Fate model
    # ======================================================

    def update_fate(self, cell, dt):
        """
        Placeholder differentiation model with simple commitment:
        a cell must have spent `fate_commit_threshold` time in its
        current fate before it's allowed to switch, to avoid
        frame-to-frame flicker near threshold boundaries.
        """

        cell.time_in_fate = getattr(cell, "time_in_fate", 0.0) + dt

        if cell.time_in_fate < self.fate_commit_threshold:
            return

        if cell.morphogen_a > 0.75:
            new_fate = "TYPE_A"
        elif cell.oxygen < 0.30:
            new_fate = "TYPE_B"
        else:
            new_fate = "PROGENITOR"

        if new_fate != cell.fate:
            cell.fate = new_fate
            cell.time_in_fate = 0.0

    # ======================================================
    # Housekeeping
    # ======================================================

    def _prune_dead_cells(self):
        """Drop dead (post-division) cells from the registry so
        the dict doesn't grow unboundedly. Safe because `divide()`
        already spawned replacement daughters before marking the
        parent dead."""



        dead_ids = [
            cell_id for cell_id, cell in self.cells.items()
            if not cell.alive
        ]

        for cell_id in dead_ids:
            del self.cells[cell_id]
