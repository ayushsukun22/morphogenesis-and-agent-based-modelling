import numpy as np


class TissueMechanics:
    """
    Experimental cell-cell mechanics.

    Hypothesis:
        Cell adhesion creates local contact geometry.
        Cells polarize relative to that geometry.
        Protrusion acts opposite to the local adhesion direction.

    This is an experimental computational rule, not a validated
    embryological model.
    """

    def __init__(
        self,
        interaction_radius=2.5,
        repulsion_strength=4.0,
        adhesion_strength=0.7,
        protrusion_strength=0.25,
        polarity_memory=0.9,
        damping=2.0,
    ):

        self.interaction_radius = interaction_radius
        self.repulsion_strength = repulsion_strength

        # Kept as fallback/global scaling knobs, but per-cell
        # adhesion_strength / protrusion_strength (driven by each
        # cell's GRN) now take precedence â€” see below.
        self.adhesion_strength = adhesion_strength
        self.protrusion_strength = protrusion_strength

        self.polarity_memory = polarity_memory
        self.damping = damping

    # ==================================================
    # Main update
    # ==================================================

    def update(self, cells, dt):

        if len(cells) < 2:
            return

        forces = {
            cell.cell_id: np.zeros(3, dtype=np.float64)
            for cell in cells
        }

        # ----------------------------------------------
        # Cell-cell interactions
        # ----------------------------------------------

        for i in range(len(cells)):

            cell_a = cells[i]

            for j in range(i + 1, len(cells)):

                cell_b = cells[j]

                delta = (
                    cell_b.position -
                    cell_a.position
                )

                distance = np.linalg.norm(delta)

                if distance < 1e-8:
                    continue

                direction = delta / distance

                contact_distance = (
                    cell_a.radius +
                    cell_b.radius
                )

                if distance > self.interaction_radius:
                    continue

                # ======================================
                # REPULSION
                # ======================================

                overlap = (
                    contact_distance -
                    distance
                )

                if overlap > 0.0:

                    repulsion = (
                        self.repulsion_strength *
                        overlap
                    )

                    forces[cell_a.cell_id] -= (
                        direction * repulsion
                    )

                    forces[cell_b.cell_id] += (
                        direction * repulsion
                    )

                # ======================================
                # ADHESION (now per-cell, from GRN readout)
                # ======================================

                adhesion = (
                    cell_a.adhesion_strength
                    * cell_b.adhesion_strength
                    * min(
                        1.0,
                        self.interaction_radius /
                        max(distance, 1e-8),
                    )
                )

                forces[cell_a.cell_id] += (
                    direction * adhesion
                )

                forces[cell_b.cell_id] -= (
                    direction * adhesion
                )

        # ----------------------------------------------
        # Polarity and protrusion
        # ----------------------------------------------

        for cell in cells:

            neighbors = []

            for other in cells:

                if other.cell_id == cell.cell_id:
                    continue

                delta = (
                    other.position -
                    cell.position
                )

                distance = np.linalg.norm(delta)

                if (
                    distance > 1e-8
                    and
                    distance <= self.interaction_radius
                ):

                    neighbors.append(
                        delta / distance
                    )

            if not neighbors:
                # No neighbors this step â€” zero out pressure so a
                # stray cell's YAP doesn't stay stuck on stale stress.
                cell.pressure = 0.0
                continue

            # ------------------------------------------
            # Mean direction toward neighbors
            # ------------------------------------------

            adhesion_direction = np.sum(
                neighbors,
                axis=0,
            )

            norm = np.linalg.norm(
                adhesion_direction
            )

            if norm < 1e-8:
                cell.pressure = 0.0
                continue

            adhesion_direction /= norm

            # ------------------------------------------
            # Proposed polarity is opposite to adhesion,
            # scaled by the GRN's polarity_bias (SNAIL-driven:
            # progenitors barely polarize, EMT cells strongly do)
            # ------------------------------------------

            proposed_polarity = (
                -adhesion_direction * cell.grn.polarity_bias
            )

            memory = np.clip(
                self.polarity_memory,
                0.0,
                1.0,
            )

            cell.polarity = (
                memory * cell.polarity
                +
                (1.0 - memory)
                * proposed_polarity
            )

            polarity_norm = np.linalg.norm(
                cell.polarity
            )

            if polarity_norm > 1e-8:

                cell.polarity /= polarity_norm

            # ------------------------------------------
            # Protrusion (per-cell strength from GRN)
            # ------------------------------------------

            protrusion = (
                cell.protrusion_strength *
                cell.polarity
            )

            forces[cell.cell_id] += protrusion

        # ----------------------------------------------
        # Integrate motion, and record pressure for YAP
        # ----------------------------------------------

        for cell in cells:

            force = forces[cell.cell_id]

            # Mechanotransduction hook: net force magnitude feeds
            # the GRN's YAP node next step, via embryo.py.
            cell.pressure = float(np.linalg.norm(force))

            # Simple overdamped mechanics.
            velocity = (
                force /
                max(self.damping, 1e-8)
            )

            cell.velocity = velocity

            cell.position += (
                velocity * dt
            )