import numpy as np


class TissueMechanics:
    """
    Experimental cell-cell mechanics.

    Updated to apply global physics multipliers (adhesion/protrusion)
    and track scalar pressure to fix the vector cancellation bug.
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

        forces = {cell.cell_id: np.zeros(3, dtype=np.float64) for cell in cells}

        # Track scalar pressure so core cells realize they are being squished
        scalar_pressures = {cell.cell_id: 0.0 for cell in cells}

        # ----------------------------------------------
        # Cell-cell interactions
        # ----------------------------------------------

        for i in range(len(cells)):
            cell_a = cells[i]

            for j in range(i + 1, len(cells)):
                cell_b = cells[j]

                delta = cell_b.position - cell_a.position
                distance = np.linalg.norm(delta)

                if distance < 1e-8 or distance > self.interaction_radius:
                    continue

                direction = delta / distance
                contact_distance = cell_a.radius + cell_b.radius
                overlap = contact_distance - distance

                # ======================================
                # REPULSION
                # ======================================
                if overlap > 0.0:
                    repulsion = self.repulsion_strength * overlap
                    repulsion_vector = direction * repulsion

                    forces[cell_a.cell_id] -= repulsion_vector
                    forces[cell_b.cell_id] += repulsion_vector

                    # Accumulate scalar magnitude for YAP mechanosensing
                    scalar_pressures[cell_a.cell_id] += repulsion
                    scalar_pressures[cell_b.cell_id] += repulsion

                # ======================================
                # ADHESION
                # ======================================
                # Multiplied by self.adhesion_strength to prevent GRN values
                # from being completely overpowered by protrusion.
                adhesion = (
                    self.adhesion_strength
                    * cell_a.adhesion_strength
                    * cell_b.adhesion_strength
                    * min(1.0, self.interaction_radius / max(distance, 1e-8))
                )

                forces[cell_a.cell_id] += direction * adhesion
                forces[cell_b.cell_id] -= direction * adhesion

        # ----------------------------------------------
        # Polarity and protrusion
        # ----------------------------------------------

        for cell in cells:
            neighbors = []
            for other in cells:
                if other.cell_id == cell.cell_id:
                    continue

                delta = other.position - cell.position
                distance = np.linalg.norm(delta)

                if 1e-8 < distance <= self.interaction_radius:
                    neighbors.append(delta / distance)

            if not neighbors:
                cell.pressure = 0.0
                continue

            # ------------------------------------------
            # Mean direction toward neighbors
            # ------------------------------------------
            adhesion_direction = np.sum(neighbors, axis=0)
            norm = np.linalg.norm(adhesion_direction)

            if norm < 1e-8:
                # Still pass scalar pressure even if perfectly surrounded
                cell.pressure = scalar_pressures[cell.cell_id]
                continue

            adhesion_direction /= norm

            # ------------------------------------------
            # Proposed polarity is opposite to adhesion
            # ------------------------------------------
            proposed_polarity = -adhesion_direction * cell.grn.polarity_bias
            memory = np.clip(self.polarity_memory, 0.0, 1.0)

            cell.polarity = (
                memory * cell.polarity
                + (1.0 - memory) * proposed_polarity
            )

            polarity_norm = np.linalg.norm(cell.polarity)
            if polarity_norm > 1e-8:
                cell.polarity /= polarity_norm

            # ------------------------------------------
            # Protrusion
            # ------------------------------------------
            # Multiplied by self.protrusion_strength so core cells can tear out
            protrusion = (
                self.protrusion_strength
                * cell.protrusion_strength
                * cell.polarity
            )

            forces[cell.cell_id] += protrusion

        # ----------------------------------------------
        # Integrate motion, and record pressure for YAP
        # ----------------------------------------------

        for cell in cells:
            force = forces[cell.cell_id]

            # Assign scalar sum of compressive forces for YAP sensing
            cell.pressure = float(scalar_pressures[cell.cell_id])

            velocity = force / max(self.damping, 1e-8)
            cell.velocity = velocity
            cell.position += velocity * dt