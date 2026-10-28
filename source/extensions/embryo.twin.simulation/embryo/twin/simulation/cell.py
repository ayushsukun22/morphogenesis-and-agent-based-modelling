from dataclasses import dataclass, field
from typing import Optional

import copy

import numpy as np

from .grn import GeneRegulatoryNetwork


@dataclass
class Cell:
    """Cell state for the embryonic digital-twin prototype."""

    cell_id: int
    position: np.ndarray

    radius: float = 0.5
    age: float = 0.0

    generation: int = 0
    parent_id: Optional[int] = None

    cell_cycle: str = "G1"
    fate: str = "PROGENITOR"

    oxygen: float = 1.0
    morphogen_a: float = 0.0

    lineage_id: str = ""

    alive: bool = True

    # ----------------------------------------
    # Gene-regulatory state
    # ----------------------------------------

    grn: GeneRegulatoryNetwork = field(default_factory=GeneRegulatoryNetwork)

    # ----------------------------------------
    # Mechanical state
    # ----------------------------------------

    velocity: np.ndarray = field(
        default_factory=lambda: np.zeros(3, dtype=np.float64)
    )

    polarity: np.ndarray = field(
        default_factory=lambda: np.array(
            [1.0, 0.0, 0.0],
            dtype=np.float64,
        )
    )

    adhesion_strength: float =   0.4 #0.7
    protrusion_strength: float =  0.8 #0.25

    pressure: float = 0.0
    time_in_fate: float = 0.0

    def __post_init__(self):

        self.position = np.asarray(
            self.position,
            dtype=np.float64,
        )

        self.velocity = np.asarray(
            self.velocity,
            dtype=np.float64,
        )

        self.polarity = np.asarray(
            self.polarity,
            dtype=np.float64,
        )

        if not self.lineage_id:
            self.lineage_id = str(self.cell_id)

    # ========================================
    # GRN sync
    # ========================================

    def sync_phenotype(self):
        """Copy GRN readout into the mechanical fields TissueMechanics reads."""
        self.adhesion_strength = self.grn.adhesion
        self.protrusion_strength = self.grn.protrusion

    # ========================================
    # Growth
    # ========================================

    def grow(
        self,
        dt: float,
        growth_rate: float,
    ):
        """Grow the cell and advance its biological age."""

        self.radius += growth_rate * dt * self.grn.division_rate
        self.age += dt

    # ========================================
    # Division
    # ========================================

    def ready_for_division(
        self,
        division_radius: float,
    ) -> bool:

        return self.radius >= division_radius

    def divide(
        self,
        daughter_id_start: int,
    ):
        """
        Divide one cell into two daughter cells.
        """

                # Division axis: mostly along the cell's current polarity
        # (so a protrusive/tip cell extends a filament instead of
        # scattering daughters randomly), with a little noise so
        # stalk cells (near-zero polarity_bias) still spread out.
        jitter = np.random.normal(0.0, 1.0, size=3)
        jitter_norm = np.linalg.norm(jitter)
        if jitter_norm > 1e-8:
            jitter /= jitter_norm

        bias = self.grn.polarity_bias   # 0 (stalk) .. 1 (tip)
        direction = bias * self.polarity + (1.0 - bias) * jitter

        length = np.linalg.norm(direction)
        if length < 1e-8:
            direction = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        else:
            direction /= length

        # Equal-volume approximation.

        daughter_radius = (
            self.radius /
            (2.0 ** (1.0 / 3.0))
        )

        # Separate daughters.

        offset = (
            direction
            * daughter_radius
            * 1.15
        )

        position_a = (
            self.position + offset
        )

        position_b = (
            self.position - offset
        )

        daughter_a = Cell(
            cell_id=daughter_id_start,
            position=position_a,
            radius=daughter_radius,
            generation=self.generation + 1,
            parent_id=self.cell_id,
            oxygen=self.oxygen,
            morphogen_a=self.morphogen_a,
            fate=self.fate,
            lineage_id=f"{self.lineage_id}.A",

            grn=copy.deepcopy(self.grn),

            adhesion_strength=self.adhesion_strength,
            protrusion_strength=self.protrusion_strength,
        )

        daughter_b = Cell(
            cell_id=daughter_id_start + 1,
            position=position_b,
            radius=daughter_radius,
            generation=self.generation + 1,
            parent_id=self.cell_id,
            oxygen=self.oxygen,
            morphogen_a=self.morphogen_a,
            fate=self.fate,
            lineage_id=f"{self.lineage_id}.B",

            grn=copy.deepcopy(self.grn),

            adhesion_strength=self.adhesion_strength,
            protrusion_strength=self.protrusion_strength,
        )

        self.alive = False

        return daughter_a, daughter_b