"""
grn.py

Per-cell gene-regulatory network.

Philosophy: mechanics never reads genes directly. Genes regulate
each other; their combined state is read out into a small set of
mechanical phenotype values every timestep. Only those phenotype
values are visible to TissueMechanics.
"""

import numpy as np


def _sigmoid(x, k=6.0, x0=0.5):
    return 1.0 / (1.0 + np.exp(-k * (x - x0)))


class GeneRegulatoryNetwork:

    def __init__(self):
        # ------------------------------------------------
        # Inputs (written by EmbryoSimulation each step)
        # ------------------------------------------------
        self.oxygen = 1.0
        self.morphogen = 0.0
        self.pressure = 0.0          # local mechanical stress proxy
        self.neighbor_signal = 0.0   # e.g. mean Notch of neighbors

        # ------------------------------------------------
        # Regulatory nodes (0..1-ish, not hard-clamped)
        # ------------------------------------------------
        self.sox2 = 1.0        # progenitor maintenance
        self.snail = 0.0       # EMT / migration
        self.yap = 0.0         # mechanosensing
        self.notch = 0.0       # lateral signaling
        self.cell_cycle = 0.5  # proliferative drive
        self.hif1a = 0.0       # hypoxia response
        self.noise=0.05

        # ------------------------------------------------
        # Mechanical phenotype (readout only â€” written here,
        # read by Cell/TissueMechanics, never fed back as an
        # "input" to update())
        # ------------------------------------------------
        self.adhesion =  0.4     #0.7
        self.protrusion = 0.8    #0.25
        self.polarity_bias = 0.8   # how strongly polarity opposes adhesion
        self.division_rate = 1.0   # multiplier on growth_rate

    # ======================================================
    # Integrate one timestep
    # ======================================================

    def update(self, dt):

        # --- HIF1a: hypoxia sensor, rises when oxygen is low
        self.hif1a += dt * (
            (1.0 - self.oxygen) - 0.4 * self.hif1a
        )

        # --- CellCycle: driven by oxygen, suppressed by HIF1a
        # (chronic hypoxia should eventually slow cycling even
        # though acute HIF1a is itself an oxygen response)
        self.cell_cycle += dt * (
            0.8 * self.oxygen
            - 0.3 * self.hif1a
            - 0.4 * self.cell_cycle
        )

        # --- SOX2: maintained by morphogen, repressed by SNAIL
        # Increase the penalty SNAIL applies to SOX2 from 0.9 to 1.5
        self.sox2 += dt * (
            0.6 * self.morphogen
            - 1.5 * self.snail
            - 0.3 * self.sox2
        )

       # --- SNAIL: driven by YAP, mutually represses SOX2,
        # laterally inhibited by neighbors' Notch signal, plus a
        # small noise term to break initial symmetry
        self.snail += dt * (
            0.9 * self.yap
            - 0.5 * self.sox2
            - 0.6 * self.notch          # lateral inhibition
            - 0.1 * self.snail
            + 0.05 * (self.noise - 0.5)   # symmetry-breaking noise
        )

        # --- YAP: mechanosensor, driven directly by pressure input
        self.yap += dt * (
            self.pressure
            - 0.5 * self.yap
        )

        # --- Notch: driven by neighbors' SNAIL/"Delta" signal
        # (classic lateral inhibition: a cell that gets ahead on
        # SNAIL suppresses its neighbors' SNAIL via their Notch)
        self.notch += dt * (
            self.neighbor_signal
            - 0.4 * self.notch
        )

        # --- Clamp regulatory nodes to a sane biological range
        for name in ("sox2", "snail", "yap", "notch", "cell_cycle", "hif1a"):
            setattr(self, name, float(np.clip(getattr(self, name), 0.0, 2.0)))

        self._read_out_phenotype()

    # ======================================================
    # Phenotype readout (mechanics-facing)
    # ======================================================

    def _read_out_phenotype(self):
        # High SOX2 -> adherent, epithelial-like progenitor state.
        self.adhesion = _sigmoid(self.sox2, k=5.0, x0=0.5)

        # High SNAIL / YAP -> motile, protrusive, low-adhesion state.
        self.protrusion = _sigmoid(
            0.6 * self.snail + 0.4 * self.yap, k=5.0, x0=0.4
        )

        # Polarity bias: how strongly the cell polarizes opposite to
        # local adhesion pull. High-adhesion progenitors barely
        # polarize; high-SNAIL cells polarize strongly.
        self.polarity_bias = _sigmoid(self.snail, k=5.0, x0=0.3)

        # Division rate multiplier from CellCycle node.
        self.division_rate = float(np.clip(self.cell_cycle, 0.0, 2.0))

    # ======================================================
    # Fate classification (discrete label from continuous state)
    # ======================================================

    def classify_fate(self):
        if self.snail > 0.8:
            return "MESENCHYMAL"
        if self.sox2 > 0.8:
            return "PROGENITOR"
        return "TRANSITIONAL"