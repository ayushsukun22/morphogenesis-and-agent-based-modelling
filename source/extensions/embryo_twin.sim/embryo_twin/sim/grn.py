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

        # ------------------------------------------------
        # Mechanical phenotype (readout only)
        # ------------------------------------------------
        self.adhesion = 0.7
        self.protrusion = 0.25
        self.polarity_bias = 0.8
        self.division_rate = 1.0

    # ======================================================
    # Integrate one timestep
    # ======================================================

    def update(self, dt):

        # HIF1a: low-oxygen alarm
        self.hif1a += dt * (
            (1.0 - self.oxygen)
            - 0.4 * self.hif1a
        )

        # CellCycle: growth drive (settles at 1.0 with full oxygen)
        self.cell_cycle += dt * (
            0.4 * self.oxygen
            - 0.3 * self.hif1a
            - 0.4 * self.cell_cycle
        )

        # SOX2: holds itself up (0.3 baseline), boosted by morphogen,
        # crushed by SNAIL
        self.sox2 += dt * (
            0.3
            + 0.6 * self.morphogen
            - 0.9 * self.snail
            - 0.3 * self.sox2
        )

        # SNAIL: driven by YAP, braked by SOX2 and by neighbours' Notch
        self.snail += dt * (
            1.2 * self.yap
            - 0.4 * self.sox2
            - 0.6 * self.notch
            - 0.4 * self.snail
            + 0.05 * (np.random.random() - 0.5)
        )

        # YAP: only squeeze above normal contact counts, and it is capped
        squeeze = min(max(self.pressure - 1.0, 0.0), 1.0)
        self.yap += dt * (0.5 * squeeze - 0.5 * self.yap)

        # Notch: hears how loud the neighbours' SNAIL is
        self.notch += dt * (
            self.neighbor_signal
            - 0.4 * self.notch
        )

        for name in ("sox2", "snail", "yap", "notch", "cell_cycle", "hif1a"):
            setattr(self, name, float(np.clip(getattr(self, name), 0.0, 2.0)))

        self._read_out_phenotype()

    # ======================================================
    # Phenotype readout (mechanics-facing)
    # ======================================================

    def _read_out_phenotype(self):
        self.adhesion = _sigmoid(self.sox2, k=5.0, x0=0.5)

        self.protrusion = _sigmoid(
            0.6 * self.snail + 0.4 * self.yap, k=5.0, x0=0.4
        )

        self.polarity_bias = _sigmoid(self.snail, k=5.0, x0=0.3)

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