"""
grn_debug.py -- test the GRN outside Omniverse.

Put this file in the same folder as grn.py, open a normal terminal there:

    python grn_debug.py                  every single-cell scenario + checks
    python grn_debug.py rest squeezed_at_10h      only the ones you name
    python grn_debug.py --pair           two squeezed neighbours (tip vs stalk test)
    python grn_debug.py --dt 1.0         the bigger time step the sim can use
    python grn_debug.py --no-plot        tables and checks only
"""

import sys
import numpy as np
from grn import GeneRegulatoryNetwork

NODES = ["sox2", "snail", "yap", "notch", "cell_cycle", "hif1a"]
PHENO = ["adhesion", "protrusion", "polarity_bias", "division_rate"]
COLS = ["t"] + NODES + PHENO


# ----------------------------------------------------------------------
# What does the cell sense over time?   (oxygen, morphogen, pressure,
# neighbour_signal).  Pressure 1.0 or less is "normal contact" and should
# be ignored; above 1.0 counts as being squeezed.
# ----------------------------------------------------------------------
SCENARIOS = {
    "rest": lambda t: dict(
        oxygen=1.0, morphogen=0.0, pressure=0.0, neighbor_signal=0.0),

    "morphogen_on": lambda t: dict(
        oxygen=1.0, morphogen=1.0, pressure=0.0, neighbor_signal=0.0),

    "mild_contact": lambda t: dict(
        oxygen=1.0, morphogen=0.5, pressure=0.5, neighbor_signal=0.0),

    "squeezed_at_10h": lambda t: dict(
        oxygen=1.0, morphogen=0.5,
        pressure=3.0 if t >= 10 else 0.0, neighbor_signal=0.0),

    "low_oxygen": lambda t: dict(
        oxygen=0.2, morphogen=0.5, pressure=0.0, neighbor_signal=0.0),

    "loud_neighbors_at_10h": lambda t: dict(
        oxygen=1.0, morphogen=0.5, pressure=3.0,
        neighbor_signal=1.0 if t >= 10 else 0.0),
}


# ----------------------------------------------------------------------
# Running
# ----------------------------------------------------------------------
def run(name, hours=30.0, dt=0.1, seed=0):
    np.random.seed(seed)                  # same noise every run
    grn = GeneRegulatoryNetwork()
    sense = SCENARIOS[name]
    rows = []
    for i in range(int(hours / dt) + 1):
        t = i * dt
        for key, value in sense(t).items():
            setattr(grn, key, value)
        grn.update(dt)
        rows.append([t] + [getattr(grn, n) for n in NODES + PHENO])
    return np.array(rows)


def v(data, col, t=None):
    """Value of a column at time t (last row if t is None)."""
    c = COLS.index(col)
    if t is None:
        return data[-1, c]
    row = int(np.argmin(np.abs(data[:, 0] - t)))
    return data[row, c]


# ----------------------------------------------------------------------
# Checks: plain-language "what should happen" for each scenario
# ----------------------------------------------------------------------
CHECKS = {
    "rest": [
        ("SOX2 settles near 1 (does not drain away)",
         lambda d: 0.8 <= v(d, "sox2") <= 1.2),
        ("SNAIL stays off",
         lambda d: v(d, "snail") < 0.1),
        ("cell stays sticky (adhesion > 0.8)",
         lambda d: v(d, "adhesion") > 0.8),
        ("division_rate settles near 1 (normal growth speed)",
         lambda d: 0.8 <= v(d, "division_rate") <= 1.2),
    ],
    "morphogen_on": [
        ("SOX2 rises above its resting level",
         lambda d: v(d, "sox2") > 1.5),
        ("SNAIL stays off",
         lambda d: v(d, "snail") < 0.1),
    ],
    "mild_contact": [
        ("pressure below the threshold does not move YAP",
         lambda d: d[:, COLS.index("yap")].max() < 0.05),
        ("SNAIL stays off",
         lambda d: d[:, COLS.index("snail")].max() < 0.1),
    ],
    "squeezed_at_10h": [
        ("YAP is quiet before the squeeze",
         lambda d: v(d, "yap", 9.0) < 0.05),
        ("YAP rises after the squeeze",
         lambda d: v(d, "yap") > 0.8),
        ("SNAIL switches on after the squeeze",
         lambda d: v(d, "snail") > 0.5),
        ("protrusion grows after the squeeze",
         lambda d: v(d, "protrusion") > v(d, "protrusion", 9.0) + 0.2),
    ],
    "low_oxygen": [
        ("HIF1a (low-oxygen alarm) rises",
         lambda d: v(d, "hif1a") > 1.0),
        ("growth slows right down",
         lambda d: v(d, "division_rate") < 0.2),
    ],
    "loud_neighbors_at_10h": [
        ("SNAIL is on before neighbours get loud",
         lambda d: v(d, "snail", 9.0) > 0.5),
        ("SNAIL is shut off once neighbours get loud (lateral inhibition)",
         lambda d: v(d, "snail") < 0.1),
    ],
}


def run_checks(name, data):
    ok = True
    for text, fn in CHECKS.get(name, []):
        passed = bool(fn(data))
        ok &= passed
        print(f"   [{'PASS' if passed else 'FAIL'}] {text}")
    return ok


# ----------------------------------------------------------------------
# Output helpers
# ----------------------------------------------------------------------
def print_table(name, data, every=5.0):
    print(f"\n=== {name} ===")
    print(" ".join(f"{c[:8]:>8}" for c in COLS))
    dt = data[1, 0] - data[0, 0]
    step = max(1, int(round(every / dt)))
    for row in data[::step]:
        print(" ".join(f"{x:8.3f}" for x in row))


def plot(results, title_suffix=""):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed: tables only. "
              "Install with: pip install matplotlib)")
        return
    n = len(results)
    fig, axes = plt.subplots(2, n, figsize=(4 * n, 6),
                             squeeze=False, sharex=True)
    for col, (name, data) in enumerate(results.items()):
        t = data[:, 0]
        for k, node in enumerate(NODES):
            axes[0][col].plot(t, data[:, 1 + k], label=node)
        for k, node in enumerate(PHENO):
            axes[1][col].plot(t, data[:, 1 + len(NODES) + k], label=node)
        axes[0][col].set_title(name + title_suffix, fontsize=9)
        axes[1][col].set_xlabel("hours")
    axes[0][0].set_ylabel("gene bars")
    axes[1][0].set_ylabel("what mechanics sees")
    axes[0][0].legend(fontsize=7)
    axes[1][0].legend(fontsize=7)
    plt.tight_layout()
    plt.show()


# ----------------------------------------------------------------------
# Two-cell test: do two equally squeezed neighbours split into one tip
# cell and one stalk cell?  (each cell hears the other's SNAIL)
# ----------------------------------------------------------------------
def run_pair(hours=40.0, dt=0.1, seed=1, pressure=3.0, morphogen=0.5):
    np.random.seed(seed)
    a, b = GeneRegulatoryNetwork(), GeneRegulatoryNetwork()
    for g in (a, b):
        g.morphogen = morphogen
        g.pressure = pressure
    rows = []
    for i in range(int(hours / dt) + 1):
        a.neighbor_signal = b.snail
        b.neighbor_signal = a.snail
        a.update(dt)
        b.update(dt)
        rows.append([i * dt, a.snail, b.snail, a.notch, b.notch,
                     a.protrusion, b.protrusion])
    return np.array(rows)


def pair_report(data, every=5.0):
    print("\n=== pair: two equally squeezed neighbours ===")
    print(f"{'t':>6} {'snailA':>8} {'snailB':>8} {'notchA':>8} "
          f"{'notchB':>8} {'protA':>8} {'protB':>8}")
    dt = data[1, 0] - data[0, 0]
    step = max(1, int(round(every / dt)))
    for row in data[::step]:
        print(" ".join(f"{x:8.3f}" for x in row))
    sa, sb = data[-1, 1], data[-1, 2]
    split = abs(sa - sb) > 1.0
    print(f"   [{'PASS' if split else 'FAIL'}] one cell became the tip, "
          f"the other stayed stalk (|snailA - snailB| = {abs(sa - sb):.2f})")
    return split


def plot_pair(data):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
    ax[0].plot(data[:, 0], data[:, 1], label="snail A")
    ax[0].plot(data[:, 0], data[:, 2], label="snail B")
    ax[0].set_title("SNAIL (tip identity)")
    ax[1].plot(data[:, 0], data[:, 3], label="notch A")
    ax[1].plot(data[:, 0], data[:, 4], label="notch B")
    ax[1].set_title("Notch (being held back)")
    for a in ax:
        a.set_xlabel("hours")
        a.legend(fontsize=8)
    plt.tight_layout()
    plt.show()


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    args = sys.argv[1:]

    dt = 0.1
    if "--dt" in args:
        i = args.index("--dt")
        dt = float(args[i + 1])
        del args[i:i + 2]

    show_plots = "--no-plot" not in args
    args = [a for a in args if a != "--no-plot"]

    if "--pair" in args:
        data = run_pair(dt=dt)
        pair_report(data)
        if show_plots:
            plot_pair(data)
        sys.exit(0)

    names = args or list(SCENARIOS)
    unknown = [n for n in names if n not in SCENARIOS]
    if unknown:
        print("Unknown scenario(s):", unknown)
        print("Available:", list(SCENARIOS))
        sys.exit(1)

    results = {}
    all_ok = True
    for n in names:
        results[n] = run(n, dt=dt)
        print_table(n, results[n])
        all_ok &= run_checks(n, results[n])

    print("\nOVERALL:", "ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    if show_plots:
        plot(results, title_suffix=f"  (dt={dt})")
