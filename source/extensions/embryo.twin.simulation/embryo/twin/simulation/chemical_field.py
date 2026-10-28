import numpy as np


class ChemicalField:
    """
    A 3D reaction-diffusion field on a regular grid.
    Cells sample it and deplete it locally as a sink.
    """

    def __init__(
        self,
        bounds_min,
        bounds_max,
        resolution=32,
        diffusion_rate=0.5,
        decay_rate=0.02,
        baseline=0.0,
    ):
        self.bounds_min = np.asarray(bounds_min, dtype=np.float64)
        self.bounds_max = np.asarray(bounds_max, dtype=np.float64)
        self.resolution = resolution

        self.spacing = (self.bounds_max - self.bounds_min) / (resolution - 1)

        # Concentration grid, initialized to baseline.
        self.grid = np.full(
            (resolution, resolution, resolution),
            baseline,
            dtype=np.float64,
        )

        self.diffusion_rate = diffusion_rate
        self.decay_rate = decay_rate

    # ------------------------------------------------
    # Coordinate mapping
    # ------------------------------------------------

    def _to_grid_coords(self, position):
        """Continuous world position -> fractional grid index (float)."""
        rel = (position - self.bounds_min) / self.spacing
        return np.clip(rel, 0.0, self.resolution - 1.001)

    # ------------------------------------------------
    # Sampling (world -> field), trilinear interpolation
    # ------------------------------------------------

    def sample(self, position):
        x, y, z = self._to_grid_coords(position)

        x0, y0, z0 = int(x), int(y), int(z)
        x1, y1, z1 = x0 + 1, y0 + 1, z0 + 1

        fx, fy, fz = x - x0, y - y0, z - z0

        g = self.grid

        c000 = g[x0, y0, z0]; c100 = g[x1, y0, z0]
        c010 = g[x0, y1, z0]; c110 = g[x1, y1, z0]
        c001 = g[x0, y0, z1]; c101 = g[x1, y0, z1]
        c011 = g[x0, y1, z1]; c111 = g[x1, y1, z1]

        c00 = c000 * (1 - fx) + c100 * fx
        c10 = c010 * (1 - fx) + c110 * fx
        c01 = c001 * (1 - fx) + c101 * fx
        c11 = c011 * (1 - fx) + c111 * fx

        c0 = c00 * (1 - fy) + c10 * fy
        c1 = c01 * (1 - fy) + c11 * fy

        return float(c0 * (1 - fz) + c1 * fz)

    # ------------------------------------------------
    # Consumption (sink), trilinear-weighted deposit
    # ------------------------------------------------

    def consume(self, position, amount):
        """Remove `amount` of concentration, spread across
        the 8 nearest voxels by the same weights used for sampling."""
        x, y, z = self._to_grid_coords(position)

        x0, y0, z0 = int(x), int(y), int(z)
        x1, y1, z1 = x0 + 1, y0 + 1, z0 + 1

        fx, fy, fz = x - x0, y - y0, z - z0

        weights = {
            (x0, y0, z0): (1 - fx) * (1 - fy) * (1 - fz),
            (x1, y0, z0): fx * (1 - fy) * (1 - fz),
            (x0, y1, z0): (1 - fx) * fy * (1 - fz),
            (x1, y1, z0): fx * fy * (1 - fz),
            (x0, y0, z1): (1 - fx) * (1 - fy) * fz,
            (x1, y0, z1): fx * (1 - fy) * fz,
            (x0, y1, z1): (1 - fx) * fy * fz,
            (x1, y1, z1): fx * fy * fz,
        }

        for (ix, iy, iz), w in weights.items():
            self.grid[ix, iy, iz] = max(
                0.0, self.grid[ix, iy, iz] - amount * w
            )



    # ------------------------------------------------
    # Production (source), trilinear-weighted deposit
    # ------------------------------------------------

    def deposit(self, position, amount):
        """Add `amount` of concentration, spread across
        the 8 nearest voxels using the same weights used
        for sampling and consumption.
        """

        x, y, z = self._to_grid_coords(position)

        x0, y0, z0 = int(x), int(y), int(z)
        x1, y1, z1 = x0 + 1, y0 + 1, z0 + 1

        fx, fy, fz = x - x0, y - y0, z - z0

        weights = {
            (x0, y0, z0): (1 - fx) * (1 - fy) * (1 - fz),
            (x1, y0, z0): fx * (1 - fy) * (1 - fz),
            (x0, y1, z0): (1 - fx) * fy * (1 - fz),
            (x1, y1, z0): fx * fy * (1 - fz),
            (x0, y0, z1): (1 - fx) * (1 - fy) * fz,
            (x1, y0, z1): fx * (1 - fy) * fz,
            (x0, y1, z1): (1 - fx) * fy * fz,
            (x1, y1, z1): fx * fy * fz,
        }

        for (ix, iy, iz), w in weights.items():
            self.grid[ix, iy, iz] += amount * w

        np.clip(
            self.grid,
            0.0,
            None,
            out=self.grid,
        )



    # ------------------------------------------------
    # Field evolution: diffusion + decay (+ optional source)
    # ------------------------------------------------

    def step(self, dt, source_mask=None, source_rate=0.0):
        g = self.grid

        # 6-neighbor discrete Laplacian, Neumann (zero-flux) boundaries
        # via edge padding.
        padded = np.pad(g, 1, mode="edge")

        laplacian = (
            padded[2:, 1:-1, 1:-1] + padded[:-2, 1:-1, 1:-1] +
            padded[1:-1, 2:, 1:-1] + padded[1:-1, :-2, 1:-1] +
            padded[1:-1, 1:-1, 2:] + padded[1:-1, 1:-1, :-2] -
            6.0 * g
        )

        g += dt * (self.diffusion_rate * laplacian - self.decay_rate * g)

        # Optional production, e.g. a vascular boundary supplying oxygen.
        if source_mask is not None:
            g += dt * source_rate * source_mask

        np.clip(g, 0.0, None, out=g)