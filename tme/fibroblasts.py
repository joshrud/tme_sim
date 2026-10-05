"""Fibroblasts and the ECM they deposit.

Fibroblasts are slow, long-lived agents that sit in the stroma and lay down matrix. Unlike the
diffusing fields, **ECM is stored, not recomputed**: it has memory, so it lives on a voxel grid
that persists across steps (this is the first field of that kind in the model).

States, following the PDAC scheme (Ohlund 2017, Biffi 2019; see docs/stromal_cells.md):
  MYCAF - juxtatumoral, TGF-beta driven, high matrix output
  ICAF  - distal, IL-1/JAK-STAT driven, low matrix output
Proximity to tumor cells is used as a proxy for TGF-beta exposure (an assumption; a real
TGF-beta field is the planned replacement).

See docs/fibroblasts.md for evidence and for the parameters that are assumptions.
"""
import numpy as np
from scipy.spatial import cKDTree

from .params import P

MYCAF, ICAF = 0, 1
STATE_NAMES = {MYCAF: "myCAF", ICAF: "iCAF"}

FIBRO = {
    "diameter": P(15.0, "um", "assumption", "modeled as a sphere of roughly tumor-cell volume"),
    "speed": P(10.0, "um/h", "proxy",
               "human fibroblasts in 3D collagen I, Hakkinen et al. 2011 doi:10.1089/ten.tea.2010.0273; "
               "normal fibroblasts, not CAFs"),
    "division_h": P(120.0, "h", "assumption", "CAF turnover in tumors; no direct measurement found"),
    "myCAF_radius": P(50.0, "um", "assumption",
                      "distance to tumor within which a fibroblast is myCAF (TGF-beta proxy)"),
    "margin": P(250.0, "um", "assumption",
                "fibroblasts stay within this distance of the tumor edge; CAFs are tissue-resident "
                "and do not disperse like circulating cells"),
    # ECM
    "ecm_voxel": P(20.0, "um", "numerical", "ECM grid spacing, ~1.3 cell diameters"),
    "ecm_radius": P(30.0, "um", "assumption", "matrix is deposited locally around the cell"),
    "ecm_deposit": P(0.0012, "1/h", "calibrated",
                     "density added per voxel per hour by a myCAF; set so the PDAC tumor edge reaches "
                     "0.72 density at 7 days (>0.7 target; desmoplasia, Whatcott 2015 / Erkan 2008). "
                     "The same rate for every indication: ordering comes from fibroblast number alone"),
    "ecm_icaf_scale": P(0.3, "-", "assumption", "iCAFs deposit less matrix than myCAFs"),
    "ecm_decay": P(0.002, "1/h", "assumption", "slow turnover; MMP activity not modeled"),
    "ecm_barrier": P(0.6, "-", "assumption",
                     "T cells cannot move into matrix denser than this. Salmon et al. 2012 report that "
                     "aligned dense fibers restrict T cells from entering tumor islets; a speed penalty "
                     "alone reproduces the opposite (slower cells dwell near tumor and kill more)"),
    "ecm_block": P(0.8, "-", "assumption",
                   "fraction of T-cell speed lost in fully dense matrix; the effect is measured "
                   "(Salmon et al. 2012 doi:10.1172/JCI45817), this functional form is not"),
}

# Fraction of all cells that are fibroblasts, per indication. ASSUMPTIONS - the ordering
# (PDAC >> OV ~ BRCA > LUAD) is well supported, the absolute values are not. See docs/fibroblasts.md.
FRACTION = {"PDAC": 0.35, "OV": 0.15, "BRCA": 0.15, "LUAD": 0.10, "CESC": 0.10}


class EcmGrid:
    """Scalar ECM density in [0, 1] on a fixed voxel grid. Persistent (stored) state."""

    def __init__(self, extent_um, voxel=None):
        self.h = voxel or FIBRO["ecm_voxel"].value
        n = int(2 * extent_um / self.h) + 4
        self.n = n
        self.origin = -np.array([extent_um, extent_um, extent_um]) - 2 * self.h
        self.rho = np.zeros((n, n, n), np.float32)

    def idx(self, pos):
        i = np.floor((pos - self.origin) / self.h).astype(int)
        return np.clip(i, 0, self.n - 1)

    def at(self, pos):
        i = self.idx(np.atleast_2d(pos))
        return self.rho[i[:, 0], i[:, 1], i[:, 2]]

    def _stencil(self):
        """Normalized spherical kernel: matrix is laid down locally around the cell.
        Built once; deposition adds this stencil, so stored matrix is never re-blurred."""
        if getattr(self, "_sten", None) is None:
            r = max(int(round(FIBRO["ecm_radius"].value / self.h)), 1)
            g = np.arange(-r, r + 1) * self.h
            d2 = g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2
            # 1 inside the deposition radius, so `amount` is read directly as
            # density added per voxel per hour (not divided across the neighborhood)
            k = (d2 <= FIBRO["ecm_radius"].value ** 2).astype(np.float32)
            self._sten, self._r = k, r
        return self._sten, self._r

    def deposit(self, pos, amount):
        """Add matrix in a sphere around each position. Amount may be scalar or per-cell."""
        pos = np.atleast_2d(pos)
        if len(pos) == 0:
            return
        sten, r = self._stencil()
        amount = np.broadcast_to(np.asarray(amount, np.float32), (len(pos),))
        i = self.idx(pos)
        lo, hi = i - r, i + r + 1
        inside = np.all((lo >= 0) & (hi <= self.n), axis=1)
        for k in np.where(inside)[0]:
            a, b = lo[k], hi[k]
            self.rho[a[0]:b[0], a[1]:b[1], a[2]:b[2]] += sten * amount[k]
        np.clip(self.rho, 0, 1, out=self.rho)

    def decay(self, dt_h):
        self.rho *= np.exp(-FIBRO["ecm_decay"].value * dt_h)

    def stats(self):
        occupied = self.rho[self.rho > 0.01]
        return dict(mean=float(self.rho.mean()), max=float(self.rho.max()),
                    mean_occupied=float(occupied.mean()) if occupied.size else 0.0,
                    frac_dense=float((self.rho > 0.7).mean()))


class Fibroblasts:
    def __init__(self, world):
        self.w = world
        self.id = np.zeros(0, int)
        self.pos = np.zeros((0, 3))
        self.dir = np.zeros((0, 3))
        self.state = np.zeros(0, int)
        self.age = np.zeros(0)
        self.ecm = None
        self.tracks = []

    @property
    def n(self):
        return len(self.id)

    def _ensure_grid(self):
        c = self.w.pos.mean(0)
        R = np.percentile(np.linalg.norm(self.w.pos - c, axis=1), 99)
        want = max(R * 2.5, 400.0)
        if self.ecm is None or want > (self.ecm.n * self.ecm.h / 2) * 0.9:
            old = self.ecm
            self.ecm = EcmGrid(want)
            if old is not None:  # carry existing matrix into the larger grid
                off = ((old.origin - self.ecm.origin) / self.ecm.h).round().astype(int)
                s = tuple(slice(o, o + old.n) for o in off)
                self.ecm.rho[s] = np.maximum(self.ecm.rho[s], old.rho)

    def seed(self, n_fib, shell=1.15):
        """Place fibroblasts around and through the tumor, denser just outside the edge."""
        w, rng = self.w, self.w.rng
        c = w.pos.mean(0)
        R = np.percentile(np.linalg.norm(w.pos - c, axis=1), 99)
        u = rng.normal(size=(n_fib, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        # radial positions: concentrated near the tumor margin, some infiltrating
        rad = R * (0.3 + shell * rng.beta(4, 2, n_fib))
        self._add(c + u * rad[:, None])

    def _add(self, pos):
        w, rng = self.w, self.w.rng
        k = len(pos)
        ids = np.arange(w.next_id, w.next_id + k)
        w.next_id += k
        d = rng.normal(size=(k, 3))
        d /= np.linalg.norm(d, axis=1)[:, None]
        self.id = np.r_[self.id, ids]
        self.pos = np.vstack([self.pos, pos])
        self.dir = np.vstack([self.dir, d])
        self.state = np.r_[self.state, np.full(k, ICAF)]
        self.age = np.r_[self.age, rng.uniform(0, FIBRO["division_h"].value, k)]
        w.events += [(w.t, "fibroblast", int(i), -1) for i in ids]
        return ids

    def step(self, hours):
        self._ensure_grid()
        if self.n == 0:
            self.ecm.decay(hours)
            return
        w, rng = self.w, self.w.rng
        # state: myCAF if close to tumor cells (proxy for TGF-beta exposure)
        tree = cKDTree(w.pos)
        d_tumor, _ = tree.query(self.pos)
        self.state = np.where(d_tumor <= FIBRO["myCAF_radius"].value, MYCAF, ICAF)

        # slow random walk, slowed further by the matrix the fibroblasts themselves deposit
        # (dense matrix traps its depositor, which is what stabilizes the stroma)
        step_um = FIBRO["speed"].value * hours * self.speed_factor(self.pos)
        turn = rng.normal(0, 0.5, (self.n, 3))
        self.dir = self.dir + turn
        self.dir /= np.linalg.norm(self.dir, axis=1)[:, None]
        self.pos = self.pos + self.dir * step_um[:, None]

        # tissue retention: CAFs are resident, so reflect any that drift past the stromal margin
        c = w.pos.mean(0)
        R_out = np.percentile(np.linalg.norm(w.pos - c, axis=1), 99) + FIBRO["margin"].value
        off = self.pos - c
        d = np.linalg.norm(off, axis=1)
        out = d > R_out
        if out.any():
            u = off[out] / d[out][:, None]
            self.pos[out] = c + u * R_out
            self.dir[out] = -u  # turn back inward

        # deposit matrix: myCAFs lay down more than iCAFs
        amount = np.where(self.state == MYCAF, FIBRO["ecm_deposit"].value,
                          FIBRO["ecm_deposit"].value * FIBRO["ecm_icaf_scale"].value) * hours
        self.ecm.deposit(self.pos, amount)
        self.ecm.decay(hours)

        # slow proliferation
        self.age += hours
        div = np.where(self.age >= FIBRO["division_h"].value)[0]
        if len(div):
            self.age[div] = 0.0
            jitter = rng.normal(0, FIBRO["diameter"].value, (len(div), 3))
            self._add(self.pos[div] + jitter)

    def ecm_at(self, pos):
        return self.ecm.at(pos) if self.ecm is not None else np.zeros(len(np.atleast_2d(pos)))

    def speed_factor(self, pos):
        """Multiplier on a migrating cell's speed from local matrix density."""
        return 1.0 - FIBRO["ecm_block"].value * self.ecm_at(pos)

    def passable(self, pos):
        """False where matrix is too dense to migrate into (barrier, not just drag)."""
        return self.ecm_at(pos) < FIBRO["ecm_barrier"].value


def table():
    rows = [("fibroblast", k, p.value, p.unit, p.status, p.source) for k, p in FIBRO.items()]
    rows += [("fibroblast", f"fraction.{k}", v, "-", "assumption",
              "fraction of all cells that are fibroblasts; ordering supported, value is not")
             for k, v in FRACTION.items()]
    return rows
