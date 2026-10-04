"""Graph-based 3D+time tumor world.

Nodes are cells (position, volume, cycle state). Edges join cells that touch.
Chemical fields live on the nodes and diffuse along edges (graph Laplacian) and
exchange with the surrounding medium through each cell's free surface.
"""
import numpy as np
import scipy.sparse as sp
from scipy import ndimage
from scipy.sparse.linalg import cg
from scipy.spatial import cKDTree

from . import params as P
from .fibroblasts import FRACTION as FIBRO_FRACTION, Fibroblasts
from .indications import DEFAULT as DEFAULT_INDICATION, INDICATIONS, MutationModel
from .tcells import TCells

ALIVE, NECROTIC, APOPTOTIC = 0, 1, 2
UM3_TO_L = 1e-15


def radius(volume):
    return np.cbrt(3.0 * volume / (4.0 * np.pi))


class Field:
    """A diffusing species defined on cell nodes.

    vmax: max per-cell uptake in (field units * L)/s, Michaelis-Menten with Km;
    k1:   alternatively a first-order uptake rate (1/s).
    """

    def __init__(self, name, D, bath, vmax=0.0, Km=0.0, k1=0.0):
        self.name, self.D, self.bath = name, D, bath
        self.vmax, self.Km, self.k1 = vmax, Km, k1
        self.c = None
        self.uptake = None  # last per-cell uptake, amount/s in field units * L


class World:
    APOPTOTIC = APOPTOTIC

    def __init__(self, seed=0, indication=DEFAULT_INDICATION):
        self.rng = np.random.default_rng(seed)
        self.indication = INDICATIONS[indication]
        self.mutations = MutationModel(indication, self.rng)
        self.t = 0.0  # hours
        self.next_id = 0
        n = 0
        self.id = np.zeros(n, int)
        self.ctype = np.zeros(n, int)  # 0 = HeLa; more identities later
        self.pos = np.zeros((n, 3))
        self.vol = np.zeros(n)
        self.age = np.zeros(n)        # hours into current cycle
        self.T = np.zeros(n)          # this cell's cycle length (h)
        self.state = np.zeros(n, int)
        self.damage = np.zeros(n)       # accumulated sublethal CTL hits (decays)
        self.t_dead = np.zeros(n)       # time of apoptosis (APOPTOTIC cells)
        self.n_mut = np.zeros(n, np.int64)      # acquired somatic mutations
        self.driver_mask = np.zeros(n, np.int64)  # bitmask of mutated driver genes
        self.events = []              # (t, kind, id, other) - birth/death/enter/exit
        self.edges = np.zeros((0, 2), int)
        self.dist = np.zeros(0)

        self.tcells = TCells(self)
        self.fibroblasts = Fibroblasts(self)
        o, g, l, m = P.OXYGEN, P.GLUCOSE, P.LACTATE, P.MITOGEN
        # oxygen in mmHg: amount/s -> mmHg*L/s via solubility
        self.fields = {
            "oxygen": Field("oxygen", o["D"].value, o["bath"].value,
                            vmax=o["uptake"].value / o["solubility"].value, Km=o["Km"].value),
            # glucose/lactate in mM: mol/s -> mM*L/s
            "glucose": Field("glucose", g["D"].value, g["bath"].value,
                             vmax=g["uptake"].value * 1e3, Km=g["Km"].value),
            "lactate": Field("lactate", l["D"].value, l["bath"].value),
            "mitogen": Field("mitogen", m["D"].value, m["bath"].value, k1=m["uptake_rate"].value),
        }

    # ------------------------------------------------------------ population
    @property
    def n(self):
        return len(self.id)

    def _draw_cycle(self, k, driver_mask=None):
        """Cycle lengths (h) for k cells: indication doubling time, lognormal spread,
        shortened slightly per acquired driver (near-neutral; see indications.py)."""
        mean, cv = self.indication.doubling_h, P.CYCLE["cv"].value
        s = np.sqrt(np.log(1 + cv**2))
        T = self.rng.lognormal(np.log(mean) - s**2 / 2, s, k)
        if driver_mask is not None and k:
            T = T * self.mutations.cycle_scale(driver_mask)
        return T

    def add_cells(self, pos, kind="seed", parent=-1, ctype=0, vol=None, age=None,
                  n_mut=None, driver_mask=None):
        """Add cells to the world (seeding, division, or immune-cell entry later)."""
        pos = np.atleast_2d(pos)
        k = len(pos)
        ids = np.arange(self.next_id, self.next_id + k)
        self.next_id += k
        T = self._draw_cycle(k, driver_mask)
        if age is None:
            age = self.rng.uniform(0, 1, k) * T
        if vol is None:
            vol = P.VOLUME["birth"].value * (1 + age / T)
        self.id = np.concatenate([self.id, ids])
        self.ctype = np.concatenate([self.ctype, np.full(k, ctype)])
        self.pos = np.vstack([self.pos, pos])
        self.vol = np.concatenate([self.vol, vol])
        self.age = np.concatenate([self.age, age])
        self.T = np.concatenate([self.T, T])
        self.state = np.concatenate([self.state, np.full(k, ALIVE)])
        self.damage = np.concatenate([self.damage, np.zeros(k)])
        self.n_mut = np.concatenate([self.n_mut, np.zeros(k, np.int64) if n_mut is None else n_mut])
        self.driver_mask = np.concatenate(
            [self.driver_mask, np.zeros(k, np.int64) if driver_mask is None else driver_mask])
        self.t_dead = np.concatenate([self.t_dead, np.zeros(k)])
        parent = np.broadcast_to(parent, k)
        self.events += [(self.t, kind, int(i), int(p)) for i, p in zip(ids, parent)]
        for f in self.fields.values():
            if f.c is not None:
                f.c = np.concatenate([f.c, np.full(k, f.bath)])
        return ids

    def remove_cells(self, mask, kind="exit"):
        """Remove cells (metastatic exit, clearance). Logged."""
        for i in self.id[mask]:
            self.events.append((self.t, kind, int(i), -1))
        keep = ~mask
        for a in ("id", "ctype", "pos", "vol", "age", "T", "state", "damage", "t_dead",
                  "n_mut", "driver_mask"):
            setattr(self, a, getattr(self, a)[keep])
        for f in self.fields.values():
            if f.c is not None:
                f.c = f.c[keep]

    def seed_fibroblasts(self, fraction=None):
        """Seed fibroblasts at the indication's fraction of current tumor cells."""
        f = FIBRO_FRACTION[self.indication.key] if fraction is None else fraction
        k = int(round(self.n * f))
        if k:
            self.fibroblasts.seed(k)
        return k

    def seed_ball(self, n_cells):
        """Random cluster of n cells, like a freshly aggregated spheroid."""
        r_cell = radius(P.VOLUME["mean"].value)
        R = r_cell * n_cells ** (1 / 3)
        pts = []
        while len(pts) < n_cells:
            p = self.rng.uniform(-R, R, 3)
            if p @ p <= R * R:
                pts.append(p)
        self.add_cells(np.array(pts))
        for _ in range(60):
            self.relax()

    # ------------------------------------------------------------ graph
    def build_graph(self):
        r = radius(self.vol)
        tol = P.NUMERICS["contact_tol"].value * P.NUMERICS["packing"].value
        tree = cKDTree(self.pos)
        pairs = tree.query_pairs(2 * tol * r.max(), output_type="ndarray")
        if len(pairs) == 0:
            self.edges, self.dist = np.zeros((0, 2), int), np.zeros(0)
            return
        i, j = pairs[:, 0], pairs[:, 1]
        d = np.linalg.norm(self.pos[i] - self.pos[j], axis=1)
        touch = d < tol * (r[i] + r[j])
        self.edges, self.dist = pairs[touch], d[touch]

    def degree(self):
        return np.bincount(self.edges.ravel(), minlength=self.n)

    # ------------------------------------------------------------ mechanics
    def relax(self):
        """One overlap-relaxation pass: push overlapping cells apart, pull near ones in."""
        self.build_graph()
        r = radius(self.vol)
        i, j = self.edges[:, 0], self.edges[:, 1]
        rest = P.NUMERICS["packing"].value * (r[i] + r[j])
        d = np.maximum(self.dist, 1e-6)
        u = (self.pos[j] - self.pos[i]) / d[:, None]
        gap = d - rest  # <0 overlap
        frac = np.where(gap < 0, P.NUMERICS["repulsion"].value, P.NUMERICS["adhesion"].value)
        disp = (0.5 * frac * gap)[:, None] * u  # move i toward j by this (away if overlap)
        dx = np.zeros_like(self.pos)
        for ax in range(3):
            dx[:, ax] += np.bincount(i, disp[:, ax], self.n) - np.bincount(j, disp[:, ax], self.n)
        self.pos += dx

    # ------------------------------------------------------------ fields
    def exposed(self):
        """Cells touching the medium, found by flood-filling empty space from outside
        on a voxel grid (one voxel ~ one cell spacing). Works for any tumor shape."""
        h = 2 * P.NUMERICS["packing"].value * radius(self.vol.mean())
        lo = self.pos.min(0) - 3 * h
        idx = np.floor((self.pos - lo) / h).astype(int)
        shape = idx.max(0) + 4
        occ = np.zeros(shape, bool)
        occ[tuple(idx.T)] = True
        occ = ndimage.binary_closing(occ, iterations=1)  # fill packing voids
        lab, _ = ndimage.label(~occ)
        medium = lab == lab[0, 0, 0]
        near_medium = ndimage.binary_dilation(medium)
        return near_medium[tuple(idx.T)] | (self.degree() == 0)  # strays sit in medium

    def solve_fields(self):
        """Quasi-steady diffusion-reaction on the contact graph (one sparse solve per field).

        Edge conductance D*kappa_ij/d_ij^2 with kappa_ij = 12/(z_i+z_j) makes the sum
        over neighbors approximate D*laplacian(c) for an isotropic neighborhood.
        Surface cells exchange with the medium through one face at half a diameter.
        """
        self.build_graph()
        n = self.n
        deg = self.degree()
        i, j = self.edges[:, 0], self.edges[:, 1]
        kappa = 12.0 / np.maximum(deg[i] + deg[j], 1)
        # space each cell occupies (cell + its share of gaps), from its mean contact
        # distance at random close packing: n = 0.64 / (pi/6 * d^3)
        dsum = np.bincount(i, self.dist, n) + np.bincount(j, self.dist, n)
        dbar = np.where(deg > 0, dsum / np.maximum(deg, 1), 2 * radius(self.vol))
        V = np.maximum(np.pi / 6 * dbar**3 / 0.64, self.vol) * UM3_TO_L
        h = 2 * P.NUMERICS["packing"].value * radius(self.vol)
        surface = self.exposed()
        self.surface = surface
        alive = self.state == ALIVE

        for f in self.fields.values():
            g = f.D * kappa / self.dist**2                           # 1/s per edge
            gb = surface * 2 * f.D / h**2                            # 1/s to medium
            L = sp.coo_matrix((np.r_[-g, -g], (np.r_[i, j], np.r_[j, i])), shape=(n, n)).tocsr()
            diag0 = np.bincount(i, g, n) + np.bincount(j, g, n) + gb
            src = np.zeros(n)
            if f.name == "lactate":
                src = np.where(alive, P.LACTATE["per_glucose"].value
                               * self.fields["glucose"].uptake / V, 0.0)
            c = np.full(n, f.bath) if f.c is None or len(f.c) != n else f.c.copy()
            for _ in range(6 if f.vmax else 1):  # Picard iterations for Michaelis-Menten
                k = np.zeros(n)
                if f.vmax:
                    k = f.vmax / V / (f.Km + np.maximum(c, 0))
                if f.k1:
                    k = np.full(n, f.k1)
                k = np.where(alive, k, 0.0)
                d = diag0 + k
                A = L + sp.diags(d)  # symmetric positive definite -> CG, warm-started
                c, _ = cg(A, gb * f.bath + src, x0=c, rtol=1e-8, M=sp.diags(1 / d))
            f.c = c
            if f.vmax:
                f.uptake = np.where(alive, f.vmax * c / (f.Km + np.maximum(c, 0)), 0.0)

    # ------------------------------------------------------------ biology
    def phase(self):
        """FUCCI-style label: 0 = G1 (red), 1 = S/G2/M (green), -1 = necrotic, -2 = apoptotic."""
        g1 = self.T * P.CYCLE["G1"].value / P.CYCLE_MEAN_H
        ph = (self.age >= g1).astype(int)
        ph[self.state == NECROTIC] = -1
        ph[self.state == APOPTOTIC] = -2
        return ph

    def step(self):
        dt = P.NUMERICS["dt"].value
        self.solve_fields()
        o2 = self.fields["oxygen"].c
        mit = self.fields["mitogen"].c
        alive = self.state == ALIVE

        # necrosis under severe hypoxia
        p_die = 1 - np.exp(-dt / P.OXYGEN["necrosis_time"].value)
        die = alive & (o2 < P.OXYGEN["necrosis"].value) & (self.rng.random(self.n) < p_die)
        self.state[die] = NECROTIC
        for i in self.id[die]:
            self.events.append((self.t, "necrosis", int(i), -1))

        # cycle progression; G1 cells arrest if O2 or mitogen is low (G1 checkpoint)
        alive = self.state == ALIVE
        in_g1 = self.phase() == 0
        arrest = in_g1 & ((o2 < P.OXYGEN["arrest"].value) | (mit < P.MITOGEN["arrest"].value))
        grow = alive & ~arrest
        self.age[grow] += dt
        self.vol[grow] = P.VOLUME["birth"].value * (1 + np.minimum(self.age[grow] / self.T[grow], 1))
        self.quiescent = alive & arrest

        # division
        div = np.where(alive & (self.age >= self.T))[0]
        if len(div):
            u = self.rng.normal(size=(len(div), 3))
            u /= np.linalg.norm(u, axis=1)[:, None]
            off = 0.5 * radius(self.vol[div] / 2)[:, None] * u
            parents = self.id[div].copy()
            child_vol = self.vol[div] / 2
            self.pos[div] -= off
            self.vol[div] = child_vol
            self.age[div] = 0.0
            # both daughters replicate the genome independently, so each draws its own mutations
            keep_mut, keep_drv = self.mutations.on_division(
                self.n_mut[div].copy(), self.driver_mask[div].copy())
            new_mut, new_drv = self.mutations.on_division(
                self.n_mut[div].copy(), self.driver_mask[div].copy())
            gained = (keep_drv != self.driver_mask[div]) | (new_drv != self.driver_mask[div])
            self.n_mut[div], self.driver_mask[div] = keep_mut, keep_drv
            self.T[div] = self._draw_cycle(len(div), keep_drv)
            for p in parents:
                self.events.append((self.t, "division", int(p), -1))
            for i in np.where(gained)[0]:
                self.events.append((self.t, "driver", int(parents[i]), -1))
            self.add_cells(self.pos[div] + 2 * off, kind="birth", parent=parents,
                           vol=child_vol, age=np.zeros(len(div)),
                           n_mut=new_mut, driver_mask=new_drv)

        for _ in range(int(P.NUMERICS["mech_iters"].value)):
            self.relax()

        # fibroblasts deposit matrix, then T cells move through it
        self.fibroblasts.step(dt)
        self.tcells.step(dt)
        gone = (self.state == APOPTOTIC) & (self.t + dt >= self.t_dead + P.TCELL["apoptotic_clearance"].value)
        if gone.any():
            self.t += dt  # log clearance at the end of the step
            self.remove_cells(gone, kind="cleared")
            self.t -= dt
        self.t += dt

    # ------------------------------------------------------------ readouts
    def diameter(self):
        c = self.pos.mean(0)
        r = np.linalg.norm(self.pos - c, axis=1)
        return 2 * (np.percentile(r, 99) + radius(P.VOLUME["mean"].value))

    def snapshot(self):
        return dict(t=self.t, id=self.id.copy(), pos=self.pos.copy(), vol=self.vol.copy(),
                    phase=self.phase(), state=self.state.copy(),
                    n_mut=self.n_mut.copy(), driver_mask=self.driver_mask.copy(),
                    **{k: f.c.copy() for k, f in self.fields.items() if f.c is not None})
