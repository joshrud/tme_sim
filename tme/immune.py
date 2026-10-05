"""Non-CD8 immune cells: macrophages, dendritic cells, NK cells, B cells, Tregs, neutrophils.

One generic agent system covers all of them: each type differs only by its parameters
(speed, lifespan, spatial preference, colour) and by two behaviours it may switch on -
suppression of CD8 killing, and antigen presentation.

CD8 T cells stay in tcells.py because they have their own multi-hit killing machinery.

Spatial placement follows documented organization (see docs/immune_cells.md):
  TAM        biased to low oxygen   (macrophages accumulate in hypoxic areas)
  B cell     tight clusters at the margin (tertiary lymphoid structures)
  Treg       co-located with CD8 T cells
  neutrophil / NK   margin-enriched (circulation-derived)
  DC         sparse, stromal

Every composition number is an ASSUMPTION reproducing the published qualitative ordering.
"""
import numpy as np
from scipy.spatial import cKDTree

from .params import P

MACROPHAGE, DC, NK, BCELL, TREG, NEUTROPHIL = range(6)
NAMES = {MACROPHAGE: "macrophage", DC: "dendritic", NK: "NK", BCELL: "B cell",
         TREG: "Treg", NEUTROPHIL: "neutrophil"}

# Viewer colours, chosen to stay distinct from tumour phase colours
# (red G1 / green S-G2 / yellow mitosis / purple apoptotic) and from white T cells.
COLORS = {MACROPHAGE: "#8a8f3a", DC: "#e08a2e", NK: "#9b6bd6",
          BCELL: "#2bb8c4", TREG: "#dfe3a8", NEUTROPHIL: "#f2e9a0"}

# radius um, speed um/min, lifespan h, spatial rule
TYPES = {
    MACROPHAGE: dict(
        radius=P(8.0, "um", "assumption", "macrophages are large; modelled as a sphere"),
        speed=P(2.0, "um/min", "assumption", "slower than T cells, faster than fibroblasts"),
        lifespan=P(20 * 24.0, "h", "assumption", "TAMs are long-lived; no clean human tumour figure"),
        rule="hypoxia", suppress=0.5, presents=False),
    DC: dict(
        radius=P(7.0, "um", "assumption", "dendritic cell body"),
        speed=P(3.0, "um/min", "assumption", "placeholder"),
        lifespan=P(7 * 24.0, "h", "assumption", "placeholder"),
        rule="stroma", suppress=0.0, presents=True),
    NK: dict(
        radius=P(5.0, "um", "assumption", "lymphocyte-sized"),
        speed=P(6.0, "um/min", "assumption", "lymphocyte-like, slower than T cells"),
        lifespan=P(7 * 24.0, "h", "assumption", "placeholder"),
        rule="margin", suppress=0.0, presents=False),
    BCELL: dict(
        radius=P(5.0, "um", "measured", "resting lymphocyte 8-11 um diameter, BioNumbers 108368"),
        speed=P(6.0, "um/min", "proxy",
                "B-cell motility coefficient ~1/5 of T cells in lymph node, Miller et al. 2002 "
                "doi:10.1126/science.1070051"),
        lifespan=P(7 * 24.0, "h", "assumption", "placeholder"),
        rule="tls", suppress=0.0, presents=True),
    TREG: dict(
        radius=P(4.0, "um", "measured", "resting lymphocyte, BioNumbers 108368"),
        speed=P(8.0, "um/min", "proxy", "T-cell-like motility"),
        lifespan=P(14 * 24.0, "h", "assumption", "placeholder"),
        rule="with_cd8", suppress=0.8, presents=False),
    NEUTROPHIL: dict(
        radius=P(5.0, "um", "assumption", "granulocyte"),
        speed=P(12.0, "um/min", "assumption", "fast-migrating myeloid cell"),
        lifespan=P(5.4 * 24, "h", "measured",
                   "human neutrophil circulatory lifespan 5.4 days by in vivo 2H2O labelling, "
                   "Pillay et al. 2010 doi:10.1182/blood-2010-01-259028 - about 10x longer than "
                   "the <1 day figure from ex vivo labelling that is often quoted"),
        rule="margin", suppress=0.3, presents=False),
}

# Monocyte -> TAM recruitment, the one well-measured kinetic here
MONOCYTE = {
    "marrow_delay": P(1.6 * 24, "h", "measured",
                      "postmitotic interval before classical monocytes leave marrow, "
                      "Patel et al. 2017 doi:10.1084/jem.20170355"),
    "circulating": P(24.0, "h", "measured", "classical monocytes circulate ~1 day, Patel et al. 2017"),
    "to_tam": P(24.0, "h", "assumption", "time in tissue before becoming a TAM"),
}

SUPPRESSION_RADIUS = P(40.0, "um", "assumption",
                       "distance over which a Treg or TAM suppresses CD8 killing")

DC_EGRESS = P(0.02, "1/h", "assumption",
              "fraction of intratumoural DCs that pick up antigen and leave for the draining node "
              "each hour; CCR7-dependent trafficking is established (Roberts et al. 2016) but the "
              "rate is not measured")

REPLENISH = P(0.5, "-", "assumption",
              "fraction of the shortfall from the target composition recruited per hour. Immune "
              "cells are continuously replaced from circulation; without this, short-lived types "
              "(neutrophils, 24 h) vanish and the compartment decays away over a few days")

# Per-indication immune composition. ALL ASSUMPTIONS - see docs/immune_cells.md.
# leukocyte = immune cells as a fraction of tumour cells; the rest are fractions of
# the immune compartment and sum to 1.
COMPOSITION = {
    "PDAC": dict(phenotype="excluded", leukocyte=0.12,
                 mix={MACROPHAGE: 0.50, NEUTROPHIL: 0.20, BCELL: 0.10, TREG: 0.07,
                      NK: 0.03, DC: 0.02}),
    "LUAD": dict(phenotype="inflamed", leukocyte=0.30,
                 mix={MACROPHAGE: 0.35, BCELL: 0.12, NEUTROPHIL: 0.10, TREG: 0.08,
                      NK: 0.06, DC: 0.04}),
    "OV": dict(phenotype="moderate", leukocyte=0.20,
               mix={MACROPHAGE: 0.40, NEUTROPHIL: 0.15, BCELL: 0.12, TREG: 0.10,
                    NK: 0.05, DC: 0.03}),
    "BRCA": dict(phenotype="variable", leukocyte=0.22,
                 mix={MACROPHAGE: 0.38, BCELL: 0.14, TREG: 0.12, NEUTROPHIL: 0.10,
                      NK: 0.05, DC: 0.03}),
    "CESC": dict(phenotype="moderate", leukocyte=0.20,
                 mix={MACROPHAGE: 0.40, NEUTROPHIL: 0.15, BCELL: 0.11, TREG: 0.08,
                      NK: 0.05, DC: 0.03}),
}
# CD8 fraction of the immune compartment, kept separate because CD8s live in tcells.py
CD8_FRACTION = {"PDAC": 0.08, "LUAD": 0.25, "OV": 0.15, "BRCA": 0.18, "CESC": 0.18}


class ImmuneCells:
    """Macrophages, DCs, NK, B cells, Tregs and neutrophils as one agent population."""

    def __init__(self, world):
        self.w = world
        self.id = np.zeros(0, int)
        self.pos = np.zeros((0, 3))
        self.dir = np.zeros((0, 3))
        self.kind = np.zeros(0, int)
        self.born = np.zeros(0)
        self.tracks = []
        self._next_track = 0.0

    @property
    def n(self):
        return len(self.id)

    # ------------------------------------------------------------------ setup
    def _place(self, rule, k, rng):
        """Positions for k cells of one type, following its documented spatial rule."""
        w = self.w
        c = w.pos.mean(0)
        r_cells = np.linalg.norm(w.pos - c, axis=1)
        R = np.percentile(r_cells, 99)
        u = rng.normal(size=(k, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]

        if rule == "hypoxia" and w.fields["oxygen"].c is not None:
            # macrophages accumulate where oxygen is low (Lewis & Pollard 2006)
            o2 = w.fields["oxygen"].c
            p = np.clip(o2.max() - o2, 1e-6, None) ** 2
            idx = rng.choice(len(o2), size=k, p=p / p.sum())
            return w.pos[idx] + rng.normal(0, 8.0, (k, 3))
        if rule == "margin":
            return c + u * (R * rng.uniform(0.85, 1.25, k))[:, None]
        if rule == "stroma":
            return c + u * (R * rng.uniform(0.9, 1.4, k))[:, None]
        if rule == "tls":
            # tertiary lymphoid structures: a few tight clusters just outside the margin
            n_clusters = max(1, k // 25)
            cu = rng.normal(size=(n_clusters, 3))
            cu /= np.linalg.norm(cu, axis=1)[:, None]
            centers = c + cu * (R * rng.uniform(1.05, 1.3, n_clusters))[:, None]
            pick = rng.integers(0, n_clusters, k)
            return centers[pick] + rng.normal(0, 25.0, (k, 3))
        if rule == "with_cd8" and w.tcells.n:
            base = w.tcells.pos[rng.integers(0, w.tcells.n, k)]
            return base + rng.normal(0, 30.0, (k, 3))
        return c + u * (R * rng.uniform(0.5, 1.2, k))[:, None]

    def seed(self, composition=None, cd8_cells=True):
        """Procedurally generate the immune compartment for this indication."""
        w, rng = self.w, self.w.rng
        comp = composition or COMPOSITION[w.indication.key]
        total = int(round(w.n * comp["leukocyte"]))
        counts = {}
        for kind, frac in comp["mix"].items():
            if kind is None:
                continue
            k = int(round(total * frac))
            if k:
                counts[kind] = k
                self._add(self._place(TYPES[kind]["rule"], k, rng), kind)
        if cd8_cells:
            k = int(round(total * CD8_FRACTION[w.indication.key]))
            if k:
                from .tcells import NAIVE  # noqa: F401 - imported lazily to avoid a cycle
                from .params import TCELL
                w.tcells.enter(k, NAIVE, TCELL["cognate_fraction_naive"].value)
                counts["CD8"] = k
        return counts

    def _add(self, pos, kind):
        w, rng = self.w, self.w.rng
        k = len(pos)
        ids = np.arange(w.next_id, w.next_id + k)
        w.next_id += k
        d = rng.normal(size=(k, 3))
        d /= np.linalg.norm(d, axis=1)[:, None]
        self.id = np.r_[self.id, ids]
        self.pos = np.vstack([self.pos, pos])
        self.dir = np.vstack([self.dir, d])
        self.kind = np.r_[self.kind, np.full(k, kind)]
        # stagger birth times so the population does not die synchronously
        life = TYPES[kind]["lifespan"].value
        self.born = np.r_[self.born, w.t - rng.uniform(0, life, k)]
        w.events += [(w.t, f"enter_{NAMES[kind]}", int(i), -1) for i in ids]
        return ids

    def recruit(self, kind, k):
        """Bring k new cells of one type in from circulation at the margin."""
        if k:
            self._add(self._place(TYPES[kind]["rule"], k, self.w.rng), kind)

    # ------------------------------------------------------------------ dynamics
    def replenish(self, hours):
        """Recruit from circulation toward the indication's target composition.

        Keeps short-lived populations (neutrophils live 24 h) present instead of letting the
        compartment decay away. The target scales with tumour size, so recruitment tracks growth.
        """
        w, rng = self.w, self.w.rng
        comp = COMPOSITION[w.indication.key]
        total = w.n * comp["leukocyte"]
        have = {k: int((self.kind == k).sum()) for k in TYPES}
        frac = min(REPLENISH.value * hours, 1.0)
        for kind, f in comp["mix"].items():
            want = total * f
            short = want - have.get(kind, 0)
            if short <= 0:
                continue
            k = rng.poisson(short * frac)
            if k:
                self._add(self._place(TYPES[kind]["rule"], int(k), rng), kind)

    def step(self, hours):
        w, rng = self.w, self.w.rng
        self.replenish(hours)
        if self.n == 0:
            return
        # death by lifespan
        life = np.array([TYPES[k]["lifespan"].value for k in self.kind])
        dead = (w.t - self.born) > life
        if dead.any():
            w.events += [(w.t, "immune_death", int(i), -1) for i in self.id[dead]]
            self._keep(~dead)
            if self.n == 0:
                return

        speed = np.array([TYPES[k]["speed"].value for k in self.kind])
        step_um = speed * (hours * 60) * w.fibroblasts.speed_factor(self.pos)
        turn = rng.normal(0, 0.6, (self.n, 3))
        self.dir = self.dir + turn
        self.dir /= np.linalg.norm(self.dir, axis=1)[:, None]
        prop = self.pos + self.dir * step_um[:, None]
        ok = w.fibroblasts.passable(prop)      # dense matrix excludes immune cells too
        self.pos[ok] = prop[ok]
        self.dir[~ok] = rng.normal(size=((~ok).sum(), 3))
        self.dir /= np.linalg.norm(self.dir, axis=1)[:, None]

        # tissue retention, same rule as fibroblasts
        c = w.pos.mean(0)
        R_out = np.percentile(np.linalg.norm(w.pos - c, axis=1), 99) + 300.0
        off = self.pos - c
        d = np.linalg.norm(off, axis=1)
        out = d > R_out
        if out.any():
            u = off[out] / d[out][:, None]
            self.pos[out] = c + u * R_out
            self.dir[out] = -u
        self._record()

    def _keep(self, mask):
        for a in ("id", "pos", "dir", "kind", "born"):
            setattr(self, a, getattr(self, a)[mask])

    # ------------------------------------------------------------------ effects
    def suppression_at(self, pos):
        """Multiplier in [0, 1] on CD8 killing from nearby suppressive cells.

        Phenomenological: Tregs and TAMs suppress CD8 function, but the model has no
        cytokines or checkpoint ligands, so suppression is a local distance-weighted factor.
        """
        pos = np.atleast_2d(pos)
        if self.n == 0 or len(pos) == 0:
            return np.ones(len(pos))
        strength = np.array([TYPES[k]["suppress"] for k in self.kind])
        sup = strength > 0
        if not sup.any():
            return np.ones(len(pos))
        tree = cKDTree(self.pos[sup])
        s = strength[sup]
        near = tree.query_ball_point(pos, SUPPRESSION_RADIUS.value)
        return np.array([float(np.clip(1.0 - s[nb].sum() * 0.25, 0.0, 1.0)) if nb else 1.0
                         for nb in near])

    def presenting_near(self, pos, radius=30.0):
        """True where an antigen-presenting cell (DC or B cell) is within `radius`."""
        pos = np.atleast_2d(pos)
        apc = np.array([TYPES[k]["presents"] for k in self.kind]) if self.n else np.zeros(0, bool)
        if not apc.any() or len(pos) == 0:
            return np.zeros(len(pos), bool)
        tree = cKDTree(self.pos[apc])
        return np.array([len(nb) > 0 for nb in tree.query_ball_point(pos, radius)])

    def dcs_leaving(self, hours):
        """Antigen-loaded DCs that depart for the draining lymph node this step."""
        n_dc = int((self.kind == DC).sum())
        if not n_dc:
            return 0
        k = self.w.rng.poisson(n_dc * DC_EGRESS.value * hours)
        k = min(int(k), n_dc)
        if k:  # they physically leave the tumour
            idx = np.where(self.kind == DC)[0][:k]
            keep = np.ones(self.n, bool)
            keep[idx] = False
            self._keep(keep)
        return k

    def counts(self):
        return {NAMES[k]: int((self.kind == k).sum()) for k in TYPES}

    def _record(self):
        w = self.w
        if w.t >= self._next_track - 1e-9:
            self.tracks.append((w.t, self.id.copy(), self.pos.copy(), self.kind.copy()))
            self._next_track = w.t + 2.0 / 60


def table():
    rows = []
    for kind, d in TYPES.items():
        for key in ("radius", "speed", "lifespan"):
            p = d[key]
            rows.append(("immune", f"{NAMES[kind]}.{key}", p.value, p.unit, p.status, p.source))
        rows.append(("immune", f"{NAMES[kind]}.suppress", d["suppress"], "-", "assumption",
                     "local reduction of CD8 killing; phenomenological"))
    for k, p in MONOCYTE.items():
        rows.append(("immune", f"monocyte.{k}", p.value, p.unit, p.status, p.source))
    rows.append(("immune", "dc_egress", DC_EGRESS.value, DC_EGRESS.unit,
                 DC_EGRESS.status, DC_EGRESS.source))
    rows.append(("immune", "suppression_radius", SUPPRESSION_RADIUS.value, "um",
                 SUPPRESSION_RADIUS.status, SUPPRESSION_RADIUS.source))
    for ind, c in COMPOSITION.items():
        rows.append(("immune", f"{ind}.leukocyte_fraction", c["leukocyte"], "-", "assumption",
                     f"immune cells per tumour cell; phenotype '{c['phenotype']}'; "
                     "reproduces published ordering, not a measurement"))
    return rows
