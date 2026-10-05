"""T cells: motile agents that enter the world, search tumor tissue and (only when activated)
kill tumor cells by perforin/granzyme hits.

States
  NAIVE     - migrates as a random walker; cognate contacts give signal 1 only.
  ANERGIC   - signal 1 without costimulation for long enough -> hyporesponsive; keeps moving.
  EFFECTOR  - activated CTL; arrests on cognate tumor cells and delivers hits.
  EXITED    - wandered beyond the tissue boundary and left the world (logged as "exit").
(Priming of naive cells by dendritic cells / macrophages comes with those cell types.)

T cells move in 1-min substeps between the tumor's 15-min steps. Tumor cells are
obstacles (T cells may squeeze 20% into gaps); T cells don't push tumor cells, consume
nutrients or divide yet.
"""
import numpy as np
from scipy.spatial import cKDTree

from . import params as P

NAIVE, ANERGIC, EFFECTOR, EXITED = 0, 1, 3, 4
STATE_NAMES = {NAIVE: "naive", ANERGIC: "anergic", EFFECTOR: "effector", EXITED: "exited"}
TC = P.TCELL


def _radius(state):
    r = TC["diameter_naive"].value / 2
    return np.where(state == EFFECTOR, r * TC["blast_volume_factor"].value ** (1 / 3), r)


class TCells:
    def __init__(self, world):
        self.w = world
        self.id = np.zeros(0, int)
        self.pos = np.zeros((0, 3))
        self.dir = np.zeros((0, 3))
        self.state = np.zeros(0, int)
        self.cognate = np.zeros(0, bool)
        self.stim_h = np.zeros(0)       # cumulative cognate contact without costimulation
        self.target = np.zeros(0, int)  # tumor-cell id in contact, -1 = none
        self.contact_end = np.zeros(0)  # hours
        self.kills = np.zeros(0, int)
        self.tracks = []                # (t, ids, pos, state), every track_every minutes
        self._next_track = 0.0

    @property
    def n(self):
        return len(self.id)

    def enter(self, k, state, cognate_p):
        """k T cells arrive at the tumor margin (extravasation / arrival from the medium)."""
        w, rng = self.w, self.w.rng
        c = w.pos.mean(0)
        R = np.percentile(np.linalg.norm(w.pos - c, axis=1), 99) + 15.0
        u = rng.normal(size=(k, 3))
        u /= np.linalg.norm(u, axis=1)[:, None]
        ids = np.arange(w.next_id, w.next_id + k)
        w.next_id += k
        self.id = np.r_[self.id, ids]
        self.pos = np.vstack([self.pos, c + u * R])
        self.dir = np.vstack([self.dir, -u])  # start heading inward
        self.state = np.r_[self.state, np.full(k, state)]
        self.cognate = np.r_[self.cognate, rng.random(k) < cognate_p]
        self.stim_h = np.r_[self.stim_h, np.zeros(k)]
        self.target = np.r_[self.target, np.full(k, -1)]
        self.contact_end = np.r_[self.contact_end, np.zeros(k)]
        self.kills = np.r_[self.kills, np.zeros(k, int)]
        w.events += [(w.t, "enter", int(i), -1) for i in ids]

    # ------------------------------------------------------------------ dynamics
    def step(self, hours):
        """Advance T cells by `hours` in 1-min substeps (tumor cells held fixed meanwhile)."""
        if self.n == 0:
            return
        w, rng = self.w, self.w.rng
        dt_min = TC["substep"].value
        r_tumor = P.NUMERICS["packing"].value * np.cbrt(3 * w.vol / (4 * np.pi))
        tree, rmax = cKDTree(w.pos), r_tumor.max()
        sigma = np.deg2rad(TC["turn_angle_median"].value) / 0.6745  # |N(0,s)| has median 0.6745 s
        tau_damage = TC["damage_recovery_median"].value / np.log(2)  # min
        center = w.pos.mean(0)
        R_out = np.percentile(np.linalg.norm(w.pos - center, axis=1), 99) + TC["world_margin"].value
        t0 = w.t
        for sub in range(int(round(hours * 60 / dt_min))):
            t = t0 + sub * dt_min / 60
            w.damage *= np.exp(-dt_min / tau_damage)
            r_T = _radius(self.state)
            speed = np.where(self.state == EFFECTOR, TC["speed_effector"].value, TC["speed_naive"].value)
            # dense matrix slows T cells (Salmon et al. 2012 doi:10.1172/JCI45817)
            speed = speed * w.fibroblasts.speed_factor(self.pos)

            # end finished contacts
            done = (self.target >= 0) & (t >= self.contact_end)
            self.target[done] = -1
            free = (self.target < 0) & (self.state != EXITED)

            # persistent random walk: turn by a random angle about a random axis
            turn = np.abs(rng.normal(0, sigma, self.n))
            axis = np.cross(self.dir, rng.normal(size=(self.n, 3)))
            axis /= np.maximum(np.linalg.norm(axis, axis=1), 1e-9)[:, None]
            self.dir = (self.dir * np.cos(turn)[:, None] + np.cross(axis, self.dir) * np.sin(turn)[:, None])
            self.dir /= np.linalg.norm(self.dir, axis=1)[:, None]
            prop = self.pos + self.dir * (speed * dt_min)[:, None]

            # dense matrix is a barrier: T cells are kept out of it entirely
            blocked = ~w.fibroblasts.passable(prop)
            for k in np.where(free)[0]:
                if blocked[k]:
                    self.dir[k] = rng.normal(size=3)
                    self.dir[k] /= np.linalg.norm(self.dir[k])
                    continue
                nb = tree.query_ball_point(prop[k], r_T[k] + rmax)
                if nb:
                    d = np.linalg.norm(w.pos[nb] - prop[k], axis=1)
                    if np.any(d < TC["squeeze"].value * (r_T[k] + r_tumor[nb])):
                        self.dir[k] = rng.normal(size=3)  # blocked: pick a new heading
                        self.dir[k] /= np.linalg.norm(self.dir[k])
                        continue
                self.pos[k] = prop[k]

            # leaving the tissue
            out = free & (np.linalg.norm(self.pos - center, axis=1) > R_out)
            self.state[out] = EXITED
            w.events += [(t, "exit", int(i), -1) for i in self.id[out]]
            free &= ~out

            # contacts with live tumor cells
            for k in np.where(free)[0]:
                nb = tree.query_ball_point(self.pos[k], r_T[k] + rmax + 1.0)
                live = [j for j in nb if w.state[j] == 0 and
                        np.linalg.norm(w.pos[j] - self.pos[k]) < 1.1 * (r_T[k] + r_tumor[j])]
                if not live or not self.cognate[k]:
                    continue
                j = live[rng.integers(len(live))]
                if self.state[k] == NAIVE:
                    self.stim_h[k] += dt_min / 60  # signal 1 from the tumour cell
                    # signal 2 requires an antigen-presenting cell nearby (DC or B cell);
                    # with one, the naive cell is primed instead of going anergic
                    if w.immune.presenting_near(self.pos[k])[0]:
                        self.state[k] = EFFECTOR
                        w.events.append((t, "primed", int(self.id[k]), int(w.id[j])))
                    elif not TC["tumor_costimulation"].value and self.stim_h[k] >= TC["anergy_signal1_h"].value:
                        self.state[k] = ANERGIC
                        w.events.append((t, "anergic", int(self.id[k]), int(w.id[j])))
                elif self.state[k] == EFFECTOR:
                    self._hit(k, j, t, rng)
            self._record(t)
        self._record(t0 + hours, force=False)

    def _hit(self, k, j, t, rng):
        """Effector CTL engages tumor cell j: arrest, perforin pore + granzyme delivery (one hit)."""
        w = self.w
        dur = TC["contact_median"].value * np.exp(rng.normal(0, 0.5))  # lognormal, median 15 min
        self.target[k] = w.id[j]
        self.contact_end[k] = t + dur / 60
        # Tregs and TAMs nearby suppress the hit (phenomenological; see immune.py)
        w.damage[j] += float(w.immune.suppression_at(self.pos[k])[0])
        w.events.append((t, "hit", int(self.id[k]), int(w.id[j])))
        if w.damage[j] >= TC["hits_to_kill"].value or rng.random() < TC["single_hit_lethal"].value:
            w.state[j] = w.APOPTOTIC
            w.t_dead[j] = t + TC["apoptosis_onset"].value / 60
            self.kills[k] += 1
            w.events.append((t, "killed", int(w.id[j]), int(self.id[k])))

    def _record(self, t=None, force=False):
        t = self.w.t if t is None else t
        if force or t >= self._next_track - 1e-9:
            self.tracks.append((t, self.id.copy(), self.pos.copy(), self.state.copy()))
            self._next_track = t + TC["track_every"].value / 60
