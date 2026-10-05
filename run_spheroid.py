"""Grow a HeLa spheroid and record a 3D+time trajectory, optionally with T cells entering.

  python run_spheroid.py                                   # tumor only, to 700 um
  python run_spheroid.py --indication PDAC --out out/pdac  # another indication
  python run_spheroid.py --tcells naive --out out/naive    # naive T cells arrive from day 3
  python run_spheroid.py --tcells effector --out out/ctl   # activated CTLs (e.g. adoptive transfer)

Output (--out, default out/):
  snapshots.npz  - positions, phase, state and every field per tumor cell, every --snap-every h
  tcells.npz     - T-cell tracks (every 2 min) and per-T-cell info
  events.csv     - every birth, division, necrosis, T-cell entry, hit, kill, anergy, clearance
  log.csv        - population summary every 6 h
"""
import argparse
import csv
import json
import os
import time

import numpy as np

from tme import params as P
from tme.tcells import EFFECTOR, NAIVE
from tme.fibroblasts import FRACTION as FIBRO_FRACTION
from tme.indications import DEFAULT as DEFAULT_INDICATION, INDICATIONS
from tme.organs import NODE_BASIN, NODE_DISTANCE_CM
from tme.world import World

ap = argparse.ArgumentParser()
ap.add_argument("diameter", nargs="?", type=float, default=700.0, help="stop at this diameter (um)")
ap.add_argument("--tcells", choices=["none", "naive", "effector"], default="none")
ap.add_argument("--t-start", type=float, default=72.0, help="hour T cells start arriving")
ap.add_argument("--t-rate", type=float, default=10.0, help="T cells arriving per hour")
ap.add_argument("--t-end", type=float, default=None, help="stop at this hour (overrides diameter)")
ap.add_argument("--snap-every", type=float, default=12.0, help="hours between tumor snapshots")
ap.add_argument("--indication", default=DEFAULT_INDICATION, choices=sorted(INDICATIONS),
                help="tumor type: sets reference cell line, doubling time, drivers, mutation rate")
ap.add_argument("--fibroblasts", action="store_true",
                help="seed fibroblasts at the indication's fraction and let them deposit ECM")
ap.add_argument("--immune", action="store_true",
                help="procedurally generate this indication's immune compartment")
ap.add_argument("--out", default="out")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

dt = P.NUMERICS["dt"].value
snap_steps = max(1, round(args.snap_every / dt))
w = World(seed=42, indication=args.indication)
print(f"indication: {w.indication.name} (line {w.indication.line}, "
      f"doubling {w.indication.doubling_h:g} h, {w.indication.mut_per_division:.1f} mutations/division)")
w.seed_ball(1000)
if args.fibroblasts:
    print(f"fibroblasts: seeded {w.seed_fibroblasts()} "
          f"({100 * FIBRO_FRACTION[w.indication.key]:.0f}% of tumor cells)")
if args.immune:
    print(f"immune: {w.seed_immune()}")
snaps, rows, organ_log = [], [], []
t0 = time.time()
step = 0


def done():
    if args.t_end is not None:
        return w.t >= args.t_end
    return w.diameter() >= args.diameter or w.t >= 40 * 24


while not done():
    if args.tcells != "none" and w.t >= args.t_start:
        k = w.rng.poisson(args.t_rate * dt)
        if k:
            naive = args.tcells == "naive"
            w.tcells.enter(k, NAIVE if naive else EFFECTOR,
                           P.TCELL["cognate_fraction_naive"].value if naive else 1.0)
    w.step()
    step += 1
    if step % 24 == 0:  # 6 h
        ph, tc = w.phase(), w.tcells
        ev = [e[1] for e in w.events]
        rows.append(dict(t_h=w.t, n=w.n, diameter_um=round(w.diameter(), 1),
                         frac_G1=np.mean(ph == 0), frac_SG2M=np.mean(ph == 1),
                         frac_necrotic=np.mean(ph == -1), frac_apoptotic=np.mean(ph == -2),
                         frac_quiescent=w.quiescent.mean(),
                         median_mut=float(np.median(w.n_mut)),
                         cells_with_driver=int((w.mutations.n_drivers(w.driver_mask) > 0).sum()),
                         fibroblasts=w.fibroblasts.n, immune=w.immune.n,
                         **{f"n_{k.replace(' ', '_')}": v for k, v in w.immune.counts().items()},
                         ecm_mean=round(float(np.mean(w.fibroblasts.ecm_at(w.pos))), 4) if w.fibroblasts.ecm is not None else 0.0,
                         tcells=tc.n, t_anergic=int(np.sum(tc.state == 1)),
                         hits=ev.count("hit"), kills=ev.count("killed"),
                         min_o2_mmHg=w.fields["oxygen"].c.min(),
                         max_lactate_mM=w.fields["lactate"].c.max(),
                         wall_s=round(time.time() - t0, 1)))
        print(rows[-1], flush=True)
    if step % snap_steps == 0:
        snaps.append(w.snapshot())
        organ_log.append((w.t, w.organs.table()))

w.solve_fields()
snaps.append(w.snapshot())

np.savez_compressed(f"{args.out}/snapshots.npz", **{f"{k}_{i:03d}": v for i, s in enumerate(snaps)
                                                    for k, v in s.items()})
tr = w.tcells.tracks
np.savez_compressed(
    f"{args.out}/tcells.npz",
    t=np.array([x[0] for x in tr]), count=np.array([len(x[1]) for x in tr]),
    ids=np.concatenate([x[1] for x in tr]) if tr else np.zeros(0, int),
    pos=np.concatenate([x[2] for x in tr]).astype(np.float32) if tr else np.zeros((0, 3), np.float32),
    state=np.concatenate([x[3] for x in tr]).astype(np.int8) if tr else np.zeros(0, np.int8),
    info_id=w.tcells.id, info_cognate=w.tcells.cognate, info_kills=w.tcells.kills)
fb = w.fibroblasts
np.savez_compressed(f"{args.out}/fibroblasts.npz",
                    pos=fb.pos.astype(np.float32), state=fb.state.astype(np.int8), id=fb.id,
                    ecm=fb.ecm.rho if fb.ecm is not None else np.zeros((1, 1, 1), np.float32),
                    ecm_origin=fb.ecm.origin if fb.ecm is not None else np.zeros(3),
                    ecm_h=np.array([fb.ecm.h if fb.ecm is not None else 0.0]))
np.savez_compressed(f"{args.out}/immune.npz", pos=w.immune.pos.astype(np.float32),
                    kind=w.immune.kind.astype(np.int8), id=w.immune.id)
organ_log.append((w.t, w.organs.table()))
with open(f"{args.out}/organs.json", "w") as fh:
    json.dump({"age": w.organs.age,
               "indication": w.indication.key,
               "node_basin": NODE_BASIN.get(w.indication.key, ""),
               "node_distance_cm": NODE_DISTANCE_CM.get(w.indication.key, 0),
               "dc_transit_h": round(w.organs.dc_transit_time(), 2),
               "precursor_frequency": w.organs.precursor_frequency,
               "t": [t for t, _ in organ_log],
               "rows": [[[o, st, float(n)] for o, st, n in tbl] for _, tbl in organ_log]}, fh)
with open(f"{args.out}/indication.txt", "w") as fh:
    fh.write(f"{w.indication.key}\t{w.indication.name}\t{w.indication.line}\n"
             f"doubling_h\t{w.indication.doubling_h}\n"
             f"mut_per_division\t{w.indication.mut_per_division:.2f}\n"
             f"truncal\t{','.join(w.indication.truncal)}\n"
             f"sources\t{w.indication.sources}\n")
with open(f"{args.out}/log.csv", "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
    wr.writeheader()
    wr.writerows(rows)
with open(f"{args.out}/events.csv", "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["t_h", "event", "cell_id", "parent_id"])  # parent_id = other cell (parent, T cell, target)
    wr.writerows(w.events)
print(f"done: t={w.t:.1f} h, n={w.n}, T cells={w.tcells.n}, diameter={w.diameter():.0f} um, "
      f"{len(w.events)} events, wall {time.time() - t0:.0f} s")
