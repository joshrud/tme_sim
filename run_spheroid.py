"""Grow a HeLa spheroid to ~700 um and record a 3D+time trajectory.

Output (out/):
  snapshots.npz  - positions, phase, state and every field per cell, every 12 h
  events.csv     - every birth, division, necrosis (later: entry / exit)
  log.csv        - population summary every 6 h
"""
import csv
import os
import sys
import time

import numpy as np

from tme.world import World

TARGET_DIAMETER = float(sys.argv[1]) if len(sys.argv) > 1 else 700.0
os.makedirs("out", exist_ok=True)

w = World(seed=42)
w.seed_ball(1000)
snaps, rows = [], []
t0 = time.time()
step = 0
while w.diameter() < TARGET_DIAMETER and w.t < 40 * 24:
    w.step()
    step += 1
    if step % 24 == 0:  # 6 h
        ph = w.phase()
        rows.append(dict(t_h=w.t, n=w.n, diameter_um=round(w.diameter(), 1),
                         frac_G1=np.mean(ph == 0), frac_SG2M=np.mean(ph == 1),
                         frac_necrotic=np.mean(ph == -1), frac_quiescent=w.quiescent.mean(),
                         min_o2_mmHg=w.fields["oxygen"].c.min(),
                         max_lactate_mM=w.fields["lactate"].c.max(),
                         wall_s=round(time.time() - t0, 1)))
        print(rows[-1], flush=True)
    if step % 48 == 0:  # 12 h
        snaps.append(w.snapshot())

w.solve_fields()
snaps.append(w.snapshot())

np.savez_compressed("out/snapshots.npz", **{f"{k}_{i:03d}": v for i, s in enumerate(snaps)
                                             for k, v in s.items()})
with open("out/log.csv", "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=rows[0].keys())
    wr.writeheader()
    wr.writerows(rows)
with open("out/events.csv", "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["t_h", "event", "cell_id", "parent_id"])
    wr.writerows(w.events)
print(f"done: t={w.t:.1f} h, n={w.n}, diameter={w.diameter():.0f} um, "
      f"{len(w.events)} events, wall {time.time() - t0:.0f} s")
