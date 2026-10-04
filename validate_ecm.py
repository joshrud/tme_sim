"""Check the fibroblast/ECM layer against the biology it is built from.

1. ECM is *stored*, not recomputed: with no fibroblasts it only decays, it does not spread.
2. One calibrated deposition rate reproduces the measured ordering of stromal density
   across indications, purely from differing fibroblast numbers.
3. The functional consequence: dense PDAC stroma slows T cells and reduces killing,
   the model's analogue of matrix-mediated immune exclusion
   (Salmon et al. 2012, J Clin Invest, doi:10.1172/JCI45817).
"""
import numpy as np

from tme.fibroblasts import FIBRO, FRACTION
from tme.tcells import EFFECTOR
from tme.world import World

DAYS = 7

print("1. ECM is stored state: with no fibroblasts it decays but does not spread")
w = World(seed=3, indication="PDAC")
w.seed_ball(1500)
w.seed_fibroblasts()
for _ in range(48):
    w.fibroblasts.step(1.0)
before = w.fibroblasts.ecm.rho.copy()
for a in ("id", "pos", "dir", "state", "age"):
    setattr(w.fibroblasts, a, getattr(w.fibroblasts, a)[:0])
for _ in range(24):
    w.fibroblasts.step(1.0)
expected = before * np.exp(-FIBRO["ecm_decay"].value * 24)
print(f"   max before {before.max():.3f} -> after 24 h {w.fibroblasts.ecm.rho.max():.3f}; "
      f"pure decay: {np.allclose(w.fibroblasts.ecm.rho, expected, atol=1e-6)}")

print(f"\n2. Stromal density after {DAYS} days, one calibrated rate for all indications")
print(f"   {'ind':5s} {'fibroblasts':>11s} {'ECM@edge':>9s} {'tumor in dense ECM':>19s}")
for ind in ("PDAC", "OV", "BRCA", "LUAD"):
    w = World(seed=3, indication=ind)
    w.seed_ball(2000)
    nf = w.seed_fibroblasts()
    for _ in range(DAYS * 24):
        w.fibroblasts.step(1.0)
    c = w.pos.mean(0)
    r = np.linalg.norm(w.pos - c, axis=1)
    R = np.percentile(r, 99)
    e = w.fibroblasts.ecm_at(w.pos)
    print(f"   {ind:5s} {nf:11d} {np.median(w.fibroblasts.ecm_at(w.pos[r > 0.8 * R])):9.3f} "
          f"{100 * (e > 0.7).mean():18.1f}%")

print(f"\n3. Does dense stroma exclude T cells? PDAC, {DAYS} d of matrix build-up, then CTLs enter")
print("   (matrix acts as a barrier, not just drag - see docs/fibroblasts.md)")
res = {}
for label, with_fib in [("with fibroblasts", True), ("no fibroblasts (control)", False)]:
    w = World(seed=9, indication="PDAC")
    w.seed_ball(2000)
    if with_fib:
        w.seed_fibroblasts()
    for _ in range(DAYS * 24):      # matrix build-up, tumor held fixed
        w.fibroblasts.step(1.0)
    w.tcells.enter(150, EFFECTOR, 1.0)
    for _ in range(96):             # 24 h of killing
        w.step()
    ev = [x[1] for x in w.events]
    res[label] = (ev.count("hit"), ev.count("killed"))
    e = w.fibroblasts.ecm_at(w.tcells.pos)
    print(f"   {label:26s} ECM at T cells {np.median(e):.2f}  "
          f"hits {ev.count('hit'):5d}  kills {ev.count('killed'):4d}")
h0, k0 = res["no fibroblasts (control)"]
h1, k1 = res["with fibroblasts"]
print(f"   -> stroma reduces contacts by {100 * (1 - h1 / h0):.0f}% and kills by {100 * (1 - k1 / k0):.0f}%")
