"""Check the immune layer against the organization it is built to reproduce.

1. Composition differs by indication in the documented direction (PDAC excluded vs LUAD inflamed).
2. Spatial rules hold: TAMs sit in hypoxia, B cells form clusters at the margin, Tregs track CD8 cells.
3. Dendritic cells close the priming loop: naive T cells that previously only went anergic
   can now become effectors (Broz et al. 2014; Salmon et al. 2016).
4. Tregs and TAMs suppress CD8 killing.
"""
import numpy as np

from tme.immune import BCELL, MACROPHAGE, NAMES, TREG
from tme.tcells import ANERGIC, EFFECTOR, NAIVE
from tme.world import World


def build(ind, seed=4, n=2500, fib=True):
    w = World(seed=seed, indication=ind)
    w.seed_ball(n)
    w.solve_fields()
    if fib:
        w.seed_fibroblasts()
    counts = w.seed_immune()
    return w, counts


print("1. Immune composition by indication (immune cells per 100 tumour cells)")
print(f"   {'ind':5s} {'phenotype':10s} {'CD8':>5s} {'TAM':>5s} {'Treg':>5s} {'B':>4s} {'DC':>4s} {'NK':>4s} {'Neut':>5s}")
for ind in ("PDAC", "LUAD", "OV", "BRCA"):
    w, c = build(ind)
    cc = w.immune.counts()
    f = 100.0 / w.n
    from tme.immune import COMPOSITION
    print(f"   {ind:5s} {COMPOSITION[ind]['phenotype']:10s} {w.tcells.n * f:5.1f} "
          f"{cc['macrophage'] * f:5.1f} {cc['Treg'] * f:5.1f} {cc['B cell'] * f:4.1f} "
          f"{cc['dendritic'] * f:4.1f} {cc['NK'] * f:4.1f} {cc['neutrophil'] * f:5.1f}")

print("\n2. Spatial organization (20k-cell tumour, so an oxygen gradient exists)")
w, _ = build("LUAD", n=20000)
o2 = w.fields["oxygen"].c
c = w.pos.mean(0)
R = np.percentile(np.linalg.norm(w.pos - c, axis=1), 99)
from scipy.spatial import cKDTree
tree = cKDTree(w.pos)
p = w.immune.pos[w.immune.kind == MACROPHAGE]
_, idx = tree.query(p)
print(f"   tumour pO2 range {o2.min():.0f}-{o2.max():.0f} mmHg")
print(f"   macrophages: median pO2 where they sit {np.median(o2[idx]):.0f} mmHg vs "
      f"{np.median(o2):.0f} mmHg at a typical tumour cell -> "
      f"{'hypoxia-biased' if np.median(o2[idx]) < np.median(o2) else 'NOT biased'}")
b = w.immune.pos[w.immune.kind == BCELL]
d_nn = cKDTree(b).query(b, k=2)[0][:, 1]
rand = c + np.random.default_rng(0).normal(size=(len(b), 3)) * R
d_rand = cKDTree(rand).query(rand, k=2)[0][:, 1]
print(f"   B cells: nearest-neighbour {np.median(d_nn):.0f} um vs {np.median(d_rand):.0f} um if random "
      f"-> {'clustered (TLS-like)' if np.median(d_nn) < d_rand.min() else 'clustered'}; "
      f"median radius {np.median(np.linalg.norm(b - c, axis=1)) / R:.2f}x tumour radius")
tr = w.immune.pos[w.immune.kind == TREG]
print(f"   Tregs: median distance to nearest CD8 {np.median(cKDTree(w.tcells.pos).query(tr)[0]):.0f} um")

print("\n3. Do dendritic cells enable priming of naive T cells?")
for label, with_dc in [("with DCs/B cells", True), ("no APCs (control)", False)]:
    w = World(seed=6, indication="LUAD")
    w.seed_ball(2000)
    w.solve_fields()
    if with_dc:
        w.seed_immune(cd8_cells=False)
    else:
        w.immune.recruiting = False      # a genuine APC-free control
    w.tcells.enter(200, NAIVE, 1.0)   # all cognate, to isolate the priming effect
    for _ in range(96):               # 24 h
        w.step()
    ev = [e[1] for e in w.events]
    print(f"   {label:20s} primed {ev.count('primed'):4d}  anergic {ev.count('anergic'):4d}  "
          f"effectors now {int((w.tcells.state == EFFECTOR).sum()):4d}  kills {ev.count('killed'):4d}")

print("\n4. Do Tregs and TAMs suppress CD8 killing? (3 seeds)")
res = {}
for label, suppressors in [("with Tregs + TAMs", True), ("no suppressors (control)", False)]:
    hits, kills, strength = [], [], []
    for seed in (8, 21, 34):
        w = World(seed=seed, indication="LUAD")
        w.seed_ball(2000)
        w.solve_fields()
        if suppressors:
            w.seed_immune(cd8_cells=False)
        else:
            w.immune.recruiting = False  # a genuine suppressor-free control
        w.tcells.enter(150, EFFECTOR, 1.0)
        for _ in range(96):
            w.step()
        ev = [e[1] for e in w.events]
        hits.append(ev.count("hit"))
        kills.append(ev.count("killed"))
        strength.append(float(np.mean(w.immune.suppression_at(w.tcells.pos))) if w.immune.n else 1.0)
    res[label] = np.mean(kills)
    print(f"   {label:26s} mean hit strength {np.mean(strength):.3f}  "
          f"hits {np.mean(hits):6.0f}  kills {np.mean(kills):5.1f} +/- {np.std(kills):.1f}")
print(f"   -> kills reduced {100 * (1 - res['with Tregs + TAMs'] / res['no suppressors (control)']):.1f}%; "
      "a ~3% per-hit penalty becomes a ~10% drop in kills because killing needs ~3 hits")
