"""Check the per-division mutation model against the literature it is built from.

1. Mutations per division match the Werner et al. 2020 rate for each indication.
2. Driver-hit frequency matches the analytic expectation mu * p_driver * n_genes.
3. Mutation counts per cell grow linearly with generation number (a division-coupled process).
4. Driver-carrying subclones expand (more cells carry drivers than there were driver events).
"""
import numpy as np

from tme.indications import INDICATIONS, MUT_PER_DIVISION_HEALTHY, MutationModel
from tme.world import World

print("1. Mutations per division (Werner et al. 2020: healthy 1.14, tumors 4-100x)")
rng = np.random.default_rng(0)
for k, ind in INDICATIONS.items():
    m = MutationModel(k, rng)
    n_mut, drv = m.founder(20000)
    n_mut, drv = m.on_division(n_mut, drv)
    print(f"   {k:5s} expected {ind.mut_per_division:6.2f}  simulated {n_mut.mean():6.2f} "
          f"(factor {ind.tumor_factor:g}x healthy {MUT_PER_DIVISION_HEALTHY.value})")

print("\n2. Driver-hit rate per division vs analytic expectation")
for k, ind in INDICATIONS.items():
    m = MutationModel(k, rng)
    N = 400000
    n_mut, drv = m.on_division(*m.founder(N))
    obs = (drv != 0).mean()
    exp = ind.mut_per_division * m.p_driver * len(m.genes)
    print(f"   {k:5s} expected {exp:.3e}  simulated {obs:.3e}  ratio {obs / exp:.2f}")

print("\n3. Mutation count scales with generation (division-coupled)")
m = MutationModel("LUAD", rng)
n_mut, drv = m.founder(5000)
for gen in range(1, 6):
    n_mut, drv = m.on_division(n_mut, drv)
    print(f"   generation {gen}: mean {n_mut.mean():6.1f}  "
          f"(expected {gen * m.ind.mut_per_division:6.1f})")

print("\n4. Driver subclones expand in a growing tumor (LUAD, 96 h)")
w = World(seed=11, indication="LUAD")
w.seed_ball(1000)
while w.t < 96 and w.n < 60000:
    w.step()
nd = w.mutations.n_drivers(w.driver_mask)
events = sum(1 for e in w.events if e[1] == "driver")
divisions = sum(1 for e in w.events if e[1] == "division")
genes = {}
for mask in w.driver_mask[nd > 0]:
    for g in w.mutations.gene_names(mask):
        genes[g] = genes.get(g, 0) + 1
print(f"   {divisions} divisions -> {events} driver events -> {(nd > 0).sum()} cells carry a driver")
print(f"   expected driver events ~{2 * divisions * w.indication.mut_per_division * w.mutations.p_driver * len(w.mutations.genes):.1f}"
      " (2 daughters per division)")
print(f"   median mutations/cell {np.median(w.n_mut):.0f}; genes hit {genes}")
print(f"   subclonal expansion: {(nd > 0).sum() / max(events, 1):.1f} cells per driver event")
