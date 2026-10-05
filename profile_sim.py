"""Profile a full TME step to find where the time actually goes.

  python profile_sim.py [n_seed] [steps]
"""
import cProfile
import io
import pstats
import sys
import time

from tme.world import World

n_seed = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
steps = int(sys.argv[2]) if len(sys.argv) > 2 else 8

w = World(seed=5, indication="LUAD")
w.seed_ball(n_seed)
w.solve_fields()
w.seed_fibroblasts()
w.seed_immune()
print(f"tumour {w.n}, fibroblasts {w.fibroblasts.n}, immune {w.immune.n}, T cells {w.tcells.n}")

# per-component wall time
comp = {}
for name, fn in [("solve_fields", w.solve_fields), ("build_graph", w.build_graph),
                 ("relax", w.relax)]:
    t0 = time.perf_counter()
    for _ in range(3):
        fn()
    comp[name] = (time.perf_counter() - t0) / 3
for name, fn in [("fibroblasts.step", lambda: w.fibroblasts.step(0.25)),
                 ("immune.step", lambda: w.immune.step(0.25)),
                 ("tcells.step", lambda: w.tcells.step(0.25)),
                 ("organs.step", lambda: w.organs.step(0.25))]:
    t0 = time.perf_counter()
    for _ in range(3):
        fn()
    comp[name] = (time.perf_counter() - t0) / 3
print("\nper-call wall time (s):")
for k, v in sorted(comp.items(), key=lambda x: -x[1]):
    print(f"  {k:20s} {v:8.3f}")

pr = cProfile.Profile()
pr.enable()
t0 = time.perf_counter()
for _ in range(steps):
    w.step()
wall = time.perf_counter() - t0
pr.disable()
print(f"\n{steps} full steps: {wall:.1f} s ({wall / steps:.2f} s/step) at {w.n} tumour cells")

s = io.StringIO()
pstats.Stats(pr, stream=s).sort_stats("cumulative").print_stats(22)
print("\n".join(s.getvalue().splitlines()[4:32]))
