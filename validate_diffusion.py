"""Check graph diffusion against the analytic O2 profile in a uniformly consuming sphere.

Zero-order uptake in a sphere of radius R with surface value c_s:
    c(r) = c_s - Q / (6 D) * (R^2 - r^2)
"""
import numpy as np

from tme import params as P
from tme.world import World, radius, UM3_TO_L

for n_cells in (2000, 8000, 30000):
    w = World(seed=1)
    w.seed_ball(n_cells)
    for k in ("glucose", "lactate", "mitogen"):
        del w.fields[k]
    f = w.fields["oxygen"]
    f.Km = 1e-6  # ~zero-order uptake
    w.solve_fields()

    r = np.linalg.norm(w.pos - w.pos.mean(0), axis=1)
    R = np.percentile(r, 99) + P.NUMERICS["packing"].value * radius(w.vol.mean())
    n_per_L = w.n / (4 / 3 * np.pi * R**3 * UM3_TO_L)
    Q = f.vmax * n_per_L  # mmHg/s
    analytic = f.bath - Q / (6 * f.D) * (R**2 - r**2)
    inner = r < 0.8 * R  # skip the ragged surface layer
    err = f.c[inner] - analytic[inner]
    drop = f.bath - analytic.min()
    print(f"N={w.n:6d} R={R:5.0f} um vol.frac={w.vol.sum() / (4/3*np.pi*R**3):.2f} "
          f"mean degree={w.degree()[inner].mean():.1f} interior cells flagged surface="
          f"{w.surface[inner].mean():.3f}")
    print(f"   center drop: analytic {drop:5.1f}  graph {f.bath - f.c[inner].min():5.1f} mmHg;"
          f"  interior RMSE {np.sqrt(np.mean(err**2)):4.1f} mmHg "
          f"({100 * np.sqrt(np.mean(err**2)) / drop:.0f}% of drop)")
