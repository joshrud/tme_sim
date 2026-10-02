"""Compare the final simulated spheroid with HeLa-Fucci spheroid observations
(Onozato et al. 2017, Cancer Sci, doi:10.1111/cas.13178; ~700 um spheroids):
  - S/G2/M (Fucci green) cells confined to an outer rim ~70 um thick
  - interior G1 (red), HIF-1a positive but pimonidazole negative (mildly hypoxic)
"""
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

z = np.load("out/snapshots.npz")
last = max(int(k.split("_")[-1]) for k in z.files)
s = {k.rsplit("_", 1)[0]: z[k] for k in z.files if k.endswith(f"_{last:03d}")}
log = list(csv.DictReader(open("out/log.csv")))

pos, ph = s["pos"], s["phase"]
r = np.linalg.norm(pos - pos.mean(0), axis=1)
R = np.percentile(r, 99) + 7.0
depth = R - r

bins = np.arange(0, R + 20, 20)
mid, green, o2, mit, lac = [], [], [], [], []
for a, b in zip(bins[:-1], bins[1:]):
    m = (depth >= a) & (depth < b)
    if m.sum() < 20:
        continue
    mid.append((a + b) / 2)
    green.append(np.mean(ph[m] == 1))
    o2.append(s["oxygen"][m].mean())
    mit.append(s["mitogen"][m].mean())
    lac.append(s["lactate"][m].mean())
mid, green = np.array(mid), np.array(green)
half = green[0] / 2
rim = mid[np.argmax(green < half)] if np.any(green < half) else np.nan
core_o2 = s["oxygen"][depth > R - 30]

print(f"time {s['t']:.0f} h ({s['t'] / 24:.1f} d), cells {len(ph)}, diameter {2 * R:.0f} um")
print(f"S/G2/M rim (green fraction falls to half its surface value): {rim:.0f} um  [HeLa-Fucci: ~70 um]")
print(f"core pO2: {core_o2.mean():.1f} mmHg (min {s['oxygen'].min():.1f})  "
      f"[HIF-1a+ / pimonidazole- implies roughly 10-40 mmHg]")
print(f"necrotic fraction: {np.mean(ph == -1):.3f}  [no necrotic core reported]")
print(f"max lactate: {s['lactate'].max():.1f} mM (proxy rates)")

fig, ax = plt.subplots(2, 2, figsize=(11, 9.5))
sl = np.abs(pos[:, 2] - pos[:, 2].mean()) < 8
col = np.where(ph == 1, "#2ca02c", np.where(ph == 0, "#d62728", "#7f7f7f"))
ax[0, 0].scatter(pos[sl, 0], pos[sl, 1], c=col[sl], s=9, lw=0)
ax[0, 0].set_title("Mid-plane slice: Fucci-style phase\nred = G1, green = S/G2/M, grey = necrotic")
sc = ax[0, 1].scatter(pos[sl, 0], pos[sl, 1], c=s["oxygen"][sl], s=9, lw=0, cmap="viridis")
plt.colorbar(sc, ax=ax[0, 1], label="pO2 (mmHg)")
ax[0, 1].set_title("Mid-plane slice: oxygen")
for a in ax[0]:
    a.set_aspect("equal")
    a.set_xlabel("x (um)")
    a.set_ylabel("y (um)")

a = ax[1, 0]
a.plot(mid, green, "g-o", ms=3, label="S/G2/M fraction")
a.plot(mid, np.array(mit), "b-", label="mitogen (rel.)")
a.axvline(70, color="k", ls="--", lw=1, label="HeLa-Fucci rim ~70 um")
a.set_xlabel("depth from surface (um)")
a.set_ylabel("fraction / relative level")
a2 = a.twinx()
a2.plot(mid, o2, "m-", label="pO2")
a2.axhspan(10, 40, color="m", alpha=0.08)
a2.set_ylabel("pO2 (mmHg); band = HIF+/pimo- range")
a.legend(loc="center right", fontsize=8)
a.set_title("Radial profiles at final size")

t = np.array([float(x["t_h"]) for x in log]) / 24
ax[1, 1].plot(t, [float(x["diameter_um"]) for x in log], "k-")
ax[1, 1].set_xlabel("days")
ax[1, 1].set_ylabel("diameter (um)")
b2 = ax[1, 1].twinx()
b2.plot(t, [float(x["frac_SG2M"]) for x in log], "g-", label="S/G2/M")
b2.plot(t, [float(x["frac_quiescent"]) for x in log], "r-", label="G1-arrested")
b2.set_ylabel("fraction of cells")
b2.legend(fontsize=8)
ax[1, 1].set_title("Growth from 1000 seeded cells")
fig.tight_layout()
fig.savefig("out/spheroid_summary.png", dpi=130)
print("saved out/spheroid_summary.png")
