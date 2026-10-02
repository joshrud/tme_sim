# TME world — v0: the base cell

A graph-based 3D+time tumor model. **Nodes are cells; edges join cells that touch.**
Chemical fields live on the nodes and move along edges. v0 has one cell identity (HeLa),
grown as a spheroid and checked against HeLa-Fucci spheroid data.

## Why HeLa
- The best-documented cell-cycle timings of any human line (Puck & Steffen 1963; live biosensors, Hahn et al. 2009).
- Measured volume and O2 uptake.
- A direct 3D benchmark: HeLa-Fucci spheroids (Onozato et al. 2017, *Cancer Sci*, doi:10.1111/cas.13178). At ~700 µm,
  S/G2/M cells sit in a **~70 µm outer rim**. The interior is G1-arrested and only **mildly hypoxic**
  (HIF-1α positive, pimonidazole negative).

## What one cell is (the whole state)
| field | meaning |
|---|---|
| `pos` (x,y,z) | center, µm |
| `vol` | volume; grows linearly from V_birth to 2·V_birth over the cycle |
| `age`, `T` | hours into the cycle / this cell's cycle length (lognormal) |
| `state` | alive / necrotic |
| `ctype` | identity (0 = HeLa; T cell, monocyte, fibroblast... later) |

Rules per 0.25 h step:
1. Rebuild the contact graph, then solve every field to quasi-steady state on the graph.
   Diffusion takes seconds to minutes, while the cell cycle takes hours.
2. **Nutrients:** O2 and glucose uptake follow Michaelis–Menten kinetics. Lactate is secreted at a fixed ratio to glucose uptake.
   A serum mitogen is consumed first-order.
3. **G1 checkpoint:** a G1 cell arrests if pO2 < 5 mmHg or mitogen < 0.5. S/G2/M cells finish their cycle.
4. **Division** when `age ≥ T`: two daughters, half the volume each, new lognormal `T`.
5. **Necrosis** below 0.8 mmHg O2 (hazard, mean 6 h).
6. Overlap relaxation pushes cells apart and gives weak adhesion.

Every birth, division, necrosis (and later immune entry and metastatic exit) goes to `World.events`.
Snapshots every 12 h make up the 4D record.

## Field solver (graph "flow fields")
- Edge conductance is `D·κ_ij/d_ij²` with `κ_ij = 12/(z_i+z_j)`. Summed over neighbors, this approximates `D∇²c`.
- Cells touching the medium are found by flood-filling empty space on a voxel grid, so it works for any shape.
  They exchange with the bath through one face.
- Each cell's uptake is spread over the space it occupies, computed from its contact distances at random close packing.
- Solved with Jacobi-preconditioned conjugate gradient, warm-started from the previous step.
  This takes ~0.4 s per step at 20k cells.
- **Validation** (`validate_diffusion.py`) against the analytic O2 profile in a uniformly consuming sphere:
  the interior error is ~9–12% of the total drop. The bias is systematic: the graph slightly over-depletes,
  by ~15% at 30k cells.

## Evidence for each parameter
See `out/param_table.md`. Each value is tagged **measured / derived / proxy / calibrated / assumption / numerical**.
Gaps flagged so far:
- **Cycle-time CV**: no HeLa-specific value verified. This is the first candidate to fit from FUCCI lineage data.
- **Glucose uptake and lactate:glucose ratio**: these are *HEK293* values (Noguchi et al. 2020). A search result had
  attributed them to HeLa, which was wrong. HeLa flux data are needed.
- **What drives interior quiescence**: O2 alone can't produce a 70 µm rim. With the measured HeLa O2 uptake,
  pO2 at 70 µm depth is still ~90 mmHg. So arrest is attributed to a consumed mitogen, supported by EGF-withdrawal
  quiescence (Laurent et al. 2013). Its uptake rate is **calibrated** to the rim. It could also be contact/mechanical
  inhibition (Montel 2011; Delarue 2014). This is the best place for a small learned rule once spatial data are in.
- Necrosis timing, hypoxic-arrest threshold (PhysiCell convention), glucose Km: all assumptions.

## Run
```bash
pip install -r requirements.txt
python validate_diffusion.py   # graph diffusion vs analytic (~1 min)
python run_spheroid.py 700     # grow to 700 um (writes out/)
python analyze.py              # compare to HeLa-Fucci spheroid
```

## Files
- `tme/params.py` — every parameter with its source
- `tme/world.py` — the world: cells, contact graph, fields, cycle, events
- `validate_diffusion.py` — graph-diffusion vs. analytic check
- `run_spheroid.py` — grow to 700 µm; writes `out/snapshots.npz`, `out/events.csv`, `out/log.csv`
- `analyze.py` — comparison with the HeLa-Fucci spheroid; writes `out/spheroid_summary.png`
