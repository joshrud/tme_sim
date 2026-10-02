"""Pack out/snapshots.npz into viewer/data/frames.bin.gz for the three.js viewer.

Layout (little-endian), header, one block per frame (each padded to 4 bytes), lineage:
  header: b"TME4", u32 version (=2), u32 n_frames,
          f32 mitosis stage durations (h): NEBD->metaphase, metaphase->anaphase,
              anaphase->cytokinesis (a division event = end of cytokinesis)
  frame:  f32 t_hours, u32 n,
          u32 id[n]          (sorted, so the viewer can match cells between frames)
          i16 pos[3n]        (x,y,z interleaved, 0.1 um, centered on the final spheroid)
          u8  radius[n]      (display radius, 0.1 um)
          i8  phase[n]       (0 = G1, 1 = S/G2/M, -1 = necrotic)
          u8  o2[n]          (pO2, mmHg, clipped to 255)
  lineage: u32 n_ids, i32 parent[n_ids] (-1 = seeded cell), f32 birth_h[n_ids]
           (indexed by cell id; a mother's division times = its daughters' birth times)
"""
import gzip
import os
import struct

import numpy as np

from tme import params as P
from tme.world import radius

z = np.load("out/snapshots.npz")
frames = sorted({int(k.rsplit("_", 1)[1]) for k in z.files})
last = frames[-1]
center = z[f"pos_{last:03d}"].mean(0)

M = P.MITOSIS
out = bytearray(b"TME4" + struct.pack("<II", 2, len(frames)))
out += struct.pack("<3f", M["nebd_to_metaphase"].value, M["metaphase_to_anaphase"].value,
                   M["anaphase_to_cytokinesis"].value)
for f in frames:
    g = lambda k: z[f"{k}_{f:03d}"]
    order = np.argsort(g("id"))
    pos = np.round((g("pos")[order] - center) * 10).astype(np.int16)
    # spheres drawn at the packing distance so neighbors touch instead of overlapping
    r = np.clip(np.round(P.NUMERICS["packing"].value * radius(g("vol")[order]) * 10), 0, 255)
    out += struct.pack("<fI", float(g("t")), len(order))
    out += g("id")[order].astype(np.uint32).tobytes()
    out += pos.tobytes()
    out += r.astype(np.uint8).tobytes()
    out += g("phase")[order].astype(np.int8).tobytes()
    out += np.clip(np.round(g("oxygen")[order]), 0, 255).astype(np.uint8).tobytes()
    out += b"\0" * (-len(out) % 4)

# lineage from the event log: who each cell came from, and when
ev = np.genfromtxt("out/events.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
born = np.isin(ev["event"], ["seed", "birth"])
n_ids = int(ev["cell_id"][born].max()) + 1
parent = np.full(n_ids, -1, np.int32)
birth = np.zeros(n_ids, np.float32)
parent[ev["cell_id"][born]] = ev["parent_id"][born]
birth[ev["cell_id"][born]] = ev["t_h"][born]
out += struct.pack("<I", n_ids) + parent.tobytes() + birth.tobytes()

os.makedirs("viewer/data", exist_ok=True)
with gzip.open("viewer/data/frames.bin.gz", "wb") as fh:
    fh.write(out)
print(f"{len(frames)} frames, {len(out) / 1e6:.1f} MB raw -> "
      f"{os.path.getsize('viewer/data/frames.bin.gz') / 1e6:.1f} MB gzipped")
