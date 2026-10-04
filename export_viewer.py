"""Pack a run (snapshots, event log, T-cell tracks) into viewer/data/frames.bin.gz.

  python export_viewer.py                 # run in out/
  python export_viewer.py out/effector    # -> viewer/data/effector.bin.gz, open with ?run=effector

Layout (little-endian), header, one block per frame (each padded to 4 bytes), lineage:
  header: b"TME4", u32 version (=3), u32 n_frames,
          f32 mitosis stage durations (h): NEBD->metaphase, metaphase->anaphase,
              anaphase->cytokinesis (a division event = end of cytokinesis)
  frame:  f32 t_hours, u32 n,
          u32 id[n]          (sorted, so the viewer can match cells between frames)
          i16 pos[3n]        (x,y,z interleaved, 0.1 um, centered on the final spheroid)
          u8  radius[n]      (display radius, 0.1 um)
          i8  phase[n]       (0 = G1, 1 = S/G2/M, -1 = necrotic, -2 = apoptotic)
          u8  o2[n]          (pO2, mmHg, clipped to 255)
  lineage: u32 n_ids, i32 parent[n_ids] (-1 = seeded cell), f32 birth_h[n_ids]
           (indexed by cell id; a mother's division times = its daughters' birth times)
  tcells:  u32 n_records; per record: f32 t, u32 count, i16 pos[3*count], u8 state[count], pad
           (record k holds T cells 0..count-1 in entry order; positions like tumor cells)
           u32 n_t; u32 id[n_t], f32 enter_h[n_t], f32 exit_h[n_t] (NaN = stayed), u8 cognate[n_t], pad
  hits:    u32 n_hits; f32 t[], u32 tcell_index[], u32 target_id[]
  deaths:  u32 n_dead; u32 target_id[], f32 t_killed[], f32 t_cleared[] (NaN = not cleared),
           u32 killer_index[]
"""
import sys
import gzip
import os
import struct

import numpy as np

from tme import params as P
from tme.world import radius

RUN = sys.argv[1] if len(sys.argv) > 1 else "out"
NAME = "frames" if os.path.normpath(RUN) == "out" else os.path.basename(os.path.normpath(RUN))
z = np.load(f"{RUN}/snapshots.npz")
frames = sorted({int(k.rsplit("_", 1)[1]) for k in z.files})
last = frames[-1]
center = z[f"pos_{last:03d}"].mean(0)

M = P.MITOSIS
out = bytearray(b"TME4" + struct.pack("<II", 3, len(frames)))
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
ev = np.genfromtxt(f"{RUN}/events.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
born = np.isin(ev["event"], ["seed", "birth"])
n_ids = int(ev["cell_id"][born].max()) + 1
parent = np.full(n_ids, -1, np.int32)
birth = np.zeros(n_ids, np.float32)
parent[ev["cell_id"][born]] = ev["parent_id"][born]
birth[ev["cell_id"][born]] = ev["t_h"][born]
out += struct.pack("<I", n_ids) + parent.tobytes() + birth.tobytes()

# T cells
pad = lambda: b"\0" * (-len(out) % 4)
try:
    tz = np.load(f"{RUN}/tcells.npz")
    t_rec, counts, tpos, tstate, tid = tz["t"], tz["count"], tz["pos"], tz["state"], tz["info_id"]
    tcog = tz["info_cognate"]
except FileNotFoundError:
    t_rec = counts = tid = tcog = np.zeros(0)
    tpos, tstate = np.zeros((0, 3)), np.zeros(0)
out += struct.pack("<I", len(t_rec))
o = 0
for t, c in zip(t_rec, counts):
    c = int(c)
    out += struct.pack("<fI", float(t), c)
    out += np.round((tpos[o:o + c] - center) * 10).astype(np.int16).tobytes()
    out += tstate[o:o + c].astype(np.uint8).tobytes()
    out += pad()
    o += c
tindex = {int(i): k for k, i in enumerate(tid)}
enter_t = np.full(len(tid), np.nan, np.float32)
for row in ev[ev["event"] == "enter"]:
    enter_t[tindex[int(row["cell_id"])]] = row["t_h"]
exit_t = np.full(len(tid), np.nan, np.float32)
for row in ev[ev["event"] == "exit"]:
    exit_t[tindex[int(row["cell_id"])]] = row["t_h"]
out += struct.pack("<I", len(tid)) + np.asarray(tid, np.uint32).tobytes() + enter_t.tobytes() + exit_t.tobytes()
out += np.asarray(tcog, np.uint8).tobytes()
out += pad()

hits = ev[ev["event"] == "hit"]
out += struct.pack("<I", len(hits)) + hits["t_h"].astype(np.float32).tobytes()
out += np.array([tindex[int(i)] for i in hits["cell_id"]], np.uint32).tobytes()
out += hits["parent_id"].astype(np.uint32).tobytes()

killed = ev[ev["event"] == "killed"]
cleared = {int(r["cell_id"]): r["t_h"] for r in ev[ev["event"] == "cleared"]}
out += struct.pack("<I", len(killed)) + killed["cell_id"].astype(np.uint32).tobytes()
out += killed["t_h"].astype(np.float32).tobytes()
out += np.array([cleared.get(int(i), np.nan) for i in killed["cell_id"]], np.float32).tobytes()
out += np.array([tindex[int(i)] for i in killed["parent_id"]], np.uint32).tobytes()

os.makedirs("viewer/data", exist_ok=True)
with gzip.open(f"viewer/data/{NAME}.bin.gz", "wb") as fh:
    fh.write(out)
print(f"{len(frames)} frames, {len(tid)} T cells, {len(hits)} hits, {len(killed)} kills, "
      f"{len(out) / 1e6:.1f} MB raw -> "
      f"{os.path.getsize(f'viewer/data/{NAME}.bin.gz') / 1e6:.1f} MB gzipped -> viewer/data/{NAME}.bin.gz")
