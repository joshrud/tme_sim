import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { makeCellMaterial, makeCellMesh, cellUniforms, SmoothTissueRenderer } from "./cellshader.js";

const $ = (id) => document.getElementById(id);

// ------------------------------------------------------------------ data
async function loadData(url) {
  const res = await fetch(url);
  const buf = await new Response(res.body.pipeThrough(new DecompressionStream("gzip"))).arrayBuffer();
  const dv = new DataView(buf);
  if (new TextDecoder().decode(new Uint8Array(buf, 0, 4)) !== "TME4") throw new Error("bad file");
  const version = dv.getUint32(4, true);
  if (version < 2) throw new Error("re-run export_viewer.py");
  const readBlock = (o) => {  // u32 n, i16 pos[3n], u8 label[n], pad
    const n = dv.getUint32(o, true); o += 4;
    const pos = new Int16Array(buf, o, 3 * n); o += 6 * n;
    const lab = new Uint8Array(buf, o, n); o += n;
    o += (4 - (o % 4)) % 4;
    return [{ n, pos, lab }, o];
  };
  const nFrames = dv.getUint32(8, true);
  const mitosis = { nebd: dv.getFloat32(12, true), meta: dv.getFloat32(16, true), ana: dv.getFloat32(20, true) };
  let o = 24;
  const frames = [];
  for (let k = 0; k < nFrames; k++) {
    const t = dv.getFloat32(o, true), n = dv.getUint32(o + 4, true);
    o += 8;
    const id = new Uint32Array(buf, o, n); o += 4 * n;
    const pos = new Int16Array(buf, o, 3 * n); o += 6 * n;
    const r = new Uint8Array(buf, o, n); o += n;
    const phase = new Int8Array(buf, o, n); o += n;
    const o2 = new Uint8Array(buf, o, n); o += n;
    o += (4 - (o % 4)) % 4;
    frames.push({ t, n, id, pos, r, phase, o2 });
  }
  const nIds = dv.getUint32(o, true); o += 4;
  const parent = new Int32Array(buf, o, nIds); o += 4 * nIds;
  const birth = new Float32Array(buf, o, nIds); o += 4 * nIds;
  const tc = { recs: [], n: 0, enter: new Float32Array(0), cognate: new Uint8Array(0),
               hits: { n: 0 }, deaths: { n: 0, target: new Uint32Array(0) } };
  if (version >= 3) {
    const nRec = dv.getUint32(o, true); o += 4;
    for (let k = 0; k < nRec; k++) {
      const t = dv.getFloat32(o, true), c = dv.getUint32(o + 4, true);
      o += 8;
      const pos = new Int16Array(buf, o, 3 * c); o += 6 * c;
      const state = new Uint8Array(buf, o, c); o += c;
      o += (4 - (o % 4)) % 4;
      tc.recs.push({ t, c, pos, state });
    }
    tc.n = dv.getUint32(o, true); o += 4;
    tc.ids = new Uint32Array(buf, o, tc.n); o += 4 * tc.n;
    tc.enter = new Float32Array(buf, o, tc.n); o += 4 * tc.n;
    tc.exit = new Float32Array(buf, o, tc.n); o += 4 * tc.n;
    tc.cognate = new Uint8Array(buf, o, tc.n); o += tc.n;
    o += (4 - (o % 4)) % 4;
    const nh = dv.getUint32(o, true); o += 4;
    tc.hits = { n: nh, t: new Float32Array(buf, o, nh), tcell: new Uint32Array(buf, o + 4 * nh, nh),
                target: new Uint32Array(buf, o + 8 * nh, nh) };
    o += 12 * nh;
    const nd = dv.getUint32(o, true); o += 4;
    tc.deaths = { n: nd, target: new Uint32Array(buf, o, nd), tk: new Float32Array(buf, o + 4 * nd, nd),
                  tc: new Float32Array(buf, o + 8 * nd, nd), killer: new Uint32Array(buf, o + 12 * nd, nd) };
    o += 16 * nd;   // the deaths block must be stepped over before anything that follows
  }
  let fib = { n: 0 }, imm = { n: 0 };
  if (version >= 4) {
    [fib, o] = readBlock(o);
    [imm, o] = readBlock(o);
  }
  return { frames, mitosis, parent, birth, tc, fib, imm };
}

// index into frame a for every cell of frame b (-1 = born since a); ids are sorted
function matchIds(a, b) {
  const m = new Int32Array(b.n).fill(-1);
  let i = 0;
  for (let j = 0; j < b.n; j++) {
    while (i < a.n && a.id[i] < b.id[j]) i++;
    if (i < a.n && a.id[i] === b.id[j]) m[j] = i;
  }
  return m;
}

function indexOf(sortedIds, n, id) {
  let lo = 0, hi = n - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1, v = sortedIds[mid];
    if (v === id) return mid;
    if (v < id) lo = mid + 1; else hi = mid - 1;
  }
  return -1;
}

// ------------------------------------------------------------------ lineage
// A mother keeps her id when she divides; each daughter gets a new id. So:
//   founder = the seeded cell a cell descends from, depth = daughter-births on that path,
//   a cell's division times = the birth times of the cells whose parent it is.
function buildLineage(parent, birth) {
  const n = parent.length;
  const founder = new Int32Array(n), depth = new Uint16Array(n), nKids = new Int32Array(n + 1);
  for (let id = 0; id < n; id++) {  // ids are assigned in time order, so parents come first
    const p = parent[id];
    founder[id] = p < 0 ? id : founder[p];
    depth[id] = p < 0 ? 0 : depth[p] + 1;
    if (p >= 0) nKids[p + 1]++;
  }
  for (let i = 0; i < n; i++) nKids[i + 1] += nKids[i];
  const divStart = nKids.slice(), divTimes = new Float32Array(nKids[n]), divChild = new Int32Array(nKids[n]);
  const fill = nKids.slice(0, n);
  for (let id = 0; id < n; id++) {
    if (parent[id] < 0) continue;
    const q = fill[parent[id]]++;
    divTimes[q] = birth[id];
    divChild[q] = id;
  }
  return { founder, depth, divStart, divTimes, divChild };
}

// ------------------------------------------------------------------ colors
const PHASE = { 0: [0.84, 0.15, 0.16], 1: [0.17, 0.63, 0.17], 2: [0.98, 0.8, 0.1], [-1]: [0.5, 0.5, 0.5],
                [-2]: [0.5, 0.22, 0.62] };  // -2 = apoptotic (killed by a CTL)
const DIM = [0.16, 0.17, 0.2];
const VIRIDIS = [[0.267, 0.005, 0.329], [0.229, 0.322, 0.546], [0.128, 0.567, 0.551],
                 [0.369, 0.789, 0.383], [0.993, 0.906, 0.144]];
function viridis(x) {
  x = Math.min(Math.max(x, 0), 1) * (VIRIDIS.length - 1);
  const i = Math.min(Math.floor(x), VIRIDIS.length - 2), f = x - i, a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
}
function cloneColor(f) {  // golden-ratio hues so neighboring founder ids look distinct
  const h = (f * 0.6180339887) % 1, s = 0.65, l = 0.55, k = (n) => (n + h * 12) % 12;
  const a = s * Math.min(l, 1 - l), g = (n) => l - a * Math.max(-1, Math.min(k(n) - 3, 9 - k(n), 1));
  return [g(0), g(8), g(4)];
}
const O2_MAX = 142;
const css = (c) => `rgb(${c.map((v) => Math.round(v * 255)).join(",")})`;

// ------------------------------------------------------------------ main
const RUN = new URLSearchParams(location.search).get("run") || "frames";  // ?run=effector
// Off-screen lymphoid organs: counts by developmental stage, updated with the timeline.
let organs = null;
try {
  organs = await (await fetch(`data/${RUN}.organs.json`)).json();
} catch (e) { organs = null; }
const { frames, mitosis, parent, birth, tc, fib, imm } = await loadData(`data/${RUN}.bin.gz`);
const lin = buildLineage(parent, birth);
$("loading").remove();
const M_TOTAL = mitosis.nebd + mitosis.meta + mitosis.ana;  // yellow window before cytokinesis
const T0 = frames[0].t, T1 = frames[frames.length - 1].t;
const matches = frames.slice(1).map((b, k) => matchIds(frames[k], b));
// tumor cells present at A but gone by B (cleared after being killed): keep drawing them until cleared
const removed = frames.slice(1).map((b, k) => {
  const seen = new Uint8Array(frames[k].n);
  for (const i of matches[k]) if (i >= 0) seen[i] = 1;
  return [...seen.keys()].filter((i) => !seen[i]);
});
const maxN = Math.max(...frames.slice(1).map((f, k) => f.n + removed[k].length), frames[0].n);
const tKilled = new Float32Array(parent.length + 1).fill(Infinity);
const tCleared = new Float32Array(parent.length + 1).fill(Infinity);
for (let q = 0; q < tc.deaths.n; q++) {
  tKilled[tc.deaths.target[q]] = tc.deaths.tk[q];
  if (!Number.isNaN(tc.deaths.tc[q])) tCleared[tc.deaths.target[q]] = tc.deaths.tc[q];
}
const APOPTOTIC_SHRINK_H = 1.0;  // visual: apoptotic cells shrink to 70% radius over ~1 h
let extent = 0;  // half-width of the final spheroid, um
{ const f = frames[frames.length - 1]; for (let i = 0; i < 3 * f.n; i++) extent = Math.max(extent, Math.abs(f.pos[i])); }
extent = Math.ceil(extent / 10 + 20);

function nextDivisionIndex(id, t) {  // index (into divTimes/divChild) of the first division after t
  for (let q = lin.divStart[id]; q < lin.divStart[id + 1]; q++) if (lin.divTimes[q] > t) return q;
  return -1;
}
function nextDivision(id, t) {
  const q = nextDivisionIndex(id, t);
  return q < 0 ? Infinity : lin.divTimes[q];
}

// current interpolated state, shared by the 3D view, the cross-section and picking
const cur = { n: 0, id: new Uint32Array(maxN), x: new Float32Array(maxN), y: new Float32Array(maxN),
              z: new Float32Array(maxN), r: new Float32Array(maxN), phase: new Int8Array(maxN),
              o2: new Float32Array(maxN),
              div: new Float32Array(4 * maxN) };  // division axis xyz + split progress (0 = none)

function computeState(t) {
  let k = 0;
  while (k < frames.length - 2 && frames[k + 1].t <= t) k++;
  const A = frames[k], B = frames[k + 1], m = matches[k];
  const f = Math.min(Math.max((t - A.t) / (B.t - A.t), 0), 1);
  const lerp = (F0, i0, F1, i1, c, w) => (1 - w) * F0.pos[3 * i0 + c] / 10 + w * F1.pos[3 * i1 + c] / 10;
  const motherAt = (jm, c, tt) => {  // mother's position at time tt within this window
    const am = m[jm], w = Math.min(Math.max((tt - A.t) / (B.t - A.t), 0), 1);
    return am < 0 ? B.pos[3 * jm + c] / 10 : lerp(A, am, B, jm, c, w);
  };
  const ANA = mitosis.ana;

  for (let j = 0; j < B.n; j++) {
    const id = B.id[j], i = m[j], rB = B.r[j] / 10;
    let x, y, z, r, ph, o2, splitAxis = null;
    if (i >= 0) {  // existed at A: interpolate; if it divides in this window, shrink during the split
      x = lerp(A, i, B, j, 0, f); y = lerp(A, i, B, j, 1, f); z = lerp(A, i, B, j, 2, f);
      const rA = A.r[i] / 10, td = nextDivision(id, A.t);
      o2 = (1 - f) * A.o2[i] + f * B.o2[j];
      if (td <= B.t) {
        r = t < td ? rA : rB;
        ph = t >= td ? (f >= 1 ? B.phase[j] : 0) : A.phase[i];
        const jd = t >= td - ANA && t < td ? indexOf(B.id, B.n, lin.divChild[nextDivisionIndex(id, A.t)]) : -1;
        if (jd >= 0) {  // anaphase -> cytokinesis: draw mother + daughter as one pinching dumbbell
          let dx = B.pos[3 * jd] - B.pos[3 * j], dy = B.pos[3 * jd + 1] - B.pos[3 * j + 1], dz = B.pos[3 * jd + 2] - B.pos[3 * j + 2];
          const len = Math.hypot(dx, dy, dz) || 1, s = (t - (td - ANA)) / ANA, half = s * (rB + B.r[jd] / 10) / 2;
          dx /= len; dy /= len; dz /= len;
          x += dx * half; y += dy * half; z += dz * half;
          splitAxis = [dx, dy, dz, Math.max(s, 1e-3)];
        }
      } else {
        r = (1 - f) * rA + f * rB;
        ph = f < 0.5 ? A.phase[i] : B.phase[j];
      }
    } else {  // daughter born in this window: emerge from the mother, then drift to its place
      const tb = birth[id], jm = indexOf(B.id, B.n, parent[id]);
      o2 = B.o2[j];
      if (jm < 0) {
        r = t >= tb ? rB : 0; x = B.pos[3 * j] / 10; y = B.pos[3 * j + 1] / 10; z = B.pos[3 * j + 2] / 10; ph = 0;
      } else if (t < tb) {
        r = 0; x = y = z = 0; ph = 0;  // not separate yet: drawn as part of the mother's dumbbell
      } else {
        let dx = B.pos[3 * j] - B.pos[3 * jm], dy = B.pos[3 * j + 1] - B.pos[3 * jm + 1], dz = B.pos[3 * j + 2] - B.pos[3 * jm + 2];
        const len = Math.hypot(dx, dy, dz) || 1, sep = rB + B.r[jm] / 10;
        dx /= len; dy /= len; dz /= len;
        r = rB;
        const u = (t - tb) / Math.max(B.t - tb, 1e-6);  // drift from the split point to its place
        x = (1 - u) * (motherAt(jm, 0, tb) + sep * dx) + u * B.pos[3 * j] / 10;
        y = (1 - u) * (motherAt(jm, 1, tb) + sep * dy) + u * B.pos[3 * j + 1] / 10;
        z = (1 - u) * (motherAt(jm, 2, tb) + sep * dz) + u * B.pos[3 * j + 2] / 10;
        ph = f >= 1 ? B.phase[j] : 0;
      }
    }
    if (t >= tKilled[id]) {
      ph = -2;
      r *= 1 - 0.3 * Math.min((t - tKilled[id]) / APOPTOTIC_SHRINK_H, 1);
    } else if (ph === -2) {
      ph = i >= 0 ? Math.max(A.phase[i], 0) : 0;  // snapshot says dead, but not yet at time t
    }
    if (ph >= 0 && nextDivision(id, t) - t <= M_TOTAL) ph = 2;  // in mitosis (NEBD -> cytokinesis)
    cur.id[j] = id; cur.x[j] = x; cur.y[j] = y; cur.z[j] = z; cur.r[j] = r; cur.phase[j] = ph; cur.o2[j] = o2;
    const D = splitAxis || [0, 0, 0, 0];
    cur.div[4 * j] = D[0]; cur.div[4 * j + 1] = D[1]; cur.div[4 * j + 2] = D[2]; cur.div[4 * j + 3] = D[3];
  }
  let n = B.n;
  for (const i of removed[k]) {  // killed and cleared before B: shrink, then vanish at clearance
    const id = A.id[i], gone = t >= tCleared[id];
    cur.id[n] = id; cur.x[n] = A.pos[3 * i] / 10; cur.y[n] = A.pos[3 * i + 1] / 10; cur.z[n] = A.pos[3 * i + 2] / 10;
    cur.r[n] = gone ? 0 : (A.r[i] / 10) * (t >= tKilled[id] ? 1 - 0.3 * Math.min((t - tKilled[id]) / APOPTOTIC_SHRINK_H, 1) : 1);
    cur.phase[n] = t >= tKilled[id] ? -2 : Math.max(A.phase[i], 0);
    cur.o2[n] = A.o2[i];
    cur.div[4 * n] = cur.div[4 * n + 1] = cur.div[4 * n + 2] = cur.div[4 * n + 3] = 0;
    n++;
  }
  cur.n = n;
  cur.nSorted = B.n;  // cur.id[0..nSorted) is sorted; removed cells follow
}

function findCur(id) {  // index of tumor cell `id` in the current state, or -1
  const j = indexOf(cur.id, cur.nSorted, id);
  if (j >= 0) return j;
  for (let i = cur.nSorted; i < cur.n; i++) if (cur.id[i] === id) return i;
  return -1;
}

function colorOf(i) {
  const id = cur.id[i];
  if (ui.clone >= 0 && lin.founder[id] !== ui.clone) return DIM;
  if (ui.color === "o2") return viridis(cur.o2[i] / O2_MAX);
  if (ui.color === "lineage") return cloneColor(lin.founder[id]);
  return PHASE[cur.phase[i]];
}

function drawLegend() {
  const l = $("legend");
  if (ui.color === "phase") {
    l.innerHTML = [[0, "G1"], [1, "S/G2"], [2, "M (mitosis)"], [-1, "necrotic"], [-2, "apoptotic (killed)"]]
      .map(([p, s]) => `<span><i style="background:${css(PHASE[p])}"></i>${s}</span>`).join("");
  } else if (ui.color === "o2") {
    l.innerHTML = `<span style="flex:1">0 <i style="width:120px;height:8px;border-radius:2px;background:linear-gradient(90deg,
      ${[0, 0.25, 0.5, 0.75, 1].map((x) => css(viridis(x))).join(",")})"></i> ${O2_MAX} mmHg</span>`;
  } else {
    l.innerHTML = `<span>One color per founding (seeded) cell. Click a cell to trace its clone.</span>`;
  }
  if (imm.n || fib.n) {
    const seen = new Set([...imm.lab]);
    const parts = IMMUNE_NAMES.map((nm, k) => seen.has(k)
      ? `<span><i style="background:${IMMUNE_COLORS[k]}"></i>${nm}</span>` : "").join("");
    const fibParts = fib.n ? FIBRO_NAMES.map((nm, k) =>
      `<span><i style="background:${FIBRO_COLORS[k]}"></i>${nm}</span>`).join("") : "";
    l.innerHTML += `<span style="flex-basis:100%;margin-top:6px"></span>${fibParts}${parts}`;
  }
}

// ------------------------------------------------------------------ 3D view
const canvas = $("view");
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
const BG = new THREE.Color(0x0e1116);
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(40, 1, 1, 10000);
camera.position.set(extent * 2.2, extent * 1.4, extent * 2.2);
camera.layers.enableAll();  // cells are on layer 1, overlays on layer 0
const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
const ambient = new THREE.AmbientLight(0xffffff, 0.55), sun = new THREE.DirectionalLight(0xffffff, 1.6);
sun.position.set(1, 2, 1.5);
scene.add(ambient, sun);

// most cells use a light mesh; the few currently dividing get a finer one so the furrow looks smooth
const material = makeCellMaterial();
const cells = makeCellMesh(1, maxN, material);
const dividing = makeCellMesh(3, 5000, material);
scene.add(cells, dividing);
const smooth = new SmoothTissueRenderer(renderer);

// ------------------------------------------------------------------ T cells
// White cells with a big "T" turned toward the viewer; the letter's color shows the state.
const T_STATE = { 0: { name: "naive", letter: "#1f2328" }, 1: { name: "anergic", letter: "#8c959f" },
                  3: { name: "effector CTL", letter: "#cf222e" } };
const T_RADIUS = { 0: 4.0, 1: 4.0, 3: 4.0 * Math.cbrt(3), 4: 4.0 };  // 4 = exited  // 8 um naive; activated blasts ~3x volume
const ENTER_MS = 1400, HIT_MS = 1000;  // wall-clock length of the entry and perforin/granzyme animations

function letterTexture(color) {
  const c = document.createElement("canvas");
  c.width = 512; c.height = 256;
  const g = c.getContext("2d");
  g.fillStyle = "#ffffff";
  g.fillRect(0, 0, 512, 256);
  g.fillStyle = color;
  g.font = "bold 170px Helvetica, Arial, sans-serif";
  g.textAlign = "center";
  g.textBaseline = "middle";
  g.fillText("T", 128, 138);  // u = 0.25 is the sphere's +z side, which is turned toward the camera
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.anisotropy = 4;
  return tex;
}
const tGeo = new THREE.SphereGeometry(1, 24, 16), tMeshes = {};
for (const [st, def] of Object.entries(T_STATE)) {
  const m = new THREE.InstancedMesh(tGeo, new THREE.MeshLambertMaterial({ map: letterTexture(def.letter) }), Math.max(tc.n, 1));
  m.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
  m.count = 0;
  m.frustumCulled = false;
  scene.add(m);  // layer 0: drawn on top of the smooth-tissue surface too
  tMeshes[st] = m;
}
// simulated T-cell state at the current time, plus where each is drawn (differs while entering)
const tcur = { x: new Float32Array(tc.n), y: new Float32Array(tc.n), z: new Float32Array(tc.n),
               dx: new Float32Array(tc.n), dy: new Float32Array(tc.n), dz: new Float32Array(tc.n),
               state: new Uint8Array(tc.n), vis: new Uint8Array(tc.n) };

function computeT(t) {
  const R = tc.recs;
  tcur.vis.fill(0);
  if (!R.length) return;
  let lo = 0, hi = R.length - 1;
  while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (R[mid].t <= t) lo = mid; else hi = mid - 1; }
  const A = R[lo], B = R[Math.min(lo + 1, R.length - 1)];
  const f = B === A ? 0 : Math.min(Math.max((t - A.t) / (B.t - A.t), 0), 1);
  for (let i = 0; i < tc.n; i++) {
    if (i >= B.c || t < tc.enter[i]) continue;
    tcur.vis[i] = !(t >= tc.exit[i]);  // exited cells are hidden (an exit animation may still draw them)
    const P = i < A.c ? A : B, w = i < A.c ? f : 1;
    tcur.x[i] = ((1 - w) * P.pos[3 * i] + w * B.pos[3 * i]) / 10;
    tcur.y[i] = ((1 - w) * P.pos[3 * i + 1] + w * B.pos[3 * i + 1]) / 10;
    tcur.z[i] = ((1 - w) * P.pos[3 * i + 2] + w * B.pos[3 * i + 2]) / 10;
    tcur.state[i] = f < 0.5 && i < A.c ? A.state[i] : B.state[i];
  }
}

// event-triggered animations (wall clock, so they are visible at any playback speed)
const anims = { enter: [], exit: [], hit: [] };
function triggerEvents(t0, t1) {  // events with time in (t0, t1]
  if (t1 <= t0) return;
  const now = performance.now(), H = tc.hits;
  for (let i = 0; i < tc.n; i++) {
    if (tc.enter[i] > t0 && tc.enter[i] <= t1) anims.enter.push({ i, start: now });
    if (tc.exit[i] > t0 && tc.exit[i] <= t1) anims.exit.push({ i, start: now });
  }
  let lo = 0, hi = H.n;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (H.t[mid] <= t0) lo = mid + 1; else hi = mid; }
  for (let q = lo; q < H.n && H.t[q] <= t1; q++) anims.hit.push({ q, start: now });
  if (anims.hit.length > 400) anims.hit.splice(0, anims.hit.length - 400);
}

const dummy = new THREE.Object3D();
function updateT(now) {
  anims.enter = anims.enter.filter((a) => now - a.start < ENTER_MS);
  anims.exit = anims.exit.filter((a) => now - a.start < ENTER_MS);
  const entering = new Map(anims.enter.map((a) => [a.i, (now - a.start) / ENTER_MS]));
  const exiting = new Map(anims.exit.map((a) => [a.i, (now - a.start) / ENTER_MS]));
  const k = { 0: 0, 1: 0, 3: 0 }, side = extent * 1.9;
  for (let i = 0; i < tc.n; i++) {
    const pe = exiting.get(i);
    if (!tcur.vis[i] && pe === undefined) continue;
    let x = tcur.x[i];
    const p = entering.get(i);
    if (p !== undefined) x = side + (x - side) * (1 - (1 - p) ** 3);  // fly in from the +x side
    if (pe !== undefined) x = x + (side - x) * pe ** 2;               // and back out when leaving
    tcur.dx[i] = x; tcur.dy[i] = tcur.y[i]; tcur.dz[i] = tcur.z[i];
    const st = tcur.state[i] === 4 ? 0 : tcur.state[i], mesh = tMeshes[st];
    dummy.position.set(x, tcur.y[i], tcur.z[i]);
    dummy.scale.setScalar(T_RADIUS[st]);
    dummy.lookAt(camera.position);  // keep the "T" facing the viewer
    dummy.updateMatrix();
    mesh.setMatrixAt(k[st]++, dummy.matrix);
  }
  for (const st in tMeshes) { tMeshes[st].count = k[st]; tMeshes[st].instanceMatrix.needsUpdate = true; }
}

// perforin + granzyme: a burst of granules from the CTL to its target on every hit
const GRAN = 8, MAX_BURSTS = 400;
const granGeo = new THREE.BufferGeometry();
granGeo.setAttribute("position", new THREE.BufferAttribute(new Float32Array(3 * GRAN * MAX_BURSTS), 3));
const dot = document.createElement("canvas");  // round granules
dot.width = dot.height = 32;
{ const g = dot.getContext("2d"); g.fillStyle = "#fff"; g.beginPath(); g.arc(16, 16, 15, 0, 2 * Math.PI); g.fill(); }
const granules = new THREE.Points(granGeo, new THREE.PointsMaterial({
  color: 0xff3fb4, size: 2.4, sizeAttenuation: true, transparent: true, depthTest: false,
  map: new THREE.CanvasTexture(dot), alphaTest: 0.5 }));
granules.renderOrder = 998;
granules.frustumCulled = false;
scene.add(granules);

function updateGranules(now) {
  anims.hit = anims.hit.filter((a) => now - a.start < HIT_MS);
  const bursts = anims.hit.map((a) => [a.q, (now - a.start) / HIT_MS]);
  if (!ui.playing) {  // paused: show hits within ~3 min of the current time, mid-flight
    const H = tc.hits;
    for (let q = 0; q < H.n && bursts.length < MAX_BURSTS; q++) if (Math.abs(H.t[q] - ui.t) < 0.05) bursts.push([q, 0.5]);
  }
  const P = granGeo.attributes.position.array;
  let n = 0;
  for (const [q, p] of bursts.slice(-MAX_BURSTS)) {
    const ti = tc.hits.tcell[q], j = findCur(tc.hits.target[q]);
    if (!tcur.vis[ti] || j < 0) continue;
    const fx = tcur.dx[ti], fy = tcur.dy[ti], fz = tcur.dz[ti], tx = cur.x[j], ty = cur.y[j], tz = cur.z[j];
    for (let g = 0; g < GRAN; g++) {
      const s = Math.min(Math.max((p - g * 0.04) / 0.55, 0), 1);       // travel to the target
      const spread = 1.2 + 4 * Math.max(p - 0.55, 0);                    // then disperse into it
      const jx = Math.sin(g * 12.99 + q) * spread, jy = Math.sin(g * 78.23 + q) * spread, jz = Math.sin(g * 39.42 + q) * spread;
      P[3 * n] = fx + (tx - fx) * s * 0.85 + jx;
      P[3 * n + 1] = fy + (ty - fy) * s * 0.85 + jy;
      P[3 * n + 2] = fz + (tz - fz) * s * 0.85 + jz;
      n++;
    }
  }
  granGeo.setDrawRange(0, n);
  granGeo.attributes.position.needsUpdate = true;
}

// Plane showing where the cross-section is taken. Deliberately much larger than the tumour
// so it reads as an infinite cutting plane running off screen, and in a warm colour that no
// cell type uses (cells are red/green/yellow/purple/white/olive) so it stays legible.
const PLANE_SPAN = 12;                       // multiples of the tumour radius
const PLANE_COLOR = 0xff8c42;
const planeSize = PLANE_SPAN * 2 * extent;
const planeMesh = new THREE.Mesh(
  new THREE.PlaneGeometry(planeSize, planeSize),
  new THREE.MeshBasicMaterial({ color: PLANE_COLOR, transparent: true, opacity: 0.16,
                                side: THREE.DoubleSide, depthWrite: false }));
planeMesh.renderOrder = 997;
// a bright grid on the plane makes its orientation and position readable at a glance
const planeGrid = new THREE.GridHelper(planeSize, PLANE_SPAN * 8, PLANE_COLOR, PLANE_COLOR);
planeGrid.material.transparent = true;
planeGrid.material.opacity = 0.35;
planeGrid.material.depthWrite = false;
planeGrid.rotation.x = Math.PI / 2;          // GridHelper lies in XZ; the plane is XY
planeMesh.add(planeGrid);
scene.add(planeMesh);

// ring marking the selected cell
const marker = new THREE.Mesh(new THREE.IcosahedronGeometry(1, 2),
  new THREE.MeshBasicMaterial({ color: 0xffffff, wireframe: true, depthTest: false, transparent: true, opacity: 0.9 }));
marker.renderOrder = 999;  // draw after the (transparent) cells so it is never covered
marker.visible = false;
scene.add(marker);

function updateInstances() {
  const meshes = [cells, dividing], k = [0, 0];
  for (let i = 0; i < cur.n; i++) {
    // when a clone is highlighted, hide everything else in 3D (the cross-section keeps it dimmed)
    const r = ui.clone >= 0 && lin.founder[cur.id[i]] !== ui.clone ? 0 : cur.r[i];
    if (r <= 0) continue;
    const which = cur.div[4 * i + 3] > 0 && k[1] < dividing.instanceMatrix.count ? 1 : 0;
    const mesh = meshes[which], n = k[which]++;
    const M = mesh.instanceMatrix.array, o = 16 * n;
    M[o] = r; M[o + 1] = 0; M[o + 2] = 0; M[o + 3] = 0;
    M[o + 4] = 0; M[o + 5] = r; M[o + 6] = 0; M[o + 7] = 0;
    M[o + 8] = 0; M[o + 9] = 0; M[o + 10] = r; M[o + 11] = 0;
    M[o + 12] = cur.x[i]; M[o + 13] = cur.y[i]; M[o + 14] = cur.z[i]; M[o + 15] = 1;
    const c = colorOf(i), C = mesh.instanceColor.array;
    C[3 * n] = c[0]; C[3 * n + 1] = c[1]; C[3 * n + 2] = c[2];
    const dv = mesh.geometry.attributes.aDiv.array;
    for (let q = 0; q < 4; q++) dv[4 * n + q] = cur.div[4 * i + q];
    mesh.geometry.attributes.aSeed.array[n] = (cur.id[i] * 0.6180339887 % 1) * 6.283;
  }
  meshes.forEach((mesh, w) => {
    mesh.count = k[w];
    mesh.instanceMatrix.needsUpdate = true;
    mesh.instanceColor.needsUpdate = true;
    mesh.geometry.attributes.aDiv.needsUpdate = true;
    mesh.geometry.attributes.aSeed.needsUpdate = true;
  });
}

function updateMarker() {  // follows the selected tumor cell or T cell
  const ti = ui.selectedT;
  if (ti >= 0) {
    marker.visible = !!tcur.vis[ti];
    marker.position.set(tcur.dx[ti], tcur.dy[ti], tcur.dz[ti]);
    marker.scale.setScalar(T_RADIUS[tcur.state[ti]] * 1.5);
    return;
  }
  const s = ui.selected >= 0 ? findCur(ui.selected) : -1;
  marker.visible = s >= 0 && cur.r[s] > 0;
  if (marker.visible) {
    marker.position.set(cur.x[s], cur.y[s], cur.z[s]);
    marker.scale.setScalar(cur.r[s] * 1.4);
  }
}

// ------------------------------------------------------------------ stroma + immune cells
// Colours match tme/immune.py; kept distinct from tumour phase colours and white T cells.
const IMMUNE_COLORS = ["#8a8f3a", "#e08a2e", "#9b6bd6", "#2bb8c4", "#dfe3a8", "#f2e9a0"];
const IMMUNE_NAMES = ["macrophage", "dendritic", "NK", "B cell", "Treg", "neutrophil"];
const IMMUNE_RADIUS = [8, 7, 5, 5, 4, 5];
const FIBRO_COLORS = ["#c58fd0", "#7a5c86"];   // myCAF (near tumour), iCAF (distal)
const FIBRO_NAMES = ["myCAF", "iCAF"];

function staticPopulation(block, colors, radii, group) {
  if (!block.n) return null;
  const geo = new THREE.IcosahedronGeometry(1, 1);
  const mat = new THREE.MeshLambertMaterial({ transparent: true, opacity: 1 });
  const mesh = new THREE.InstancedMesh(geo, mat, block.n);
  mesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(3 * block.n), 3);
  mesh.frustumCulled = false;
  const d = new THREE.Object3D(), col = new THREE.Color();
  for (let i = 0; i < block.n; i++) {
    const k = block.lab[i];
    d.position.set(block.pos[3 * i] / 10, block.pos[3 * i + 1] / 10, block.pos[3 * i + 2] / 10);
    d.scale.setScalar(typeof radii === "number" ? radii : radii[k] || 5);
    d.updateMatrix();
    mesh.setMatrixAt(i, d.matrix);
    col.set(colors[k] || "#888888");
    mesh.instanceColor.setXYZ(i, col.r, col.g, col.b);
  }
  mesh.instanceMatrix.needsUpdate = true;
  mesh.instanceColor.needsUpdate = true;
  scene.add(mesh);
  return mesh;
}
const fibMesh = staticPopulation(fib, FIBRO_COLORS, 7.5);
const immMesh = staticPopulation(imm, IMMUNE_COLORS, IMMUNE_RADIUS);

// ------------------------------------------------------------------ cross-section
const xs = $("xs"), xctx = xs.getContext("2d");
const AXES = { 0: ["y", "z", "x"], 1: ["x", "z", "y"], 2: ["x", "y", "z"] };  // [u, v, normal]

function drawSection() {
  const W = xs.width, s = W / (2 * extent), [u, v, w] = AXES[ui.axis].map((k) => cur[k]);
  xctx.fillStyle = "#000";
  xctx.fillRect(0, 0, W, W);
  for (let i = 0; i < cur.n; i++) {
    const d = w[i] - ui.slice, r = cur.r[i];
    if (d > r || d < -r) continue;
    xctx.fillStyle = css(colorOf(i));
    xctx.beginPath();
    xctx.arc(W / 2 + u[i] * s, W / 2 - v[i] * s, Math.max(Math.sqrt(r * r - d * d) * s, 0.5), 0, 2 * Math.PI);
    xctx.fill();
  }
  const T = { 0: tcur.x, 1: tcur.y, 2: tcur.z }, tu = T[{ 0: 1, 1: 0, 2: 0 }[ui.axis]],
        tv = T[{ 0: 2, 1: 2, 2: 1 }[ui.axis]], tw = T[ui.axis];
  xctx.lineWidth = 2;
  for (let i = 0; i < tc.n; i++) {  // T cells: white discs, outlined in their letter color
    if (!tcur.vis[i]) continue;
    const d = tw[i] - ui.slice, r = T_RADIUS[tcur.state[i]];
    if (d > r || d < -r) continue;
    xctx.beginPath();
    xctx.arc(W / 2 + tu[i] * s, W / 2 - tv[i] * s, Math.max(Math.sqrt(r * r - d * d) * s, 1.5), 0, 2 * Math.PI);
    xctx.fillStyle = "#ffffff";
    xctx.fill();
    xctx.strokeStyle = T_STATE[tcur.state[i]].letter;
    xctx.stroke();
  }
  xctx.fillStyle = "#e6edf3";
  xctx.fillRect(16, W - 24, 100 * s, 4);  // 100 um scale bar
  xctx.font = `${Math.round(W / 26)}px sans-serif`;
  xctx.fillText("100 µm", 16, W - 32);
  const name = { 0: "Sagittal", 1: "Coronal", 2: "Transverse" }[ui.axis];
  xctx.fillText(`${name}  ${["x", "y", "z"][ui.axis]} = ${ui.slice.toFixed(0)} µm`, 16, 16 + W / 26);
}

function placePlane() {
  planeMesh.rotation.set(0, 0, 0);
  planeMesh.position.set(0, 0, 0);
  if (ui.axis === 0) { planeMesh.rotation.y = Math.PI / 2; planeMesh.position.x = ui.slice; }
  if (ui.axis === 1) { planeMesh.rotation.x = Math.PI / 2; planeMesh.position.y = ui.slice; }
  if (ui.axis === 2) { planeMesh.position.z = ui.slice; }
}

// ------------------------------------------------------------------ picking + lineage panel
const raycaster = new THREE.Raycaster();
function pick(clientX, clientY) {  // -> { kind: "tumor", id } | { kind: "t", i } | null
  const rect = canvas.getBoundingClientRect();
  raycaster.setFromCamera(new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1,
                                            -((clientY - rect.top) / rect.height) * 2 + 1), camera);
  const { origin: o, direction: d } = raycaster.ray;
  let best = null, bestT = Infinity;
  const test = (x, y, z, r, hit) => {
    const ox = x - o.x, oy = y - o.y, oz = z - o.z;
    const tca = ox * d.x + oy * d.y + oz * d.z, d2 = ox * ox + oy * oy + oz * oz - tca * tca;
    if (tca <= 0 || d2 > r * r) return;
    const tHit = tca - Math.sqrt(r * r - d2);
    if (tHit < bestT) { bestT = tHit; best = hit; }
  };
  for (let i = 0; i < cur.n; i++) {
    const r = cur.r[i];
    if (r <= 0 || (ui.clone >= 0 && lin.founder[cur.id[i]] !== ui.clone)) continue;
    test(cur.x[i], cur.y[i], cur.z[i], r, { kind: "tumor", id: cur.id[i] });
  }
  for (let i = 0; i < tc.n; i++) {
    if (tcur.vis[i]) test(tcur.dx[i], tcur.dy[i], tcur.dz[i], T_RADIUS[tcur.state[i]], { kind: "t", i });
  }
  return best;
}

function lineagePath(id) {  // founder -> ... -> id, following parents
  const path = [];
  for (let c = id; c >= 0; c = parent[c]) path.push(c);
  return path.reverse();
}

function showSelection() {
  const box = $("lineage");
  if (ui.selectedT >= 0) {
    const i = ui.selectedT;
    box.innerHTML = `
      <div><b>T cell #${tc.ids[i]}</b> · <span id="tstate"></span></div>
      <div class="muted">Entered at ${tc.enter[i].toFixed(1)} h · recognizes HeLa: ${tc.cognate[i] ? "yes" : "no"}</div>
      <div class="muted">Hits delivered: <span id="thits"></span> · kills: <span id="tkills"></span> (so far)</div>
      <div class="row"><button id="unsel">Clear</button></div>`;
    $("unsel").onclick = () => { ui.selectedT = -1; ui.dirty = true; showSelection(); };
    updateCloneCount();
    return;
  }
  if (ui.selected < 0) { box.innerHTML = `<span class="muted">Click a cell to see where it came from.</span>`; return; }
  const id = ui.selected, f = lin.founder[id], path = lineagePath(id);
  const fmt = (c) => `#${c}` + (parent[c] < 0 ? " (seeded)" : ` (born ${birth[c].toFixed(1)} h)`);
  box.innerHTML = `
    <div><b>Cell #${id}</b> · ${parent[id] < 0 ? "seeded at 0 h" : `born ${birth[id].toFixed(1)} h from #${parent[id]}`}</div>
    <div class="muted">Founder #${f} · lineage depth ${lin.depth[id]} · clone now <span id="clonecount"></span> cells</div>
    <div class="muted path">${path.map(fmt).join(" → ")}</div>
    <div class="row"><button id="hl">${ui.clone === f ? "Show all cells" : "Highlight clone"}</button>
    <button id="unsel">Clear</button></div>`;
  $("hl").onclick = () => { ui.clone = ui.clone === f ? -1 : f; ui.dirty = true; showSelection(); };
  $("unsel").onclick = () => { ui.selected = -1; ui.clone = -1; ui.dirty = true; showSelection(); };
  updateCloneCount();
}

function updateCloneCount() {  // live numbers in the selection panel
  if (ui.selectedT >= 0) {
    const i = ui.selectedT, H = tc.hits, D = tc.deaths;
    let h = 0, kills = 0;
    for (let q = 0; q < H.n && H.t[q] <= ui.t; q++) h += H.tcell[q] === i;
    for (let q = 0; q < D.n; q++) kills += D.killer[q] === i && D.tk[q] <= ui.t;
    if ($("thits")) { $("thits").textContent = h; $("tkills").textContent = kills; }
    if ($("tstate")) $("tstate").textContent = tcur.vis[i] ? T_STATE[tcur.state[i]].name
      : ui.t >= tc.exit[i] ? `left the tissue at ${tc.exit[i].toFixed(1)} h` : "not yet entered";
    return;
  }
  const el = $("clonecount");
  if (!el || ui.selected < 0) return;
  const f = lin.founder[ui.selected];
  let k = 0;
  for (let i = 0; i < cur.n; i++) k += cur.r[i] > 0 && lin.founder[cur.id[i]] === f;
  el.textContent = k.toLocaleString();
}

let down = null;
canvas.addEventListener("pointerdown", (e) => { down = [e.clientX, e.clientY]; });
canvas.addEventListener("pointerup", (e) => {
  if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 4) return;
  const hit = pick(e.clientX, e.clientY);
  if (!hit) return;
  if (hit.kind === "t") { ui.selectedT = hit.i; ui.selected = -1; } else { ui.selected = hit.id; ui.selectedT = -1; }
  ui.dirty = true;
  showSelection();
});

// ------------------------------------------------------------------ UI state
const ui = { t: T0, playing: false, speed: 12, color: "phase", axis: 2, slice: 0,
             selected: -1, selectedT: -1, clone: -1, style: "spheres", smoothing: 10, opacity: 1, dirty: true };
const time = $("time");
time.min = T0; time.max = T1;
$("slice").min = -extent; $("slice").max = extent;

function setTime(t) {
  ui.t = Math.min(Math.max(t, T0), T1);
  time.value = ui.t;
  ui.dirty = true;
}
time.oninput = () => setTime(+time.value);
$("play").onclick = () => {
  if (ui.t >= T1) setTime(T0);
  ui.playing = !ui.playing;
  $("play").textContent = ui.playing ? "❚❚ Pause" : "▶ Play";
};
// log-scale speed: 0.05 h/s (3 min of biology per second, to watch mitosis) to 48 h/s
const SPEED_MIN = 0.05, SPEED_MAX = 48;
$("speed").oninput = (e) => {
  ui.speed = SPEED_MIN * (SPEED_MAX / SPEED_MIN) ** (+e.target.value / 100);
  $("sval").textContent = ui.speed < 1 ? `${(ui.speed * 60).toFixed(0)} min/s` : `${ui.speed.toFixed(1)} h/s`;
};
$("speed").value = (100 * Math.log(12 / SPEED_MIN) / Math.log(SPEED_MAX / SPEED_MIN)).toFixed(0);
$("speed").oninput({ target: $("speed") });
$("opacity").oninput = (e) => {
  ui.opacity = material.opacity = +e.target.value;
  material.depthWrite = material.opacity > 0.99;  // see inner cells through outer ones
  $("oval").textContent = material.opacity.toFixed(2);
};
$("style").onchange = (e) => { ui.style = e.target.value; $("smoothRow").hidden = ui.style !== "smooth"; };
$("smoothing").oninput = (e) => { ui.smoothing = +e.target.value; $("smval").textContent = `${ui.smoothing} µm`; };
$("wobble").onchange = (e) => { cellUniforms.uWobble.value = e.target.checked ? 0.035 : 0; };
$("showPlane").onchange = (e) => { planeMesh.visible = e.target.checked; };
$("showStroma").onchange = (e) => {
  if (fibMesh) fibMesh.visible = e.target.checked;
  if (immMesh) immMesh.visible = e.target.checked;
};
$("color").onchange = (e) => { ui.color = e.target.value; drawLegend(); ui.dirty = true; };
$("slice").oninput = (e) => { ui.slice = +e.target.value; $("pval").textContent = `${ui.slice} µm`; placePlane(); ui.dirty = true; };
for (const b of document.querySelectorAll("#plane button")) {
  b.onclick = () => {
    document.querySelectorAll("#plane button").forEach((x) => x.classList.toggle("on", x === b));
    ui.axis = +b.dataset.axis;
    placePlane();
    ui.dirty = true;
  };
}
drawLegend();
placePlane();
showSelection();

// ------------------------------------------------------------------ video export
// Composite the 3D view and the cross-section into one canvas and record it.
const rec = { recorder: null, canvas: document.createElement("canvas") };
rec.ctx = rec.canvas.getContext("2d");
function compositeFrame() {
  const c = rec.canvas, W = renderer.domElement.width, H = renderer.domElement.height;
  if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }
  rec.ctx.drawImage(renderer.domElement, 0, 0);
  const k = devicePixelRatio > 1 ? 2 : 1, iw = 300 * k;
  rec.ctx.drawImage(xs, W - iw - 16 * k, 16 * k, iw, iw);
  rec.ctx.fillStyle = "#e6edf3";
  rec.ctx.font = `${16 * k}px sans-serif`;
  rec.ctx.fillText(`t = ${ui.t.toFixed(1)} h   ${$("status").textContent}`, 16 * k, H - 16 * k);
}
$("record").onclick = () => {
  if (rec.recorder) { rec.recorder.stop(); return; }
  const type = ["video/webm;codecs=vp9", "video/webm", "video/mp4"].find((t) => MediaRecorder.isTypeSupported(t));
  const chunks = [];
  rec.recorder = new MediaRecorder(rec.canvas.captureStream(30), { mimeType: type, videoBitsPerSecond: 8e6 });
  rec.recorder.ondataavailable = (e) => chunks.push(e.data);
  rec.recorder.onstop = () => {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob(chunks, { type }));
    a.download = `tme_4d.${type.startsWith("video/mp4") ? "mp4" : "webm"}`;
    a.click();
    rec.recorder = null;
    $("record").textContent = "● Record video";
  };
  compositeFrame();
  rec.recorder.start();
  $("record").textContent = "■ Stop & save";
  setTime(T0);
  ui.playing = true;
  $("play").textContent = "❚❚ Pause";
};

// ------------------------------------------------------------------ loop
function placeOrganPanel() {
  // sit just below the cross-section panel, whose height depends on the pane width
  const sec = $("section"), org = $("organs");
  if (!org || org.hidden) return;
  const top = sec.getBoundingClientRect().bottom + 12;
  org.style.top = `${top}px`;
  org.style.maxHeight = `${Math.max(innerHeight - top - 16, 120)}px`;
}

function resize() {
  renderer.setSize(innerWidth, innerHeight, false);
  placeOrganPanel();
  const db = renderer.getDrawingBufferSize(new THREE.Vector2());
  smooth.setSize(db.x, db.y);
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
}
addEventListener("resize", resize);
resize();

// handle for console inspection / scripting
window.tme = { cur, lin, parent, birth, camera, controls, ui, setTime, renderer, scene, cells, dividing,
  computeState, updateInstances, drawSection, tc, tcur, computeT, anims, triggerEvents, updateT, updateGranules,
  select: (id) => { ui.selected = id; ui.dirty = true; showSelection(); },
  refresh: () => { computeState(ui.t); computeT(ui.t); updateInstances(); drawSection(); updateCloneCount(); } };  // sync, for scripts

function drawOrgans(t) {
  if (!organs) return;
  let k = 0;
  while (k < organs.t.length - 1 && organs.t[k + 1] <= t) k++;
  const rows = organs.rows[k];
  const fmt = (n) => n >= 1e9 ? (n / 1e9).toFixed(1) + "B"
    : n >= 1e6 ? (n / 1e6).toFixed(1) + "M"
    : n >= 1e3 ? (n / 1e3).toFixed(1) + "k" : Math.round(n).toLocaleString();
  let html = "", lastOrgan = null;
  for (const [org, stage, n] of rows) {
    if (org !== lastOrgan) { html += `<tr class="hdr"><td colspan="2">${org}</td></tr>`; lastOrgan = org; }
    html += `<tr><td class="s">${stage}</td><td class="n">${fmt(n)}</td></tr>`;
  }
  $("organTable").innerHTML = html;
}
if (organs) {
  $("organs").hidden = false;
  requestAnimationFrame(placeOrganPanel);
  $("organAge").textContent = `· age ${organs.age}`;
  $("organMeta").innerHTML =
    `Draining basin: ${organs.node_basin} (~${organs.node_distance_cm} cm)<br>` +
    `DC transit to node: ${organs.dc_transit_h} h · cognate precursors ` +
    `${organs.precursor_frequency.toExponential(0)} of repertoire`;
}

let last = performance.now();
function loop(now) {
  const dt = Math.min((now - last) / 1000, 0.1);
  last = now;
  if (ui.playing) {
    const prev = ui.t;
    setTime(ui.t + dt * ui.speed);
    triggerEvents(prev, ui.t);
    if (ui.t >= T1) {
      ui.playing = false;
      $("play").textContent = "▶ Play";
      if (rec.recorder) setTimeout(() => rec.recorder && rec.recorder.stop(), 300);
    }
  }
  if (ui.dirty) {
    computeState(ui.t);
    computeT(ui.t);
    updateInstances();
    drawSection();
    updateCloneCount();
    drawOrgans(ui.t);
    let nM = 0, nVis = 0;
    for (let i = 0; i < cur.n; i++) { nM += cur.phase[i] === 2 && cur.r[i] > 0; nVis += cur.r[i] > 0; }
    $("tval").textContent = `${ui.t.toFixed(2)} h (day ${(ui.t / 24).toFixed(1)})`;
    const nT = { 0: 0, 1: 0, 3: 0 };
    for (let i = 0; i < tc.n; i++) if (tcur.vis[i]) nT[tcur.state[i]]++;
    let killed = 0;
    for (let q = 0; q < tc.deaths.n; q++) killed += tc.deaths.tk[q] <= ui.t;
    $("status").textContent = `${nVis.toLocaleString()} cells · ${nM.toLocaleString()} in mitosis`
      + (tc.n ? ` · T cells: ${nT[0]} naive, ${nT[1]} anergic, ${nT[3]} effector · ${killed} killed` : "");
    ui.dirty = false;
  }
  controls.update();
  updateT(now);
  updateGranules(now);
  updateMarker();
  cellUniforms.uTime.value = now / 1000;
  if (ui.style === "smooth") {
    smooth.render(scene, camera, [cells, dividing], { smoothing: ui.smoothing, opacity: ui.opacity,
                  lightDir: sun.position.clone().normalize(), clearColor: BG });
  } else {
    renderer.setClearColor(BG, 1);
    renderer.render(scene, camera);
  }
  if (rec.recorder) compositeFrame();
  requestAnimationFrame(loop);
}
requestAnimationFrame(loop);
