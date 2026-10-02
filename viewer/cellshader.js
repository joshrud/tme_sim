// Cell shading: per-cell shape deformation (membrane wobble, cleavage furrow during division)
// and an optional screen-space "smooth tissue" renderer that turns the cells into one soft,
// continuous surface. Both are visual only - nothing here changes the simulation.
import * as THREE from "three";

// Shared vertex code. Geometry is a unit sphere, so position == normal.
//   aDiv.xyz = division axis (unit), aDiv.w = split progress 0..1 (0 = not dividing)
//   aSeed    = per-cell phase offset for the wobble
const DEFORM_GLSL = /* glsl */ `
attribute vec4 aDiv;
attribute float aSeed;
uniform float uTime;
uniform float uWobble;

vec3 deformCell(vec3 p) {
  // membrane wobble: product of two slow travelling waves over the cell surface
  float w = sin(dot(p, vec3(2.9, 2.1, 1.3)) + uTime * 1.1 + aSeed)
          * sin(dot(p, vec3(-1.6, 2.4, 2.8)) + uTime * 0.8 + aSeed * 1.7);
  p *= 1.0 + uWobble * w;
  // dividing: elongate along the axis and pinch a furrow at the equator as cytokinesis nears
  float s = aDiv.w;
  if (s > 0.0) {
    vec3 ax = aDiv.xyz;
    float a = dot(p, ax);
    vec3 lat = (p - a * ax) * (1.0 - 0.92 * pow(s, 1.5) * exp(-a * a / 0.06));
    p = ax * a * (1.0 + 0.65 * s) + lat;
  }
  return p;
}

vec3 deformedNormal(vec3 n) {  // normal of the deformed surface by finite differences
  vec3 t1 = normalize(cross(n, abs(n.y) < 0.99 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0)));
  vec3 t2 = cross(n, t1);
  vec3 p0 = deformCell(n);
  vec3 N = normalize(cross(deformCell(normalize(n + 0.02 * t1)) - p0,
                           deformCell(normalize(n + 0.02 * t2)) - p0));
  return dot(N, p0) < 0.0 ? -N : N;
}
`;

export const cellUniforms = { uTime: { value: 0 }, uWobble: { value: 0.035 } };

// Lit material for the normal sphere view, with the deformation patched in.
export function makeCellMaterial() {
  const mat = new THREE.MeshLambertMaterial({ transparent: true, opacity: 1 });
  mat.onBeforeCompile = (shader) => {
    Object.assign(shader.uniforms, cellUniforms);
    shader.vertexShader = DEFORM_GLSL + shader.vertexShader
      .replace("#include <beginnormal_vertex>", "vec3 objectNormal = deformedNormal(normalize(position));")
      .replace("#include <begin_vertex>", "vec3 transformed = deformCell(position);");
  };
  return mat;
}

// Instanced sphere mesh with the per-cell attributes the shader expects.
export function makeCellMesh(detail, capacity, material) {
  const mesh = new THREE.InstancedMesh(new THREE.IcosahedronGeometry(1, detail), material, capacity);
  mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
  mesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(3 * capacity), 3);
  mesh.geometry.setAttribute("aDiv", new THREE.InstancedBufferAttribute(new Float32Array(4 * capacity), 4));
  mesh.geometry.setAttribute("aSeed", new THREE.InstancedBufferAttribute(new Float32Array(capacity), 1));
  mesh.frustumCulled = false;
  mesh.layers.set(1);  // cells live on layer 1 so the smooth renderer can draw them separately
  return mesh;
}

// ------------------------------------------------------------------ smooth tissue renderer
// Screen-space surface smoothing (as used for particle fluids):
//   1. draw cells as (color, view depth) into a float target
//   2. blur color + depth with a depth-aware (bilateral) separable filter
//   3. rebuild normals from the smoothed depth and light the resulting surface
const QUAD_VS = /* glsl */ `varying vec2 vUv; void main() { vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }`;

const PREPASS = {
  vertexShader: DEFORM_GLSL + /* glsl */ `
    varying vec3 vColor; varying float vViewZ;
    void main() {
      vec4 mv = modelViewMatrix * instanceMatrix * vec4(deformCell(position), 1.0);
      vViewZ = -mv.z; vColor = instanceColor;
      gl_Position = projectionMatrix * mv;
    }`,
  fragmentShader: /* glsl */ `
    varying vec3 vColor; varying float vViewZ;
    void main() { gl_FragColor = vec4(vColor, vViewZ); }`,
};

const BLUR_FS = /* glsl */ `
  uniform sampler2D tSrc; uniform vec2 uDir; uniform float uProj; uniform float uSmooth;
  varying vec2 vUv;
  void main() {
    vec4 c = texture2D(tSrc, vUv);
    if (c.a <= 0.0) { gl_FragColor = vec4(0.0); return; }
    float rpx = clamp(uSmooth * uProj / c.a, 1.0, 40.0);  // world radius -> pixels at this depth
    float sigma = rpx * 0.5;
    vec3 col = vec3(0.0); float z = 0.0, wsum = 0.0;
    for (int i = -12; i <= 12; i++) {
      float o = float(i) * rpx / 12.0;
      vec4 s = texture2D(tSrc, vUv + uDir * o);
      if (s.a <= 0.0) continue;
      float dz = (s.a - c.a) / (1.5 * uSmooth);  // don't blend across depth jumps (silhouettes)
      float w = exp(-o * o / (2.0 * sigma * sigma)) * exp(-dz * dz);
      col += s.rgb * w; z += s.a * w; wsum += w;
    }
    gl_FragColor = vec4(col / wsum, z / wsum);
  }`;

const COMPOSITE_FS = /* glsl */ `
  uniform sampler2D tSmooth; uniform vec2 uTexel; uniform vec2 uTanHalf;
  uniform vec3 uLight; uniform float uOpacity; uniform float uProjA; uniform float uProjB;
  varying vec2 vUv;
  vec3 viewPos(vec2 uv, float z) { return vec3((uv * 2.0 - 1.0) * uTanHalf * z, -z); }
  float depthAt(vec2 uv) { return texture2D(tSmooth, uv).a; }
  vec3 edge(vec3 P, vec2 uv, vec2 d) {  // smaller one-sided difference avoids silhouette artifacts
    float zp = depthAt(uv + d), zm = depthAt(uv - d);
    vec3 a = zp > 0.0 ? viewPos(uv + d, zp) - P : vec3(1e9);
    vec3 b = zm > 0.0 ? P - viewPos(uv - d, zm) : vec3(1e9);
    return abs(a.z) < abs(b.z) ? a : b;
  }
  void main() {
    vec4 c = texture2D(tSmooth, vUv);
    if (c.a <= 0.0) discard;
    vec3 P = viewPos(vUv, c.a);
    vec3 N = normalize(cross(edge(P, vUv, vec2(uTexel.x, 0.0)), edge(P, vUv, vec2(0.0, uTexel.y))));
    if (N.z < 0.0) N = -N;
    vec3 V = normalize(-P), H = normalize(uLight + V);
    float diff = max(dot(N, uLight), 0.0);
    float spec = pow(max(dot(N, H), 0.0), 48.0) * 0.35;
    float rim = pow(1.0 - max(dot(N, V), 0.0), 3.0) * 0.25;
    gl_FragColor = vec4(c.rgb * (0.45 + 0.75 * diff) + spec + rim * vec3(0.8, 0.9, 1.0), uOpacity);
    float ndc = (uProjA * -c.a + uProjB) / c.a;  // write real depth so overlays sort correctly
    gl_FragDepth = ndc * 0.5 + 0.5;
  }`;

export class SmoothTissueRenderer {
  constructor(renderer) {
    this.renderer = renderer;
    const opt = { type: THREE.FloatType, minFilter: THREE.NearestFilter, magFilter: THREE.NearestFilter };
    this.rtScene = new THREE.WebGLRenderTarget(1, 1, { ...opt, depthBuffer: true });
    this.rtA = new THREE.WebGLRenderTarget(1, 1, { ...opt, depthBuffer: false });
    this.rtB = new THREE.WebGLRenderTarget(1, 1, { ...opt, depthBuffer: false });
    this.prepass = new THREE.ShaderMaterial({ ...PREPASS, uniforms: cellUniforms });
    this.blur = new THREE.ShaderMaterial({
      vertexShader: QUAD_VS, fragmentShader: BLUR_FS, depthTest: false, depthWrite: false,
      uniforms: { tSrc: { value: null }, uDir: { value: new THREE.Vector2() }, uProj: { value: 1 }, uSmooth: { value: 10 } },
    });
    this.composite = new THREE.ShaderMaterial({
      vertexShader: QUAD_VS, fragmentShader: COMPOSITE_FS, transparent: true,
      depthTest: true, depthWrite: true, depthFunc: THREE.AlwaysDepth,
      uniforms: { tSmooth: { value: this.rtB.texture }, uTexel: { value: new THREE.Vector2() },
                  uTanHalf: { value: new THREE.Vector2() }, uLight: { value: new THREE.Vector3() },
                  uOpacity: { value: 1 }, uProjA: { value: 0 }, uProjB: { value: 0 } },
    });
    this.quad = new THREE.Mesh(new THREE.PlaneGeometry(2, 2));
    this.quad.frustumCulled = false;
    this.quadScene = new THREE.Scene().add(this.quad);
    this.quadCam = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
    this.size = new THREE.Vector2();
  }

  setSize(w, h) {
    for (const rt of [this.rtScene, this.rtA, this.rtB]) rt.setSize(w, h);
    this.size.set(w, h);
  }

  // smoothing: world radius (um); lightDir: world-space direction toward the light
  render(scene, camera, cellMeshes, { smoothing, opacity, lightDir, clearColor }) {
    const r = this.renderer, { x: w, y: h } = this.size;
    const tanHalf = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);

    // 1. color + view depth of the cells
    const saved = cellMeshes.map((m) => m.material);
    cellMeshes.forEach((m) => { m.material = this.prepass; });
    camera.layers.set(1);
    r.setRenderTarget(this.rtScene);
    r.setClearColor(0x000000, 0);
    r.clear();
    r.render(scene, camera);
    cellMeshes.forEach((m, k) => { m.material = saved[k]; });

    // 2. bilateral blur, horizontal then vertical, twice
    this.quad.material = this.blur;
    this.blur.uniforms.uProj.value = h / (2 * tanHalf);
    this.blur.uniforms.uSmooth.value = smoothing;
    const H = [1 / w, 0], V = [0, 1 / h];
    for (const [src, dst, dir] of [[this.rtScene, this.rtA, H], [this.rtA, this.rtB, V],
                                   [this.rtB, this.rtA, H], [this.rtA, this.rtB, V]]) {
      this.blur.uniforms.tSrc.value = src.texture;
      this.blur.uniforms.uDir.value.set(...dir);
      r.setRenderTarget(dst);
      r.render(this.quadScene, this.quadCam);
    }

    // 3. lit surface to the screen, then the non-cell overlays (plane, marker) on top
    const u = this.composite.uniforms, P = camera.projectionMatrix.elements;
    u.uTexel.value.set(1 / w, 1 / h);
    u.uTanHalf.value.set(tanHalf * camera.aspect, tanHalf);
    u.uLight.value.copy(lightDir).transformDirection(camera.matrixWorldInverse);
    u.uOpacity.value = opacity;
    u.uProjA.value = P[10];
    u.uProjB.value = P[14];
    this.quad.material = this.composite;
    r.setRenderTarget(null);
    r.setClearColor(clearColor, 1);
    r.clear();
    r.render(this.quadScene, this.quadCam);
    camera.layers.set(0);
    r.autoClear = false;
    r.render(scene, camera);
    r.autoClear = true;
    camera.layers.enableAll();
  }
}
