"""Compiled inner loops.

Most of the heavy numerics in this model already run as compiled C inside SciPy and NumPy
(KD-tree queries, sparse matrix-vector products, BLAS). Those are not worth reimplementing.
What is *not* covered by a library call is the per-edge force accumulation in the mechanical
relaxation: a scatter-add over the contact graph. In NumPy that needs several full-size
temporaries plus six `bincount` passes; as a compiled kernel it is one pass with no allocation.

numba compiles this to machine code via LLVM at first call, so there is no separate build step
and no C toolchain requirement. If numba is unavailable the NumPy path is used instead, so the
model runs either way - `HAVE_NUMBA` reports which.
"""
import numpy as np

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except ImportError:  # pragma: no cover - exercised only where numba is absent
    HAVE_NUMBA = False

    def njit(*a, **k):
        def wrap(fn):
            return fn
        return wrap

    prange = range


@njit(cache=True, fastmath=True)
def relax_forces(pos, r, edges, packing, repulsion, adhesion, out):
    """Accumulate pairwise overlap-relaxation displacements into `out` (n, 3).

    For each contact edge, cells that overlap are pushed apart by `repulsion` of the overlap
    and cells merely near each other are pulled together by `adhesion` of the gap. Equal and
    opposite, so the centre of mass is preserved.
    """
    out[:] = 0.0
    for e in range(edges.shape[0]):
        a = edges[e, 0]
        b = edges[e, 1]
        dx = pos[b, 0] - pos[a, 0]
        dy = pos[b, 1] - pos[a, 1]
        dz = pos[b, 2] - pos[a, 2]
        d = np.sqrt(dx * dx + dy * dy + dz * dz)
        if d < 1e-6:
            d = 1e-6
        rest = packing * (r[a] + r[b])
        gap = d - rest
        frac = repulsion if gap < 0.0 else adhesion
        s = 0.5 * frac * gap / d
        ux = dx * s
        uy = dy * s
        uz = dz * s
        out[a, 0] += ux
        out[a, 1] += uy
        out[a, 2] += uz
        out[b, 0] -= ux
        out[b, 1] -= uy
        out[b, 2] -= uz
    return out


@njit(cache=True, fastmath=True)
def edge_lengths(pos, edges, out):
    """Euclidean length of every contact edge, without building intermediate arrays."""
    for e in range(edges.shape[0]):
        a = edges[e, 0]
        b = edges[e, 1]
        dx = pos[b, 0] - pos[a, 0]
        dy = pos[b, 1] - pos[a, 1]
        dz = pos[b, 2] - pos[a, 2]
        out[e] = np.sqrt(dx * dx + dy * dy + dz * dz)
    return out


def warmup():
    """Trigger compilation once so the cost does not land inside a timed step."""
    if not HAVE_NUMBA:
        return False
    pos = np.zeros((2, 3))
    pos[1, 0] = 1.0
    edges = np.zeros((1, 2), np.int64)
    edges[0, 1] = 1
    relax_forces(pos, np.ones(2), edges, 0.86, 0.25, 0.1, np.zeros((2, 3)))
    edge_lengths(pos, edges, np.zeros(1))
    return True
