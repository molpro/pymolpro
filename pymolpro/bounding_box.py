"""
Axis-aligned bounding box search for an expensive, batch-vectorized 3D
inside/outside oracle.

Given:
  - inside(points) -> bool array, where `points` is an (N, 3) numpy array
    and the oracle is much more efficient when called on a big batch of
    points at once than when called once per point -- and, crucially,
    each call may carry a meaningful *fixed* cost independent of how many
    points are in the batch (e.g. parsing/setup work redone every call)
  - a rough estimate of the object's center and size

finds an axis-aligned box that contains the whole object, including any
disjoint pieces, using as few oracle *calls* as practical -- total point
count is a secondary concern. The box is deliberately allowed to be
larger than the minimal bounding box.

Used by :meth:`pymolpro.orbital.Orbital.cube_data` to locate the region
where an orbital's amplitude exceeds a threshold. There, the "inside"
oracle is ``|orbital value| >= threshold`` evaluated with
:meth:`Orbital.evaluate`, which is already vectorized over an (N, 3)
point array. evaluateOrbitals() also redoes a non-trivial amount of XML
basis-set bookkeeping on every call, independent of how many points are
passed in, which is exactly the "meaningful fixed per-call cost" this
module is designed around -- hence the emphasis on batching everything
possible into as few calls as practical, rather than minimizing total
point count.

Algorithm
---------
1. Seed: sample a grid over the rough region (one batched call) to find
   points known to be inside the object. This is what lets disjoint
   pieces get picked up, as long as they lie roughly within the initial
   estimate.
2. Grow, per pass:
   a. Ladder: for every face, precompute a whole exponentially-spaced
      ladder of candidate offsets (step, step*growth_factor,
      step*growth_factor^2, ...) up front, and probe all of them, for
      all 6 faces, in a single batched call. From the hit/miss pattern
      along each face's ladder we get a bracket -- last hit, first miss
      -- in one round trip instead of one call per doubling.
   b. Bisect: a handful of bisection rounds, each batched across all 6
      faces at once, tighten each bracket to "first confirmed empty"
      (bisection always keeps the box on the confirmed-empty side, so it
      can only overshoot the true boundary, never cut into the object).
   c. Verify: one more batched call, denser grid, re-probing each face's
      final bound in case the coarser grid stepped over a narrow
      protrusion; any face that finds something there reopens (rare) for
      another ladder+bisect cycle at the denser resolution.
3. Pad: add a final fixed-fraction safety margin, since the box doesn't
   need to be minimal.

A handful of calls per pass (roughly 2 + refine_steps, typically 6-8),
independent of the object's size or distance from the initial guess --
that distance is absorbed into the ladder's exponential spacing rather
than into the number of round trips.
"""

from dataclasses import dataclass
from typing import Callable
import numpy as np

# points: (N, 3) float array -> (N,) bool array
InsideFn = Callable[[np.ndarray], np.ndarray]


@dataclass
class BBoxResult:
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    zmin: float
    zmax: float
    points_evaluated: int
    oracle_calls: int
    seeds_found: int

    def as_tuple(self):
        return (self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax)


class CountingOracle:
    """Wraps the user's batched inside() function and counts usage."""

    def __init__(self, fn: InsideFn):
        self.fn = fn
        self.points_evaluated = 0
        self.calls = 0

    def __call__(self, points: np.ndarray) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError(f"expected an (N, 3) array, got shape {points.shape}")
        self.calls += 1
        self.points_evaluated += points.shape[0]
        result = np.asarray(self.fn(points), dtype=bool)
        if result.shape != (points.shape[0],):
            raise ValueError(
                f"inside() must return a ({points.shape[0]},) bool array, "
                f"got shape {result.shape}"
            )
        return result


AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
OTHER_AXES = {"x": ("y", "z"), "y": ("x", "z"), "z": ("x", "y")}
FACES = [("x", -1), ("x", +1), ("y", -1), ("y", +1), ("z", -1), ("z", +1)]


def _linspace(a, b, n):
    if n == 1:
        return np.array([(a + b) / 2])
    return np.linspace(a, b, n)


def _grid_3d(cx, hx, cy, hy, cz, hz, n):
    """Full (n^3, 3) grid of points spanning a box."""
    xs = _linspace(cx - hx, cx + hx, n)
    ys = _linspace(cy - hy, cy + hy, n)
    zs = _linspace(cz - hz, cz + hz, n)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    return np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)


def _face_points(axis, u_axis, v_axis, coord, u_center, u_half, v_center, v_half, n):
    """(n^2, 3) grid of points spanning one axis-aligned face plane."""
    us = _linspace(u_center - u_half, u_center + u_half, n)
    vs = _linspace(v_center - v_half, v_center + v_half, n)
    U, V = np.meshgrid(us, vs, indexing="ij")
    pts = np.empty((U.size, 3))
    pts[:, AXIS_INDEX[axis]] = coord
    pts[:, AXIS_INDEX[u_axis]] = U.ravel()
    pts[:, AXIS_INDEX[v_axis]] = V.ravel()
    return pts


def _probe_batch(faces_offsets, grid_n, probe_points_fn, oracle):
    """Evaluate several candidate offsets for several faces in a single
    oracle call. `faces_offsets` is {face: [coord, ...]}. Returns
    {face: [bool, ...]} aligned with the input offsets."""
    parts = []
    slices = {}
    pos = 0
    for face, offsets in faces_offsets.items():
        face_slices = []
        for coord in offsets:
            pts = probe_points_fn(face, coord, grid_n)
            start, end = pos, pos + pts.shape[0]
            parts.append(pts)
            face_slices.append((start, end))
            pos = end
        slices[face] = face_slices
    if not parts:
        return {}
    all_hits = oracle(np.concatenate(parts, axis=0))
    return {face: [bool(all_hits[s:e].any()) for s, e in idx] for face, idx in slices.items()}


def _grow_pass(
    oracle, bounds, faces, step0, face_grid_n, dense_n,
    growth_factor, refine_steps, ladder_size, probe_points_fn,
    max_ladder_extensions=6,
):
    """Grow and tighten `faces` (a subset of the 6 box faces) by one
    pass, mutating `bounds` in place. Every stage is batched across all
    of `faces` simultaneously, so the number of oracle calls is roughly
    constant (about 2 + refine_steps) regardless of how far the faces
    have to travel or how many faces are being grown at once."""

    # ---- Stage A: exponential ladder, one batched call ----
    ladders = {}
    for face in faces:
        axis, sign = face
        step = step0[axis]
        cur = bounds[face]
        offsets = []
        for _ in range(ladder_size):
            cur = cur + sign * step
            offsets.append(cur)
            step *= growth_factor
        ladders[face] = offsets

    hits = _probe_batch(ladders, face_grid_n, probe_points_fn, oracle)

    lo, hi, unresolved, last_step = {}, {}, [], {}
    for face in faces:
        axis, sign = face
        prev = bounds[face]
        bracket_found = False
        for coord, hit in zip(ladders[face], hits[face]):
            if hit:
                prev = coord
            else:
                lo[face], hi[face] = prev, coord
                bracket_found = True
                break
        if not bracket_found:
            lo[face] = prev
            last_step[face] = step0[axis] * (growth_factor ** ladder_size)
            unresolved.append(face)

    # ---- Stage A2 (rare): the object is bigger than the whole ladder
    # reached. Keep doubling, one more offset per still-unresolved face
    # per round, batched across all of them.
    extensions = 0
    while unresolved and extensions < max_ladder_extensions:
        extensions += 1
        next_offsets = {face: [lo[face] + face[1] * last_step[face]] for face in unresolved}
        hits2 = _probe_batch(next_offsets, face_grid_n, probe_points_fn, oracle)
        still_unresolved = []
        for face in unresolved:
            coord = next_offsets[face][0]
            if hits2[face][0]:
                lo[face] = coord
                last_step[face] *= growth_factor
                still_unresolved.append(face)
            else:
                hi[face] = coord
        unresolved = still_unresolved
    for face in unresolved:  # gave up -- treat current reach as the boundary
        hi[face] = lo[face] + face[1] * last_step[face]

    # ---- Stage B: bisection, batched across all faces per round ----
    for _ in range(refine_steps):
        mids = {face: (lo[face] + hi[face]) / 2 for face in faces}
        hits3 = _probe_batch({f: [mids[f]] for f in faces}, face_grid_n, probe_points_fn, oracle)
        for face in faces:
            if hits3[face][0]:
                lo[face] = mids[face]
            else:
                hi[face] = mids[face]

    for face in faces:
        bounds[face] = hi[face]

    # ---- Stage C: one denser verification call; reopen any face where
    # it catches something the coarser grid stepped over ----
    hits4 = _probe_batch({f: [hi[f]] for f in faces}, dense_n, probe_points_fn, oracle)
    reopen = [face for face in faces if hits4[face][0]]
    if reopen:
        _grow_pass(
            oracle, bounds, reopen, step0, dense_n, dense_n,
            growth_factor, refine_steps, ladder_size, probe_points_fn,
            max_ladder_extensions,
        )


def find_bounding_box(
    inside: InsideFn,
    center: tuple,
    size,
    seed_grid_n: int = 5,
    face_grid_n: int = 4,
    initial_step_frac: float = 0.1,
    growth_factor: float = 2.0,
    ladder_size: int = 10,
    refine_steps: int = 4,
    verify_multiplier: int = 2,
    margin_frac: float = 0.05,
    max_seed_expansions: int = 6,
    growth_passes: int = 2,
) -> BBoxResult:
    """
    Parameters
    ----------
    inside : function(points: (N, 3) ndarray) -> (N,) bool ndarray
        Expensive, batch-vectorized membership test. Assume each call has
        a meaningful fixed cost on top of its per-point cost -- that's
        what this function is designed to economize on.
    center : (x, y, z)
        Rough center of the object.
    size : float or (sx, sy, sz)
        Rough overall size (diameter) of the object along each axis.
    seed_grid_n : points per axis for the initial seed grid. Raise this
        if the object could have several small, separated pieces within
        the rough region.
    face_grid_n : points per axis for the grid sampled on each growing
        face. Raise this if pieces could be narrow/off-center -- a
        protrusion much thinner than the face grid spacing can be missed.
    initial_step_frac : first ladder step, as a fraction of `size`.
    growth_factor : ladder step multiplier (rung k is at distance
        initial_step_frac*size * growth_factor**k from the starting
        bound).
    ladder_size : number of rungs precomputed and probed in the single
        Stage-A call. With the defaults this reaches roughly
        0.1 * 2**10 ~= 100x the initial size guess, which is normally
        far more than enough; on the rare object that's still bigger,
        a slower (but still batched) fallback extends the ladder.
    refine_steps : number of bisection rounds used to pull each face back
        in from "first confirmed empty" to a tight-but-still-safe bound.
        Each round is one oracle call, batched across every active face.
    verify_multiplier : how much denser the one-time verification grid is
        relative to face_grid_n.
    margin_frac : final safety padding, as a fraction of `size`, added to
        every face after growth (cheap insurance since the box need not
        be minimal).
    max_seed_expansions : how many times to widen the seed search if the
        rough estimate turns out to miss the object entirely.
    growth_passes : number of full passes over all 6 faces. A second pass
        lets each face's bound benefit from the other faces having
        already grown once (their cross-section is only known accurately
        after the other axes have grown). Each pass costs roughly
        2 + refine_steps oracle calls.

    Returns
    -------
    BBoxResult
    """
    cx, cy, cz = center
    if isinstance(size, (int, float)):
        sx = sy = sz = float(size)
    else:
        sx, sy, sz = size

    oracle = CountingOracle(inside)

    # ---- 1. Seed (one batched call, widened/densified if needed) ----
    seeds = np.empty((0, 3))
    radius = (sx, sy, sz)
    n = seed_grid_n
    expansions = 0
    while seeds.shape[0] == 0 and expansions <= max_seed_expansions:
        rx, ry, rz = radius
        pts = _grid_3d(cx, rx / 2, cy, ry / 2, cz, rz / 2, n)
        hits = oracle(pts)
        seeds = pts[hits]
        if seeds.shape[0] == 0:
            # Widen and densify -- a coarser grid over a larger volume
            # would only get less likely to hit a small object.
            radius = (rx * 1.8, ry * 1.8, rz * 1.8)
            n += 2
            expansions += 1

    if seeds.shape[0] == 0:
        raise RuntimeError(
            f"No point inside the object found after {max_seed_expansions} "
            f"seed expansions (final search radius {radius}, grid {n}^3). "
            f"Check the rough center/size estimate."
        )

    # ---- 2. Initial box from seeds ----
    bounds = {
        ("x", -1): float(seeds[:, 0].min()), ("x", +1): float(seeds[:, 0].max()),
        ("y", -1): float(seeds[:, 1].min()), ("y", +1): float(seeds[:, 1].max()),
        ("z", -1): float(seeds[:, 2].min()), ("z", +1): float(seeds[:, 2].max()),
    }
    nominal = {"x": sx, "y": sy, "z": sz}
    for axis in ("x", "y", "z"):
        if bounds[(axis, +1)] - bounds[(axis, -1)] < 1e-9:
            pad = 0.05 * nominal[axis]
            bounds[(axis, -1)] -= pad
            bounds[(axis, +1)] += pad

    # ---- 3. Grow all 6 faces, a handful of batched calls per pass ----
    step0 = {"x": sx * initial_step_frac, "y": sy * initial_step_frac, "z": sz * initial_step_frac}
    dense_n = face_grid_n * verify_multiplier

    def face_probe_points(face, coord, grid_n):
        axis, _ = face
        u_axis, v_axis = OTHER_AXES[axis]
        u_center = (bounds[(u_axis, -1)] + bounds[(u_axis, +1)]) / 2
        u_half = max((bounds[(u_axis, +1)] - bounds[(u_axis, -1)]) / 2, 1e-9)
        v_center = (bounds[(v_axis, -1)] + bounds[(v_axis, +1)]) / 2
        v_half = max((bounds[(v_axis, +1)] - bounds[(v_axis, -1)]) / 2, 1e-9)
        return _face_points(axis, u_axis, v_axis, coord, u_center, u_half, v_center, v_half, grid_n)

    for _ in range(growth_passes):
        _grow_pass(
            oracle, bounds, FACES, step0, face_grid_n, dense_n,
            growth_factor, refine_steps, ladder_size, face_probe_points,
        )

    # ---- 4. Pad with safety margin ----
    xmin = bounds[("x", -1)] - margin_frac * sx
    xmax = bounds[("x", +1)] + margin_frac * sx
    ymin = bounds[("y", -1)] - margin_frac * sy
    ymax = bounds[("y", +1)] + margin_frac * sy
    zmin = bounds[("z", -1)] - margin_frac * sz
    zmax = bounds[("z", +1)] + margin_frac * sz

    return BBoxResult(
        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
        points_evaluated=oracle.points_evaluated, oracle_calls=oracle.calls,
        seeds_found=seeds.shape[0],
    )


if __name__ == "__main__":
    # Demo: two disjoint spheres, oracle vectorized with plain numpy ops
    # (this is the shape of speedup a real batched oracle -- e.g. a GPU
    # simulation or a compiled/parallel routine -- would get).
    def two_spheres(points: np.ndarray) -> np.ndarray:
        x, y, z = points[:, 0], points[:, 1], points[:, 2]
        in_a = x * x + y * y + z * z <= 4  # sphere A: center (0,0,0), r=2
        dx, dy, dz = x - 6, y - 3, z + 2
        in_b = dx * dx + dy * dy + dz * dz <= 1.5 ** 2  # sphere B: disjoint piece
        return in_a | in_b

    # Rough estimate: object sits around (3, 1.5, -1) with overall size ~9
    # (covers both blobs plus slack).
    result = find_bounding_box(two_spheres, center=(3, 1.5, -1), size=9)

    print("Computed box:")
    print(f"  x: [{result.xmin:.3f}, {result.xmax:.3f}]")
    print(f"  y: [{result.ymin:.3f}, {result.ymax:.3f}]")
    print(f"  z: [{result.zmin:.3f}, {result.zmax:.3f}]")
    print(f"Points evaluated: {result.points_evaluated}")
    print(f"Oracle calls (batches): {result.oracle_calls}")
    print(f"Seed points found inside: {result.seeds_found}")

    print("\nTrue extent (analytic, for comparison):")
    print("  x: [-2.000, 7.500]")
    print("  y: [-2.000, 4.500]")
    print("  z: [-3.500, 2.000]")
