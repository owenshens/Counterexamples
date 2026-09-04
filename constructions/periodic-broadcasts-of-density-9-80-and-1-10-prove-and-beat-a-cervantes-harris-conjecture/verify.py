#!/usr/bin/env python3
"""verify.py -- machine check of the exhibited broadcasts and the hand arithmetic of the note
"Periodic Broadcasts of Density 1/10, 1/8 and 5/36: Confirming and Strengthening Three Conjectured
Density Bounds of Cervantes and Harris".

What it bears on: the upper-bound theorem, i.e. that the lattice-periodic sets W5, W6, W7
(the note's T_5, T_6, T_7) are (4,2), (4,3), (4,4) broadcasts of densities 1/10, 1/8, 5/36,
together with the presentation of the graph and the shell sizes. The lower bound r/C(r) and the
general (3,1)-implies-(4,2) implication are NOT verified here; only the integers C(2), C(3), C(4),
the fractions r/C(r), and one (3,1) witness re-checked as a (4,2) broadcast are computed.

Subject: the infinite truncated square tiling graph H = H_{oo,oo} (the Archimedean 4.8.8 net) and
(t,r) broadcast domination on it, as defined by Cervantes and Harris.

Contract:
  * Python 3.9+, STANDARD LIBRARY ONLY. No numpy, no sympy, no networkx, no external data file.
  * Exact arithmetic only: integers and fractions.Fraction. No float is ever compared.
  * Every object it checks is defined by a congruence written out in this file. W5, W6 and W7 are
    transcribed from the note's own statements of T_5, T_6 and T_7; the remaining predicates
    (in_W1, in_W2, in_W3, in_W4, in_T31) are further patterns written out here and are not
    statements of the note.
  * One `PASS <name> [detail]` line per check; closes with `VERDICT: ALL <n> CHECKS PASS`;
    exits 0 if and only if every check passed.

What it does NOT re-run is printed at the end under `NOT RE-RUN:`.
"""

import sys
from collections import deque
from fractions import Fraction

# ----------------------------------------------------------------------------------------------
# check bookkeeping
# ----------------------------------------------------------------------------------------------
_N = 0
_F = 0


def check(name, cond, detail=''):
    global _N, _F
    _N += 1
    if cond:
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _F += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))
    return bool(cond)


# ----------------------------------------------------------------------------------------------
# 1. THE GRAPH, from the note's normal form and its six adjacency rules
#
# A vertex is a triple (a, x, y) with a in Z_4 and (x, y) in Z^2, written a@(x,y).  The six rules
# printed in the note are
#     0@(x,y) ~ 1@(x,y)   1@(x,y) ~ 2@(x,y)   2@(x,y) ~ 3@(x,y)   3@(x,y) ~ 0@(x,y)
#     0@(x,y) ~ 2@(x,y+1)                     1@(x,y) ~ 3@(x+1,y)
# ----------------------------------------------------------------------------------------------
def nbrs(v):
    a, x, y = v
    if a == 0:
        return ((1, x, y), (3, x, y), (2, x, y + 1))
    if a == 1:
        return ((0, x, y), (2, x, y), (3, x + 1, y))
    if a == 2:
        return ((1, x, y), (3, x, y), (0, x, y - 1))
    return ((2, x, y), (0, x, y), (1, x - 1, y))


def ball(v, radius):
    """{vertex: distance} for every vertex at distance <= radius of v (plain BFS, no lattice)."""
    dist = {v: 0}
    q = deque([v])
    while q:
        u = q.popleft()
        d = dist[u]
        if d == radius:
            continue
        for w in nbrs(u):
            if w not in dist:
                dist[w] = d + 1
                q.append(w)
    return dist


def shells(v, radius):
    d = ball(v, radius)
    out = []
    for k in range(radius + 1):
        out.append(sum(1 for u in d if d[u] == k))
    return tuple(out)


# --- 1a. the adjacency relation is symmetric and 3-regular ------------------------------------
WINDOW = [(a, x, y) for a in range(4) for x in range(-4, 5) for y in range(-4, 5)]
check('graph.3-regular',
      all(len(set(nbrs(v))) == 3 for v in WINDOW),
      'every vertex of a %d-vertex window has exactly 3 distinct neighbours' % len(WINDOW))
check('graph.symmetric',
      all(v in nbrs(w) for v in WINDOW for w in nbrs(v)),
      'w in N(v) implies v in N(w) throughout the window')
check('graph.no-loop', all(v not in nbrs(v) for v in WINDOW), 'no vertex is its own neighbour')

# --- 1b. translation by Z^2 is an automorphism (this is what makes a period lattice legal) -----
def shift(v, dx, dy):
    return (v[0], v[1] + dx, v[2] + dy)


check('graph.Z2-translation-is-an-automorphism',
      all(set(shift(w, dx, dy) for w in nbrs(v)) == set(nbrs(shift(v, dx, dy)))
          for v in WINDOW for (dx, dy) in ((1, 0), (0, 1), (1, 1), (-1, 2))),
      'N(v+(dx,dy)) = N(v)+(dx,dy) for (dx,dy) in {(1,0),(0,1),(1,1),(-1,2)}')

# --- 1c. the coordination data of this coordinate model -----------------------------------------
SH = {a: shells((a, 0, 0), 4) for a in range(4)}
check('shells.1-3-5-8-11-from-all-four-vertex-types',
      all(SH[a] == (1, 3, 5, 8, 11) for a in range(4)),
      'shells (d=0..4) = (1,3,5,8,11) for a=0,1,2,3'
      + '; ' + ' '.join('a=%d:%s' % (a, SH[a]) for a in range(4)))
B2 = 1 + 3 + 5
B3 = 1 + 3 + 5 + 8
check('ball-sizes', B2 == 9 and B3 == 17, '|B_2| = 1+3+5 = 9, |B_3| = 1+3+5+8 = 17')
check('shells.prefix-3-5-8-11',
      tuple(SH[0][1:]) == (3, 5, 8, 11),
      'the computed coordination prefix from a=0 is 3,5,8,11')


# --- 1d. the coordinate model really is the 4.8.8 net, not merely some cubic graph --------------
# Everything in the paper is a statement about H_{oo,oo}, so the coordinate rules must define THAT
# graph.  In the truncated square tiling every vertex lies on exactly one 4-face and two 8-faces
# and there is no triangle, no pentagon, no hexagon and no heptagon through any vertex.  Counting
# simple cycles by length through a vertex is a face census here because the girth is 4 and the
# only short cycles of a planar 3-regular tiling through a vertex are its faces.
def simple_cycles_through(v, length):
    """number of simple cycles of exactly `length` vertices through v, counted once per cycle."""
    found = set()

    def walk(path):
        u = path[-1]
        if len(path) == length:
            if v in nbrs(u):
                found.add(min(tuple(path), (path[0],) + tuple(path[:0:-1])))
            return
        for w in nbrs(u):
            if w not in path:
                walk(path + [w])

    walk([v])
    return len(found)


FACE_SAMPLE = [(a, x, y) for a in range(4) for x in range(-1, 2) for y in range(-1, 2)]
FACE_PROFILE = {v: tuple(simple_cycles_through(v, L) for L in (3, 4, 5, 6, 7, 8))
                for v in FACE_SAMPLE}
check('graph.short-cycle-profile-at-sampled-vertices',
      all(FACE_PROFILE[v] == (0, 1, 0, 0, 0, 2) for v in FACE_SAMPLE),
      'through each of %d sampled vertices: no 3-, 5-, 6- or 7-cycle, exactly one 4-cycle and '
      'exactly two 8-cycles -- i.e. one square and two octagons through each sampled vertex, the '
      'vertex figure of the truncated square tiling; no classification over Archimedean nets and '
      'no global isomorphism is checked here' % len(FACE_SAMPLE))


# ----------------------------------------------------------------------------------------------
# 2. RECEPTION, from the paper's definition
#    f(u) = sum over towers v with d(u,v) < t of (t - d(u,v));  the cutoff d < t is STRICT.
# ----------------------------------------------------------------------------------------------
def reception(u, in_T, t):
    d = ball(u, t - 1)
    return sum(t - dd for w, dd in d.items() if in_T(w))


def towers_within(u, in_T, radius):
    d = ball(u, radius)
    return sorted(w for w, dd in d.items() if in_T(w))


# ----------------------------------------------------------------------------------------------
# 3. THE PATTERNS. W5, W6, W7 are transcribed from the note; W1-W4 and T31 are defined here only.
#
# Each is (name, t, r, membership predicate, HNF lattice generators <(A,0),(B,D)>, expected
# density).  The coset representatives of Z^2 / L for L = <(A,0),(B,D)> are {(x,y) : 0<=x<A,
# 0<=y<D}, so the class set is {0,1,2,3} x those, of size 4*A*D.
# ----------------------------------------------------------------------------------------------
W1_ZERO = frozenset({0, 5, 10, 12, 14, 16, 18})
W1_ONE = frozenset({2, 7})


def in_W1(v):                                    # (4,2), density 9/80
    a, x, y = v
    c = (x - y) % 20
    return (a == 0 and c in W1_ZERO) or (a == 1 and c in W1_ONE)


def in_W5(v):                                    # (4,2), density 1/10
    a, x, y = v
    i = (x + y) % 5
    return (a == 0 and i == 0) or (a == 2 and i == 3)


def in_W2(v):                                    # (4,3), density 1/7
    a, x, y = v
    c = (x - y) % 7
    return (a == 0 and c in (0, 1, 3)) or (a == 3 and c == 5)


def in_W3(v):                                    # (4,4), density 1/6
    a, x, y = v
    return (a == 1 and x % 3 == 0) or (a == 3 and x % 3 == 2)


def in_W6(v):                                    # (4,3), density 1/8
    a, x, y = v
    c = (x - y) % 4
    return (a == 0 and c == 0) or (a == 3 and c == 2)


W7_T = frozenset({(0, 0, 0), (0, 2, 0), (0, 4, 1), (1, 0, 2), (1, 3, 0),
                  (2, 1, 1), (2, 3, 2), (2, 5, 2), (3, 2, 2), (3, 5, 0)})


def in_W7(v):                                    # (4,4), density 5/36
    a, x, y = v
    return (a, x % 6, y % 3) in W7_T


W4_T = frozenset({(0, 0, 0), (0, 1, 1), (1, 1, 0), (1, 2, 2)})


def in_W4(v):                                    # (4,2), density 1/9
    a, x, y = v
    return (a, x % 3, y % 3) in W4_T


def in_T31(v):
    """The broadcast T = {2@(x, x+4y)} u {1@(x, x+4y+2)}, written out here.

    It is checked below BOTH as a (3,1) broadcast of density 1/8 and as a (4,2) broadcast of
    density 1/8.
    """
    a, x, y = v
    c = (y - x) % 4
    return (a == 2 and c == 0) or (a == 1 and c == 2)


def reps(A, B, D):
    return [(a, x, y) for a in range(4) for x in range(A) for y in range(D)]


def audit_witness(name, t, r, in_T, gens, want_density, want_min, want_max):
    """Class count, lattice invariance, density and reception for one pattern, from its
    congruence alone."""
    (A, Bc, D) = gens
    index = A * D
    R = reps(A, Bc, D)
    nclass = len(R)
    check('%s.class-count' % name, nclass == 4 * index,
          '4 vertex types x index %d = %d classes' % (index, nclass))

    # (i) T is invariant under the stated period lattice -> its density needs no boundary term.
    gv = ((A, 0), (Bc, D))
    inv = all(in_T(v) == in_T(shift(v, dx, dy)) for v in WINDOW for (dx, dy) in gv)
    check('%s.lattice-invariant' % name, inv,
          'membership of T is unchanged by both generators (%d,%d) and (%d,%d)' % (A, 0, Bc, D))

    # (ii) the density is exactly (tower classes)/(4 * index).
    k = sum(1 for v in R if in_T(v))
    dens = Fraction(k, nclass)
    check('%s.density' % name, dens == want_density,
          '%d tower classes of %d -> %s = %s' % (k, nclass, dens, want_density))

    # (iii) reception, by fresh BFS in the INFINITE graph at every class representative.
    fs = {}
    for v in R:
        fs[v] = reception(v, in_T, t)
    lo, hi = min(fs.values()), max(fs.values())
    bad = sorted(v for v in R if fs[v] < r)
    check('%s.is-a-(%d,%d)-broadcast' % (name, t, r), not bad,
          'min f = %d >= r = %d over all %d classes; violations = %d' % (lo, r, nclass, len(bad)))
    check('%s.f-range' % name, (lo, hi) == (want_min, want_max),
          'min f = %d, max f = %d' % (lo, hi))
    return fs


print('')
print('--- W1: a pattern of density 9/80 defined in this file, giving delta_{4,2}(H) <= 9/80 '
      '(upper bound only; no matching lower bound is checked) ---')
F1 = audit_witness('W1', 4, 2, in_W1, (20, 1, 1), Fraction(9, 80), 2, 5)

# W1's full 80-entry reception table, as listed in this program (rows a, columns c). This table
# is not printed in the paper; it is transcribed here only from the congruence in_W1 above.
W1_TABLE = {
    0: [4, 2, 3, 4, 3, 4, 2, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    1: [3, 2, 4, 3, 4, 3, 2, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    2: [2, 3, 4, 2, 4, 2, 3, 4, 2, 4, 2, 5, 2, 5, 2, 5, 2, 5, 2, 5],
    3: [4, 3, 2, 4, 2, 4, 3, 2, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
}
check('W1.matches-transcribed-table',
      all(F1[(a, c, 0)] == W1_TABLE[a][c] for a in range(4) for c in range(20)),
      'all 80 entries of the W1 table listed in this program reproduced by fresh BFS')

W1_TOWER_CLASSES = sorted((a, c) for a in range(4) for c in range(20) if in_W1((a, c, 0)))
check('W1.nine-tower-classes', len(W1_TOWER_CLASSES) == 9,
      '9 of 80: ' + ' '.join('%d@%d' % p for p in W1_TOWER_CLASSES))

W1_TIGHT_LISTED_HERE = [(0, 1), (0, 6), (1, 1), (1, 6), (2, 0), (2, 3), (2, 5), (2, 8), (2, 10),
                        (2, 12), (2, 14), (2, 16), (2, 18), (3, 2), (3, 4), (3, 7), (3, 9)]
W1_TIGHT = sorted((a, c) for a in range(4) for c in range(20) if F1[(a, c, 0)] == 2)
check('W1.tight-classes', W1_TIGHT == sorted(W1_TIGHT_LISTED_HERE),
      '%d classes with f = 2 exactly, and they are the ones listed in this program' % len(W1_TIGHT))

# Four W1 classes re-derived with the tower certificates transcribed in this program.
for (u, want_f, note) in ((( 3, 0, 0), 4, '0@(0,0) at d=1 and 0@(-1,-1) at d=3'),
                          (( 0, 3, 0), 4, '1@(2,0) and 1@(3,1), both at d=2'),
                          (( 0, 4, 0), 3, 'three towers at d=3 exactly, none within d=2'),
                          (( 2, 4, 0), 4, '0@(4,-1) at d=1 and 0@(5,0) at d=3')):
    check('W1.hand-check.%d@(%d,%d)' % u, reception(u, in_W1, 4) == want_f,
          'f = %d; %s' % (want_f, note))

check('W1.outer-shell-mechanism-at-0@(4,0)',
      towers_within((0, 4, 0), in_W1, 2) == []
      and sorted(towers_within((0, 4, 0), in_W1, 3)) == [(0, 4, -1), (0, 5, 0), (1, 3, 1)]
      and reception((0, 4, 0), in_W1, 4) == 3,
      'no tower within distance 2; f = 3 from 0@(4,-1), 0@(5,0), 1@(3,1) at distance exactly 3')

print('')
print('--- W5: the note\'s T_5.  delta_{4,2}(H) <= 1/10 < 9/80 ---')
F5 = audit_witness('W5', 4, 2, in_W5, (5, 4, 1), Fraction(1, 10), 2, 4)
W5_TABLE = {0: [4, 2, 4, 2, 3], 1: [4, 3, 4, 3, 2], 2: [2, 4, 2, 4, 3], 3: [3, 4, 3, 4, 2]}
check('W5.printed-table',
      all(F5[(a, i, 0)] == W5_TABLE[a][i] for a in range(4) for i in range(5)),
      'all 20 entries of the T_5 table printed in the note reproduced by fresh BFS')
check('W5.beats-the-conjectured-value', Fraction(1, 10) < Fraction(9, 80),
      '1/10 = 0.100000 < 9/80 = 0.112500, an improvement of exactly 1/80, i.e. 11.1% of 9/80')
check('W5.below-the-radius-2-packing-floor', Fraction(1, 10) < Fraction(1, B2),
      '1/10 < 1/|B_2| = 1/9, so W5 cannot rely on within-distance-2 reception alone')
check('W5.outer-shell-mechanism-at-0@(1,0)',
      towers_within((0, 1, 0), in_W5, 2) == []
      and sorted(towers_within((0, 1, 0), in_W5, 3)) == [(0, 0, 0), (0, 1, -1)]
      and reception((0, 1, 0), in_W5, 4) == 2,
      'no tower within distance 2; f = 1+1 = 2 from 0@(1,-1) and 0@(0,0) at distance exactly 3')

print('')
print('--- W2, W3: patterns defined in this file at densities 1/7 and 1/6 ---')
F2 = audit_witness('W2', 4, 3, in_W2, (7, 1, 1), Fraction(1, 7), 3, 7)
W2_TABLE = {0: [6, 6, 4, 4, 4, 4, 4], 1: [7, 4, 4, 3, 4, 3, 6],
            2: [6, 3, 5, 3, 5, 3, 4], 3: [6, 6, 4, 4, 4, 4, 4]}
check('W2.matches-transcribed-table',
      all(F2[(a, c, 0)] == W2_TABLE[a][c] for a in range(4) for c in range(7)),
      'all 28 entries of the W2 table listed in this program reproduced by fresh BFS; '
      'min 3 attained at %d classes'
      % sum(1 for a in range(4) for c in range(7) if F2[(a, c, 0)] == 3))
F3 = audit_witness('W3', 4, 4, in_W3, (3, 0, 1), Fraction(1, 6), 4, 6)
W3_TABLE = {0: [5, 6, 5], 1: [6, 4, 5], 2: [5, 6, 5], 3: [5, 4, 6]}
check('W3.matches-transcribed-table',
      all(F3[(a, x, 0)] == W3_TABLE[a][x] for a in range(4) for x in range(3)),
      'all 12 entries of the W3 table listed in this program reproduced by fresh BFS')

print('')
print('--- W6, W7: the note\'s T_6 and T_7, at densities 1/8 and 5/36 ---')
F6 = audit_witness('W6', 4, 3, in_W6, (4, 1, 1), Fraction(1, 8), 3, 5)
W6_TABLE = {0: [4, 4, 4, 4], 1: [3, 4, 3, 5], 2: [3, 5, 3, 4], 3: [4, 4, 4, 4]}
check('W6.printed-table',
      all(F6[(a, i, 0)] == W6_TABLE[a][i] for a in range(4) for i in range(4)),
      'all 16 entries of the T_6 table printed in the note reproduced by fresh BFS')
check('W6.beats-the-conjectured-value', Fraction(1, 8) < Fraction(1, 7),
      '1/8 = 0.125000 < 1/7 = 0.142857')
check('W6.tight-class-1@(0,0)',
      towers_within((1, 0, 0), in_W6, 1) == [(0, 0, 0)] and reception((1, 0, 0), in_W6, 4) == 3,
      'f = 3 = r from the single tower 0@(0,0) at d = 1')
check('W6.tight-class-2@(0,0)',
      towers_within((2, 0, 0), in_W6, 1) == []
      and sorted(towers_within((2, 0, 0), in_W6, 3)) == [(0, 0, 0), (3, 1, -1)]
      and reception((2, 0, 0), in_W6, 4) == 3,
      'f = 3 = 2+1 from 0@(0,0) at d = 2 and 3@(1,-1) at d = 3')

F7 = audit_witness('W7', 4, 4, in_W7, (6, 0, 3), Fraction(5, 36), 4, 5)
W7_TABLE = {0: [[4, 4, 5], [5, 4, 4], [4, 4, 5], [4, 4, 4], [5, 4, 5], [4, 4, 4]],
            1: [[5, 5, 4], [4, 4, 4], [5, 4, 4], [4, 5, 5], [4, 4, 4], [4, 4, 5]],
            2: [[4, 4, 4], [5, 4, 5], [4, 4, 4], [5, 4, 4], [4, 4, 5], [5, 4, 4]],
            3: [[5, 4, 4], [4, 4, 4], [5, 5, 4], [4, 4, 5], [4, 4, 4], [4, 5, 5]]}
check('W7.printed-table',
      all(F7[(a, x, y)] == W7_TABLE[a][x][y]
          for a in range(4) for x in range(6) for y in range(3)),
      'all 72 entries of the T_7 table printed in the note reproduced by fresh BFS')
check('W7.beats-the-conjectured-value', Fraction(5, 36) < Fraction(1, 6),
      '5/36 = 0.138889 < 1/6 = 0.166667')

print('')
print('--- W4: delta_{4,2} <= 1/9, and it is NOT a perfect radius-2 code ---')
F4 = audit_witness('W4', 4, 2, in_W4, (3, 0, 3), Fraction(1, 9), 2, 5)
W4_TABLE = {0: [[4, 3, 3], [5, 4, 3], [3, 3, 4]],
            1: [[4, 4, 2], [4, 4, 4], [4, 2, 4]],
            2: [[3, 4, 2], [4, 4, 3], [5, 2, 3]],
            3: [[3, 2, 3], [4, 5, 3], [4, 3, 3]]}
check('W4.matches-transcribed-table',
      all(F4[(a, x, y)] == W4_TABLE[a][x][y]
          for a in range(4) for x in range(3) for y in range(3)),
      'all 36 entries of the W4 table listed in this program reproduced by fresh BFS')
check('W4.beats-the-conjectured-value', Fraction(1, 9) < Fraction(9, 80),
      '1/9 = 0.111111 < 9/80 = 0.112500')

# The four radius-2 balls of the four tower classes, reduced mod 3.  Arithmetic says 4 x |B_2| =
# 4 x 9 = 36 = the number of classes; the check below shows this is NOT a partition.
def cls4(v):
    return (v[0], v[1] % 3, v[2] % 3)


W4_TOWER_REPS = [(0, 0, 0), (0, 1, 1), (1, 1, 0), (1, 2, 2)]
entries = []
for tv in W4_TOWER_REPS:
    b = ball(tv, 2)
    check('W4.ball-of-%d@(%d,%d)-has-9-vertices' % tv, len(b) == 9 and len(set(map(cls4, b))) == 9,
          '|B_2| = 9 and its 9 vertices lie in 9 distinct classes mod 3')
    entries.extend(cls4(w) for w in b)
check('W4.entry-count', len(entries) == 36, '4 balls x 9 = 36 entries, i.e. exactly the class count')
distinct = set(entries)
check('W4.only-32-distinct-classes-covered', len(distinct) == 32,
      '36 entries cover only 32 distinct classes, so the balls do NOT partition the domain')

collisions = sorted(c for c in distinct if entries.count(c) > 1)
COLL_PRINTED = sorted([(3, 1, 0), (1, 0, 1), (2, 1, 1), (2, 2, 0)])
check('W4.four-collisions', collisions == COLL_PRINTED and all(entries.count(c) == 2 for c in collisions),
      'exactly the four doubly-covered classes listed in this program: '
      + ' '.join('%d@(%d,%d)' % c for c in collisions))

uncovered = sorted(cls4(v) for v in reps(3, 0, 3) if not towers_within(v, in_W4, 2))
UNCOV_PRINTED = {(0, 0, 1): 3, (1, 0, 2): 2, (1, 2, 1): 2, (2, 2, 1): 2}
check('W4.four-uncovered-classes', uncovered == sorted(UNCOV_PRINTED),
      'exactly four classes have NO tower within distance 2: '
      + ' '.join('%d@(%d,%d)' % c for c in uncovered))
check('W4.uncovered-classes-are-served-from-shell-3',
      all(reception(u, in_W4, 4) == want and towers_within(u, in_W4, 2) == []
          for u, want in UNCOV_PRINTED.items()),
      'their f values are ' + ', '.join('%d@(%d,%d):%d' % (c + (v,))
                                        for c, v in sorted(UNCOV_PRINTED.items())))
check('W4.is-not-a-distance-2-dominating-set', bool(uncovered),
      'four classes lie at distance >= 3 from every tower, so W4 does NOT settle delta_{3,1}')

print('')
print('--- the baseline: the (3,1) broadcast T31 written out in this file, of density 1/8, '
      'checked here to be a (4,2) broadcast as well (this one witness only; the general '
      'implication is not verified) ---')
# The congruence below is written out in this file.  It is checked as a (3,1) broadcast of density
# 1/8 and then as a (4,2) broadcast, so that nothing here claims a first upper bound on
# delta_{4,2}.
audit_witness('T31.as-(3,1)', 3, 1, in_T31, (4, 1, 1), Fraction(1, 8), 1, 3)
audit_witness('T31.as-(4,2)', 4, 2, in_T31, (4, 1, 1), Fraction(1, 8), 3, 5)
check('baseline.ordering',
      Fraction(9, 80) < Fraction(1, 8),
      '9/80 = 0.112500 < 1/8 = 0.125000')
check('baseline.novelty-interval',
      Fraction(1, 10) < Fraction(9, 80) < Fraction(1, 8)
      and Fraction(1, 8) - Fraction(1, 10) == Fraction(1, 40),
      '1/10 < 9/80 < 1/8 and 1/8 - 1/10 = 1/40, i.e. 20% of 1/8')

print('')
print('--- the lower bound, from the published shell sizes alone ---')
# Capping lemma: g(u) = sum min(4-d, r) over towers within distance 3 still satisfies g >= r, so
# double counting over a large window gives r <= density * C(r) with
#     C(r) = 1*min(4,r) + 3*min(3,r) + 5*min(2,r) + 8*min(1,r).
def C(r):
    return sum(n * min(4 - d, r) for d, n in enumerate(SH[0][:4]))


check('lower.C(2)', C(2) == 26, 'C(2) = 1*2 + 3*2 + 5*2 + 8*1 = 26')
check('lower.C(3)', C(3) == 30, 'C(3) = 1*3 + 3*3 + 5*2 + 8*1 = 30')
check('lower.C(4)', C(4) == 31, 'C(4) = 1*4 + 3*3 + 5*2 + 8*1 = 31 (the cap does not bite at r = 4)')
check('lower.C(2)-arithmetic', Fraction(2, C(2)) == Fraction(1, 13),
      'NOTE (capping/double-counting lemma, comment above, not verified here): given '
      'r <= density * C(r), C(2) = 26 yields the bound 2/26 = 1/13 = 0.076923; '
      'only the identity 2/26 = 1/13 is checked')
check('lower.C(3)-arithmetic', Fraction(3, C(3)) == Fraction(1, 10),
      'NOTE (capping/double-counting lemma, not verified here): given r <= density * C(r), '
      'C(3) = 30 yields 3/30 = 1/10 = 0.100000; only the identity 3/30 = 1/10 is checked')
check('lower.C(4)-arithmetic', Fraction(4, C(4)) == Fraction(4, 31),
      'NOTE (capping/double-counting lemma, not verified here): given r <= density * C(r), '
      'C(4) = 31 yields 4/31 = 0.129032; only the value 4/C(4) = 4/31 is checked')
check('lower.beats-the-naive-weight-count',
      Fraction(2, 31) < Fraction(1, 13) and C(4) == 31,
      'the uncapped per-tower emission is 4+9+10+8 = 31, giving only 2/31 = 0.064516; '
      'the capped count 1/13 is strictly stronger')

print('')
print('--- the constants of the note, and their consistency ---')
check('record.sandwich-4-2', Fraction(1, 13) <= Fraction(1, 10) and Fraction(1, 10) < Fraction(9, 80),
      'constants ordered: 1/13 <= 1/10 < 9/80; the upper bound delta_{4,2} <= 1/10 is '
      'witnessed by W5, while the lower end 1/13 rests on the capping lemma NOT verified here')
check('record.sandwich-4-3', Fraction(1, 10) < Fraction(1, 8) < Fraction(1, 7),
      'constants ordered: 1/10 < 1/8 < 1/7; the upper bound delta_{4,3} <= 1/8 is '
      'witnessed by W6, while the lower end 1/10 rests on the capping lemma NOT verified here')
check('record.sandwich-4-4', Fraction(4, 31) < Fraction(5, 36) < Fraction(1, 6),
      'constants ordered: 4/31 < 5/36 < 1/6; the upper bound delta_{4,4} <= 5/36 is '
      'witnessed by W7, while the lower end 4/31 rests on the capping lemma NOT verified here')
check('record.four-constants-are-ordered',
      Fraction(1, 12) < Fraction(1, 10) < Fraction(1, 8) < Fraction(5, 36),
      'the four constants 1/12, 1/10, 1/8, 5/36 increase in this order; 1/10, 1/8, 5/36 are the '
      'upper bounds witnessed above by W5, W6, W7, and 1/12 is a literal written here and not '
      'checked as a bound. No claim that these are bounds of record and no monotonicity of '
      'delta_{4,r} is verified')
check('record.two-(4,4)-witnesses-are-also-(4,3)-broadcasts',
      all(reception(v, in_W7, 4) >= 3 for v in reps(6, 0, 3))
      and all(reception(v, in_W3, 4) >= 3 for v in reps(3, 0, 1)),
      'the two (4,4) witnesses W3 and W7 are checked to be (4,3) broadcasts as well, which is the '
      'monotonicity the ordering above relies on')

print('')
print('NOT RE-RUN: this program checks the exhibited patterns and the arithmetic printed above, and')
print('NOT RE-RUN: nothing else. It does not search over period lattices and so establishes no lower')
print('NOT RE-RUN: bound and no optimality: no claim that any density above is least is checked.')
print('NOT RE-RUN: The capping/double-counting step behind r/C(r) is not verified; only the integers')
print('NOT RE-RUN: C(2), C(3), C(4) and the fractions r/C(r) are computed.')
print('NOT RE-RUN: It reads no external file and no figure, and it makes no claim about the figures')
print('NOT RE-RUN: of any cited paper.')
print('')

if _F:
    print('VERDICT: %d of %d CHECKS FAILED' % (_F, _N))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _N)
sys.exit(0)
