#!/usr/bin/env python3
"""Verification program for

    "A 24-vertex tree refuting the Marchessault--Mynhardt uniquely-radial equality"

Python 3.9+, standard library only: no third-party package and no external data file.

The tree T and the broadcast f are printed in the paper and are read from it.  The other
objects this program uses -- the second example T5, the further broadcasts g and h, and the
small control trees -- are transcribed HERE and are not printed in the paper, which says so in
its Scope section.  Every quantity the paper asserts about T and f is re-derived here from the
definitions of

    S. Marchessault and C. M. Mynhardt, "Lower boundary independent broadcasts in
    trees", Discuss. Math. Graph Theory 44 (2024) 75--99 (= arXiv:2105.04611v2),

namely: broadcast, sigma(f), hearing, dominating, "overlaps only in boundaries"
(bn-independence), N_f(u), B_f(u), PB_f(u), covered/uncovered edges, U_f^E, split-set,
radial subtree, uniquely radial, and the two quoted maximality criteria (the
boundary/private-boundary criterion and the uncovered-edge-component criterion).
All arithmetic is exact integer arithmetic;
there is no floating point anywhere in the file.

One `PASS <name>` line per check, then

    VERDICT: ALL <n> CHECKS PASS

and exit status 0 iff every check passed.

WHAT IS DELIBERATELY NOT DONE HERE: the program performs no exhaustive search over the
broadcasts of T, so it establishes no exact value of i_bn(T); the note claims no exact
value of i_bn(T) either.  What is established below is the note's theorem: i_bn(T) in
{8,9} and hence i_bn(T) != gamma_b(T) + ceil((|M|+1)/3) = 10.
"""

import itertools
import sys
from collections import deque

# ---------------------------------------------------------------------------
# check bookkeeping
# ---------------------------------------------------------------------------
_N = 0
_BAD = []


def CHECK(cond, name, detail=""):
    global _N
    if cond:
        _N += 1
        print("PASS %s%s" % (name, ("  --  " + detail) if detail else ""))
    else:
        _BAD.append(name)
        print("FAIL %s%s" % (name, ("  --  " + detail) if detail else ""))


# ---------------------------------------------------------------------------
# 1. THE OBJECTS, EXACTLY AS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------
# T: the chain of four copies of the spider S(2,2,1).  Spine v0..v19; pendants u1..u4
# attached at the block centres, which are the spine vertices v2, v7, v12, v17.
EDGES_T = [
    ("v0", "v1"), ("v1", "v2"), ("v2", "v3"), ("v3", "v4"), ("v4", "v5"),
    ("v5", "v6"), ("v6", "v7"), ("v7", "v8"), ("v8", "v9"), ("v9", "v10"),
    ("v10", "v11"), ("v11", "v12"), ("v12", "v13"), ("v13", "v14"), ("v14", "v15"),
    ("v15", "v16"), ("v16", "v17"), ("v17", "v18"), ("v18", "v19"),
    ("v2", "u1"), ("v7", "u2"), ("v12", "u3"), ("v17", "u4"),
]

# T5: the chain of five copies of S(2,2,1).  Spine w0..w24; pendants x1..x5 at the
# block centres w2, w7, w12, w17, w22.
EDGES_T5 = ([("w%d" % i, "w%d" % (i + 1)) for i in range(24)]
            + [("w%d" % (5 * i + 2), "x%d" % (i + 1)) for i in range(5)])

# The unique maximum split-set of T, as printed in the paper.
M_T = [("v4", "v5"), ("v9", "v10"), ("v14", "v15")]
# The unique maximum split-set of T5, as printed in the paper.
M_T5 = [("w4", "w5"), ("w9", "w10"), ("w14", "w15"), ("w19", "w20")]

# Witness A: f = 1 on these nine vertices, 0 elsewhere.
WITNESS_A = ["v0", "v2", "v4", "v7", "v9", "v12", "v14", "v17", "v19"]
# Witness B, the structurally different cost-9 broadcast.
WITNESS_B = {"v4": 4, "v9": 1, "v13": 2, "v17": 2}
# The gamma_b-broadcast of T exhibited in the paper: strength 2 at each block centre.
GB_WITNESS_T = {"v2": 2, "v7": 2, "v12": 2, "v17": 2}
# The gamma_b-broadcast of T5 exhibited in the paper.
GB_WITNESS_T5 = {"w2": 2, "w7": 2, "w12": 2, "w17": 2, "w22": 2}
# Witness H on T5: h = 1 on these eleven spine vertices.
WITNESS_H = ["w0", "w2", "w4", "w7", "w9", "w12", "w14", "w17", "w19", "w22", "w24"]
# The nine pairwise-disjoint closed neighbourhoods of the paper's gamma(T) >= 9 count.
PACKING_CENTRES = ["u1", "u2", "u3", "u4", "v0", "v19", "v5", "v10", "v15"]


# ---------------------------------------------------------------------------
# 2. GRAPH MACHINERY
# ---------------------------------------------------------------------------
def _e(a, b):
    return (a, b) if a <= b else (b, a)


class Graph:
    """An undirected simple graph, given by its edge list.  All distances are exact."""

    def __init__(self, edges, verts=None):
        self.E = sorted({_e(a, b) for a, b in edges})
        vs = set(verts or [])
        for a, b in self.E:
            vs.add(a)
            vs.add(b)
        self.V = sorted(vs)
        self.adj = {v: [] for v in self.V}
        for a, b in self.E:
            self.adj[a].append(b)
            self.adj[b].append(a)
        self.D = {v: self.bfs(v) for v in self.V}
        self.ecc = {v: max(self.D[v].values()) for v in self.V}
        self.diam = max(self.ecc.values())
        self.rad = min(self.ecc.values())

    def bfs(self, s, cut=frozenset()):
        d = {s: 0}
        q = deque([s])
        while q:
            x = q.popleft()
            for y in self.adj[x]:
                if y in d or _e(x, y) in cut:
                    continue
                d[y] = d[x] + 1
                q.append(y)
        return d

    def connected(self):
        return len(self.bfs(self.V[0])) == len(self.V)

    def is_tree(self):
        return self.connected() and len(self.E) == len(self.V) - 1

    def components(self, cut):
        cut = frozenset(_e(a, b) for a, b in cut)
        seen, out = set(), []
        for v in self.V:
            if v in seen:
                continue
            c = set(self.bfs(v, cut))
            seen |= c
            out.append(sorted(c))
        return out

    def induced(self, verts):
        vs = set(verts)
        return Graph([e for e in self.E if e[0] in vs and e[1] in vs], verts=vs)

    def path_between(self, a, b):
        p, cur = [a], a
        while cur != b:
            for y in self.adj[cur]:
                if self.D[y][b] == self.D[cur][b] - 1:
                    cur = y
                    break
            p.append(cur)
        return p


# ---------------------------------------------------------------------------
# 3. THE PAPER'S BROADCAST DEFINITIONS, TRANSCRIBED
# ---------------------------------------------------------------------------
def vplus(g, f):
    return [v for v in g.V if f.get(v, 0) > 0]


def is_broadcast(g, f):
    """f : V -> {0,...,diam} with f(v) <= e(v) for every v."""
    return all(0 <= f.get(v, 0) <= g.ecc[v] for v in g.V)


def sigma(f):
    return sum(f.values())


def heard_by(g, f, x):
    """The broadcasters x hears: those u in V_f^+ with d(x,u) <= f(u)."""
    return [u for u in vplus(g, f) if g.D[x][u] <= f[u]]


def dominating(g, f):
    return all(heard_by(g, f, x) for x in g.V)


def multiply_heard(g, f):
    return [x for x in g.V if len(heard_by(g, f, x)) > 1]


def bn_independent(g, f):
    """"Overlaps only in boundaries": every vertex x that hears more than one
    broadcasting vertex also satisfies d(x,u) >= f(u) for all u in V_f^+."""
    vp = vplus(g, f)
    for x in multiply_heard(g, f):
        for u in vp:
            if g.D[x][u] < f[u]:
                return False
    return True


def Nf(g, f, v):
    return {u for u in g.V if g.D[u][v] <= f[v]}


def Bf(g, f, v):
    return {u for u in g.V if g.D[u][v] == f[v]}


def PBf(g, f, v):
    """The private boundary of v: the vertices that stop hearing anything when v's
    strength is lowered by one."""
    f2 = dict(f)
    f2[v] = f[v] - 1
    return {u for u in Nf(g, f, v) if not heard_by(g, f2, u)}


def uncovered_edges(g, f):
    """U_f^E.  Edge xy is covered by u when x,y in N_f(u) and at least one of x,y is
    not in B_f(u)."""
    vp = vplus(g, f)
    out = []
    for (x, y) in g.E:
        cov = False
        for u in vp:
            dx, dy = g.D[x][u], g.D[y][u]
            if dx <= f[u] and dy <= f[u] and (dx != f[u] or dy != f[u]):
                cov = True
                break
        if not cov:
            out.append((x, y))
    return out


def maximal_thm_minimal(g, f):
    """The quoted boundary criterion: f dominating and B_f(v) - PB_f(v) nonempty for every
    v in V_f^+.  Applied only when |V_f^+| >= 2, where the overlap it asks for can
    exist at all."""
    if not (bn_independent(g, f) and dominating(g, f)):
        return False
    vp = vplus(g, f)
    if len(vp) < 2:
        return no_single_increment(g, f)
    return all(Bf(g, f, v) - PBf(g, f, v) for v in vp)


def maximal_prop_2bv(g, f):
    """The quoted uncovered-edge criterion: |V_f^+| >= 2 and every component of
    T - U_f^E contains at
    least two broadcasting vertices."""
    if not bn_independent(g, f):
        return False
    vp = set(vplus(g, f))
    if len(vp) < 2:
        return False
    U = uncovered_edges(g, f)
    return all(len(vp & set(c)) >= 2 for c in g.components(U))


def no_single_increment(g, f):
    """The direct test: no unit increment of f, at any vertex, stays bn-independent."""
    for v in g.V:
        if f.get(v, 0) >= g.ecc[v]:
            continue
        f2 = dict(f)
        f2[v] = f.get(v, 0) + 1
        if is_broadcast(g, f2) and bn_independent(g, f2):
            return False
    return True


# ---------------------------------------------------------------------------
# 4. EXHAUSTIVE ENUMERATION (used only on the two 12-vertex controls)
# ---------------------------------------------------------------------------
def enum_broadcasts(g, budget):
    """Every broadcast f with sigma(f) <= budget, as a dict."""
    V = g.V
    n = len(V)
    f = {v: 0 for v in V}

    def rec(i, rem):
        if i == n:
            yield dict(f)
            return
        v = V[i]
        for k in range(0, min(g.ecc[v], rem) + 1):
            f[v] = k
            yield from rec(i + 1, rem - k)
        f[v] = 0

    yield from rec(0, budget)


def gamma_b_exact(g, cap):
    """(gamma_b, all gamma_b-broadcasts) by scanning cost upwards; None if > cap."""
    for c in range(1, cap + 1):
        wits = [f for f in enum_broadcasts(g, c)
                if sigma(f) == c and dominating(g, f)]
        if wits:
            return c, wits
    return None, []


def min_maximal_bn(g, budget):
    """(cost, witness) of a cheapest maximal bn-independent broadcast of cost <= budget;
    (None, None) if there is none.  Both maximality criteria are demanded to agree."""
    best, wit = None, None
    for f in enum_broadcasts(g, budget):
        c = sigma(f)
        if c == 0 or (best is not None and c >= best):
            continue
        if not (dominating(g, f) and bn_independent(g, f)):
            continue
        m1 = maximal_thm_minimal(g, f)
        m2 = (maximal_prop_2bv(g, f) if len(vplus(g, f)) >= 2
              else no_single_increment(g, f))
        if m1 != m2:
            return ("CRITERIA-DISAGREE", f)
        if m1:
            best, wit = c, dict(f)
    return best, wit


# ---------------------------------------------------------------------------
# 5. SPLIT-SETS, FROM THE DEFINITION RESTATED IN THE NOTE
# ---------------------------------------------------------------------------
def is_split_set(g, p, cuts):
    """Definition of a split-set, applied to the cut positions `cuts` (indices into
    the edge list of the path p): every component-subpath has even positive length and
    is diametrical in its own component of T - M."""
    L = len(p) - 1
    M = frozenset(_e(p[j], p[j + 1]) for j in cuts)
    if len(M) != len(cuts):
        return False, M
    prev = -1
    for j in list(cuts) + [L]:
        lo, hi = prev + 1, j
        ln = hi - lo
        if ln <= 0 or ln % 2 != 0:
            return False, M
        comp = set(g.bfs(p[lo], M))
        if set(p[lo:hi + 1]) - comp or (comp & set(p)) != set(p[lo:hi + 1]):
            return False, M
        sub = g.induced(comp)
        if sub.diam != ln:
            return False, M
        prev = j
    return True, M


def cut_candidates(L):
    """Every set of cut positions whose parts all have even positive length.  Complete:
    the condition is a condition on each part separately, so a depth-first search that
    closes one part at a time and prunes as soon as a closed part is odd or empty
    reaches every admissible set exactly once."""
    out = []

    def rec(prev, cuts):
        for j in range(prev + 1, L):
            ln = j if prev < 0 else j - prev - 1
            if ln <= 0 or ln % 2 != 0:
                continue
            tail = L - j - 1
            if tail > 0 and tail % 2 == 0:
                out.append(tuple(cuts + [j]))
            rec(j, cuts + [j])

    rec(-1, [])
    return out


def split_sets(g):
    """Every split-set of g, taken over every diametrical path."""
    out = {}
    for a, b in itertools.combinations(g.V, 2):
        if g.D[a][b] != g.diam:
            continue
        p = g.path_between(a, b)
        for cuts in cut_candidates(len(p) - 1):
            ok, M = is_split_set(g, p, cuts)
            if ok:
                out[M] = True
    return sorted(out, key=lambda s: (len(s), sorted(s)))


def split_sets_bruteforce(g):
    """The same thing by exhaustive enumeration of ALL subsets of every diametrical
    path's edge set.  Exponential, so it is used only to audit `split_sets` on small
    controls."""
    out = {}
    for a, b in itertools.combinations(g.V, 2):
        if g.D[a][b] != g.diam:
            continue
        p = g.path_between(a, b)
        L = len(p) - 1
        for r in range(1, L + 1):
            for cuts in itertools.combinations(range(L), r):
                ok, M = is_split_set(g, p, cuts)
                if ok:
                    out[M] = True
    return sorted(out, key=lambda s: (len(s), sorted(s)))


# ---------------------------------------------------------------------------
# 6. SMALL BUILDERS FOR THE CONTROLS
# ---------------------------------------------------------------------------
def path_graph(n, tag="p"):
    return Graph([("%s%02d" % (tag, i), "%s%02d" % (tag, i + 1)) for i in range(n - 1)])


def spider_chain(k, tag="s"):
    """The chain of k copies of S(2,2,1), built the same way the paper builds T."""
    sp = ["%sv%02d" % (tag, i) for i in range(5 * k)]
    edges = [(sp[i], sp[i + 1]) for i in range(5 * k - 1)]
    edges += [(sp[5 * i + 2], "%su%d" % (tag, i + 1)) for i in range(k)]
    return Graph(edges), sp


def ceil_div(a, b):
    return -((-a) // b)


def uniquely_radial(sub):
    """(is uniquely radial, gamma_b, number of gamma_b-broadcasts).  A radial subtree is
    uniquely radial when gamma_b equals its radius AND it has exactly one
    gamma_b-broadcast (the paper fixes the phrase by that use)."""
    gb, wits = gamma_b_exact(sub, sub.rad)
    return (gb == sub.rad and len(wits) == 1), gb, len(wits)


def vkey(v):
    """Sort key that orders v9 before v10 -- for display only; set membership always
    goes through the canonical _e()."""
    i = 0
    while i < len(v) and not v[i].isdigit():
        i += 1
    return (v[:i], int(v[i:]) if i < len(v) else -1)


def fmt_edges(es):
    return " ".join("%s%s" % tuple(sorted(e, key=vkey))
                    for e in sorted(es, key=lambda e: sorted(map(vkey, e))))


def fmt_verts(vs):
    return " ".join(sorted(vs, key=vkey))


# ---------------------------------------------------------------------------
# 7. THE CHECKS
# ---------------------------------------------------------------------------
def main():
    T = Graph(EDGES_T)
    spine = ["v%d" % i for i in range(20)]
    pend = ["u1", "u2", "u3", "u4"]

    # -- 7.1 the object -----------------------------------------------------
    CHECK(T.is_tree() and len(T.V) == 24 and len(T.E) == 23,
          "T is a tree on 24 vertices with 23 edges",
          "n=%d m=%d connected=%s" % (len(T.V), len(T.E), T.connected()))
    CHECK(sorted(T.V) == sorted(spine + pend),
          "the vertex set of T is exactly the printed spine plus the four pendants",
          "|V|=%d" % len(T.V))
    deg = {v: len(T.adj[v]) for v in T.V}
    CHECK(sorted(v for v in T.V if deg[v] == 1) == sorted(["v0", "v19"] + pend)
          and sorted(v for v in T.V if deg[v] == 3) == ["v12", "v17", "v2", "v7"]
          and all(deg[v] == 2 for v in T.V
                  if deg[v] not in (1, 3)),
          "the degrees of T are as printed: six leaves, four vertices of degree 3, "
          "the rest of degree 2",
          "leaves=%s | deg-3=%s" % (fmt_verts(v for v in T.V if deg[v] == 1),
                                    fmt_verts(v for v in T.V if deg[v] == 3)))
    CHECK(T.diam == 19 and T.D["v0"]["v19"] == 19,
          "diam(T) = 19, realised by the pair (v0,v19)", "diam=%d" % T.diam)
    CHECK(T.rad == 10 and sorted(v for v in T.V if T.ecc[v] == 10) == ["v10", "v9"],
          "rad(T) = 10, attained exactly at v9 and v10",
          "rad=%d centre=%s" % (T.rad, fmt_verts(v for v in T.V if T.ecc[v] == T.rad)))
    diapairs = [(a, b) for a, b in itertools.combinations(T.V, 2)
                if T.D[a][b] == T.diam]
    CHECK(diapairs == [("v0", "v19")],
          "the diametrical path of T is unique, so every split-set lies in the spine",
          "diametrical pairs = %s" % (diapairs,))

    # -- 7.2 the split-set enumerator is audited against brute force --------
    for name, g in (("P12", path_graph(12)),
                    ("the two-block chain", spider_chain(2, "c")[0])):
        a = split_sets(g)
        b = split_sets_bruteforce(g)
        CHECK(a == b and len(a) > 0,
              "the pruned split-set enumerator agrees with exhaustive enumeration of "
              "all edge subsets on %s" % name,
              "pruned=%d subsets brute=%d subsets, sizes=%s"
              % (len(a), len(b), sorted({len(s) for s in a})))

    # -- 7.3 m = 3 and the maximum split-set is unique ----------------------
    ss = split_sets(T)
    m = max(len(s) for s in ss)
    maxss = [s for s in ss if len(s) == m]
    CHECK(m == 3 and len(maxss) == 1,
          "T has maximum split-set cardinality |M| = 3, and exactly one maximum "
          "split-set",
          "m=%d #maximum split-sets=%d all cardinalities present=%s"
          % (m, len(maxss), sorted({len(s) for s in ss})))
    CHECK(maxss[0] == frozenset(_e(a, b) for a, b in M_T),
          "that unique maximum split-set is the one printed in the paper",
          "M = %s" % fmt_edges(sorted(maxss[0])))
    M = sorted(maxss[0])

    # -- 7.4 the radial subtrees --------------------------------------------
    comps = T.components(M)
    CHECK(len(comps) == 4 and all(len(c) == 6 for c in comps),
          "T - M has four components, each on six vertices",
          "sizes=%s" % [len(c) for c in comps])
    ok_ur = True
    detail = []
    for c in comps:
        sub = T.induced(c)
        ur, gb, ngb = uniquely_radial(sub)
        degs = sorted(len(sub.adj[v]) for v in sub.V)
        ok_ur &= (ur and sub.rad == 2 and sub.diam == 4 and gb == 2
                  and degs == [1, 1, 1, 2, 2, 3])
        detail.append("rad=%d diam=%d gamma_b=%d #gamma_b-broadcasts=%d" %
                      (sub.rad, sub.diam, gb, ngb))
    CHECK(ok_ur,
          "every radial subtree of T is a copy of S(2,2,1): radius exactly 2, and "
          "uniquely radial (gamma_b = rad = 2 attained by exactly one broadcast)",
          " | ".join(detail))
    # the anti-control: radius 2 alone is not enough
    p5 = path_graph(5)
    ur5, gb5, n5 = uniquely_radial(p5)
    CHECK((not ur5) and p5.rad == 2 and gb5 == 2 and n5 == 4,
          "anti-control: P5 also has radius 2 and gamma_b = 2, but it has four "
          "gamma_b-broadcasts and so is NOT uniquely radial -- the hypothesis really "
          "does need the six-vertex spider",
          "rad=%d gamma_b=%d #gamma_b-broadcasts=%d" % (p5.rad, gb5, n5))

    # -- 7.5 gamma_b(T) = 8 and the value the problem predicts --------------
    CHECK(is_broadcast(T, GB_WITNESS_T) and dominating(T, GB_WITNESS_T)
          and sigma(GB_WITNESS_T) == 8,
          "the printed broadcast of strength 2 at each block centre is a dominating "
          "broadcast of cost 8, so gamma_b(T) <= 8",
          "sigma=%d" % sigma(GB_WITNESS_T))
    CHECK((T.diam - m) % 2 == 0 and (T.diam - m) // 2 == 8,
          "the cited Herke--Mynhardt formula gives gamma_b(T) = (diam T - |M|)/2 = 8, "
          "agreeing with that witness",
          "(19 - 3)/2 = %d" % ((T.diam - m) // 2))
    gb_T = (T.diam - m) // 2
    CHECK(gb_T < T.rad and len(ss) > 0,
          "T is non-radial (it has a split-set, and gamma_b = 8 < 10 = rad), so the "
          "paper's bound is in force",
          "gamma_b=%d rad=%d #split-sets=%d" % (gb_T, T.rad, len(ss)))
    predicted = gb_T + ceil_div(m + 1, 3)
    CHECK(predicted == 10,
          "the equality of the problem therefore demands "
          "i_bn(T) = gamma_b + ceil((|M|+1)/3) = 10",
          "8 + ceil(4/3) = %d" % predicted)

    # -- 7.6 witness A ------------------------------------------------------
    f = {v: (1 if v in WITNESS_A else 0) for v in T.V}
    CHECK(is_broadcast(T, f) and sigma(f) == 9 and len(WITNESS_A) == 9,
          "witness A is a legal broadcast of cost 9 (f = 1 on nine vertices)",
          "sigma=%d" % sigma(f))
    CHECK(all(_e(a, b) not in set(T.E)
              for a, b in itertools.combinations(WITNESS_A, 2)),
          "the support of witness A is an independent set of T",
          "no two of the nine are adjacent")
    CHECK(dominating(T, f),
          "witness A is dominating: every one of the 24 vertices hears a broadcaster",
          "unheard vertices = %s" % ([x for x in T.V if not heard_by(T, f, x)] or "none"))
    mh = sorted(multiply_heard(T, f))
    CHECK(bn_independent(T, f) and mh == ["v1", "v13", "v18", "v3", "v8"],
          "witness A is bn-independent, and the vertices hearing more than one "
          "broadcaster are exactly v1, v3, v8, v13, v18",
          "multiply heard = %s" % fmt_verts(mh))
    U = uncovered_edges(T, f)
    CHECK(U == [_e("v10", "v11"), _e("v15", "v16"), _e("v5", "v6")],
          "U_f^E for witness A is exactly the three edges v5v6, v10v11, v15v16",
          "U_f^E = %s" % fmt_edges(U))
    cU = T.components(U)
    counts = [len(set(c) & set(WITNESS_A)) for c in cU]
    CHECK(maximal_prop_2bv(T, f) and len(cU) == 4 and min(counts) >= 2,
          "witness A is maximal bn-independent by the quoted uncovered-edge "
          "criterion: T - U_f^E has four components and each holds at least two "
          "broadcasters",
          "broadcasters per component = %s" % counts)
    CHECK(maximal_thm_minimal(T, f),
          "witness A is maximal bn-independent by the quoted boundary criterion: "
          "B_f(v) - PB_f(v) is non-empty at all nine broadcasters",
          "|B_f - PB_f| = %s" % [len(Bf(T, f, v) - PBf(T, f, v)) for v in WITNESS_A])
    CHECK(no_single_increment(T, f),
          "witness A is maximal by the direct test as well: no unit increment of f, at "
          "any vertex of T, stays bn-independent",
          "24 vertices tested")

    # -- 7.7 witness B ------------------------------------------------------
    g = {v: WITNESS_B.get(v, 0) for v in T.V}
    CHECK(is_broadcast(T, g) and sigma(g) == 9,
          "witness B is a legal broadcast of cost 9 with a different shape "
          "(g(v4)=4, g(v9)=1, g(v13)=2, g(v17)=2)",
          "sigma=%d, e(v4)=%d >= 4" % (sigma(g), T.ecc["v4"]))
    CHECK(dominating(T, g) and bn_independent(T, g),
          "witness B is dominating and bn-independent",
          "multiply heard = %s" % fmt_verts(multiply_heard(T, g)))
    Ug = uncovered_edges(T, g)
    CHECK(Ug == [_e("v10", "v11")],
          "U_g^E for witness B is the single edge v10v11",
          "U_g^E = %s" % fmt_edges(Ug))
    CHECK(maximal_prop_2bv(T, g) and maximal_thm_minimal(T, g),
          "witness B is maximal bn-independent by both of the paper's criteria",
          "components of T - U_g^E carry %s broadcasters"
          % [len(set(c) & set(WITNESS_B)) for c in T.components(Ug)])
    inside = [x for x in T.V
              if any(T.D[x][u] < g[u] for u in vplus(T, g))]
    mhg = multiply_heard(T, g)
    CHECK(sorted(mhg) == ["v15", "v8"]
          and all(T.D[x][u] == g[u] for x in mhg for u in heard_by(T, g, x))
          and not (set(mhg) & set(inside)),
          "under witness B exactly two vertices hear more than one broadcaster -- v8 "
          "and v15 -- and each lies on the boundary of both broadcasters it hears; the "
          "many vertices strictly inside some ball all hear exactly one broadcaster, "
          "which \"overlaps only in boundaries\" permits (reading the condition as a "
          "ban on strict interiors would give a false negative here)",
          "multiply heard = %s ; strictly inside some ball = %d vertices, none of them "
          "multiply heard" % (fmt_verts(mhg), len(inside)))

    # -- 7.8 the conclusion for T -------------------------------------------
    CHECK(sigma(f) == 9 and predicted == 10 and gb_T == 8,
          "i_bn(T) <= 9 from witness A, and i_bn(T) >= gamma_b(T) = 8 from the cited "
          "inequality, so i_bn(T) is 8 or 9",
          "i_bn(T) in {8,9}, predicted %d" % predicted)
    CHECK(8 != predicted and 9 != predicted,
          "both admissible values differ from the predicted 10, so the equality of the "
          "problem fails on T",
          "8 != 10 and 9 != 10")
    CHECK(sigma(f) <= predicted,
          "the paper's inequality i_bn(T) <= gamma_b(T) + ceil((|M|+1)/3) is NOT "
          "contradicted: 9 <= 10",
          "9 <= 10")

    # -- 7.9 the mechanism: gamma(T) = i(T) = 9 -----------------------------
    nbhds = [set(T.adj[v]) | {v} for v in PACKING_CENTRES]
    disjoint = all(not (nbhds[i] & nbhds[j])
                   for i, j in itertools.combinations(range(len(nbhds)), 2))
    CHECK(disjoint and len(nbhds) == 9,
          "the nine closed neighbourhoods N[u1..u4], N[v0], N[v19], N[v5], N[v10], "
          "N[v15] are pairwise disjoint, so gamma(T) >= 9",
          "centres = %s" % fmt_verts(PACKING_CENTRES))
    dom_ok = all(set(T.adj[x]) & set(WITNESS_A) or x in WITNESS_A for x in T.V)
    CHECK(dom_ok,
          "the support of witness A is an independent dominating set of T of size 9, so "
          "gamma(T) = i(T) = 9 = gamma_b(T) + 1",
          "|S| = 9")

    # -- 7.10 both-polarity controls, exhaustive on the cheap side ----------
    # CTRL1: the two-block chain must return EQUALITY (the quoted |M| in {1,2} theorem).
    c1, sp1 = spider_chain(2, "c")
    ss1 = split_sets(c1)
    m1 = max(len(s) for s in ss1)
    ur1 = all(uniquely_radial(c1.induced(c))[0] and c1.induced(c).rad == 2
              for c in c1.components(sorted([s for s in ss1 if len(s) == m1][0])))
    gb1, _w1 = gamma_b_exact(c1, 6)
    CHECK(len(c1.V) == 12 and c1.diam == 9 and m1 == 1
          and len([s for s in ss1 if len(s) == m1]) == 1 and ur1 and gb1 == 4,
          "control 1 (chain of two S(2,2,1)) satisfies the hypothesis with |M| = 1: "
          "n=12, diam=9, one maximum split-set, both radial subtrees uniquely radial "
          "of radius 2, gamma_b = 4",
          "m=%d gamma_b=%d" % (m1, gb1))
    fw1 = {v: 0 for v in c1.V}
    fw1[sp1[5]] = 5
    b1, _ = min_maximal_bn(c1, 4)
    CHECK(is_broadcast(c1, fw1) and maximal_thm_minimal(c1, fw1) and sigma(fw1) == 5
          and b1 is None,
          "control 1 returns EQUALITY: a maximal bn-independent broadcast of cost 5 "
          "exists and an exhaustive search finds none of cost <= 4, so i_bn = 5 = "
          "gamma_b + 1 exactly as the quoted |M| in {1,2} theorem requires",
          "cost 5 witness maximal=True; cheapest of cost <=4 = %s" % (b1,))
    # CTRL2: P12 must return a STRICT inequality (Neilson's i_bn(P_n) = ceil(2n/5)).
    p12 = path_graph(12)
    ss2 = split_sets(p12)
    m2 = max(len(s) for s in ss2)
    gb2, _w2 = gamma_b_exact(p12, 6)
    fw2 = {v: 0 for v in p12.V}
    fw2.update({"p02": 2, "p05": 1, "p08": 1, "p10": 1})
    b2, _ = min_maximal_bn(p12, 4)
    CHECK(gb2 == 4 and m2 == 3 and is_broadcast(p12, fw2)
          and maximal_thm_minimal(p12, fw2) and sigma(fw2) == 5 and b2 is None,
          "control 2 returns a STRICT inequality: on P12, gamma_b = 4 and |M| = 3, a "
          "maximal bn-independent broadcast of cost 5 exists and none of cost <= 4 "
          "does, so i_bn = 5 < 6 = gamma_b + ceil(4/3), matching Neilson's ceil(2n/5)",
          "gamma_b=%d m=%d i_bn=5 predicted=%d ceil(2*12/5)=%d"
          % (gb2, m2, gb2 + ceil_div(m2 + 1, 3), ceil_div(24, 5)))
    CHECK(all(uniquely_radial(p12.induced(c))[1] == 1
              for c in p12.components(sorted([s for s in ss2 if len(s) == m2][0]))),
          "and control 2 is outside the hypothesis for the right reason: its four "
          "radial subtrees are copies of P3, of radius 1, not >= 2",
          "component gamma_b values are all 1")

    # -- 7.11 a second example, with |M| even -------------------------------
    T5 = Graph(EDGES_T5)
    CHECK(T5.is_tree() and len(T5.V) == 30 and len(T5.E) == 29 and T5.diam == 24,
          "the second exhibited tree T5 (chain of five S(2,2,1)) is a tree on 30 "
          "vertices with diameter 24",
          "n=%d m=%d diam=%d" % (len(T5.V), len(T5.E), T5.diam))
    ss5 = split_sets(T5)
    m5 = max(len(s) for s in ss5)
    max5 = [s for s in ss5 if len(s) == m5]
    CHECK(m5 == 4 and len(max5) == 1
          and max5[0] == frozenset(_e(a, b) for a, b in M_T5),
          "T5 has maximum split-set cardinality |M| = 4 -- an EVEN one -- and exactly "
          "one maximum split-set, the printed one",
          "m=%d #maximum=%d M=%s" % (m5, len(max5), fmt_edges(sorted(max5[0]))))
    M5 = sorted(max5[0])
    c5 = T5.components(M5)
    CHECK(len(c5) == 5 and all(uniquely_radial(T5.induced(c))[0]
                               and T5.induced(c).rad == 2 for c in c5),
          "all five radial subtrees of T5 are uniquely radial of radius exactly 2, so "
          "T5 satisfies the hypothesis too",
          "components=%d" % len(c5))
    CHECK(is_broadcast(T5, GB_WITNESS_T5) and dominating(T5, GB_WITNESS_T5)
          and sigma(GB_WITNESS_T5) == 10 and (T5.diam - m5) // 2 == 10,
          "gamma_b(T5) = (24 - 4)/2 = 10, matched by the printed dominating broadcast "
          "of cost 10 at the five block centres",
          "(24-4)/2 = %d, sigma = %d" % ((T5.diam - m5) // 2, sigma(GB_WITNESS_T5)))
    pred5 = 10 + ceil_div(m5 + 1, 3)
    h = {v: (1 if v in WITNESS_H else 0) for v in T5.V}
    CHECK(pred5 == 12 and is_broadcast(T5, h) and sigma(h) == 11
          and dominating(T5, h) and bn_independent(T5, h)
          and maximal_thm_minimal(T5, h) and maximal_prop_2bv(T5, h),
          "the problem demands i_bn(T5) = 12, while the printed cost-11 broadcast is a "
          "maximal bn-independent broadcast by both criteria, so i_bn(T5) in {10,11}",
          "predicted=%d witness cost=%d U_h^E=%s broadcasters per component of "
          "T5-U_h^E=%s"
          % (pred5, sigma(h), fmt_edges(uncovered_edges(T5, h)),
             [len(set(c) & set(WITNESS_H))
              for c in T5.components(uncovered_edges(T5, h))]))
    CHECK(10 != pred5 and 11 != pred5,
          "so the equality fails on T5 as well, and there is no parity obstruction to "
          "|M| being even",
          "10 != 12 and 11 != 12")

    # -- verdict ------------------------------------------------------------
    print("")
    if _BAD:
        print("VERDICT: %d CHECK(S) FAILED: %s" % (len(_BAD), "; ".join(_BAD)))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % _N)
    print("")
    print("NOT RE-RUN HERE: (1) no exhaustive search over the broadcasts of T or of T5 "
          "is performed, so no exact value of i_bn(T) is established here, and no exact "
          "value of i_bn(T5) is established either; what is established above is "
          "i_bn(T) in {8,9} and i_bn(T5) in {10,11}, which is all the refutation uses. "
          "(2) The lower "
          "bound gamma_b(G) <= i_bn(G) and the formula gamma_b(T) = (diam T - |M|)/2 "
          "for a maximum split-set are transcribed from the source as cited results "
          "and are not proved here; only their arithmetic on these trees is checked, "
          "and the value 8 (resp. 10) is independently corroborated by the exhibited "
          "dominating broadcast of that cost, which bounds gamma_b from above only. "
          "(3) The two maximality criteria the paper quotes from the source are "
          "transcribed as criteria; they are applied and cross-checked against each "
          "other and, on witness A, against a direct no-single-increment test, but "
          "they are not proved. (4) No bibliographic or attribution claim is checked: "
          "the wording of the problem, its numbering in either version, the reading of "
          "\"uniquely radial\", and the claim that no earlier work answers it. (5) The "
          "k-fold family of which T and T5 are the members k=4 and k=5 is not treated: "
          "nothing here bears on any other k. (6) Nothing is checked about trees other "
          "than T, T5 and the three small controls, and no repair of the equality is "
          "proposed or tested.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
