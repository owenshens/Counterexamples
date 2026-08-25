#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verifier for "A Counterexample to the Hibi--Seyed Fakhari Conjecture on
Join--Meet Ideals", which refutes [HSF, Conjecture 2.4]:
"the join-meet ideal of every nonmodular finite lattice is not linearly
related".  Standard library only; exact integer / GF(2) arithmetic.

=====================================================================
TAKEN FROM THE PAPER (inputs -- not checked, they define the object or
are the claims we test against)
=====================================================================
  * The lattice L: 12 elements {0hat,1hat} u {a_i,b_i : i in Z/5}, with
    a_0..a_4 the atoms, b_0..b_4 the coatoms, and b_i covering exactly
    a_i and a_{i+1} (indices mod 5).  Only these COVER relations are
    hard-coded; the order, joins, meets, ranks and the ideal are derived.
  * The graph6 string "IhCGGC@_G" for the atom-coatom incidence graph in
    the vertex order a_0,b_0,a_1,b_1,...,a_4,b_4.
  * The modular-law witness (x,y,z) = (a_0, a_2, b_0).
  * The claimed numbers, each compared against a recomputation:
        #generators 35 = 10 atom pairs + 10 coatom pairs + 15 mixed;
        generator shapes 5 x (A_iA_{i+1} - y B_i), 5 x (B_iB_{i+1} - A_{i+1} z),
        25 x (UV - yz);
        (E_d, V_d, c_d, dim Z_d) = (420,312,41,149), (2730,1293,91,1528),
        (12740,4276,161,8625) for d = 3,4,5;  dim Z_2 = 0;
        rank_F2(mu_4) = 1528 from a 1788 x 1528 matrix;
        rank_F2(mu_5) = 8625 from a 18336 x 8625 matrix;
        beta_{1,3} = 149, beta_{1,4} = beta_{1,5} = 0.

=====================================================================
DERIVED HERE (the actual checks; every one is computed from the cover
relations above and can fail if the object is corrupted)
=====================================================================
  1. L is a poset with 12 elements, one bottom, one top; every pair has a
     unique least upper and unique greatest lower bound (so L is a lattice);
     it is graded of rank 3 and pure; 5 atoms, 5 coatoms; the cover
     relations recovered from the derived order are exactly the input ones.
  2. The atom-coatom incidence graph is the 10-cycle C_10 and equals the
     decoding of the paper's graph6 string; and, as the paper's
     lattice-hood argument asserts, two distinct atoms lie under a common
     coatom exactly when they are consecutive (and then under a unique
     one, which is their join, else their join is 1hat), dually for
     coatoms, and every incomparable atom-coatom pair meets at 0hat and
     joins at 1hat.
  3. HYPOTHESIS of Conjecture 2.4: L is nonmodular -- found by an
     exhaustive search over all triples x <= z, and the paper's witness is
     re-evaluated element by element.
  4. The join-meet ideal I_L: all 66 unordered pairs are reduced, the 35
     nonzero quadrics are formed from the derived joins/meets, their shapes
     are classified, and minimality of the generating set is checked
     (distinct leading monomials, no leading monomial is a trailing one).
  5. Multigraphs Gamma_d for d = 2,3,4,5,6 (and d = 7 with --deep) built
     from the derived generators; E_d, V_d, c_d counted by union-find and
     dim Z_d = E_d - V_d + c_d.  dim Z_d is then recomputed independently
     as E_d - rank(incidence matrix) by real elimination over F_2, over
     F_3 and over F_{2^31-1}; and field-independence is certified
     structurally, by an exact upper bound (the oriented-incidence columns
     of each component sum to zero, so rank <= V_d - c_d over every field)
     together with an exact lower bound (the tree-edge x nonroot-vertex
     minor is triangular with +-1 diagonal, so rank >= V_d - c_d over
     every field).  The fundamental cycle basis is verified cycle by
     cycle: zero integer boundary and coordinate vector e_j.
  6. CONCLUSION of Conjecture 2.4 is violated: mu_d : S_1 (x) Z_{d-1} ->
     Z_d is built explicitly in the fundamental-cycle basis and its rank
     is computed by exact elimination over F_2 (the paper's computation)
     and, from the signed integer matrix, mod 2^31-1, which forces full
     rank over Q and hence over every field of characteristic zero (and,
     being a rank over F_{2^31-1}, over that field too; no other
     characteristic p >= 3 is covered).  This
     yields beta_{1,3} = 149 (the nonzero linear strand) and
     beta_{1,4} = beta_{1,5} = beta_{1,6} = 0 -- one degree beyond the
     paper.  Every one of the 1788 / 18336 / 103500 rows is certified to
     lie in the cycle space of Gamma_d (even degree at every vertex, zero
     integer boundary) and to be rebuilt exactly from its extracted
     coordinates, over GF(2) and over Z.
  7. NOT RE-RUN, printed as a NOTE by the program.  This is NOT the
     paper's full claim: the theorem asserts beta_{1,j} = 0 for every
     j >= 4, whereas only j = 4,5,6 (7 with --deep) are recomputed here.
     beta_{1,j} = 0 for j >= 7 (j >= 8 with --deep) is not recomputed;
     the paper covers all j >= 6 by citing [HSF, Lemma 2.8] for
     reg(I_L) = 4, a literature input, and no finite computation can
     settle infinitely many degrees.  Nor is full target rank of mu_j
     over a general field of characteristic p >= 3 re-run: see item 6.

Usage:  python3 verify.py            (64 checks, about 10 s, 0.2 GB)
        python3 verify.py --deep     (72 checks, about 50 s, 1.4 GB;
                                          adds degree 7)
Every line is "PASS <name> [detail]" or "FAIL <name> [detail]", the last
line is "VERDICT: ALL <n> CHECKS PASS" or "VERDICT: <k> OF <n> CHECKS
FAILED", and the exit status is 0 if and only if every check passed.
"""

import sys
import time
from itertools import combinations

FAILURES = []
TOTAL = 0


def check(name, ok, detail=""):
    """Record and print one independent check."""
    global TOTAL
    TOTAL += 1
    tag = "PASS" if ok else "FAIL"
    if not ok:
        FAILURES.append(name)
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    sys.stdout.flush()


def note(text):
    print("NOTE %s" % text)
    sys.stdout.flush()


GRAPH6 = "IhCGGC@_G"          # paper input: incidence graph, order a0,b0,a1,b1,...
ELEMENTS = ["0hat", "1hat"] + ["a%d" % i for i in range(5)] + ["b%d" % i for i in range(5)]


def paper_cover_pairs():
    """The ONLY hard-coded structure of L: its cover relations (lo, hi).

    b_i covers exactly a_i and a_{i+1}; the a_i are atoms, the b_i coatoms.
    """
    pairs = []
    for i in range(5):
        pairs.append(("0hat", "a%d" % i))
        pairs.append(("b%d" % i, "1hat"))
        pairs.append(("a%d" % i, "b%d" % i))
        pairs.append(("a%d" % ((i + 1) % 5), "b%d" % i))
    return sorted(set(pairs))


def transitive_closure(elements, cover_pairs):
    """Reflexive-transitive closure of the cover relation -> set of (x,y) with x<=y."""
    leq = set((e, e) for e in elements)
    leq |= set(cover_pairs)
    changed = True
    while changed:
        changed = False
        for (x, y) in list(leq):
            for (u, v) in list(cover_pairs):
                if y == u and (x, v) not in leq:
                    leq.add((x, v))
                    changed = True
    return leq


def derived_covers(elements, leq):
    """Cover relations read back off the derived order."""
    out = set()
    for x in elements:
        for y in elements:
            if x != y and (x, y) in leq:
                mid = [w for w in elements if w != x and w != y
                       and (x, w) in leq and (w, y) in leq]
                if not mid:
                    out.add((x, y))
    return out


def is_partial_order(elements, leq):
    """antisymmetry + transitivity (reflexivity holds by construction)."""
    for x in elements:
        if (x, x) not in leq:
            return False, "not reflexive at %s" % x
    for (x, y) in leq:
        if x != y and (y, x) in leq:
            return False, "antisymmetry fails on %s,%s" % (x, y)
    for (x, y) in leq:
        for z in elements:
            if (y, z) in leq and (x, z) not in leq:
                return False, "transitivity fails %s<=%s<=%s" % (x, y, z)
    return True, "12 elements, %d relations" % len(leq)


def join_meet(elements, leq):
    """Least upper / greatest lower bounds; returns (join, meet, problems)."""
    join, meet, problems = {}, {}, []
    for x in elements:
        for y in elements:
            ub = [w for w in elements if (x, w) in leq and (y, w) in leq]
            lb = [w for w in elements if (w, x) in leq and (w, y) in leq]
            least = [w for w in ub if all((w, v) in leq for v in ub)]
            greatest = [w for w in lb if all((v, w) in leq for v in lb)]
            if len(least) != 1:
                problems.append("join(%s,%s) has %d candidates" % (x, y, len(least)))
            else:
                join[(x, y)] = least[0]
            if len(greatest) != 1:
                problems.append("meet(%s,%s) has %d candidates" % (x, y, len(greatest)))
            else:
                meet[(x, y)] = greatest[0]
    return join, meet, problems


def rank_function(elements, leq, covers, bottom):
    """Chain-length rank; returns (rank dict, pure?, longest maximal chain)."""
    up = {}
    for (lo, hi) in covers:
        up.setdefault(lo, []).append(hi)
    rank, lengths = {bottom: 0}, []
    stack = [(bottom, 0)]
    seen_max = {}
    while stack:
        v, d = stack.pop()
        seen_max[v] = max(seen_max.get(v, -1), d)
        if v not in up:
            lengths.append(d)
        for w in up.get(v, []):
            stack.append((w, d + 1))
    # rank is well defined iff every maximal chain to v has the same length
    graded = True
    for v in elements:
        mins = min_chain_len(up, bottom, v)
        if mins != seen_max.get(v, None):
            graded = False
        rank[v] = seen_max.get(v, None)
    pure = len(set(lengths)) == 1
    return rank, graded, pure, (max(lengths) if lengths else -1)


def min_chain_len(up, bottom, target):
    """BFS shortest cover-path length bottom -> target (None if unreachable)."""
    if bottom == target:
        return 0
    frontier, dist = [bottom], {bottom: 0}
    while frontier:
        nxt = []
        for v in frontier:
            for w in up.get(v, []):
                if w not in dist:
                    dist[w] = dist[v] + 1
                    nxt.append(w)
        frontier = nxt
    return dist.get(target)


def graph6_decode(s):
    """Decode a graph6 string (n < 63) into (n, sorted edge list)."""
    if not s:
        raise ValueError("empty graph6")
    n = ord(s[0]) - 63
    if n < 0 or n > 62:
        raise ValueError("unsupported graph6 order byte")
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        if v < 0 or v > 63:
            raise ValueError("bad graph6 byte %r" % ch)
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    need = n * (n - 1) // 2
    if len(bits) < need:
        raise ValueError("graph6 payload too short: %d < %d" % (len(bits), need))
    edges, pos = [], 0
    for j in range(1, n):
        for i in range(j):
            if bits[pos]:
                edges.append((i, j))
            pos += 1
    if any(bits[need:]):
        raise ValueError("nonzero padding bits in graph6")
    return n, sorted(edges)


def cycle_edge_set(n):
    """Edges of the n-cycle 0-1-2-...-(n-1)-0 as sorted 2-tuples."""
    return sorted(tuple(sorted((i, (i + 1) % n))) for i in range(n))


def nonmodular_violations(elements, leq, join, meet):
    """All triples with x <= z and x v (y ^ z) != (x v y) ^ z."""
    bad = []
    for x in elements:
        for z in elements:
            if (x, z) not in leq:
                continue
            for y in elements:
                lhs = join[(x, meet[(y, z)])]
                rhs = meet[(join[(x, y)], z)]
                if lhs != rhs:
                    bad.append((x, y, z, lhs, rhs))
    return bad


VAR = dict((e, k) for k, e in enumerate(ELEMENTS))   # 0hat->0 (y), 1hat->1 (z)


def mono(*names):
    """Squarefree/any degree-2 monomial as a sorted tuple of variable indices."""
    return tuple(sorted(VAR[t] for t in names))


def build_generators(elements, leq, join, meet):
    """All x_p x_q - x_{p^q} x_{pvq} over unordered pairs; returns
    (nonzero generators as (lead, trail) monomial pairs, #zero pairs)."""
    gens, zero = [], 0
    for (p, q) in combinations(elements, 2):
        lead = mono(p, q)
        trail = mono(meet[(p, q)], join[(p, q)])
        if lead == trail:
            zero += 1
        else:
            gens.append((lead, trail))
    return gens, zero


def classify_pairs(elements, leq, atoms, coatoms):
    """Counts of incomparable pairs by kind: (atom-atom, coatom-coatom, mixed)."""
    aa = cc = mixed = 0
    for (p, q) in combinations(elements, 2):
        if (p, q) in leq or (q, p) in leq:
            continue
        if p in atoms and q in atoms:
            aa += 1
        elif p in coatoms and q in coatoms:
            cc += 1
        else:
            mixed += 1
    return aa, cc, mixed


def paper_predicted_generators():
    """The 10 non-'UV-yz' quadrics exactly as displayed in the paper."""
    pred = set()
    for i in range(5):
        pred.add((mono("a%d" % i, "a%d" % ((i + 1) % 5)),
                  mono("0hat", "b%d" % i)))
        pred.add((mono("b%d" % i, "b%d" % ((i + 1) % 5)),
                  mono("a%d" % ((i + 1) % 5), "1hat")))
    return pred


NVARS = len(ELEMENTS)          # 12 variables


def monomials(k):
    """All monomials of degree k in NVARS variables, as sorted index tuples."""
    from itertools import combinations_with_replacement
    return list(combinations_with_replacement(range(NVARS), k))


def build_gamma(gens, d):
    """Multigraph Gamma_d: an edge (g,u) for every generator g and every
    monomial u of degree d-2, joining the degree-d monomials u*lead(g) and
    u*trail(g).  Edge index = gi * (#u) + ui.  Returns a dict."""
    us = monomials(d - 2)
    uidx = dict((u, i) for i, u in enumerate(us))
    nu = len(us)
    vid = {}
    edges = [None] * (len(gens) * nu)
    for gi, (lead, trail) in enumerate(gens):
        base = gi * nu
        for ui, u in enumerate(us):
            m1 = tuple(sorted(u + lead))
            m2 = tuple(sorted(u + trail))
            v1 = vid.get(m1)
            if v1 is None:
                v1 = len(vid)
                vid[m1] = v1
            v2 = vid.get(m2)
            if v2 is None:
                v2 = len(vid)
                vid[m2] = v2
            edges[base + ui] = (v1, v2)
    return {"d": d, "us": us, "uidx": uidx, "nu": nu, "vid": vid,
            "edges": edges, "E": len(edges), "V": len(vid)}


def spanning_forest(G):
    """Union-find spanning forest; adds tree/nontree edge lists, c, adjacency."""
    V, edges = G["V"], G["edges"]
    par = list(range(V))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    tree, nontree = [], []
    for e, (a, b) in enumerate(edges):
        ra, rb = find(a), find(b)
        if ra == rb:
            nontree.append(e)
        else:
            par[ra] = rb
            tree.append(e)
    G["tree"], G["nontree"] = tree, nontree
    G["c"] = V - len(tree)
    G["col_of"] = dict((e, j) for j, e in enumerate(nontree))
    adj = [[] for _ in range(V)]
    for e in tree:
        a, b = edges[e]
        adj[a].append((b, e))
        adj[b].append((a, e))
    # rooted forest: parent vertex, parent edge, depth
    pv = [-1] * V
    pe = [-1] * V
    dep = [-1] * V
    for s in range(V):
        if dep[s] != -1:
            continue
        dep[s] = 0
        stack = [s]
        while stack:
            v = stack.pop()
            for (w, e) in adj[v]:
                if dep[w] == -1:
                    dep[w] = dep[v] + 1
                    pv[w], pe[w] = v, e
                    stack.append(w)
    G["pv"], G["pe"], G["dep"] = pv, pe, dep
    return G


def tree_path_steps(G, a, b):
    """Directed steps (edge, from, to) along the forest path a -> b."""
    pv, pe, dep = G["pv"], G["pe"], G["dep"]
    x, y = a, b
    up_a, up_b = [], []
    while dep[x] > dep[y]:
        up_a.append((pe[x], x, pv[x]))
        x = pv[x]
    while dep[y] > dep[x]:
        up_b.append((pe[y], y, pv[y]))
        y = pv[y]
    while x != y:
        up_a.append((pe[x], x, pv[x]))
        x = pv[x]
        up_b.append((pe[y], y, pv[y]))
        y = pv[y]
    return up_a + [(e, to, frm) for (e, frm, to) in reversed(up_b)]


def fundamental_cycle(G, e):
    """Fundamental cycle of non-tree edge e: (edge list, integer sign list).
    Orientation: boundary of edge f is [lead vertex] - [trail vertex]."""
    edges = G["edges"]
    a, b = edges[e]
    ed, sg = [e], [1]
    for (f, frm, to) in tree_path_steps(G, a, b):
        ed.append(f)
        sg.append(1 if edges[f] == (to, frm) else -1)
    return ed, sg


def boundary_is_zero(G, ed, sg):
    """Exact integer boundary check: sum of signs * ([lead]-[trail]) == 0."""
    acc = {}
    edges = G["edges"]
    for f, s in zip(ed, sg):
        v1, v2 = edges[f]
        acc[v1] = acc.get(v1, 0) + s
        acc[v2] = acc.get(v2, 0) - s
    return all(v == 0 for v in acc.values())


def even_degree(G, ed):
    """GF(2) cycle-space membership: every vertex meets an even number of edges."""
    deg = {}
    edges = G["edges"]
    for f in ed:
        v1, v2 = edges[f]
        deg[v1] = deg.get(v1, 0) + 1
        deg[v2] = deg.get(v2, 0) + 1
    return all(v % 2 == 0 for v in deg.values())


def f2_rank(rows, ncols, early=True):
    """Exact GF(2) Gaussian elimination on rows given as column bitmasks.
    Returns (rank, #rows consumed).  Stops early once rank == ncols, which
    is legitimate because rank <= ncols always."""
    pivots = {}
    used = 0
    for row in rows:
        used += 1
        r = row
        while r:
            b = r.bit_length() - 1
            q = pivots.get(b)
            if q is None:
                pivots[b] = r
                break
            r ^= q
        if early and len(pivots) >= ncols:
            break
    return len(pivots), used


def f2_incidence_rank(G):
    """GF(2) rank of the (unsigned = signed mod 2) incidence matrix of G."""
    rows = (((1 << a) | (1 << b)) for (a, b) in G["edges"])
    rk, _ = f2_rank(rows, G["V"], early=False)
    return rk


def modp_rank(rows, p):
    """Exact rank mod an odd prime p; rows are dicts {col: value}."""
    pivots = {}
    for row in rows:
        r = dict((c, v % p) for c, v in row.items() if v % p)
        while r:
            c = max(r)
            piv = pivots.get(c)
            if piv is None:
                inv = pow(r[c], p - 2, p)
                pivots[c] = dict((k, (v * inv) % p) for k, v in r.items())
                break
            f = r[c]
            for k, v in piv.items():
                nv = (r.get(k, 0) - f * v) % p
                if nv:
                    r[k] = nv
                elif k in r:
                    del r[k]
    return len(pivots)


def oriented_incidence_rows(G):
    """Rows of the oriented incidence matrix: boundary of each edge."""
    for (a, b) in G["edges"]:
        if a == b:
            yield {}
        else:
            yield {a: 1, b: -1}


def all_fundamental_cycles(G):
    """Fundamental cycle basis of the cycle space of G, indexed by column."""
    return [fundamental_cycle(G, e) for e in G["nontree"]]


def multiply_cycle(Glow, Ghigh, ed, sg, r):
    """Image of a cycle of Gamma_{d-1} under multiplication by variable r:
    edge (g,u) |-> edge (g, x_r u).  Signs are preserved because lead and
    trail monomials both get multiplied by x_r."""
    nu_lo, us_lo = Glow["nu"], Glow["us"]
    nu_hi, uidx_hi = Ghigh["nu"], Ghigh["uidx"]
    ed2 = []
    for f in ed:
        gi, ui = divmod(f, nu_lo)
        u2 = tuple(sorted(us_lo[ui] + (r,)))
        ed2.append(gi * nu_hi + uidx_hi[u2])
    return ed2, list(sg)


def coords_f2(Ghigh, ed2):
    """GF(2) coordinates of a cycle in the fundamental basis = its incidences
    with the non-tree edges."""
    col_of = Ghigh["col_of"]
    mask = 0
    for e in ed2:
        c = col_of.get(e)
        if c is not None:
            mask ^= (1 << c)
    return mask


def coords_z(Ghigh, ed2, sg2):
    """Integer coordinates in the oriented fundamental basis."""
    col_of = Ghigh["col_of"]
    out = {}
    for e, s in zip(ed2, sg2):
        c = col_of.get(e)
        if c is not None:
            out[c] = out.get(c, 0) + s
    return dict((c, v) for c, v in out.items() if v)


def reconstruct_f2(cycles, mask):
    """Symmetric difference of the fundamental cycles selected by mask."""
    acc = set()
    while mask:
        b = mask.bit_length() - 1
        mask ^= (1 << b)
        acc ^= set(cycles[b][0])
    return acc


def reconstruct_z(cycles, coeffs):
    """Integer combination of oriented fundamental cycles -> {edge: coeff}."""
    acc = {}
    for c, a in coeffs.items():
        ed, sg = cycles[c]
        for f, s in zip(ed, sg):
            nv = acc.get(f, 0) + a * s
            if nv:
                acc[f] = nv
            elif f in acc:
                del acc[f]
    return acc


def mu_row_masks(Glow, Ghigh, Clow):
    """GF(2) rows of mu_d : S_1 (x) Z_{d-1} -> Z_d, in the fundamental basis."""
    for r in range(NVARS):
        for (ed, sg) in Clow:
            ed2, sg2 = multiply_cycle(Glow, Ghigh, ed, sg, r)
            yield coords_f2(Ghigh, ed2)


def certify_mu_rows(Glow, Ghigh, Clow, Chigh, stride=1):
    """For every (stride-th) row of mu_d check that the image really is a
    cycle of Gamma_d and that the extracted coordinates rebuild it exactly,
    over GF(2) and over Z.  Returns (tested, bad_cycle, bad_f2, bad_z)."""
    tested = bad_cycle = bad_f2 = bad_z = 0
    k = -1
    for r in range(NVARS):
        for (ed, sg) in Clow:
            k += 1
            if k % stride:
                continue
            tested += 1
            ed2, sg2 = multiply_cycle(Glow, Ghigh, ed, sg, r)
            if len(set(ed2)) != len(ed2):
                bad_cycle += 1
                continue
            if not (even_degree(Ghigh, ed2) and boundary_is_zero(Ghigh, ed2, sg2)):
                bad_cycle += 1
            mask = coords_f2(Ghigh, ed2)
            if reconstruct_f2(Chigh, mask) != set(ed2):
                bad_f2 += 1
            want = dict(zip(ed2, sg2))
            if reconstruct_z(Chigh, coords_z(Ghigh, ed2, sg2)) != want:
                bad_z += 1
    return tested, bad_cycle, bad_f2, bad_z


def component_column_relations(G):
    """Sum of the oriented-incidence columns over each component is zero, so
    rank <= V - c over EVERY field.  Verified here as an exact integer check."""
    par = {}

    def find(x):
        while par.get(x, x) != x:
            x = par[x]
        return x

    for e in G["tree"]:
        a, b = G["edges"][e]
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    root = [find(v) for v in range(G["V"])]
    ncomp = len(set(root))
    if ncomp != G["c"]:
        return False, "component count mismatch %d vs %d" % (ncomp, G["c"])
    for (a, b) in G["edges"]:
        if root[a] != root[b]:
            return False, "edge crosses two forest components"
    return True, "%d components, all column sums vanish" % ncomp


def part_lattice():
    """Checks 1-12: the exhibited object, decoded, counted and printed back."""
    covers = paper_cover_pairs()
    check("element_count_is_12",
          len(ELEMENTS) == 12 and len(set(ELEMENTS)) == 12,
          "L = " + ",".join(ELEMENTS))
    leq = transitive_closure(ELEMENTS, covers)
    ok, why = is_partial_order(ELEMENTS, leq)
    check("is_partial_order", ok, why)
    dc = derived_covers(ELEMENTS, leq)
    check("cover_relations_recovered_from_order",
          dc == set(covers) and len(dc) == 20,
          "%d covers (20 expected), input cover set reproduced: %s"
          % (len(dc), dc == set(covers)))
    below_b = dict(("b%d" % i,
                    sorted(x for (x, y) in dc if y == "b%d" % i)) for i in range(5))
    want_b = dict(("b%d" % i, sorted(["a%d" % i, "a%d" % ((i + 1) % 5)]))
                  for i in range(5))
    check("each_coatom_covers_exactly_two_consecutive_atoms",
          below_b == want_b,
          "; ".join("%s>%s" % (k, ",".join(v)) for k, v in sorted(below_b.items())))
    bottoms = [x for x in ELEMENTS if all((x, y) in leq for y in ELEMENTS)]
    tops = [x for x in ELEMENTS if all((y, x) in leq for y in ELEMENTS)]
    check("bounded_with_unique_0hat_and_1hat",
          bottoms == ["0hat"] and tops == ["1hat"],
          "min=%s max=%s" % (bottoms, tops))
    join, meet, problems = join_meet(ELEMENTS, leq)
    check("every_pair_has_unique_join_and_meet_so_L_is_a_lattice",
          not problems and len(join) == 144 and len(meet) == 144,
          "%d/144 joins, %d/144 meets well defined, %d defects%s"
          % (len(join), len(meet), len(problems),
             (": " + problems[0]) if problems else ""))
    atoms = sorted(y for (x, y) in dc if x == "0hat")
    coatoms = sorted(x for (x, y) in dc if y == "1hat")
    check("five_atoms_five_coatoms",
          atoms == ["a%d" % i for i in range(5)]
          and coatoms == ["b%d" % i for i in range(5)],
          "atoms=%s coatoms=%s" % (",".join(atoms), ",".join(coatoms)))
    rank, graded, pure, longest = rank_function(ELEMENTS, leq, covers, "0hat")
    check("graded_pure_of_rank_3",
          graded and pure and longest == 3
          and all(rank[x] == 1 for x in atoms) and all(rank[x] == 2 for x in coatoms)
          and rank["1hat"] == 3,
          "graded=%s pure=%s longest maximal chain=%d, ranks a_i=%s b_i=%s 1hat=%s"
          % (graded, pure, longest, sorted(set(rank[x] for x in atoms)),
             sorted(set(rank[x] for x in coatoms)), rank["1hat"]))
    return leq, join, meet, atoms, coatoms, covers, dc


def part_incidence_graph(leq, atoms, coatoms):
    """Checks: the atom-coatom incidence graph is C_10 and matches graph6."""
    vid = {}
    for i in range(5):
        vid["a%d" % i] = 2 * i
        vid["b%d" % i] = 2 * i + 1
    edges = sorted(tuple(sorted((vid[a], vid[b])))
                   for a in atoms for b in coatoms if (a, b) in leq)
    c10 = cycle_edge_set(10)
    check("incidence_graph_is_the_10_cycle",
          edges == c10 and len(edges) == 10,
          "edges=%s" % ("".join("(%d,%d)" % e for e in edges)))
    deg = {}
    for (u, v) in edges:
        deg[u] = deg.get(u, 0) + 1
        deg[v] = deg.get(v, 0) + 1
    connected = len(edges) == 10 and len(deg) == 10 and set(deg.values()) == set([2])
    check("incidence_graph_2_regular_on_10_vertices", connected,
          "degrees %s" % sorted(deg.values()))
    try:
        n6, e6 = graph6_decode(GRAPH6)
        ok = (n6 == 10 and e6 == edges)
        detail = "graph6 %s -> n=%d, %d edges, equal to derived graph: %s" % (
            GRAPH6, n6, len(e6), ok)
    except ValueError as exc:
        ok, detail = False, "graph6 decode error: %s" % exc
    check("graph6_string_decodes_to_the_derived_incidence_graph", ok, detail)


def part_covering_structure(leq, join, meet, atoms, coatoms):
    """The paper's lattice-hood argument, checked pair by pair."""
    def consec(i, j, n=5):
        return (j - i) % n == 1 or (i - j) % n == 1

    bad = []
    for i in range(5):
        for j in range(i + 1, 5):
            a, b = "a%d" % i, "a%d" % j
            common = [c for c in coatoms if (a, c) in leq and (b, c) in leq]
            if consec(i, j):
                if len(common) != 1 or join[(a, b)] != common[0]:
                    bad.append("%s,%s common coatoms %s join %s"
                               % (a, b, common, join[(a, b)]))
            elif common or join[(a, b)] != "1hat":
                bad.append("%s,%s common coatoms %s join %s"
                           % (a, b, common, join[(a, b)]))
    check("atom_pairs_have_a_unique_common_coatom_iff_consecutive_else_join_1hat",
          not bad, "all 10 atom pairs behave as claimed" if not bad else bad[0])
    bad = []
    for i in range(5):
        for j in range(i + 1, 5):
            a, b = "b%d" % i, "b%d" % j
            common = [c for c in atoms if (c, a) in leq and (c, b) in leq]
            if consec(i, j):
                if len(common) != 1 or meet[(a, b)] != common[0]:
                    bad.append("%s,%s common atoms %s meet %s"
                               % (a, b, common, meet[(a, b)]))
            elif common or meet[(a, b)] != "0hat":
                bad.append("%s,%s common atoms %s meet %s"
                           % (a, b, common, meet[(a, b)]))
    check("coatom_pairs_have_a_unique_common_atom_iff_consecutive_else_meet_0hat",
          not bad, "all 10 coatom pairs behave as claimed" if not bad else bad[0])
    bad = []
    n = 0
    for a in atoms:
        for b in coatoms:
            if (a, b) in leq:
                continue
            n += 1
            if meet[(a, b)] != "0hat" or join[(a, b)] != "1hat":
                bad.append("%s,%s meet %s join %s" % (a, b, meet[(a, b)], join[(a, b)]))
    check("incomparable_atom_coatom_pairs_meet_at_0hat_and_join_at_1hat",
          not bad and n == 15,
          "%d incomparable mixed pairs, %d misbehaving" % (n, len(bad)))


def part_nonmodular(leq, join, meet):
    """HYPOTHESIS of Conjecture 2.4: L must be nonmodular."""
    x, y, z = "a0", "a2", "b0"
    lhs = join[(x, meet[(y, z)])]
    rhs = meet[(join[(x, y)], z)]
    ok = ((x, z) in leq and meet[(y, z)] == "0hat" and join[(x, y)] == "1hat"
          and lhs == "a0" and rhs == "b0" and lhs != rhs)
    check("paper_witness_a0_a2_b0_breaks_the_modular_law", ok,
          "a0<=b0 is %s, a2^b0=%s, a0va2=%s, a0v(a2^b0)=%s, (a0va2)^b0=%s, differ=%s"
          % ((x, z) in leq, meet[(y, z)], join[(x, y)], lhs, rhs, lhs != rhs))
    bad = nonmodular_violations(ELEMENTS, leq, join, meet)
    check("L_is_nonmodular_exhaustive_search", len(bad) > 0,
          "%d violating triples (x<=z) out of %d; first=%s"
          % (len(bad), sum(1 for a in ELEMENTS for b in ELEMENTS
                           for c in ELEMENTS if (a, c) in leq), bad[0][:3]))
    return len(bad)


DISPLAY = dict([(VAR["0hat"], "y"), (VAR["1hat"], "z")]
               + [(VAR["a%d" % i], "A%d" % i) for i in range(5)]
               + [(VAR["b%d" % i], "B%d" % i) for i in range(5)])


def mono_str(m):
    return "*".join(DISPLAY[v] for v in m)


def gen_str(g):
    return "%s-%s" % (mono_str(g[0]), mono_str(g[1]))


def part_ideal(leq, join, meet, atoms, coatoms):
    """Checks: the 35 join-meet quadrics, their shapes, and minimality."""
    gens, zero = build_generators(ELEMENTS, leq, join, meet)
    check("join_meet_ideal_has_35_nonzero_quadrics",
          len(gens) == 35 and zero == 31 and len(gens) + zero == 66,
          "%d unordered pairs = %d nonzero quadrics (35 claimed) + %d vanishing "
          "(comparable) pairs" % (len(gens) + zero, len(gens), zero))
    note("the 35 generators: " + "  ".join(sorted(gen_str(g) for g in gens)))
    aa, cc, mixed = classify_pairs(ELEMENTS, leq, set(atoms), set(coatoms))
    check("incomparable_pairs_are_10_atom_10_coatom_15_mixed",
          (aa, cc, mixed) == (10, 10, 15) and aa + cc + mixed == 35,
          "atom-atom=%d coatom-coatom=%d mixed=%d" % (aa, cc, mixed))
    gset = set(gens)
    pred = paper_predicted_generators()
    found = sorted(gen_str(g) for g in pred if g in gset)
    check("ten_displayed_quadrics_AiAi1_yBi_and_BiBi1_Ai1z_are_exactly_reproduced",
          pred <= gset and len(pred) == 10,
          "%d of the 10 displayed quadrics are among the derived generators: %s"
          % (len(found), "  ".join(found)))
    yz = mono("0hat", "1hat")
    rest = [g for g in gens if g not in pred]
    middle = set(VAR[e] for e in atoms) | set(VAR[e] for e in coatoms)
    ok_rest = (len(rest) == 25
               and all(g[1] == yz for g in rest)
               and all(len(set(g[0])) == 2 and set(g[0]) <= middle for g in rest)
               and len(set(g[0] for g in rest)) == 25)
    check("remaining_25_quadrics_all_have_shape_UV_minus_yz", ok_rest,
          "%d remaining quadrics (25 claimed), %d with trailing term y*z, %d with "
          "leading term U*V in the ten middle variables, %d distinct leading terms"
          % (len(rest), sum(1 for g in rest if g[1] == yz),
             sum(1 for g in rest if len(set(g[0])) == 2 and set(g[0]) <= middle),
             len(set(g[0] for g in rest))))
    leads = [g[0] for g in gens]
    trails = [g[1] for g in gens]
    check("leading_monomials_distinct_and_never_trailing",
          len(set(leads)) == 35 and not (set(leads) & set(trails)),
          "%d distinct leading monomials, %d distinct trailing, overlap %d"
          % (len(set(leads)), len(set(trails)), len(set(leads) & set(trails))))
    return gens


PAPER_TABLE = {3: (420, 312, 41, 149),        # (E_d, V_d, c_d, dim Z_d), paper input
               4: (2730, 1293, 91, 1528),
               5: (12740, 4276, 161, 8625)}
PAPER_MU = {4: (1788, 1528), 5: (18336, 8625)}   # (rows, cols) = (rank), paper input
BIGPRIME = 2147483647                            # 2^31 - 1


def unimodular_forest_minor(G):
    """The tree-edge x nonroot-vertex minor of the oriented incidence matrix is
    triangular with +-1 diagonal (order rows/cols by decreasing depth), so its
    determinant is +-1 and rank >= V - c over EVERY field."""
    pv, pe, dep, edges = G["pv"], G["pe"], G["dep"], G["edges"]
    nonroot = [v for v in range(G["V"]) if pv[v] != -1]
    if len(nonroot) != G["V"] - G["c"]:
        return False, "nonroot count %d != V-c %d" % (len(nonroot), G["V"] - G["c"])
    if sorted(pe[v] for v in nonroot) != sorted(G["tree"]):
        return False, "child->tree-edge map is not a bijection"
    for v in nonroot:
        a, b = edges[pe[v]]
        if set([a, b]) != set([v, pv[v]]) or dep[pv[v]] != dep[v] - 1:
            return False, "parent edge of vertex %d is not (v,parent)" % v
    return True, "%dx%d triangular minor, determinant +-1" % (len(nonroot), len(nonroot))


def part_graphs(gens, degrees):
    """Checks: Gamma_d counts, dim Z_d by three independent computations."""
    Gs, Cs = {}, {}
    for d in degrees:
        G = spanning_forest(build_gamma(gens, d))
        Gs[d] = G
        E, V, c = G["E"], G["V"], G["c"]
        dim = E - V + c
        note("gamma_%d: E=%d edges, V=%d vertices, c=%d components, "
             "E-V+c=%d, self-loops=%d (so phi_%d has no zero rows)"
             % (d, E, V, c, dim, sum(1 for (a, b) in G["edges"] if a == b), d))
        if d in PAPER_TABLE:
            pE, pV, pc, pdim = PAPER_TABLE[d]
            check("gamma_%d_E_V_c_match_the_paper_table" % d,
                  (E, V, c) == (pE, pV, pc),
                  "computed (E,V,c)=(%d,%d,%d), paper (%d,%d,%d)" % (E, V, c, pE, pV, pc))
            check("dim_Z_%d_equals_paper_value_%d" % (d, pdim),
                  dim == pdim and len(G["nontree"]) == pdim,
                  "E-V+c=%d, non-tree edges=%d, paper=%d" % (dim, len(G["nontree"]), pdim))
        elif d == 2:
            check("dim_Z_2_is_zero_so_the_35_quadrics_are_linearly_independent",
                  dim == 0 and len(G["nontree"]) == 0,
                  "E-V+c = %d-%d+%d = %d, paper claims 0" % (E, V, c, dim))
        rk2 = f2_incidence_rank(G)
        check("dim_Z_%d_reproduced_by_GF2_elimination" % d,
              rk2 == V - c and E - rk2 == dim,
              "rank_F2(incidence)=%d, V-c=%d, E-rank=%d, E-V+c=%d"
              % (rk2, V - c, E - rk2, dim))
        rk3 = modp_rank(oriented_incidence_rows(G), 3)
        rkp = modp_rank(oriented_incidence_rows(G), BIGPRIME)
        check("dim_Z_%d_same_over_F3_and_over_F_%d" % (d, BIGPRIME),
              rk3 == V - c and rkp == V - c,
              "rank mod 3 = %d, rank mod %d = %d, V-c = %d" % (rk3, BIGPRIME, rkp, V - c))
        ok_up, why_up = component_column_relations(G)
        ok_lo, why_lo = unimodular_forest_minor(G)
        check("dim_Z_%d_is_field_independent_by_rank_bounds" % d, ok_up and ok_lo,
              "upper: %s; lower: %s" % (why_up, why_lo))
        Cs[d] = all_fundamental_cycles(G)
        bad = 0
        for jj, (ed, sg) in enumerate(Cs[d]):
            if not (boundary_is_zero(G, ed, sg) and even_degree(G, ed)
                    and coords_f2(G, ed) == (1 << jj)
                    and coords_z(G, ed, sg) == {jj: 1}):
                bad += 1
        check("fundamental_cycle_basis_of_Z_%d_is_valid_and_unitriangular" % d,
              len(Cs[d]) == dim and bad == 0,
              "%d cycles for dim Z_%d = %d, %d with nonzero integer boundary or "
              "coordinate vector != e_j" % (len(Cs[d]), d, dim, bad))
    return Gs, Cs


def mu_row_coeffs(Glow, Ghigh, Clow):
    """Integer (signed) rows of mu_d in the oriented fundamental basis."""
    for r in range(NVARS):
        for (ed, sg) in Clow:
            ed2, sg2 = multiply_cycle(Glow, Ghigh, ed, sg, r)
            yield coords_z(Ghigh, ed2, sg2)


def part_betti(Gs, Cs, degrees):
    """Checks: mu_d shapes, row certification, GF(2) and char-0 ranks, betti."""
    betti = {}
    d3 = min(d for d in degrees if d >= 3)
    if d3 == 3 and 2 in Cs:
        rk3, used3 = f2_rank(mu_row_masks(Gs[2], Gs[3], Cs[2]), len(Cs[3]))
        betti[3] = len(Cs[3]) - rk3
        check("beta_1_3_equals_149_the_nonzero_linear_strand",
              betti[3] == 149 and rk3 == 0 and used3 == 12 * len(Cs[2]),
              "dim Z_3 - rank(mu_3) = %d - %d = %d (mu_3 has %d rows since Z_2=0)"
              % (len(Cs[3]), rk3, betti[3], 12 * len(Cs[2])))
    for d in sorted(degrees):
        if d < 4 or (d - 1) not in Cs:
            continue
        Glo, Ghi, Clo, Chi = Gs[d - 1], Gs[d], Cs[d - 1], Cs[d]
        nrows, ncols = NVARS * len(Clo), len(Chi)
        if d in PAPER_MU:
            check("mu_%d_matrix_shape_is_%dx%d_as_stated" % ((d,) + PAPER_MU[d]),
                  (nrows, ncols) == PAPER_MU[d],
                  "12 * dim Z_%d = %d rows, dim Z_%d = %d columns"
                  % (d - 1, nrows, d, ncols))
        else:
            note("mu_%d is a %d x %d matrix over F_2 (no paper claim; extra degree)"
                 % (d, nrows, ncols))
        tested, bc, bf, bz = certify_mu_rows(Glo, Ghi, Clo, Chi)
        check("mu_%d_every_row_is_a_certified_cycle_with_exact_coordinates" % d,
              tested == nrows and bc == 0 and bf == 0 and bz == 0,
              "%d/%d rows: cycle-space %d bad, GF(2) rebuild %d bad, "
              "integral rebuild %d bad" % (tested, nrows, bc, bf, bz))
        rk, used = f2_rank(mu_row_masks(Glo, Ghi, Clo), ncols)
        detail = "rank_F2 = %d of %d columns (%d of %d rows consumed)" % (
            rk, ncols, used, nrows)
        if d in PAPER_MU:
            check("mu_%d_rank_over_F2_is_%d" % (d, PAPER_MU[d][1]),
                  rk == PAPER_MU[d][1] and rk == ncols, detail)
        else:
            check("mu_%d_has_full_target_rank_over_F2" % d, rk == ncols, detail)
        rkp = modp_rank(mu_row_coeffs(Glo, Ghi, Clo), BIGPRIME)
        check("mu_%d_has_full_target_rank_in_characteristic_zero" % d,
              rkp == ncols,
              "signed integer matrix has rank %d mod %d, columns = %d; equality "
              "forces full rank over Q and hence over every field of "
              "characteristic zero, and it also certifies full rank over "
              "F_%d itself, but NOT over a general field of characteristic "
              "p>=3" % (rkp, BIGPRIME, ncols, BIGPRIME))
        betti[d] = ncols - rk
        check("beta_1_%d_is_zero_over_F2_and_in_characteristic_zero" % d,
              betti[d] == 0 and ncols - rkp == 0,
              "beta_1_%d = dim Z_%d - rank(mu_%d) = %d - %d = %d"
              % (d, d, d, ncols, rk, betti[d]))
    return betti


def main(argv):
    deep = "--deep" in argv[1:]
    t0 = time.time()
    print("verify.py -- A Counterexample to the Hibi-Seyed Fakhari Conjecture")
    print("               on Join-Meet Ideals (refutes [HSF, Conjecture 2.4])")
    print("stdlib only, exact integer / GF(2) arithmetic, python %d.%d"
          % sys.version_info[:2])
    leq, join, meet, atoms, coatoms, covers, dc = part_lattice()
    if len(join) != 144 or len(meet) != 144:
        note("ABORT: the exhibited object is not a lattice, so the join-meet "
             "ideal is undefined and the remaining checks cannot be run.")
        return verdict()
    part_incidence_graph(leq, atoms, coatoms)
    part_covering_structure(leq, join, meet, atoms, coatoms)
    nviol = part_nonmodular(leq, join, meet)
    gens = part_ideal(leq, join, meet, atoms, coatoms)
    degrees = [2, 3, 4, 5, 6] + ([7] if deep else [])
    Gs, Cs = part_graphs(gens, degrees)
    betti = part_betti(Gs, Cs, degrees)
    computed = sorted(d for d in betti if d >= 4)
    check("counterexample_stands_nonmodular_lattice_with_no_nonlinear_syzygies",
          nviol > 0 and betti.get(3) == 149
          and computed == list(range(4, max(degrees) + 1))
          and all(betti[d] == 0 for d in computed),
          "12-element lattice, %d modular-law failures, beta_1,3=%s, and "
          "(beta_1,j for j=%s) = (%s) over F_2 and in characteristic 0"
          % (nviol, betti.get(3), ",".join(str(d) for d in computed),
             ",".join(str(betti[d]) for d in computed)))
    note("NOT RE-RUN: this program does NOT re-derive the paper's full claim. "
         "The published theorem asserts beta_{1,j}=0 for EVERY j>=4, and only "
         "finitely many degrees can be recomputed. beta_{1,j}=0 for j>=%d is "
         "NOT recomputed here; the paper obtains it from the cited "
         "[HSF, Lemma 2.8] regularity bound reg(I_L)=4, a literature input, "
         "and no finite computation can cover infinitely many degrees. "
         "Degrees j=4,5%s were recomputed in full; the paper itself computes "
         "only j=4,5. Also NOT re-run: full target rank of mu_j over a "
         "general field of characteristic p>=3 -- the ranks printed above are "
         "certified over F_2, over Q and hence over every field of "
         "characteristic zero, and over F_%d, and over no other "
         "characteristic."
         % (max(degrees) + 1, "" if max(degrees) < 6 else
            "," + ",".join(str(d) for d in range(6, max(degrees) + 1)),
            BIGPRIME))
    if not deep:
        note("re-run with --deep to add degree 7 (mu_7 is 429144 x 121549 over "
             "F_2; needs about 1.3 GB of RAM and a few more seconds).")
    note("elapsed %.1f s" % (time.time() - t0))
    return verdict()


def verdict():
    if FAILURES:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(FAILURES), TOTAL))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % TOTAL)
    return 0


if __name__ == "__main__":
    try:
        RC = main(sys.argv)
    except Exception as exc:                       # never exit without a verdict
        check("verifier_ran_to_completion", False, "unhandled exception: %r" % (exc,))
        RC = verdict()
    sys.exit(RC)
