#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- independent re-derivation of every computational claim in

    "An Ordered Pair of Graphs of Nonzero Twin-Width Attaining the Pettersson--Sylvester Tensor-Product Boundn-Width for the
     Pettersson-Sylvester Tensor-Product Bound"

The result under test: for G = C_8(1,2) (the square of the 8-cycle) and H = C_5,

    tww(G x H) = 8 = max{tww(G)*Delta(H) + Delta(H), tww(H)} + Delta(H)
                   = max{2*2 + 2, 2} + 2,

with tww(G) = tww(H) = 2, both NON-ZERO -- which answers the first question of
the Conclusion of Pettersson & Sylvester, DMTCS 25:1 (2023) Paper No. 18.

WHAT THIS PROGRAM READS: only what the paper prints -- the two circulant
parameter sets, the two contraction sequences, and the 131-character graph6
string of the 40-vertex product. No external data file, no third-party module.

WHAT IT DOES:
  * decodes the paper's graph6 string and checks LABEL-FOR-LABEL that it is
    G x H under the labelling k = 5u + i printed in the paper;
  * recomputes the full s- and b-profiles of both factors (the minimum over ALL
    pairs, which is what Ahn-Hendrey-Kim-Oum Lemma 3.1 bounds tww below by);
  * re-verifies the paper's two contraction sequences at width 2 in a
    from-scratch trigraph simulator, giving tww = 2 EXACTLY for both factors
    with no exhaustive search anywhere in the load-bearing path;
  * enumerates all C(40,2) = 780 pairs of the product, splits them into the
    paper's three cases and reports the per-case minima;
  * evaluates both orderings of the bound;
  * BUILDS an explicit width-8 contraction sequence of the 40-vertex product
    and re-verifies it in the same trigraph simulator, so that tww(G x H) = 8
    is certified WITHOUT invoking the theorem under test -- the upper bound is
    then a shipped, re-checkable object rather than an appeal to thm:tensor;
  * repeats both bounds for the family H = C_n, n = 5..12;
  * censuses all 121 ordered pairs of non-cographs on at most 5 vertices and
    decides every one of them NOT tight, which is the only completeness
    statement this note makes;
  * runs controls of both polarities against independently known values.

WHAT IT DOES *NOT* DO is stated in the NOT RE-RUN lines at the end of the run,
and in the `## Scope` section of REVIEW_NOTE.md.

Exact integer and set arithmetic throughout; no floating point, no randomness,
no solver. Python 3.9+, standard library only. Exits 0 iff every check passed.
"""

import sys
from itertools import combinations, permutations

# ---------------------------------------------------------------------------
# the check harness
# ---------------------------------------------------------------------------
_STATE = {'pass': 0, 'fail': 0}


def ck(name, ok, detail=''):
    tag = 'PASS' if ok else 'FAIL'
    _STATE['pass' if ok else 'fail'] += 1
    print('%s %s%s' % (tag, name, ('   [%s]' % detail) if detail else ''))


def note(*a):
    print('  ', *a)


# ---------------------------------------------------------------------------
# graphs as dict: vertex -> frozenset of neighbours
# ---------------------------------------------------------------------------
def circulant(n, diffs):
    """C_n(d1,...,dk): u ~ v iff u - v = +-di mod n."""
    g = {u: set() for u in range(n)}
    for u in range(n):
        for d in diffs:
            for v in ((u + d) % n, (u - d) % n):
                if v != u:
                    g[u].add(v)
                    g[v].add(u)
    return {u: frozenset(s) for u, s in g.items()}


def cycle(n):
    return circulant(n, [1])


def complete(n):
    return {u: frozenset(v for v in range(n) if v != u) for u in range(n)}


def path(n):
    g = {u: set() for u in range(n)}
    for u in range(n - 1):
        g[u].add(u + 1)
        g[u + 1].add(u)
    return {u: frozenset(s) for u, s in g.items()}


def edges_of(g):
    return sorted((u, v) for u in g for v in g[u] if u < v)


def tensor(g, h):
    """(u,i)(v,j) is an edge iff uv in E(G) AND ij in E(H)."""
    p = {(u, i): set() for u in g for i in h}
    for u in g:
        for i in h:
            for v in g[u]:
                for j in h[i]:
                    p[(u, i)].add((v, j))
    return {k: frozenset(s) for k, s in p.items()}


def complement(g):
    vs = set(g)
    return {u: frozenset(vs - g[u] - {u}) for u in g}


# ---------------------------------------------------------------------------
# the two Ahn et al. quantities
# ---------------------------------------------------------------------------
def s_val(g, u, v):
    """|N(u) sym-diff N(v)|, no deletion."""
    return len(g[u] ^ g[v])


def b_val(g, u, v):
    """|(N(u) sym-diff N(v)) \\ {u,v}| -- the quantity of Ahn et al. Lemma 3.1."""
    return len((g[u] ^ g[v]) - {u, v})


def b_min(g):
    """(minimum of b over ALL pairs u != v, an argmin, number of pairs)."""
    best, arg, cnt = None, None, 0
    for u, v in combinations(sorted(g), 2):
        cnt += 1
        x = b_val(g, u, v)
        if best is None or x < best:
            best, arg = x, (u, v)
    return best, arg, cnt


def s_min(g):
    return min(s_val(g, u, v) for u, v in combinations(sorted(g), 2))


# ---------------------------------------------------------------------------
# a from-scratch trigraph simulator (the contraction-sequence definition)
# ---------------------------------------------------------------------------
def sequence_width(g, merges):
    """Parts start as singletons.  `merges` is a list of (a, b): merge the parts
    currently containing a and b.  Between two distinct parts P, Q the edge is
    BLACK if every p-q pair is an edge of g, RED if some pair is an edge and
    some is not, ABSENT if no pair is an edge.  The width of the sequence is the
    maximum, over EVERY state including the initial singleton state, of the
    maximum red degree of a part.  Asserts the sequence ends in one part."""
    parts = [frozenset([v]) for v in sorted(g)]

    def red_deg_max():
        w = 0
        for i, pp in enumerate(parts):
            d = 0
            for j, qq in enumerate(parts):
                if i == j:
                    continue
                cnt = sum(1 for p in pp for q in qq if q in g[p])
                if 0 < cnt < len(pp) * len(qq):
                    d += 1
            w = max(w, d)
        return w

    width = red_deg_max()                     # 0 at the singleton state
    for a, b in merges:
        ia = next(i for i, pp in enumerate(parts) if a in pp)
        ib = next(i for i, pp in enumerate(parts) if b in pp)
        assert ia != ib, 'merge %r,%r: already in the same part' % (a, b)
        parts = ([pp for k, pp in enumerate(parts) if k not in (ia, ib)]
                 + [parts[ia] | parts[ib]])
        width = max(width, red_deg_max())
    assert len(parts) == 1, 'sequence left %d parts, not 1' % len(parts)
    return width


def is_cograph(g):
    """A graph is a cograph iff it has no induced P_4.  Exact, exhaustive over
    4-subsets; used only in the controls, never in the load-bearing path."""
    vs = sorted(g)
    for quad in combinations(vs, 4):
        sub = {u: frozenset(v for v in quad if v in g[u]) for u in quad}
        deg = sorted(len(sub[u]) for u in quad)
        ne = sum(deg) // 2
        if ne == 3 and deg == [1, 1, 2, 2]:
            return False                      # induced P_4
    return True


# ---------------------------------------------------------------------------
# a DETERMINISTIC greedy that PROPOSES a contraction sequence
# ---------------------------------------------------------------------------
# It is a heuristic and it is allowed to be: the width of whatever it returns is
# certified afterwards by sequence_width() above, which recomputes every
# trigraph relation from the ORIGINAL adjacency by the definition.  A verified
# sequence of width w is a proof that tww <= w, whatever produced it -- so the
# UPPER bounds below are sound even though the search that found them is not
# exhaustive.  No randomness and no floating point: ties break on vertex index,
# so the sequence printed in the transcript is reproducible byte for byte.
def greedy_contraction(g, shortlist=12):
    """-> list of merges (a, b) reducing g to one part, chosen greedily."""
    vs = sorted(g)
    n = len(vs)
    ix = {v: i for i, v in enumerate(vs)}
    adj = [0] * n
    for v in vs:
        for w in g[v]:
            adj[ix[v]] |= 1 << ix[w]
    pc = lambda x: bin(x).count('1')          # noqa: E731
    act = (1 << n) - 1                        # parts still alive
    BLACK, RED, rdeg = list(adj), [0] * n, [0] * n
    seq, alive = [], n
    while alive > 1:
        ids = [i for i in range(n) if (act >> i) & 1]
        cands = []
        for x, i in enumerate(ids):
            for j in ids[x + 1:]:
                # merging i and j makes RED exactly the parts they disagree on
                nr = ((RED[i] | RED[j] | (BLACK[i] ^ BLACK[j]))
                      & act & ~((1 << i) | (1 << j)))
                cands.append((pc(nr), pc(nr & ~RED[i] & ~RED[j]), i, j, nr))
        cands.sort()                          # ints only: a total, stable order
        scored = []
        for (cd, _new, i, j, nr) in cands[:shortlist]:
            mx = cd                           # red degree of the merged part
            for k in ids:                     # ... and of everybody else
                if k == i or k == j:
                    continue
                d = (rdeg[k] - ((RED[i] >> k) & 1) - ((RED[j] >> k) & 1)
                     + ((nr >> k) & 1))
                if d > mx:
                    mx = d
            scored.append((mx, cd, i, j, nr))
        scored.sort()
        _mx, _cd, i, j, nr = scored[0]
        for k in ids:                          # commit: j is absorbed into i
            if k == i or k == j:
                continue
            BLACK[k] &= ~((1 << i) | (1 << j))
            RED[k] &= ~((1 << i) | (1 << j))
            if (nr >> k) & 1:
                RED[k] |= 1 << i
            elif ((BLACK[i] >> k) & 1) and ((BLACK[j] >> k) & 1):
                BLACK[k] |= 1 << i
            rdeg[k] = pc(RED[k])
        BLACK[i] = (BLACK[i] & BLACK[j]) & act & ~nr & ~((1 << i) | (1 << j))
        RED[i], rdeg[i] = nr, pc(nr)
        act &= ~(1 << j)
        seq.append((vs[i], vs[j]))
        alive -= 1
    return seq


def verified_ub(g, shortlist=12):
    """(width, sequence) -- the greedy's proposal, re-certified from scratch."""
    seq = greedy_contraction(g, shortlist)
    return sequence_width(g, seq), seq


# ---------------------------------------------------------------------------
# small-graph enumeration up to isomorphism, and exact twin-width
# ---------------------------------------------------------------------------
def graphs_upto_iso(n):
    """Every graph on n vertices, one representative per isomorphism class.
    Brute force over all 2^C(n,2) labelled graphs and all n! relabellings; only
    ever called with n <= 5, where that is 1024 x 120."""
    pairs = list(combinations(range(n), 2))
    seen = {}
    for bits in range(1 << len(pairs)):
        a = {u: set() for u in range(n)}
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                a[i].add(j)
                a[j].add(i)
        g = {u: frozenset(s) for u, s in a.items()}
        best = None
        for p in permutations(range(n)):
            code = 0
            for k, (i, j) in enumerate(pairs):
                if p[j] in g[p[i]]:
                    code |= 1 << k
            if best is None or code < best:
                best = code
        if best not in seen:
            seen[best] = g
    return list(seen.values())


def tww_exhaustive(g):
    """EXACT twin-width, by minimising the width over ALL contraction sequences.
    Used only on the <= 5-vertex factors of the census, never on the product and
    never in the load-bearing path of the closure."""
    def red_deg_max(parts):
        w = 0
        for i, pp in enumerate(parts):
            d = 0
            for j, qq in enumerate(parts):
                if i == j:
                    continue
                cnt = sum(1 for p in pp for q in qq if q in g[p])
                if 0 < cnt < len(pp) * len(qq):
                    d += 1
            w = max(w, d)
        return w

    memo = {}

    def best_from(parts):
        if len(parts) <= 1:
            return 0
        if parts in memo:
            return memo[parts]
        best = None
        pl = list(parts)
        for a, b in combinations(range(len(pl)), 2):
            nxt = tuple(sorted(
                [pl[c] for c in range(len(pl)) if c not in (a, b)]
                + [pl[a] | pl[b]], key=sorted))
            w = max(red_deg_max(list(nxt)), best_from(nxt))
            if best is None or w < best:
                best = w
            if best == 0:
                break
        memo[parts] = best
        return best

    start = tuple(frozenset([v]) for v in sorted(g))
    return max(red_deg_max(list(start)), best_from(start))


# ---------------------------------------------------------------------------
# graph6
# ---------------------------------------------------------------------------
def decode_graph6(s):
    """(n, sorted edge list).  Handles n < 63 and the 63+2-byte form."""
    s = s.strip()
    if s.startswith('>>graph6<<'):
        s = s[len('>>graph6<<'):]
    data = [ord(c) - 63 for c in s]
    for x in data:
        if not 0 <= x < 64:
            raise ValueError('character out of graph6 range')
    if data[0] < 63:
        n, rest = data[0], data[1:]
    elif data[1] < 63:
        n = (data[2] << 12) | (data[3] << 6) | data[4]
        rest = data[5:]
    else:
        raise ValueError('graph6 n-form too large for this decoder')
    bits = []
    for x in rest:
        for k in range(5, -1, -1):
            bits.append((x >> k) & 1)
    need = n * (n - 1) // 2
    if len(bits) < need:
        raise ValueError('graph6 payload too short: %d bits for %d' % (len(bits), need))
    edges, idx = [], 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                edges.append((i, j))
            idx += 1
    return n, sorted(edges)


# ---------------------------------------------------------------------------
# the bound of Pettersson & Sylvester (thm:tensor)
# ---------------------------------------------------------------------------
def ps_bound(tww_g, delta_h, tww_h):
    return max(tww_g * delta_h + delta_h, tww_h) + delta_h


def delta(g):
    return max(len(g[u]) for u in g)


# ===========================================================================
print('=== re-verification of "A Tight Pair with Both Factors of Nonzero '
      'Twin-Width for the')
print('    Pettersson-Sylvester Tensor-Product Bound"  --  exact integer/set '
      'arithmetic ===')

# --- the objects exactly as the paper prints them -------------------------
G = circulant(8, [1, 2])          # C_8(1,2), the square of the 8-cycle
H = cycle(5)                      # C_5

# the graph6 string printed in the paper, as the concatenation of its three
# printed lines, of lengths 44, 44 and 43 -- transcribed here character for
# character from the paper's verbatim block
G6 = ('g?@LA_gcASsgSgISQc?Ch@Q_Ad?Ad?HQ??Ch?Ad??IS?'
      '?Sg?AS_??AS_?Ad???Sg??@Q_??QcAO??dL??@Q`O??S'
      'gI??AdC_?@IOHQ??ChhO??SgSg??ISDI??AdCh??AS_')

# the two contraction sequences printed in the paper
SEQ_G = [(0, 1), (2, 3), (4, 5), (6, 7), (0, 4), (2, 6), (0, 2)]
SEQ_H = [(0, 2), (3, 4), (0, 1), (0, 3)]

print('\n-- STEP 0: the two factors, as printed')
ck('G = C_8(1,2) is 4-regular on 8 vertices with 16 edges',
   len(G) == 8 and all(len(G[u]) == 4 for u in G) and len(edges_of(G)) == 16,
   'N(0) = %s' % sorted(G[0]))
ck('H = C_5 is 2-regular on 5 vertices with 5 edges',
   len(H) == 5 and all(len(H[i]) == 2 for i in H) and len(edges_of(H)) == 5)
ck('C_8(1,2) is the SQUARE of the 8-cycle (u ~ v iff dist_{C_8}(u,v) <= 2)',
   all((v in G[u]) == (0 < min((u - v) % 8, (v - u) % 8) <= 2)
       for u in range(8) for v in range(8) if u != v))
c34 = circulant(8, [3, 4])
iso = all(((3 * u) % 8 in c34[(3 * v) % 8]) == (u in circulant(8, [1, 4])[v])
          for u in range(8) for v in range(8) if u != v)
ck('complement of C_8(1,2) is C_8(3,4), which is the Wagner graph C_8(1,4) via '
   'x -> 3x on Z_8',
   complement(G) == c34 and iso)

print('\n-- STEP 1: tww(C_5) = 2 and tww(C_8(1,2)) = 2, both NON-ZERO')
note('Ahn-Hendrey-Kim-Oum Lemma 3.1 bounds tww below by the MINIMUM of b over '
     'ALL pairs,')
note('so the FULL profile is computed; one exhibited small pair would prove '
     'nothing.')

profH_b = {d: b_val(H, 0, d % 5) for d in (1, 2)}
profH_s = {d: s_val(H, 0, d % 5) for d in (1, 2)}
note('C_5      b by difference d: %s      s by difference d: %s   '
     '(d=1 adjacent, d=2 not)' % (profH_b, profH_s))
bH, argH, cntH = b_min(H)
ck('b_min(C_5) = 2 over all 10 pairs', bH == 2 and cntH == 10,
   'argmin %s' % (argH,))

profG_b = {d: b_val(G, 0, d) for d in (1, 2, 3, 4)}
profG_s = {d: s_val(G, 0, d) for d in (1, 2, 3, 4)}
note('C_8(1,2) b by difference d: %s   s by difference d: %s'
     % (profG_b, profG_s))
ck('the b-profile of C_8(1,2) is {1:2, 2:4, 3:4, 4:4}, as printed in the paper',
   profG_b == {1: 2, 2: 4, 3: 4, 4: 4})
ck('the s-profile of C_8(1,2) is {1:4, 2:6, 3:4, 4:4}, as printed in the paper',
   profG_s == {1: 4, 2: 6, 3: 4, 4: 4})
bG, argG, cntG = b_min(G)
ck('b_min(C_8(1,2)) = 2 over all 28 pairs, attained at the ADJACENT pair (0,1)',
   bG == 2 and cntG == 28 and argG == (0, 1) and 1 in G[0], 'argmin %s' % (argG,))
ck('s_min(C_8(1,2)) = 4 = tww(G) + 2, the structural condition the argument uses',
   s_min(G) == 4)

wH = sequence_width(H, SEQ_H)
ck('the paper\'s C_5 sequence (0,2)(3,4)({0,2},1)(rest) has width 2, so '
   'tww(C_5) <= 2', wH == 2, 'width %d' % wH)
wG = sequence_width(G, SEQ_G)
ck('the paper\'s C_8(1,2) sequence (01)(23)(45)(67)(A,C)(B,D)(rest) has width 2, '
   'so tww <= 2', wG == 2, 'width %d' % wG)
ck('therefore tww(C_5) = 2 EXACTLY (2 <= tww <= 2); no exhaustive search is '
   'load-bearing', bH == 2 and wH == 2)
ck('therefore tww(C_8(1,2)) = 2 EXACTLY (2 <= tww <= 2); no exhaustive search '
   'is load-bearing', bG == 2 and wG == 2)
ck('both factors have NON-ZERO twin-width, which is exactly what the open '
   'question demands', bH > 0 and bG > 0)

print('\n-- STEP 2: the value of the bound at the ORDERED pair (G,H)')
T_GH = ps_bound(2, delta(H), 2)
T_HG = ps_bound(2, delta(G), 2)
note('T(G,H) = max{2*%d + %d, 2} + %d = %d       T(H,G) = max{2*%d + %d, 2} + %d '
     '= %d' % (delta(H), delta(H), delta(H), T_GH,
               delta(G), delta(G), delta(G), T_HG))
ck('the bound at (G,H) = (C_8(1,2), C_5) is 8', T_GH == 8)
ck('the ACTIVE branch of the max is tww(G)*Delta(H)+Delta(H) = 6 > tww(H) = 2, '
   'so the term attained carries tww(G) != 0',
   2 * delta(H) + delta(H) == 6 > 2)
ck('the SWAPPED order gives 16, and the bound is therefore ASYMMETRIC', T_HG == 16)

print('\n-- STEP 3: the lower bound tww(G x H) >= 8, over ALL C(40,2) = 780 pairs')
P = tensor(G, H)
ck('G x H has 40 vertices, 160 edges and is 8-regular',
   len(P) == 40 and len(edges_of(P)) == 160
   and all(len(P[k]) == 8 for k in P))
ck('the neighbourhood identity N_P((u,i)) = N_G(u) x N_H(i) holds at every vertex',
   all(P[(u, i)] == frozenset((v, j) for v in G[u] for j in H[i])
       for u in G for i in H))

vs = sorted(P)
caseA, caseB, caseC = [], [], []
for a, c in combinations(range(len(vs)), 2):
    (u, i), (v, j) = vs[a], vs[c]
    val = b_val(P, vs[a], vs[c])
    (caseA if u == v else caseB if i == j else caseC).append((val, (vs[a], vs[c])))
minA, minB, minC = min(caseA), min(caseB), min(caseC)
note('pairs:  A (u=v) %d   B (i=j) %d   C (both differ) %d   total %d'
     % (len(caseA), len(caseB), len(caseC),
        len(caseA) + len(caseB) + len(caseC)))
note('per-case minima:  A = %d at %s   B = %d at %s   C = %d at %s'
     % (minA[0], minA[1], minB[0], minB[1], minC[0], minC[1]))
ck('all 780 pairs enumerated, and the three cases partition them',
   len(caseA) + len(caseB) + len(caseC) == 780 == 40 * 39 // 2
   and len(caseA) == 80 and len(caseB) == 140 and len(caseC) == 560)
ck('case (A) u = v: minimum is 8, attained at a distance-2 pair of C_5',
   minA[0] == 8)
ck('case (B) i = j: minimum is 8, attained at an adjacent pair of C_8(1,2)',
   minB[0] == 8)
ck('case (C) both differ: minimum is 12, so case (C) is never the minimum',
   minC[0] == 12)
ck('every case-(C) pair with uv in E(G) and ij in E(H) has b >= 14, as the '
   'proof claims',
   all(val >= 14 for val, ((u, i), (v, j)) in caseC
       if v in G[u] and j in H[i]))
bP, argP, cntP = b_min(P)
ck('b_min(G x H) = 8, so tww(G x H) >= 8 by Ahn et al. Lemma 3.1',
   bP == 8 and cntP == 780, 'argmin %s' % (argP,))
ck('*** CLOSURE: tww(C_8(1,2) x C_5) = 8 = the bound, both factors of '
   'twin-width 2 ***', bP == T_GH == 8)

print('\n-- STEP 3b: the UPPER bound tww(G x H) <= 8 WITHOUT the theorem under '
      'test')
note('thm:tensor already gives <= 8, but a tightness claim whose upper bound is '
     'the very')
note('theorem being certified is worth making independent of it.  A width-8 '
     'contraction')
note('sequence of the product, re-certified from the original adjacency, does '
     'that.')
ubP, seqP = verified_ub(P)
ck('a contraction sequence of the 40-vertex product is built here and '
   're-certified at width 8, so tww(G x H) <= 8 with NO appeal to thm:tensor',
   ubP == 8 and len(seqP) == 39, 'width %d, %d merges' % (ubP, len(seqP)))
ck('*** tww(C_8(1,2) x C_5) = 8 EXACTLY and NON-CIRCULARLY: 8 <= tww by Ahn et '
   'al. Lemma 3.1, tww <= 8 by the sequence just verified', bP == ubP == 8)
ck('thm:tensor evaluated at (G,H) AGREES with that independent value, so the '
   'theorem is a consistency check here and not an input', T_GH == ubP)
note('the sequence, as (part containing a, part containing b) merges in order:')
for _r in range(0, len(seqP), 5):
    note('   ' + '  '.join('%s+%s' % (a, b) for a, b in seqP[_r:_r + 5]))

print('\n-- STEP 4: the witness, and label-for-label agreement with the paper')
n6, e6 = decode_graph6(G6)
lab = sorted(tuple(sorted((5 * u + i, 5 * v + j)))
             for (u, i), (v, j) in edges_of(P))
note('graph6 string: %d characters, n = %d, %d edges' % (len(G6), n6, len(e6)))
ck('the paper\'s graph6 string has 131 characters (44 + 44 + 43) and decodes '
   'to 40 vertices', len(G6) == 131 and n6 == 40)
ck('the decoded graph is G x H LABEL-FOR-LABEL under k = 5u + i (not merely '
   'isomorphic)', e6 == lab, '%d edges both ways' % len(e6))
ck('every decoded edge {5u+i, 5v+j} has u-v = +-1,+-2 mod 8 and i-j = +-1 mod 5',
   all(min((a // 5 - b // 5) % 8, (b // 5 - a // 5) % 8) in (1, 2)
       and min((a % 5 - b % 5) % 5, (b % 5 - a % 5) % 5) == 1
       for a, b in e6))

# the one structural fact the note's prior-art paragraph turns on
_mu = sorted({len(P[a] & P[b]) for a, b in combinations(sorted(P), 2)
              if b not in P[a]})
_lam = sorted({len(P[a] & P[b]) for a, b in combinations(sorted(P), 2)
               if b in P[a]})
note('common neighbours in P: non-adjacent pairs %s, adjacent pairs %s'
     % (_mu, _lam))
ck('P is NOT strongly regular -- non-adjacent pairs have 0, 1, 2 or 4 common '
   'neighbours, so mu is not constant (this is what separates P from the '
   'strongly regular families of the nearest prior work)',
   _mu == [0, 1, 2, 4] and _lam == [0])

print('\n-- STEP 5: the infinite family H = C_n for every n >= 5')
fam = 0
for n in range(4, 13):
    Cn = cycle(n)
    bn = b_min(Cn)[0]
    wn = sequence_width(Cn, [(0, 1)] + [(0, k) for k in range(2, n)])
    Pn = tensor(G, Cn)
    bpn = b_min(Pn)[0]
    Tn = ps_bound(2, delta(Cn), 2)
    ubn = verified_ub(Pn)[0] if n >= 5 else None
    note('n=%-2d  b_min(C_n)=%d  round-the-cycle width=%d   '
         'b_min(C_8(1,2) x C_n)=%d   verified UB=%s   bound=%d   %s'
         % (n, bn, wn, bpn, ubn, Tn, 'TIGHT' if bpn == Tn and bn == 2 else
            'EXCLUDED (C_4 is a cograph, tww 0)' if n == 4 else 'NOT TIGHT'))
    if n == 4:
        ck('n=4: C_4 is a cograph with b_min = 0, so it is correctly EXCLUDED '
           'from the family', bn == 0 and is_cograph(Cn))
    else:
        ck('n=%d: tww(C_n) = 2 (b_min = 2 and a width-2 sequence) and '
           'tww(C_8(1,2) x C_n) = %d = the bound, the upper half again by a '
           'sequence verified here rather than by thm:tensor' % (n, Tn),
           bn == 2 and wn == 2 and bpn == Tn == ubn == 8)
        fam += 1
ck('the family is confirmed by machine for all %d values n = 5..12, both bounds '
   'independently' % fam, fam == 8)

print('\n-- STEP 6: the ONE completeness statement made anywhere here -- the 121'
      ' ordered')
print('           pairs of non-cographs on at most 5 vertices, every one NOT '
      'tight')
note('This bears on NOTHING in the closure above: the witness G = C_8(1,2) lives '
     'on 8')
note('vertices, outside the band. It is here only so that the paper\'s "what is '
     'not')
note('settled" list can point at a checked fact instead of at an assertion. A '
     'pair is')
note('decided NOT TIGHT by one of two SOUND routes and by nothing else:')
note('  (i)  T(H,G) < T(G,H): then tww(G x H) = tww(H x G) <= T(H,G) < T(G,H) by '
     'thm:tensor;')
note('  (ii) a contraction sequence of G x H re-certified here at width < '
     'T(G,H).')
note('Neither route can call a tight pair non-tight. A pair the two routes miss '
     'is')
note('reported as a SURVIVOR, not as non-tight.')

pool, iso_counts = [], {}
for _n in (4, 5):
    _reps = graphs_upto_iso(_n)
    _nc = [a for a in _reps if not is_cograph(a)]
    iso_counts[_n] = (len(_reps), len(_reps) - len(_nc), len(_nc))
    pool += _nc
note('graphs up to isomorphism / cographs / non-cographs: %s' % iso_counts)
ck('the enumeration reproduces OEIS A000088 (graphs on n unlabelled nodes) at '
   'n = 4, 5: 11 and 34',
   iso_counts[4][0] == 11 and iso_counts[5][0] == 34)
ck('and the published cograph counts 10 and 24, leaving 1 and 10 non-cographs, '
   'so the pool is 11 graphs and the census is 11 x 11 = 121 ORDERED pairs',
   iso_counts[4][1] == 10 and iso_counts[5][1] == 24
   and iso_counts[4][2] == 1 and iso_counts[5][2] == 10 and len(pool) == 11)
POOLTWW = [tww_exhaustive(a) for a in pool]
note('exact twin-width of the 11 pool graphs (exhaustive over all sequences): %s'
     % sorted(POOLTWW))
ck('every pool graph has twin-width 1 or 2, and the unique 4-vertex one is P_4 '
   'with tww 1',
   all(t in (1, 2) for t in POOLTWW)
   and [t for a, t in zip(pool, POOLTWW) if len(a) == 4] == [1])


def decide_tight(g, tg, h, th):
    """-> 'NOT_TIGHT_SWAP' | 'NOT_TIGHT' | 'SURVIVOR', with T and the evidence."""
    T = ps_bound(tg, delta(h), th)
    if ps_bound(th, delta(g), tg) < T:
        return 'NOT_TIGHT_SWAP', T, ps_bound(th, delta(g), tg)
    w = verified_ub(tensor(g, h))[0]
    return ('NOT_TIGHT' if w < T else 'SURVIVOR'), T, w


tally, survivors = {}, []
for _a, _ta in zip(pool, POOLTWW):
    for _b, _tb in zip(pool, POOLTWW):
        _st, _T, _w = decide_tight(_a, _ta, _b, _tb)
        tally[_st] = tally.get(_st, 0) + 1
        if _st == 'SURVIVOR':
            survivors.append((sorted(edges_of(_a)), sorted(edges_of(_b)), _T, _w))
note('census tally: %s   (products of order up to 5*5 = 25)'
     % sorted(tally.items()))
ck('all 121 ordered pairs are enumerated and the tally adds up',
   sum(tally.values()) == 121 == len(pool) ** 2)
ck('37 pairs fall to route (i), the asymmetry of the bound',
   tally.get('NOT_TIGHT_SWAP', 0) == 37)
ck('the other 84 fall to route (ii), a sequence verified here below the bound',
   tally.get('NOT_TIGHT', 0) == 84)
ck('*** 121 of 121 NOT TIGHT, zero survivors: no ordered pair of non-cographs on '
   'at most 5 vertices attains the bound ***',
   not survivors and tally.get('NOT_TIGHT', 0) + tally.get('NOT_TIGHT_SWAP', 0)
   == 121)
_wst, _wT, _ww = decide_tight(G, 2, H, 2)
ck('CONTROL ON THE DECIDER: run on the paper\'s own pair it does NOT say '
   'NOT_TIGHT -- T(H,G) = 16 > 8 blocks route (i) and the verified width is 8, '
   'not below 8, so route (ii) is blocked too',
   _wst == 'SURVIVOR' and _wT == 8 and _ww == 8)
note('so the 121 NOT_TIGHT verdicts are not the output of a decider that refuses '
     'everything: the same decider, on the same code path, declines to refuse '
     'the tight pair.')

print('\n-- CONTROLS, BOTH POLARITIES, EACH AGAINST AN INDEPENDENTLY KNOWN VALUE')
for k in (3, 4):
    Kk = complete(k)
    Pk = tensor(Kk, Kk)
    ck('FORCED POSITIVE (K_%d, K_%d): the clique corollary gives tww = 2(%d-1) = '
       '%d, the bound gives %d, b_min measures %d -- all three agree'
       % (k, k, k, 2 * (k - 1), ps_bound(0, delta(Kk), 0), b_min(Pk)[0]),
       2 * (k - 1) == ps_bound(0, delta(Kk), 0) == b_min(Pk)[0])
K2 = complete(2)
P22 = tensor(K2, K2)
ck('FORCED SILENT (K_2, K_2): the bound gives 2 but b_min = 0 and the product '
   'is a cograph (tww 0), so no tightness may be inferred',
   ps_bound(0, 1, 0) == 2 and b_min(P22)[0] == 0 and is_cograph(P22))
P4 = path(4)
ck('anti-control on the smallest non-cograph: b_min(P_4) = 1 and the sequence '
   '(0,2)(1,3)(rest) has width 1, so tww(P_4) = 1',
   b_min(P4)[0] == 1 and sequence_width(P4, [(0, 2), (1, 3), (0, 1)]) == 1)
P44 = tensor(P4, P4)
ck('FORCED SILENT (P_4, P_4): the bound gives 6 but b_min(P_4 x P_4) = 1 < 6, '
   'so Lemma 3.1 CANNOT certify tightness there',
   ps_bound(1, delta(P4), 1) == 6 and b_min(P44)[0] == 1)
ck('the published equivalence "tww = 0 iff cograph" spot-checked: C_4 and K_5 '
   'are cographs, P_4 and C_5 are not',
   is_cograph(cycle(4)) and is_cograph(complete(5))
   and not is_cograph(P4) and not is_cograph(H))
try:
    sequence_width(H, [(0, 2), (3, 4)])          # leaves 3 parts, must be refused
    _sim_ok = False
except AssertionError:
    _sim_ok = True
ck('sanity on the trigraph simulator: it REFUSES a sequence that does not '
   'reduce to a single part', _sim_ok)

# ---------------------------------------------------------------------------
print()
print('NOT RE-RUN: no EXHAUSTIVE search for a contraction sequence of the '
      '40-vertex product. STEP 3b builds one greedily and certifies it at width '
      '8, which proves tww <= 8; it does not prove that no width-7 sequence '
      'exists. That direction is not needed and is not claimed here -- it is '
      'Ahn et al. Lemma 3.1, via b_min = 8, that rules width 7 out.')
print('NOT RE-RUN: no MINIMALITY. Nothing here bears on whether 40 is the least '
      'order of a tight product, nor on whether a tight pair with non-regular H '
      'exists. The only completeness statement made anywhere is STEP 6, the 121 '
      'ordered pairs of non-cographs on at most 5 vertices; the witness lives on '
      '8 vertices, outside that band. Two larger censuses (all ordered pairs of '
      'non-cographs on at most 6 vertices, and a minimal-tight-pair sweep) were '
      'dispatched in the discovery run and never delivered; those channels are '
      'UNREAD, not negative, and nothing above leans on them.')
print('NOT RE-RUN: the family claim for general n > 12 is the hand argument '
      'printed in the paper. Machine coverage here is n = 5..12 only, for both '
      'bounds.')
print('NOT RE-RUN: no prior-art search. OpenAlex citation listing, Semantic '
      'Scholar keyword search and MathSciNet were all UNREAD in the '
      'accompanying prior-art pass; see the paper\'s bibliography for what was '
      'read.')
print()
if _STATE['fail']:
    print('VERDICT: %d CHECKS FAILED of %d'
          % (_STATE['fail'], _STATE['fail'] + _STATE['pass']))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _STATE['pass'])
sys.exit(0)
