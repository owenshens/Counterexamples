#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- re-derivation of every computational claim of

    "Explicit OF(6,4) and OF(8,4): three of the four even exceptions at g = 4 are closed"

Python 3.9+, STANDARD LIBRARY ONLY (this file imports nothing but `sys` and `itertools`).
No external data file, no network, no solver, no randomness, no seed.
All arithmetic is exact integer arithmetic on ints, tuples and frozensets; no float
decides any check.

THE ONLY INPUTS ARE THE TWO BASE FACTORS PRINTED IN SECTION 2 OF THE PAPER, transcribed
below as the string constants PAPER_M6 and PAPER_M8 in exactly the text the paper prints.
The program parses them, rebuilds the ambient graph and the group from the paper's own
description, and re-derives every number the paper states.  It reads nothing from the
pipeline that produced the result: no row, no database, no earlier artifact.

Layout:
    A  the printed objects parse and are perfect matchings
    B  the ambient graphs, the groups, freeness, the orbit counts
    C  the two hypotheses of the reduction lemma (transversal, distinct part-pairs)
    D  the factorizations by direct expansion, independent of the lemma
    E  controls on the decider, in both polarities
    F  OF(12,4) by the source's doubling lemma, and the g = 4 spectrum
"""

import itertools
import sys

# ---------------------------------------------------------------------------
# THE OBJECTS, EXACTLY AS PRINTED IN SECTION 2 OF THE PAPER
# One edge per line, "part value part value"; parts are integers or the literal inf.
# ---------------------------------------------------------------------------

PAPER_M6 = """
0 0  1 1
1 2  2 2
3 3  4 1
4 3  0 2
0 1  2 1
1 3  3 1
2 3  4 2
3 2  0 3
inf 0  1 0
inf 1  2 0
inf 2  3 0
inf 3  4 0
"""

PAPER_M8 = """
0 0  1 0
0 1  2 1
0 2  4 1
0 3  5 0
1 1  2 2
1 2  inf 2
1 3  5 1
2 0  3 2
2 3  6 0
3 0  4 3
3 1  5 2
3 3  6 3
4 0  6 2
4 2  inf 1
5 3  inf 0
6 1  inf 3
"""

# ---------------------------------------------------------------------------
# check harness
# ---------------------------------------------------------------------------

_PASSED = []
_FAILED = []


def check(name, ok, detail=''):
    if ok:
        _PASSED.append(name)
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILED.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def section(title):
    print('')
    print('=== %s ===' % title)


# ---------------------------------------------------------------------------
# the ambient graph K_{n x g} on parts Z_p u {inf}, n = p+1, and the group Z_p x Z_g
# ---------------------------------------------------------------------------

def parse_factor(text):
    """The paper's printed block -> a list of edges, each a frozenset of two vertices.
    A vertex is (part, value); a part is an int or the string 'inf'."""
    edges = []
    for raw in text.strip().splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        tok = line.split()
        if len(tok) != 4:
            raise ValueError('not four tokens: %r' % raw)
        u = (tok[0] if tok[0] == 'inf' else int(tok[0]), int(tok[1]))
        v = (tok[2] if tok[2] == 'inf' else int(tok[2]), int(tok[3]))
        edges.append(frozenset((u, v)))
    return edges


def parts_of(p):
    return list(range(p)) + ['inf']


def vertices(p, g):
    return set((u, a) for u in parts_of(p) for a in range(g))


def edges_of(p, g):
    """E(K_{n x g}) with n = p+1: every pair of vertices in DIFFERENT parts."""
    E = set()
    V = sorted(vertices(p, g), key=lambda v: (str(v[0]), v[1]))
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            if V[i][0] != V[j][0]:
                E.add(frozenset((V[i], V[j])))
    return E


def group(p, g):
    return [(t, s) for t in range(p) for s in range(g)]


def act_v(v, t, s, p, g):
    x, a = v
    return (x if x == 'inf' else (x + t) % p, (a + s) % g)


def act_e(e, t, s, p, g):
    u, v = tuple(e)
    return frozenset((act_v(u, t, s, p, g), act_v(v, t, s, p, g)))


def part_pair(e):
    u, v = tuple(e)
    return frozenset((u[0], v[0]))


def invariant(e, p, g):
    """A G-invariant of an edge.  For a finite-finite edge take the part-difference class
    d in 1..(p-1)/2 with the orientation that realises it, and the value shift along that
    orientation; for an edge at infinity take the value shift from infinity.  Both are
    untouched by (t,s), which shifts the two parts equally and the two values equally."""
    u, v = tuple(e)
    if u[0] == 'inf' or v[0] == 'inf':
        (_, a), (_, b) = (u, v) if u[0] == 'inf' else (v, u)
        return ('inf', None, (b - a) % g)
    (x, a), (y, b) = u, v
    d0 = (y - x) % p
    if d0 <= (p - 1) // 2:
        return ('fin', d0, (b - a) % g)
    return ('fin', p - d0, (a - b) % g)


def orbits(p, g):
    """Explicit orbit enumeration under G.  -> (list of orbits as sets of edges)."""
    E = edges_of(p, g)
    seen = set()
    out = []
    for e in sorted(E, key=lambda e: sorted((str(x[0]), x[1]) for x in e)):
        if e in seen:
            continue
        orb = set(act_e(e, t, s, p, g) for (t, s) in group(p, g))
        seen |= orb
        out.append(orb)
    return out


def expand(M, p, g):
    """The |G| translates of the base factor M."""
    return [[act_e(e, t, s, p, g) for e in M] for (t, s) in group(p, g)]


def codeword(e, p, g):
    """The weight-2 word of length n = p+1 that the edge is, in the source's coding normal
    form: coordinate u carries a+1 if (u,a) is an end of the edge, and 0 otherwise."""
    w = []
    ends = dict(tuple(e))
    for u in parts_of(p):
        w.append(ends[u] + 1 if u in ends else 0)
    return tuple(w)


def hamming(w1, w2):
    return sum(1 for a, b in zip(w1, w2) if a != b)


def factor_report(F, p, g):
    """Every OF axiom, re-derived from the factor as given."""
    V = vertices(p, g)
    ends = [x for e in F for x in e]
    pps = [part_pair(e) for e in F]
    words = [codeword(e, p, g) for e in F]
    mind = min((hamming(w1, w2) for w1, w2 in itertools.combinations(words, 2)),
               default=None)
    return {
        'size': len(F),
        'is_perfect_matching': len(ends) == len(set(ends)) and set(ends) == V,
        'part_pairs_distinct': len(pps) == len(set(pps)),
        'min_distance': mind,
        'part_pairs': set(pps),
    }


def latin_ok(sq):
    k = len(sq)
    for i in range(k):
        if sorted(sq[i]) != list(range(k)):
            return False
        if sorted(sq[r][i] for r in range(k)) != list(range(k)):
            return False
    return True


def enumerate_pm_kng(p, g, want_parallel_free=True):
    """Exhaustive enumeration of the perfect matchings of K_{n x g}, n = p+1.
    -> (total, number that are parallel-edge-free)."""
    V = sorted(vertices(p, g), key=lambda v: (str(v[0]), v[1]))
    nv = len(V)
    tot = [0]
    pf = [0]

    def rec(used, pairs):
        if len(used) == nv:
            tot[0] += 1
            if want_parallel_free and len(pairs) == nv // 2:
                pf[0] += 1
            return
        i = 0
        while V[i] in used:
            i += 1
        u = V[i]
        for j in range(nv):
            v = V[j]
            if v in used or v == u or v[0] == u[0]:
                continue
            rec(used | {u, v}, pairs | {frozenset((u[0], v[0]))})

    rec(frozenset(), frozenset())
    return tot[0], pf[0]


def double_factorial_odd(m):
    """m!! for odd m."""
    r = 1
    k = m
    while k > 1:
        r *= k
        k -= 2
    return r


def poly_pow(coeffs, e):
    out = [1]
    for _ in range(e):
        new = [0] * (len(out) + len(coeffs) - 1)
        for i, a in enumerate(out):
            for j, b in enumerate(coeffs):
                new[i + j] += a * b
        out = new
    return out


# ===========================================================================
def main():
    print('verify.py -- Explicit OF(6,4) and OF(8,4): three of the four even')
    print('exceptions at g = 4 are closed.  Exact integer arithmetic throughout.')

    G = 4
    M6 = parse_factor(PAPER_M6)
    M8 = parse_factor(PAPER_M8)

    # -------------------------------------------------------------- A
    section('A. THE OBJECTS PRINTED IN SECTION 2 OF THE PAPER')

    check('A1-M6-parses',
          len(M6) == 12 and len(set(M6)) == 12 and all(len(e) == 2 for e in M6),
          '12 lines, 12 distinct undirected edges, no loop')
    ends6 = [x for e in M6 for x in e]
    check('A2-M6-is-a-perfect-matching-of-K_6x4',
          len(ends6) == 24 and len(set(ends6)) == 24 and set(ends6) == vertices(5, G),
          '24 ends, all distinct, covering all 6*4 = 24 vertices exactly once')
    check('A3-M6-every-edge-crosses-parts',
          all(len(set(x[0] for x in e)) == 2 for e in M6),
          'no edge inside a part, so M6 lies in E(K_{6x4})')
    check('A4-M6-each-part-shows-every-value-once',
          all(sorted(a for (u, a) in ends6 if u == P) == [0, 1, 2, 3] for P in parts_of(5)),
          'each of the 6 parts contributes the values {0,1,2,3}')

    check('A5-M8-parses',
          len(M8) == 16 and len(set(M8)) == 16 and all(len(e) == 2 for e in M8),
          '16 lines, 16 distinct undirected edges, no loop')
    ends8 = [x for e in M8 for x in e]
    check('A6-M8-is-a-perfect-matching-of-K_8x4',
          len(ends8) == 32 and len(set(ends8)) == 32 and set(ends8) == vertices(7, G),
          '32 ends, all distinct, covering all 8*4 = 32 vertices exactly once')
    check('A7-M8-every-edge-crosses-parts',
          all(len(set(x[0] for x in e)) == 2 for e in M8),
          'no edge inside a part, so M8 lies in E(K_{8x4})')
    check('A8-M8-each-part-shows-every-value-once',
          all(sorted(a for (u, a) in ends8 if u == P) == [0, 1, 2, 3] for P in parts_of(7)),
          'each of the 8 parts contributes the values {0,1,2,3}')

    # -------------------------------------------------------------- B
    section('B. THE AMBIENT GRAPHS, THE GROUPS, FREENESS, THE ORBIT COUNTS')

    E6 = edges_of(5, G)
    E8 = edges_of(7, G)
    c = lambda m: m * (m - 1) // 2
    check('B1-ambient-K_6x4',
          len(vertices(5, G)) == 24 and len(E6) == 240 == c(24) - 6 * c(4),
          '|V| = 24, |E| = 240 = C(24,2) - 6*C(4,2) = 276 - 36')
    check('B2-ambient-K_8x4',
          len(vertices(7, G)) == 32 and len(E8) == 448 == c(32) - 8 * c(4),
          '|V| = 32, |E| = 448 = C(32,2) - 8*C(4,2) = 496 - 48')

    for p, tag, nE, nG, nOrb in ((5, 'p=5', 240, 20, 12), (7, 'p=7', 448, 28, 16)):
        Ep = E6 if p == 5 else E8
        Gp = group(p, G)
        # a group of permutations of E, acting so that parts go to parts
        closed = all(act_e(e, t, s, p, G) in Ep for e in Ep for (t, s) in Gp)
        idfix = all(act_e(e, 0, 0, p, G) == e for e in Ep)
        check('B3-%s-G-acts-on-E-and-preserves-the-parts' % tag,
              closed and idfix and len(Gp) == nG,
              '|G| = %d = %d*%d; every element maps E onto E and each part onto a part'
              % (nG, p, G))
        free = True
        witness = None
        for e in Ep:
            for (t, s) in Gp:
                if (t, s) != (0, 0) and act_e(e, t, s, p, G) == e:
                    free = False
                    witness = (e, t, s)
        check('B4-%s-G-acts-FREELY-on-E' % tag, free,
              'no non-identity element of G fixes an edge; %d edge x %d element tests'
              % (nE, nG - 1) if free else 'counterexample %r' % (witness,))
        orbs = orbits(p, G)
        check('B5-%s-orbit-count-and-orbit-sizes' % tag,
              len(orbs) == nOrb and set(len(o) for o in orbs) == {nG}
              and nOrb * nG == nE,
              '%d orbits, every orbit of size %d, %d*%d = %d = |E|'
              % (nOrb, nG, nOrb, nG, nE))
        inv_classes = {}
        for e in Ep:
            inv_classes.setdefault(invariant(e, p, G), set()).add(e)
        as_partition = lambda fams: set(frozenset(fam) for fam in fams)
        check('B6-%s-the-printed-invariant-is-a-COMPLETE-orbit-invariant' % tag,
              as_partition(orbs) == as_partition(inv_classes.values()),
              'the %d invariant classes are exactly the %d G-orbits, so "one edge per orbit" '
              'can be read off the invariants' % (len(inv_classes), nOrb))
        check('B7-%s-a-one-factor-has-exactly-one-edge-per-orbit' % tag,
              G * (p + 1) // 2 == nOrb,
              'factor size gn/2 = %d equals the number of orbits %d, so a base factor MUST '
              'be an orbit transversal' % (G * (p + 1) // 2, nOrb))

    # -------------------------------------------------------------- C
    section('C. THE TWO HYPOTHESES OF THE REDUCTION LEMMA')

    inv6 = [invariant(e, 5, G) for e in M6]
    check('C1-M6-is-an-orbit-transversal',
          len(set(inv6)) == 12,
          'the 12 edges of M6 lie in 12 distinct orbits, one per orbit')
    pp6 = [part_pair(e) for e in M6]
    check('C2-M6-part-pairs-pairwise-distinct',
          len(set(pp6)) == 12 and len(pp6) == 12,
          '12 distinct part-pairs out of C(6,2) = 15')
    missed6 = set(frozenset(q) for q in itertools.combinations(parts_of(5), 2)) - set(pp6)
    m6names = sorted(tuple(sorted(map(str, q))) for q in missed6)
    deg6 = {}
    for q in missed6:
        for u in q:
            deg6[u] = deg6.get(u, 0) + 1
    check('C3-M6-missed-part-pairs-form-a-perfect-matching-of-K_6',
          m6names == [('0', 'inf'), ('1', '4'), ('2', '3')]
          and all(deg6.get(u, 0) == 1 for u in parts_of(5)),
          'missed = {inf,0},{1,4},{2,3}; every part in exactly one of them')
    check('C4-M6-orbit-invariants-in-the-printed-order',
          tuple(i[2] for i in inv6[0:4]) == (1, 0, 2, 3)
          and tuple(i[2] for i in inv6[4:8]) == (0, 2, 3, 1)
          and tuple(i[2] for i in inv6[8:12]) == (0, 3, 2, 1)
          and all(i[1] == 1 for i in inv6[0:4])
          and all(i[1] == 2 for i in inv6[4:8])
          and all(i[0] == 'inf' for i in inv6[8:12]),
          'pentagon d=1 gives (1,0,2,3), pentagram d=2 gives (0,2,3,1), infinity gives '
          '(0,3,2,1) -- each a permutation of Z_4')

    inv8 = [invariant(e, 7, G) for e in M8]
    check('C5-M8-is-an-orbit-transversal',
          len(set(inv8)) == 16,
          'the 16 edges of M8 lie in 16 distinct orbits, one per orbit')
    pp8 = [part_pair(e) for e in M8]
    check('C6-M8-part-pairs-pairwise-distinct',
          len(set(pp8)) == 16 and len(pp8) == 16,
          '16 distinct part-pairs out of C(8,2) = 28')
    missed8 = set(frozenset(q) for q in itertools.combinations(parts_of(7), 2)) - set(pp8)
    deg8 = {}
    for q in missed8:
        for u in q:
            deg8[u] = deg8.get(u, 0) + 1
    check('C7-M8-misses-12-part-pairs-forming-a-3-regular-complement',
          len(missed8) == 12 and set(deg8.get(u, 0) for u in parts_of(7)) == {3},
          '28 - 16 = 12 missed part-pairs, every part in exactly 3 of them, 12*2 = 8*3')
    cls8 = {}
    for i in inv8:
        cls8.setdefault((i[0], i[1]), []).append(i[2])
    check('C8-M8-orbit-invariants-cover-every-class-once',
          sorted(cls8) == [('fin', 1), ('fin', 2), ('fin', 3), ('inf', None)]
          and all(sorted(v) == [0, 1, 2, 3] for v in cls8.values()),
          'difference classes d = 1,2,3 and infinity each carry the shifts 0,1,2,3 once')

    # -------------------------------------------------------------- D
    section('D. THE FACTORIZATIONS BY DIRECT EXPANSION')

    for p, tag, M, nfac, fsz, nE in ((5, 'OF(6,4)', M6, 20, 12, 240),
                                     (7, 'OF(8,4)', M8, 28, 16, 448)):
        Ep = E6 if p == 5 else E8
        facs = expand(M, p, G)
        reps = [factor_report(F, p, G) for F in facs]
        check('D1-%s-number-of-factors' % tag,
              len(facs) == nfac == G * p,
              '%d translates = g(n-1) = %d*%d, the degree of K_{%dx4}' % (nfac, G, p, p + 1))
        check('D2-%s-the-factors-are-pairwise-distinct' % tag,
              len(set(frozenset(F) for F in facs)) == nfac,
              'the %d translates are %d different one-factors' % (nfac, nfac))
        check('D3-%s-every-factor-is-a-perfect-matching' % tag,
              all(r['is_perfect_matching'] and r['size'] == fsz for r in reps),
              'all %d factors are perfect matchings of size %d' % (nfac, fsz))
        check('D4-%s-no-factor-has-two-edges-on-one-part-pair' % tag,
              all(r['part_pairs_distinct'] for r in reps),
              'the OF axiom: no parallel edges in any of the %d factors' % nfac)
        check('D5-%s-minimum-Hamming-distance-inside-every-factor-is-at-least-3' % tag,
              all(r['min_distance'] >= 3 for r in reps),
              'each factor read as a constant-weight code has minimum distance %d >= 3'
              % min(r['min_distance'] for r in reps))
        allE = [e for F in facs for e in F]
        check('D6-%s-the-factors-PARTITION-the-edge-set' % tag,
              len(allE) == nE and len(set(allE)) == nE and set(allE) == Ep,
              '%d*%d = %d edges, all distinct, and equal to E(K_{%dx4}) as a set'
              % (nfac, fsz, nE, p + 1))
        npairs = c(p + 1)
        missed_count = {}
        for r in reps:
            for q in itertools.combinations(parts_of(p), 2):
                if frozenset(q) not in r['part_pairs']:
                    missed_count[frozenset(q)] = missed_count.get(frozenset(q), 0) + 1
        each = nfac * (npairs - fsz) // npairs
        check('D7-%s-every-part-pair-is-missed-the-same-number-of-times' % tag,
              set(missed_count.values()) == {each}
              and sum(missed_count.values()) == nfac * (npairs - fsz) == npairs * each,
              'each of the %d part-pairs is missed by exactly %d of the %d factors; '
              '%d*%d = %d = %d*%d' % (npairs, each, nfac, nfac, npairs - fsz,
                                      nfac * (npairs - fsz), npairs, each))
        check('D8-%s-part-pair-multiplicity-accounting' % tag,
              nfac * fsz == nE and npairs * G * G == nE,
              'every part-pair carries g^2 = 16 edges and %d*16 = %d = |E| = %d*%d'
              % (npairs, nE, nfac, fsz))

    # -------------------------------------------------------------- E
    section('E. CONTROLS ON THE DECIDER, IN BOTH POLARITIES')

    # E1/E2: a NEGATIVE control -- a perfect matching that is NOT a valid base factor
    # must be rejected, and its expansion must fail to be an OF.
    bad = [e for e in M6 if e not in (frozenset(((0, 1), (2, 1))), frozenset(((1, 3), (3, 1))))]
    bad = bad + [frozenset(((0, 1), (1, 3))), frozenset(((2, 1), (3, 1)))]
    badends = [x for e in bad for x in e]
    badrep = factor_report(bad, 5, G)
    check('E1-negative-control-is-still-a-perfect-matching',
          len(bad) == 12 and len(set(badends)) == 24 and set(badends) == vertices(5, G),
          'the perturbed base factor swaps two edges and is still a perfect matching, so '
          'the control isolates the two hypotheses and nothing else')
    badfacs = expand(bad, 5, G)
    badall = [e for F in badfacs for e in F]
    check('E2-negative-control-is-REJECTED',
          (not badrep['part_pairs_distinct'])
          and badrep['min_distance'] == 2
          and len(set(badall)) != 240,
          'the perturbed factor repeats the part-pair {0,1}, its minimum distance drops to 2, '
          'and its %d translates cover only %d distinct edges of 240 -- so the decider is not '
          'one that says yes unconditionally' % (len(badfacs), len(set(badall))))

    # E3: a POSITIVE control whose answer is classical and independent of this work --
    # the same free-action + transversal code path at g = 1 must reproduce the standard
    # 1-factorization of K_6 from the starter {inf,0},{1,4},{2,3}.
    K6start = [frozenset((('inf', 0), (0, 0))),
               frozenset(((1, 0), (4, 0))),
               frozenset(((2, 0), (3, 0)))]
    k6facs = expand(K6start, 5, 1)
    k6reps = [factor_report(F, 5, 1) for F in k6facs]
    k6all = [e for F in k6facs for e in F]
    check('E3-positive-control-classical-1-factorization-of-K_6',
          len(k6facs) == 5 and all(r['is_perfect_matching'] and r['size'] == 3 for r in k6reps)
          and len(set(k6all)) == 15 == len(k6all) and set(k6all) == edges_of(5, 1),
          'g = 1, p = 5: the 5 translates of the starter are 5 perfect matchings partitioning '
          'the 15 edges of K_6 -- the textbook 1-factorization, produced by the SAME code path')

    # E4/E5/E6: a PROVED NEGATIVE at n = g = 4, where no distance-3 one-factor can exist
    # because a factor needs gn/2 = 8 edges on distinct part-pairs and C(4,2) = 6 < 8.
    tot44, pf44 = enumerate_pm_kng(3, 4)
    check('E4-proved-negative-perfect-matchings-of-K_4x4-enumerated',
          tot44 == 368064,
          'exhaustive enumeration finds 368064 perfect matchings of K_{4x4}')
    coeffs = poly_pow([1, 6, 3], 4)
    ie = sum((-1) ** k * coeffs[k] * double_factorial_odd(15 - 2 * k) for k in range(len(coeffs)))
    check('E5-proved-negative-the-same-count-by-inclusion-exclusion',
          coeffs == [1, 24, 228, 1080, 2646, 3240, 2052, 648, 81]
          and double_factorial_odd(15) == 2027025 and ie == 368064 == tot44,
          '(1+6x+3x^2)^4 = %s and sum_k (-1)^k c_k (15-2k)!! = 368064, against 15!! = 2027025 '
          'for K_16 -- a second, independent route to the same integer' % (coeffs,))
    check('E6-proved-negative-none-of-them-is-parallel-edge-free',
          pf44 == 0 and 4 * 4 // 2 > c(4),
          '0 of the 368064 matchings has 8 distinct part-pairs, as gn/2 = 8 > C(4,2) = 6 '
          'forces; the decider therefore returns NO where the answer is provably NO')

    # -------------------------------------------------------------- F
    section('F. OF(12,4) BY THE SOURCE\'S DOUBLING LEMMA, AND THE g = 4 SPECTRUM')

    check('F1-doubling-lemma-hypotheses-hold-at-(n,g)=(6,4)',
          (6 * 4) % 2 == 0 and 6 > 4,
          'ng = 24 is even and n = 6 > g = 4, which is all the quoted lemma asks; it yields '
          'OF(12,4) from the OF(6,4) exhibited here')
    L1 = [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
    L2 = [[0, 1, 2, 3], [2, 3, 0, 1], [3, 2, 1, 0], [1, 0, 3, 2]]
    supers = set((L1[i][j], L2[i][j]) for i in range(4) for j in range(4))
    check('F2-the-lemma-s-side-conditions-at-g=4',
          4 not in (2, 6) and latin_ok(L1) and latin_ok(L2) and len(supers) == 16,
          'g = 4 is neither 2 nor 6, and the exhibited pair of order-4 Latin squares is '
          'orthogonal: their superimposition realises all 16 ordered pairs')

    odd_lo, odd_hi = 4 + 3, 2 * 4 - 3
    even_lo, even_hi = 4 + 2, 4 * 4 - 4
    printed_even = [n for n in range(even_lo, even_hi + 1) if n % 2 == 0]
    check('F3-the-source-s-printed-exception-ranges-at-g=4',
          odd_lo > odd_hi and (even_lo, even_hi) == (6, 12) and printed_even == [6, 8, 10, 12],
          'odd range [g+3, 2g-3] = [7,5] is EMPTY; even range [g+2, 4g-4] = [6,12] gives '
          'exactly {6,8,10,12}')

    def source(n):
        """Which result settles OF(n,4)?  First applicable clause wins."""
        if n <= 4:
            return None                       # excluded by the necessary condition n > g
        if n == 5:
            return 'source: OF(g+1,g), g even'
        if n % 2 == 1 and n >= 7:
            return 'source: odd n >= max(2g-1, g+3) = 7'
        if n in (6, 8):
            return 'THIS PAPER: explicit base factor'
        if n == 12:
            return 'THIS PAPER: OF(6,4) + the source\'s doubling lemma'
        if n % 4 == 2 and n >= 10:
            return 'source: doubling from odd m = n/2 >= 5'
        if n % 4 == 0 and n >= 16:
            return 'source: even n = 0 mod 4, n >= 16'
        return None

    NMAX = 4000
    gaps = [n for n in range(5, NMAX + 1) if source(n) is None]
    check('F4-the-g=4-spectrum-is-completely-covered',
          gaps == [],
          'every n with 5 <= n <= %d is settled by one of the six clauses, so an OF(n,4) '
          'exists for every n > 4' % NMAX)
    ours = sorted(n for n in range(5, NMAX + 1) if source(n).startswith('THIS PAPER'))
    check('F5-exactly-{6,8,12}-comes-from-this-paper',
          ours == [6, 8, 12]
          and sorted(set(printed_even) - set(ours)) == [10]
          and source(10) == 'source: doubling from odd m = n/2 >= 5',
          'the clauses of the source leave exactly {6,8,12}, which this paper supplies; '
          'n = 10 is NOT ours -- it is the source\'s OF(5,4) doubled')
    check('F6-no-claim-below-the-necessary-condition',
          all(source(n) is None for n in range(1, 5)) and (4 * 4) % 2 == 0,
          'nothing is asserted for n <= g = 4; at g = 4 the parity clause ng even is vacuous, '
          'so n > 4 is the whole necessary condition')

    # -------------------------------------------------------------- scope
    print('')
    print('NOT RE-RUN: the source\'s doubling lemma itself.  Checks F1-F2 verify its stated '
          'hypotheses and side conditions at (n,g) = (6,4); the lemma is QUOTED, not reproved, '
          'and NO OF(12,4) is constructed here.  Only OF(6,4) and OF(8,4) are exhibited and '
          'verified edge by edge by this program.')
    print('NOT RE-RUN: the clauses labelled "source:" in F4/F5.  They are the target paper\'s '
          'own results (OF(g+1,g) for even g, odd n >= 2g-1, the doubling lemma, and the '
          'n = 0 mod 4 proposition) and are taken on trust; F4 checks only that the six '
          'clauses TOGETHER (four of the source, two of this paper) leave no n > 4 uncovered, '
          'and it checks it over the finite range 5 <= n <= 4000; beyond it the clauses are '
          'residue conditions and the coverage argument is the paper\'s, not the program\'s.')
    print('NOT RE-RUN: MINIMALITY and UNIQUENESS.  No census of the prescribed-automorphism '
          'lane was run here, neither base factor is claimed unique or smallest, and no OF is '
          'claimed to be the only one on its parameters.')
    print('NOT RE-RUN: g >= 5, the range n <= g (the source\'s ODAR objects), and Hamming '
          'distances other than 3.  Nothing in this program touches them.')
    print('NOT RE-RUN: the wording, numbering and line offsets quoted from arXiv:2602.16319v1, '
          'and the prior-art channels (arXiv API, Crossref, Semantic Scholar, zbMATH) that '
          'bound the novelty.  MathSciNet was never consulted and OpenAlex never answered.  '
          'This program checks mathematics, not provenance.')

    print('')
    if _FAILED:
        print('FAILURES: %d -- %s' % (len(_FAILED), ', '.join(_FAILED)))
        print('VERDICT: NOT ALL CHECKS PASS')
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % len(_PASSED))
    return 0


if __name__ == '__main__':
    sys.exit(main())
