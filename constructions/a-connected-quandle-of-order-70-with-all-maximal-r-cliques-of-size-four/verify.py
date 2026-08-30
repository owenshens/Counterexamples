#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verification program for

    "A connected quandle of order 70 all of whose maximal R-cliques have size 4"

It re-derives every quantity the paper states, from the objects the paper prints and from
nothing else.  Python 3.9+, standard library only (no numpy, no sympy, no networkx, no
external data file).  All arithmetic is exact integer / Fraction arithmetic; no float ever
decides anything.

The three objects it consumes are exactly the three the paper exhibits:

  (Q)  the conjugation quandle on the 70 three-cycles of A_7, x*y = y^{-1} x y  (paper, Theorem 1 and Proposition 2);
  (C)  the maximal R-clique C = {(1 2 3), (1 3 2), (4 5 6), (4 6 5)}, parsed below from the
       paper's own printed cycle strings (paper, Theorem 1);
  (E)  the five right translations of the target paper's Example 5.16, quoted in Section 1 of the
       paper as read from arXiv:2504.09368v1, main.tex L1017.

Every check prints one `PASS <name> [detail]` line.  The program closes with
`VERDICT: ALL <n> CHECKS PASS` and exits 0 if and only if every check passed.
"""

import sys
from fractions import Fraction
from itertools import combinations, permutations

# ---------------------------------------------------------------------------
# the check harness
# ---------------------------------------------------------------------------
_N = [0]
_BAD = [0]


def check(name, got, want, detail=''):
    """One check.  Prints PASS on equality, FAIL otherwise; never raises."""
    _N[0] += 1
    if got == want:
        print('PASS %-34s %s' % (name, detail if detail else '%r' % (got,)))
    else:
        _BAD[0] += 1
        print('FAIL %-34s got %r, want %r  %s' % (name, got, want, detail))


# ---------------------------------------------------------------------------
# permutations of {0,...,n-1} as tuples of images
# ---------------------------------------------------------------------------
def comp(p, q):
    """(p then q):  i |-> q[p[i]].  Right action, matching x*y = y^{-1} x y below."""
    return tuple(q[i] for i in p)


def inv(p):
    r = [0] * len(p)
    for i, v in enumerate(p):
        r[v] = i
    return tuple(r)


def conj(x, y):
    """y^{-1} x y."""
    return comp(inv(y), comp(x, y))


def cyc(n, points):
    """The cycle (points[0] points[1] ... ) as a permutation of {0,...,n-1}."""
    p = list(range(n))
    k = len(points)
    for i in range(k):
        p[points[i]] = points[(i + 1) % k]
    return tuple(p)


def parse_cycle(s, n):
    """Parse a one-line-notation cycle exactly as the paper prints it, e.g. '(1 2 3)'.
    The paper's points are 1..n; we shift to 0..n-1."""
    s = s.strip()
    assert s.startswith('(') and s.endswith(')'), s
    pts = [int(t) - 1 for t in s[1:-1].replace(',', ' ').split()]
    assert all(0 <= t < n for t in pts) and len(set(pts)) == len(pts), s
    return cyc(n, pts)


def perm_of_cycles(n, *cycles):
    p = tuple(range(n))
    for c in cycles:
        p = comp(p, cyc(n, list(c)))
    return p


def three_cycles(n):
    """The 2*C(n,3) three-cycles of A_n, in a canonical order."""
    out = []
    for s in combinations(range(n), 3):
        out.append(cyc(n, [s[0], s[1], s[2]]))
        out.append(cyc(n, [s[0], s[2], s[1]]))
    return sorted(set(out))


# ---------------------------------------------------------------------------
# quandle machinery, all from the multiplication table
# ---------------------------------------------------------------------------
def table_from_class(els):
    """T[x][y] = index of y^{-1} x y.  Returns (T, None) or (None, reason)."""
    idx = {e: i for i, e in enumerate(els)}
    n = len(els)
    T = [[0] * n for _ in range(n)]
    for i, a in enumerate(els):
        for j, b in enumerate(els):
            c = conj(a, b)
            if c not in idx:
                return None, 'not closed under conjugation'
            T[i][j] = idx[c]
    return T, None


def right_translations(T):
    n = len(T)
    return [tuple(T[x][a] for x in range(n)) for a in range(n)]


def is_right_quasigroup(T):
    n = len(T)
    return all(sorted(c) == list(range(n)) for c in right_translations(T))


def is_idempotent(T):
    return all(T[i][i] == i for i in range(len(T)))


def rack_violations(T):
    """(x*y)*z == (x*z)*(y*z) over all n^3 triples; returns the violation count."""
    n = len(T)
    bad = 0
    for x in range(n):
        Tx = T[x]
        for y in range(n):
            xy = T[x][y]
            Ty = T[y]
            Txy = T[xy]
            for z in range(n):
                if Txy[z] != T[Tx[z]][Ty[z]]:
                    bad += 1
    return bad


def mlt_r_orbit(T, start=0):
    n = len(T)
    seen = {start}
    stack = [start]
    while stack:
        u = stack.pop()
        for a in range(n):
            v = T[u][a]
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


def commuting_graph_from_table(T):
    """a ~ b  iff  R_a R_b = R_b R_a as maps on Q."""
    n = len(T)
    cols = right_translations(T)
    adj = [set() for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if comp(cols[a], cols[b]) == comp(cols[b], cols[a]):
                adj[a].add(b)
                adj[b].add(a)
    return adj


def commuting_graph_from_group(els):
    """a ~ b  iff  ab = ba in the ambient symmetric group."""
    n = len(els)
    adj = [set() for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if comp(els[a], els[b]) == comp(els[b], els[a]):
                adj[a].add(b)
                adj[b].add(a)
    return adj


def edge_set(adj):
    return set((a, b) for a in range(len(adj)) for b in adj[a] if a < b)


def maximal_cliques(adj):
    """Bron-Kerbosch with a deterministic pivot.  Yields every maximal clique."""
    n = len(adj)

    def expand(R, P, X):
        if not P and not X:
            yield tuple(sorted(R))
            return
        pivot = max(P | X, key=lambda u: len(adj[u] & P))
        for v in sorted(P - adj[pivot]):
            for c in expand(R | {v}, P & adj[v], X & adj[v]):
                yield c
            P = P - {v}
            X = X | {v}

    return expand(set(), set(range(n)), set())


def is_clique(adj, S):
    S = list(S)
    return all(b in adj[a] for a, b in combinations(S, 2))


def is_maximal_clique(adj, S):
    S = set(S)
    if not is_clique(adj, S):
        return False
    return not any(S <= adj[v] for v in range(len(adj)) if v not in S)


def multinomial_packings(n, k):
    """Number of ways to pick k pairwise disjoint unordered triples from n points:
    n! / ((n-3k)! * 6^k * k!).  Exact integer."""
    num = 1
    for i in range(n - 3 * k + 1, n + 1):
        num *= i
    den = 1
    for i in range(1, k + 1):
        den *= 6 * i
    assert num % den == 0
    return num // den


def binom(n, k):
    if k < 0 or k > n:
        return 0
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r


# ===========================================================================
banner = ('=== verification of "A connected quandle of order 70 all of whose maximal '
          'R-cliques have size 4" ===')
print(banner)
print('python %s' % sys.version.split()[0])
print('stdlib only; exact integer / Fraction arithmetic throughout')
print('')

# ---------------------------------------------------------------------------
# PART A.  The object of Theorem 1: Q = 3-cycles of A_7 under x*y = y^{-1} x y
# ---------------------------------------------------------------------------
print('--- A. the witness Q: the 70 three-cycles of A_7 under conjugation ---')

els7 = three_cycles(7)
check('A01 |Q| = 70', len(els7), 70, '|Q| = %d' % len(els7))
check('A02 |Q| = 2*C(7,3)', len(els7), 2 * binom(7, 3), '2*C(7,3) = 2*%d = %d' % (binom(7, 3), 2 * binom(7, 3)))
check('A03 every element is a 3-cycle',
      sorted(set(len([i for i in range(7) if p[i] != i]) for p in els7)), [3],
      'each moves exactly 3 of the 7 points')
check('A04 every element is even',
      all(sum(1 for i, j in combinations(range(7), 2) if p[i] > p[j]) % 2 == 0 for p in els7),
      True, 'all 70 lie in A_7')

T7, why = table_from_class(els7)
check('A05 Q closed under conjugation', why, None, 'the 70x70 table T[x][y] = index(y^{-1} x y) is defined')
check('A06 right quasigroup', is_right_quasigroup(T7), True, 'each of the 70 maps R_y is a permutation of Q')
check('A07 idempotent x*x = x', is_idempotent(T7), True, 'so Q is a quandle, not merely a rack')
rv = rack_violations(T7)
check('A08 self-distributivity', rv, 0,
      '0 violations of (x*y)*z = (x*z)*(y*z) over all 70^3 = %d triples' % (70 ** 3))
check('A09 70^3 triples were checked', 70 ** 3, 343000, 'the sweep is exhaustive, not sampled')

cols7 = right_translations(T7)
check('A10 right translations distinct', len(set(cols7)), 70,
      'all 70 R_y differ, so the Cayley kernel is trivial')
orb = mlt_r_orbit(T7, 0)
check('A11 Q is connected', len(orb), 70, 'the Mlt_r-orbit of one point is all 70 elements')

# the class does not split in A_7: the odd permutation (4 5) centralises (1 2 3)
a123 = parse_cycle('(1 2 3)', 7)
t45 = parse_cycle('(4 5)', 7)
check('A12 (4 5) centralises (1 2 3)', comp(a123, t45), comp(t45, a123),
      'an ODD centralising element, so the S_7-class does not split in A_7')
cent = [p for p in permutations(range(7)) if comp(a123, p) == comp(p, a123)]
check('A13 |C_{S_7}((1 2 3))| = 72', len(cent), 72, '5040/70 = 72 = 3*4!, so the S_7-class has 70 elements')
check('A14 5040/70 = 72', 5040 // 70, 72, 'index computation closes')

adjT = commuting_graph_from_table(T7)
adjG = commuting_graph_from_group(els7)
check('A15 R-commuting = group commuting', edge_set(adjT), edge_set(adjG),
      'the two edge sets are equal (%d edges each), which is the step Z(A_7) = 1 predicts'
      % len(edge_set(adjT)))
check('A16 commuting graph is 9-regular', sorted(set(len(s) for s in adjT)), [9],
      'every vertex has 1 + 2*C(4,3) = %d neighbours' % (1 + 2 * binom(4, 3)))
check('A17 edge count 315', len(edge_set(adjT)), 70 * 9 // 2, '70*9/2 = 315')
check('A18 neighbour formula', 1 + 2 * binom(4, 3), 9, 'inverse, plus both orders on each triple in the other 4 points')

# --- the explicit C, parsed from the paper's printed cycle strings ----------
C_AS_PRINTED = '(1 2 3), (1 3 2), (4 5 6), (4 6 5)'


def split_cycles(s):
    """Split a printed comma-separated list of cycles, e.g. '(1 2 3), (1 3 2)'."""
    out, depth, cur = [], 0, ''
    for ch in s:
        if ch == '(':
            depth += 1
        if ch == ',' and depth == 0:
            out.append(cur)
            cur = ''
            continue
        cur += ch
        if ch == ')':
            depth -= 1
    if cur.strip():
        out.append(cur)
    return [t.strip() for t in out]


idx7 = {e: i for i, e in enumerate(els7)}
Cperms = [parse_cycle(s, 7) for s in split_cycles(C_AS_PRINTED)]
check('A19 C parsed from the paper', len(Cperms), 4, 'read from the printed string %r' % C_AS_PRINTED)
check('A20 C lies in Q', all(p in idx7 for p in Cperms), True, 'all four printed cycles are elements of Q')
C = sorted(idx7[p] for p in Cperms)
check('A21 |C| = 4', len(set(C)), 4, 'the four printed elements are distinct')
check('A22 C is an R-clique', is_clique(adjT, C), True, 'pairwise R-commuting')
check('A23 C is MAXIMAL', is_maximal_clique(adjT, C), True,
      'no element of Q outside C R-commutes with all of C')
check('A24 4 does not divide 70', 70 % 4, 2, '|Q| mod |C| = 70 mod 4 = 2, so |C| does not divide |Q|')
check('A25 |C| divides |Q| is False', (70 % 4 == 0), False, 'Problem 5.24 is answered YES by this object')

# controls on the maximality decider: it must be able to say NO
check('A26 CONTROL C minus one not maximal', is_maximal_clique(adjT, C[:3]), False,
      'a proper subset of C is a clique but must not be reported maximal')
pair = sorted(idx7[parse_cycle(s, 7)] for s in ('(1 2 3)', '(1 3 2)'))
check('A27 CONTROL {(1 2 3),(1 3 2)} not maximal', is_maximal_clique(adjT, pair), False,
      'an extendable clique must not be reported maximal')
check('A28 CONTROL a non-edge exists', is_clique(adjT, [idx7[parse_cycle('(1 2 3)', 7)],
                                                       idx7[parse_cycle('(1 2 4)', 7)]]), False,
      '(1 2 3) and (1 2 4) share two points and do not commute, so the adjacency test is not vacuous')

# --- every maximal R-clique -------------------------------------------------
cl7 = list(maximal_cliques(adjT))
check('A29 number of maximal R-cliques', len(cl7), 70, 'Bron-Kerbosch enumerates 70 of them')
check('A30 every maximal R-clique has size 4', sorted(set(len(c) for c in cl7)), [4],
      'the size multiset is {4: %d}' % len(cl7))
check('A31 C is among them', tuple(sorted(C)) in set(cl7), True, 'the exhibited clique is one of the 70')
check('A32 incidence identity', binom(7, 3) * binom(4, 3) // 2, 70,
      'C(7,3)*C(4,3)/2 = 35*4/2 = 70 = the clique count, an independent route to the same number')
check('A33 all maximal cliques are non-dividing',
      sorted(set(70 % len(c) for c in cl7)), [2], 'every one of the 70 leaves remainder 2')
check('A34 maximal cliques are packings',
      all(len(set(i for p in c for i in range(7) if p[i] != i)) == 6
          for c in [[els7[i] for i in cl] for cl in cl7]),
      True, 'each clique covers exactly 6 of the 7 points, i.e. two disjoint triples')
check('A35 packing count formula', multinomial_packings(7, 2), 70,
      '7!/(1! * 6^2 * 2!) = 70 maximum triple packings of 7 points')

# --- tamper control: the rack checker must not pass vacuously ---------------
Tbad = [row[:] for row in T7]
Tbad[0][0], Tbad[1][0] = Tbad[1][0], Tbad[0][0]
check('A36 TAMPER rack checker is not vacuous', rack_violations(Tbad) > 0, True,
      'transposing two entries in one column of the table yields %d violations' % rack_violations(Tbad))
check('A37 TAMPER quandle checker is not vacuous', is_idempotent(Tbad), False,
      'the same tampered table is no longer idempotent')

# ---------------------------------------------------------------------------
# PART B.  The target paper's own near-witness, Example 5.16 (main.tex L1017)
# ---------------------------------------------------------------------------
print('')
print('--- B. the target paper\'s Example 5.16, as printed in Section 1 of our paper ---')

R_5_16 = [perm_of_cycles(5, (2, 3, 4)),
          perm_of_cycles(5, (2, 4, 3)),
          perm_of_cycles(5, (0, 1), (3, 4)),
          perm_of_cycles(5, (0, 1), (2, 4)),
          perm_of_cycles(5, (0, 1), (2, 3))]
T5 = [[R_5_16[y][x] for y in range(5)] for x in range(5)]     # x*y = R_y(x)
check('B01 the five printed maps are permutations',
      all(sorted(p) == list(range(5)) for p in R_5_16), True, 'each R_y is a bijection of a 5-element set')
check('B02 Example 5.16 is idempotent', is_idempotent(T5), True, 'R_y fixes y for each y')
check('B03 Example 5.16 is a rack', rack_violations(T5), 0, '0 violations over 5^3 = 125 triples')
check('B04 recovered R_y equal the printed ones', right_translations(T5), R_5_16,
      'the table convention T[x][y] = R_y(x) round-trips')
adj5 = commuting_graph_from_table(T5)
cl5 = sorted(list(maximal_cliques(adj5)))
check('B05 its maximal R-cliques', cl5, [(0, 1), (2,), (3,), (4,)],
      'sizes 2, 1, 1, 1 -- the paper states only that two of them are maximal, and it is right')
check('B06 a non-dividing size occurs', sorted(set(5 % len(c) for c in cl5)), [0, 1],
      '|C| = 2 does not divide |Q| = 5, so the decider genuinely can return a non-dividing hit')
check('B07 Example 5.16 is NOT connected', sorted(mlt_r_orbit(T5, 0)), [0, 1],
      'the Mlt_r-orbit of 0 is the block {0,1}, so the connectedness hypothesis fails there')
check('B08 the block {0,1} is Mlt_r-invariant',
      all(T5[x][y] in (0, 1) for x in (0, 1) for y in range(5)), True,
      'every generator preserves {0,1}, which is why the authors\' near-witness is not a witness')

# ---------------------------------------------------------------------------
# PART C.  The infinite family Q_n, n = 6m+1
# ---------------------------------------------------------------------------
print('')
print('--- C. the family Q_n = 3-cycles of A_n, n = 6m+1 (paper, Theorem 3) ---')

# C.1 the combinatorial claim, by direct enumeration on the commuting graph
FULL_N = 8          # rack axioms are re-swept in full for n <= FULL_N
GRAPH_N = 12        # the commuting graph and its maximal cliques for n <= GRAPH_N
for n in range(5, GRAPH_N + 1):
    els = three_cycles(n)
    k = n // 3
    expect_size = 2 * k
    expect_count = multinomial_packings(n, k)
    adj = commuting_graph_from_group(els)
    cls = list(maximal_cliques(adj))
    sizes = sorted(set(len(c) for c in cls))
    check('C%02d n=%-2d clique sizes' % (n, n), sizes, [expect_size],
          '|Q_%d| = %d, every maximal R-clique has size 2*floor(%d/3) = %d' % (n, len(els), n, expect_size))
    check('C%02d n=%-2d clique count' % (n + 20, n), len(cls), expect_count,
          'matches n!/((n-3k)! 6^k k!) = %d' % expect_count)
    check('C%02d n=%-2d divisibility' % (n + 40, n), (len(els) % expect_size == 0), (n % 6 != 1),
          '|C| divides |Q| is %s; a hit occurs exactly when n = 1 mod 6'
          % (len(els) % expect_size == 0))
    if n <= FULL_N:
        Tn, whyn = table_from_class(els)
        check('C%02d n=%-2d quandle axioms' % (n + 60, n),
              (whyn, is_right_quasigroup(Tn), is_idempotent(Tn), rack_violations(Tn),
               len(mlt_r_orbit(Tn, 0)) == len(els)),
              (None, True, True, 0, True),
              'closed, right quasigroup, idempotent, 0 rack violations over %d^3 triples, connected' % len(els))
        check('C%02d n=%-2d table vs group commuting' % (n + 80, n),
              edge_set(commuting_graph_from_table(Tn)), edge_set(adj),
              'R_aR_b = R_bR_a  iff  ab = ba, on labels, %d edges' % len(edge_set(adj)))

# C.2 the arithmetic of the family, exactly, for m = 1..40
MFAM = 40
bad_int = []
rows = []
for m in range(1, MFAM + 1):
    n = 6 * m + 1
    order = 2 * binom(n, 3)
    check_order = 2 * m * (36 * m * m - 1)
    size = 2 * (n // 3)
    ratio = Fraction(order, size)
    rows.append((m, n, order, size, ratio))
    if order != check_order or size != 4 * m or ratio != Fraction(36 * m * m - 1, 2) or ratio.denominator == 1:
        bad_int.append(m)
check('C101 family closed form', bad_int, [],
      'for m = 1..%d: |Q_n| = 2*C(n,3) = 2m(36m^2-1), |C| = 4m, |Q_n|/|C| = (36m^2-1)/2, never an integer'
      % MFAM)
check('C102 36m^2-1 is odd', sorted(set((36 * m * m - 1) % 2 for m in range(1, MFAM + 1))), [1],
      'which is exactly why the ratio has denominator 2')
check('C103 the four displayed ratios',
      [rows[i][4] for i in range(4)],
      [Fraction(35, 2), Fraction(143, 2), Fraction(323, 2), Fraction(575, 2)],
      '70/4 = 35/2, 572/8 = 143/2, 1938/12 = 323/2, 4600/16 = 575/2')
check('C104 the four displayed orders', [rows[i][2] for i in range(4)], [70, 572, 1938, 4600],
      'n = 7, 13, 19, 25')
check('C105 the four displayed clique sizes', [rows[i][3] for i in range(4)], [4, 8, 12, 16],
      '|C| = 4m')

# C.3 the divisibility criterion over a long range
hits = [n for n in range(4, 401) if binom(n, 3) % (n // 3) != 0]
check('C106 criterion floor(n/3) does not divide C(n,3)',
      hits, [n for n in range(7, 401) if n % 6 == 1],
      'over 4 <= n <= 400 the non-dividing n are exactly n = 1 mod 6, n >= 7 (%d values)' % len(hits))
check('C107 the named non-hits', [n for n in (5, 6, 8, 9, 10, 11, 12) if n in hits], [],
      'n = 5,6,8,9,10,11,12 all divide, as the paper states')
check('C108 n = 7 is the least hit', min(hits), 7, 'the family first produces a witness at n = 7, order 70')

# ---------------------------------------------------------------------------
print('')
print('NOT RE-RUN: the minimality question.  This program does NOT enumerate connected quandles by')
print('NOT RE-RUN: order, so it neither confirms nor contradicts the claim that no connected quandle of')
print('NOT RE-RUN: order at most 24 is a witness; that census is reported in the paper as unreproduced,')
print('NOT RE-RUN: orders 25-69 were never examined, and connected RACKS THAT ARE NOT QUANDLES were')
print('NOT RE-RUN: never examined at any order.  Nothing here asserts that 70 is the least order.')
print('NOT RE-RUN: the version of record.  The problem statement and its locator in S1 are quoted from')
print('NOT RE-RUN: arXiv:2504.09368v1; the journal text (J. Algebra 698 (2026) 493-532) was not')
print('NOT RE-RUN: accessible to us and its problem numbering may differ.  No program can check that.')
print('NOT RE-RUN: Part C sweeps the commuting graph for 5 <= n <= %d and the full rack identity only for'
      % GRAPH_N)
print('NOT RE-RUN: n <= %d; for n = 13 and beyond only the closed-form arithmetic (C101-C108) is checked.'
      % FULL_N)
print('NOT RE-RUN: no claim about Problem 5.22 (do all maximal R-cliques of a connected rack have the')
print('NOT RE-RUN: same size?) is made or tested here; in Q every maximal R-clique has size 4.')
print('')

if _BAD[0]:
    print('VERDICT: %d of %d CHECKS FAILED' % (_BAD[0], _N[0]))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _N[0])
sys.exit(0)
