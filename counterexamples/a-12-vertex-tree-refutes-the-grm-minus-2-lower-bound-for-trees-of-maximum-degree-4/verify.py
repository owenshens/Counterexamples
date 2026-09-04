#!/usr/bin/env python3
"""verify.py -- checker for the note:

    "A 12-vertex tree refuting the GRM_{-2} lower bound for trees of maximum degree four"

It reads the edge lists of W1 (n=12), W2 (n=16), W3 (n=20) and the n=8 double star D8, together
with the parametric family F(k), and re-derives numerical quantities from them.  The edge lists of
W1 and W2 and the values -13 and -18 are printed in the note; the W3 edge list and the value -22
are supplied in this file and are not printed in the note.

WHAT IS CHECKED
  * each exhibited graph really is a tree with maximum degree exactly 4 (so it lies in the class
    fixed at main.tex l.797 of arXiv:2604.06044v1);
  * its value of GRM_{-2} = sum_{uv in E} (deg u - 2)(deg v - 2) computed THREE independent ways:
    directly from the definition, from the source's pre-substitution formula at l.795, and from
    the source's equation (main2) at l.800 -- all three forced to agree;
  * its edge partition m_ij and vertex counts n_i, against the tables printed in the paper;
  * the source's OWN degree-count system (prva)-(sesta) at l.785-790, evaluated in EXACT rational
    arithmetic on the witness's own free parameters, returns the witness's actual
    n_1, n_2, n_4, m_14, m_24, m_33 -- so the counterexample is inside the source's
    parameterisation, not outside it;
  * the arithmetic identity GRM_{-2} = -(n+3) + S and a companion counting identity, on each
    exhibited graph;
  * the family F(k) built by this file has order 4k+4, maximum degree 4 and GRM_{-2} = -(n+2) for
    k = 3..12, and F(3), F(4) are isomorphic to W2, W3;
  * the source's claimed extremal parameter sets at l.843 and l.844 both evaluate to -n, so the
    claimed extremal family sits strictly above the exhibited witnesses.

CONTROLS, BOTH POLARITIES
  * FORCED POSITIVES against values written as literals in this file, not taken from the note:
    K_{1,3} -> -3, K_{1,4} -> -8, K_{1,5} -> -15, P_5 -> 0.
  * THE ANTI-CONTROL, and it is the load-bearing one: the n=8 double star -- the source's OWN
    m_44 = 1 case at k = 1 -- must be reported as SATISFYING the refuted bound, because n = 8 is
    the single order in the residue class where the source is right. A checker that flagged n = 8
    as well would be flagging everything and would decide nothing.
  * FORCED NEGATIVES: the source's earlier and weaker bounds at l.822 (>= -(n+2)) and l.830
    (>= -(n+1)) are checked to be RESPECTED by W1 at n = 12; only l.838 is refuted there.

Python 3.9+, standard library only (fractions), no external data file, no arguments, no network,
no randomness. Exact integer / Fraction arithmetic throughout; no floating point is used anywhere.
Exits 0 if and only if every check passes.
"""

from fractions import Fraction as F

CHECKS = []          # (name, ok, detail)


def ck(name, ok, detail=''):
    CHECKS.append((name, bool(ok), detail))
    print('%s %s%s' % ('PASS' if ok else 'FAIL', name, ('   ' + detail) if detail else ''))
    return bool(ok)


# ---------------------------------------------------------------------------
# The objects, exactly as printed in the paper.
# ---------------------------------------------------------------------------
W1 = (12, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (5, 6), (5, 7),
           (6, 8), (8, 9), (8, 10), (8, 11)])

W2 = (16, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (5, 6), (5, 7), (6, 8),
           (7, 9), (8, 10), (8, 11), (8, 14), (9, 12), (9, 13), (9, 15)])

W3 = (20, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (2, 6), (5, 7), (5, 8),
           (6, 9), (6, 10), (6, 14), (7, 11), (8, 12), (11, 13), (11, 15),
           (11, 16), (12, 17), (12, 18), (12, 19)])

# The source's own m_44 = 1 case at k = 1 (main.tex l.844): the n = 8 double star. ANTI-CONTROL.
D8 = (8, [(0, 1), (0, 2), (0, 3), (0, 4), (4, 5), (4, 6), (4, 7)])

STARS = [('K_{1,3}', 4, [(0, 1), (0, 2), (0, 3)], -3),
         ('K_{1,4}', 5, [(0, 1), (0, 2), (0, 3), (0, 4)], -8),
         ('K_{1,5}', 6, [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)], -15),
         ('P_5', 5, [(0, 1), (1, 2), (2, 3), (3, 4)], 0)]

# The tables printed in the paper, one per exhibited graph:
#   name -> (n, edges, claimed GRM_{-2}, claimed m_ij, claimed n_i, claimed bound status)
PRINTED = [
    ('W1', W1[0], W1[1], -13,
     {(1, 3): 1, (1, 4): 6, (2, 3): 2, (2, 4): 2},
     {1: 7, 2: 2, 3: 1, 4: 2}, 'VIOLATED', 1),
    ('W2', W2[0], W2[1], -18,
     {(1, 4): 9, (2, 3): 3, (2, 4): 3},
     {1: 9, 2: 3, 3: 1, 4: 3}, 'VIOLATED', 2),
    ('W3', W3[0], W3[1], -22,
     {(1, 4): 11, (2, 3): 3, (2, 4): 5},
     {1: 11, 2: 4, 3: 1, 4: 4}, 'VIOLATED', 2),
    ('D8', D8[0], D8[1], -8,
     {(1, 4): 6, (4, 4): 1},
     {1: 6, 4: 2}, 'HOLDS', 0),
]


# ---------------------------------------------------------------------------
# Graph primitives (integer arithmetic only)
# ---------------------------------------------------------------------------
def degrees(n, E):
    d = [0] * n
    for a, b in E:
        d[a] += 1
        d[b] += 1
    return d


def is_tree(n, E):
    if len(E) != n - 1:
        return False, 'edge count %d != n-1 = %d' % (len(E), n - 1)
    seen = set()
    for a, b in E:
        if a == b:
            return False, 'loop at %d' % a
        k = (min(a, b), max(a, b))
        if k in seen:
            return False, 'repeated edge %s' % (k,)
        seen.add(k)
    adj = dict((v, []) for v in range(n))
    for a, b in E:
        adj[a].append(b)
        adj[b].append(a)
    comp, stack = set([0]), [0]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in comp:
                comp.add(w)
                stack.append(w)
    if len(comp) != n:
        return False, 'not connected: traversal from 0 reached %d of %d vertices' % (len(comp), n)
    return True, '%d vertices, %d edges, no loop, no repeated edge, connected' % (n, len(E))


def grm_direct(n, E):
    """The definition at main.tex l.121 / l.154 with lambda = -2."""
    d = degrees(n, E)
    return sum((d[a] - 2) * (d[b] - 2) for a, b in E)


def profile(n, E):
    d = degrees(n, E)
    m, ni = {}, {}
    for a, b in E:
        k = (min(d[a], d[b]), max(d[a], d[b]))
        m[k] = m.get(k, 0) + 1
    for v in range(n):
        ni[d[v]] = ni.get(d[v], 0) + 1
    return m, ni, d


def canon(n, E):
    """AHU canonical form of a free tree, rooted at its centre(s). Used for two isomorphism tests
    only; it is O(n log n) on graphs of at most 20 vertices."""
    if n == 1:
        return '()'
    adj = dict((v, set()) for v in range(n))
    for a, b in E:
        adj[a].add(b)
        adj[b].add(a)
    deg = dict((v, len(adj[v])) for v in range(n))
    leaves = [v for v in range(n) if deg[v] <= 1]
    removed, remaining = set(), n
    while remaining > 2:
        nxt = []
        for v in leaves:
            removed.add(v)
            remaining -= 1
            for w in adj[v]:
                if w not in removed:
                    deg[w] -= 1
                    if deg[w] == 1:
                        nxt.append(w)
        leaves = nxt
    centres = sorted(v for v in range(n) if v not in removed)

    def enc(v, parent):
        kids = sorted(enc(w, v) for w in adj[v] if w != parent)
        return '(' + ''.join(kids) + ')'

    forms = []
    for c in centres:
        forms.append(enc(c, None))
    return min(forms)


def family(k):
    """The family F(k) as constructed here: order n = 4k+4, GRM_{-2} = -(n+2).

    Skeleton on k+1 nodes: node 0 (the future degree-3 vertex) joined to nodes 1, 2, 3, and nodes
    4..k hung in a path off node 3. Every skeleton edge is subdivided exactly once, then each of
    the nodes 1..k receives 4 - (its skeleton degree) pendant vertices.
    """
    assert k >= 3
    sk = [(0, 1), (0, 2), (0, 3)] + [(i, i + 1) for i in range(3, k)]
    E, nxt = [], k + 1
    for a, b in sk:                       # subdivide each skeleton edge once
        E.append((a, nxt))
        E.append((nxt, b))
        nxt += 1
    sdeg = dict((v, 0) for v in range(k + 1))
    for a, b in sk:
        sdeg[a] += 1
        sdeg[b] += 1
    for v in range(1, k + 1):             # top every non-hub-of-degree-3 node up to degree 4
        for _ in range(4 - sdeg[v]):
            E.append((v, nxt))
            nxt += 1
    return nxt, E


def paper_system(n, m, ni):
    """(prva)-(sesta), main.tex l.785-790, in EXACT rational arithmetic. Returns
    [(label, predicted, actual)] for n_1, n_2, n_4, m_14, m_24, m_33."""
    g = lambda i, j: m.get((i, j), 0)
    m12, m13, m14 = g(1, 2), g(1, 3), g(1, 4)
    m22, m23, m24 = g(2, 2), g(2, 3), g(2, 4)
    m33, m34, m44 = g(3, 3), g(3, 4), g(4, 4)
    n1, n2, n3, n4 = ni.get(1, 0), ni.get(2, 0), ni.get(3, 0), ni.get(4, 0)
    q = F
    p_n1 = (-q(m12, 2) - q(m13, 4) - q(m22, 2) - q(m23, 4) + q(m34, 4) + q(m44, 2)
            + q(n, 2) + q(n3, 4) + q(3, 2))
    p_n2 = (q(3 * m12, 4) + q(3 * m13, 8) + q(3 * m22, 4) + q(3 * m23, 8) - q(3 * m34, 8)
            - q(3 * m44, 4) + q(n, 4) - q(7 * n3, 8) - q(5, 4))
    p_n4 = (-q(m12, 4) - q(m13, 8) - q(m22, 4) - q(m23, 8) + q(m34, 8) + q(m44, 4)
            + q(n, 4) - q(3 * n3, 8) - q(1, 4))
    p_m14 = (-q(3 * m12, 2) - q(5 * m13, 4) - q(m22, 2) - q(m23, 4) + q(m34, 4) + q(m44, 2)
             + q(n, 2) + q(n3, 4) + q(3, 2))
    p_m24 = (q(m12, 2) + q(3 * m13, 4) - q(m22, 2) - q(m23, 4) - q(3 * m34, 4) - q(3 * m44, 2)
             + q(n, 2) - q(7 * n3, 4) - q(5, 2))
    p_m33 = -q(m13, 2) - q(m23, 2) - q(m34, 2) + q(3 * n3, 2)
    return [('n_1', p_n1, n1), ('n_2', p_n2, n2), ('n_4', p_n4, n4),
            ('m_14', p_m14, m14), ('m_24', p_m24, m24), ('m_33', p_m33, m33)]


# ---------------------------------------------------------------------------
print('CHECKER for "A 12-vertex tree refuting the GRM_{-2} lower bound for trees of maximum')
print('degree four".  Source under test: arXiv:2604.06044v1, main.tex l.838,')
print('  "Therefore, based on equation (\\ref{main2}), if $n \\not\\equiv \\{1, 2, 3\\} \\pmod{4}$,')
print('   it follows that $GRM_{-2}(T) \\geq -n$."')
print('GRM_{-2}(T) = sum_{uv in E(T)} (deg u - 2)(deg v - 2)   [l.121/l.154 at lambda = -2]')
print('')

print('--- FORCED POSITIVE CONTROLS: values written as literals in this file ------------------')
for nm, n, E, want in STARS:
    got = grm_direct(n, E)
    ck('control-%s' % nm.replace('_', '').replace('{', '').replace('}', '').replace(',', ''),
       got == want, 'GRM_{-2}(%s) = %d, literal in this file %d' % (nm, got, want))
print('')

for nm, n, E, want_grm, want_m, want_ni, want_status, want_gap in PRINTED:
    print('--- %s,  n = %d ---------------------------------------------------------------' % (nm, n))
    ok, why = is_tree(n, E)
    ck('%s-is-a-tree' % nm, ok, why)
    m, ni, d = profile(n, E)
    ck('%s-degree-sum-is-twice-the-edge-count' % nm, sum(d) == 2 * len(E),
       'sum of degrees %d = 2 * %d edges; degree sequence %s'
       % (sum(d), len(E), sorted(d, reverse=True)))
    ck('%s-maximum-degree-is-exactly-4' % nm, max(d) == 4,
       'Delta = %d, attained at %s -- inside the class fixed at main.tex l.797'
       % (max(d), sorted(v for v in range(n) if d[v] == 4)))
    ck('%s-edge-partition-matches-the-expected-table' % nm, m == want_m,
       'm_ij = %s' % dict(('m%d%d' % k, v) for k, v in sorted(m.items())))
    ck('%s-vertex-counts-match-the-expected-table' % nm, ni == want_ni,
       'n_i = %s' % dict(('n%d' % k, v) for k, v in sorted(ni.items())))

    g1 = grm_direct(n, E)
    ck('%s-GRM-minus-2-equals-the-expected-value' % nm, g1 == want_grm,
       'direct edge sum = %d, expected value recorded in this file %d' % (g1, want_grm))

    gg = lambda i, j: m.get((i, j), 0)
    m12, m13, m14 = gg(1, 2), gg(1, 3), gg(1, 4)
    m22, m23, m24 = gg(2, 2), gg(2, 3), gg(2, 4)
    m33, m34, m44 = gg(3, 3), gg(3, 4), gg(4, 4)
    n1, n2, n3, n4 = ni.get(1, 0), ni.get(2, 0), ni.get(3, 0), ni.get(4, 0)

    g2 = -(m13 + 2 * m14 - m33 - 2 * m34 - 4 * m44)                    # main.tex l.795
    ck('%s-source-formula-l795-agrees' % nm, g2 == g1,
       '-(m13+2m14-m33-2m34-4m44) = %d' % g2)
    g3 = -(-3 * m12 - m13 - m22 - m34 - 3 * m44 + n - n3 + 3)          # main.tex l.800 (main2)
    ck('%s-source-equation-main2-l800-agrees' % nm, g3 == g1,
       '-(-3m12-m13-m22-m34-3m44+n-n3+3) = %d' % g3)

    S = n3 + 3 * m12 + m13 + m22 + m34 + 3 * m44
    ck('%s-identity-GRM-equals-minus-n-plus-3-plus-S' % nm, -(n + 3) + S == g1,
       'S = %d, so -(n+3)+S = %d' % (S, -(n + 3) + S))

    rows = paper_system(n, m, ni)
    bad = [(lab, p, a) for lab, p, a in rows if p != a]
    ck('%s-satisfies-the-source-system-prva-to-sesta' % nm, not bad,
       'exact rational evaluation at l.785-790 returns ' +
       ', '.join('%s=%s' % (lab, p) for lab, p, a in rows) +
       ('' if not bad else '   MISMATCH: %s' % (bad,)))

    if n % 4 == 0:
        violated = g1 < -n
        gap = (-n) - g1 if violated else 0
        status = 'VIOLATED' if violated else 'HOLDS'
        ck('%s-l838-bound-status-is-%s' % (nm, want_status), status == want_status and gap == want_gap,
           'claimed bound -n = %d, actual GRM_{-2} = %d  ==>  l.838 %s%s'
           % (-n, g1, status, (' by %d' % gap) if violated else ''))

    if m12 == 0 and m44 == 0:
        ck('%s-counting-identity-n-equals-3n3-4n4-1-m22-m33-m34' % nm,
           n == 3 * n3 + 4 * n4 + 1 + m22 - m33 - m34,
           '3*%d + 4*%d + 1 + %d - %d - %d = %d = n'
           % (n3, n4, m22, m33, m34, 3 * n3 + 4 * n4 + 1 + m22 - m33 - m34))
    print('')

print('--- FORCED NEGATIVES: the source\'s weaker bounds must SURVIVE where the paper says ----')
g_w1, g_w2, g_w3 = grm_direct(*W1), grm_direct(*W2), grm_direct(*W3)
ck('W1-respects-the-source-bound-at-l822', g_w1 >= -(12 + 2),
   'l.822 asks GRM_{-2} >= -(n+2) = %d and W1 gives %d' % (-(12 + 2), g_w1))
ck('W1-meets-the-source-bound-at-l830-with-equality', g_w1 == -(12 + 1),
   'l.830 asks GRM_{-2} >= -(n+1) = %d and W1 gives exactly that, so W1 refutes l.838 ALONE'
   % (-(12 + 1),))
ck('W2-also-violates-the-source-bound-at-l830', g_w2 < -(16 + 1),
   'l.830 asks >= %d at n=16 and W2 gives %d' % (-(16 + 1), g_w2))
ck('W2-meets-the-source-bound-at-l822-with-equality', g_w2 == -(16 + 2),
   'l.822 asks >= %d and W2 gives exactly that' % (-(16 + 2),))
ck('W3-meets-the-source-bound-at-l822-with-equality', g_w3 == -(20 + 2),
   'l.822 asks >= %d and W3 gives exactly that' % (-(20 + 2),))
print('')

print('--- THE FAMILY F(k) AS BUILT HERE: order 4k+4 and GRM_{-2} = -(n+2) CHECKED FOR k = 3..12')
print('    (the general k >= 3 claim is NOT machine-checked here) ---------------------------')
for k in range(3, 13):
    n, E = family(k)
    ok, why = is_tree(n, E)
    m, ni, d = profile(n, E)
    g = grm_direct(n, E)
    good = (ok and n == 4 * k + 4 and max(d) == 4 and g == -(n + 2)
            and ni.get(3, 0) == 1 and m.get((1, 3), 0) == 0 and m.get((2, 2), 0) == 0
            and m.get((3, 4), 0) == 0 and m.get((4, 4), 0) == 0 and m.get((1, 2), 0) == 0)
    ck('family-k%d-has-order-4k+4-Delta-4-and-GRM-minus-n-minus-2' % k, good,
       'k=%d: n=%d, Delta=%d, GRM_{-2}=%d = -(n+2), n_3=%d, m_14=%d  (%s)'
       % (k, n, max(d), g, ni.get(3, 0), m.get((1, 4), 0), why))
n3f, E3f = family(3)
n4f, E4f = family(4)
ck('family-k3-is-isomorphic-to-W2', canon(n3f, E3f) == canon(W2[0], W2[1]),
   'AHU canonical forms agree on %d vertices' % n3f)
ck('family-k4-is-isomorphic-to-W3', canon(n4f, E4f) == canon(W3[0], W3[1]),
   'AHU canonical forms agree on %d vertices' % n4f)
print('')

print('--- THE SOURCE\'S OWN CLAIMED EXTREMAL PARAMETERS AT l.843 AND l.844 EVALUATE TO -n ----')
bad1, bad2 = [], []
for k in range(2, 13):
    n = 4 * k + 4
    # l.843, first case: n1=2k+2, n2=k+2, n4=k, m22=3, m14=2k+2, m24=2k-2 (m13=m33=m34=m44=0)
    v1 = -(0 + 2 * (2 * k + 2) - 0 - 0 - 0)
    if v1 != -n:
        bad1.append((k, v1))
    # l.844, second case: n1=2k+4, n2=k-1, n4=k+1, m44=1, m14=2k+4, m24=2k-2
    v2 = -(0 + 2 * (2 * k + 4) - 0 - 0 - 4 * 1)
    if v2 != -n:
        bad2.append((k, v2))
ck('l843-first-claimed-extremal-case-has-GRM-equal-minus-n', not bad1,
   'checked k = 2..12 through the source formula at l.795; every case gives -n%s'
   % ('' if not bad1 else '   MISMATCH %s' % (bad1,)))
ck('l844-second-claimed-extremal-case-has-GRM-equal-minus-n', not bad2,
   'checked k = 2..12; every case gives -n, so TT^4_opt(k) sits 1 above W1 at k=2 and 2 above '
   'F(k) for k >= 3%s' % ('' if not bad2 else '   MISMATCH %s' % (bad2,)))
print('')

failed = [nm for nm, ok, _ in CHECKS if not ok]
print('SCOPE OF THIS RUN: it recomputes GRM_{-2} and the degree and edge counts of the graphs')
print('listed below from their edge lists. It does NOT enumerate trees, so it establishes no exact')
print('minimum at any order and identifies no minimisers. The statement of the family claim for')
print('all k >= 3 is not machine-checked here: only k = 3..12 are instantiated. This run does not')
print('decide the case n = 8 beyond the single double star D8.')
print('')
print('objects checked: W1 (n=12), W2 (n=16), W3 (n=20), D8 (n=8 anti-control), F(k) for k=3..12,')
print('and the source\'s claimed extremal parameter sets at l.843 and l.844.')
if failed:
    print('FAILED CHECKS: %s' % (failed,))
    raise SystemExit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(CHECKS))
raise SystemExit(0)
