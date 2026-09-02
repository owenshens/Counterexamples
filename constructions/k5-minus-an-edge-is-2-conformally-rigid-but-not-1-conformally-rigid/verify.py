#!/usr/bin/env python3
"""verify.py -- exact re-derivation of the computational claims of

    "K_5 minus an edge is 2-conformally rigid but not 1-conformally rigid"

Python 3.9+, STANDARD LIBRARY ONLY (fractions, itertools, sys).  No numpy, no sympy,
no external data file.  Every decision is made in exact integer or Fraction
arithmetic; no floating-point number is ever compared.

The objects below -- the graph6 string, the explicit edge list on {1,...,5}, the integer
matrix M, and the weight w* -- were TRANSCRIBED BY HAND from the paper and are hard-coded
here; the program never reads paper.tex, so the faithfulness of that transcription is not
machine-checked.  From those hard-coded objects it re-derives the spectra, the certificate
properties and the sums.  Some checks below go BEYOND the paper (see the closing SCOPE
note); a check that establishes more than the paper claims is deliberate, not a claim
about the paper.

One `PASS <name> [detail]` line per check; closes with

    VERDICT: ALL <n> CHECKS PASS

and exits 0 iff every check passed.
"""

import itertools
import sys
from fractions import Fraction as F

# =============================================================================
# 0.  THE OBJECTS AS PRINTED IN THE PAPER
# =============================================================================

N = 5                                          # paper, Theorem 1 (Section 1): n = 5
GRAPH6 = 'D^{'                                 # paper, Section 3 (a relabelling; see below)

# The labelled graph of the paper: K_5 on V = {1,2,3,4,5} with the pair {4,5} deleted.
MISSING_PAIR = (4, 5)
EDGES = [(1, 2), (1, 3), (1, 4), (1, 5),
         (2, 3), (2, 4), (2, 5),
         (3, 4), (3, 5)]

# The bijection printed in the paper between the graph6 vertices 0..4 and V = 1..5.
G6_TO_PAPER = {0: 4, 1: 5, 2: 1, 3: 2, 4: 3}

# The integer matrix of the paper, Proposition 2, equation (3):   X = M/45.
M_INT = [[22, -3, -3, -8, -8],
         [-3, 22, -3, -8, -8],
         [-3, -3, 22, -8, -8],
         [-8, -8, -8, 12, 12],
         [-8, -8, -8, 12, 12]]
X_DEN = 45

# The weight of the paper, Proposition 4 (Section 2):  9/7 on the six edges meeting {4,5},
# 3/7 on the three edges inside {1,2,3}.
def w_star(e):
    a, b = e
    return F(9, 7) if (a in MISSING_PAIR or b in MISSING_PAIR) else F(3, 7)


# =============================================================================
# 1.  A TINY EXACT LINEAR-ALGEBRA / POLYNOMIAL KIT (stdlib only)
# =============================================================================

def zeros(n, m=None):
    m = n if m is None else m
    return [[F(0) for _ in range(m)] for _ in range(n)]


def eye(n):
    A = zeros(n)
    for i in range(n):
        A[i][i] = F(1)
    return A


def mat(A):
    """Coerce a list of lists of ints/Fractions to a Fraction matrix."""
    return [[F(x) for x in row] for row in A]


def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    C = zeros(n, m)
    for i in range(n):
        Ai = A[i]
        for j in range(m):
            s = F(0)
            for t in range(k):
                s += Ai[t] * B[t][j]
            C[i][j] = s
    return C


def mvec(A, v):
    return [sum((A[i][j] * v[j] for j in range(len(v))), F(0)) for i in range(len(A))]


def madd(A, B, cb=F(1)):
    return [[A[i][j] + cb * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def smul(c, A):
    return [[c * x for x in row] for row in A]


def trace(A):
    return sum((A[i][i] for i in range(len(A))), F(0))


def is_sym(A):
    n = len(A)
    return all(A[i][j] == A[j][i] for i in range(n) for j in range(n))


def is_zero_mat(A):
    return all(x == 0 for row in A for x in row)


def det(A):
    """Exact determinant by Fraction Gaussian elimination."""
    A = [row[:] for row in A]
    n = len(A)
    d = F(1)
    for c in range(n):
        p = None
        for r in range(c, n):
            if A[r][c] != 0:
                p = r
                break
        if p is None:
            return F(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        inv = F(1) / A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] * inv
            if f:
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return d


# --- polynomials as coefficient lists, LOWEST degree first --------------------

def pmul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                out[i + j] += a * b
    return out


def pprod(ps):
    out = [F(1)]
    for p in ps:
        out = pmul(out, p)
    return out


def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def charpoly(A):
    """Faddeev--LeVerrier, exact.  Returns det(xI - A) as coefficients, lowest first."""
    n = len(A)
    c = [F(1)]                      # c[0] = 1 is the coefficient of x^n
    Mk = zeros(n)
    for k in range(1, n + 1):
        Mk = mmul(A, madd(Mk, eye(n), c[k - 1]))
        c.append(-trace(Mk) / k)
    # p(x) = x^n + c1 x^{n-1} + ... + cn   -> lowest-first
    return list(reversed(c))


def linear(root):
    """(x - root) as a lowest-first coefficient list."""
    return [-F(root), F(1)]


def roots_from_factored(pairs):
    """pairs = [(root, multiplicity), ...] -> the monic polynomial, and the multiset."""
    ps = []
    spec = []
    for r, m in pairs:
        for _ in range(m):
            ps.append(linear(r))
            spec.append(F(r))
    return pprod(ps), sorted(spec)


# --- graph helpers ------------------------------------------------------------

def laplacian(edges, weight=None):
    L = zeros(N)
    for e in edges:
        a, b = e
        w = F(1) if weight is None else weight(e)
        i, j = a - 1, b - 1
        L[i][i] += w
        L[j][j] += w
        L[i][j] -= w
        L[j][i] -= w
    return L


def Lstar(X, e):
    """(L^*(X))_e = X_aa + X_bb - 2 X_ab, so that tr(L(w)X) = <w, L^*(X)>."""
    a, b = e[0] - 1, e[1] - 1
    return X[a][a] + X[b][b] - 2 * X[a][b]


def decode_graph6(s):
    """Minimal graph6 decoder (n < 63).  Returns (n, set of 0-based edges)."""
    n = ord(s[0]) - 63
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        if not (0 <= v < 64):
            raise ValueError('bad graph6 byte %r' % ch)
        bits.extend((v >> k) & 1 for k in range(5, -1, -1))
    E = set()
    t = 0
    for j in range(1, n):                      # column order: (0,j),(1,j),...,(j-1,j)
        for i in range(j):
            if bits[t]:
                E.add((i, j))
            t += 1
    return n, E


# =============================================================================
# 2.  THE CHECK HARNESS
# =============================================================================

_PASSES = []
_FAILS = []


def check(name, ok, detail=''):
    if ok:
        _PASSES.append(name)
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILS.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(text):
    print('NOTE %s' % text)


def hr(title):
    print()
    print('=== %s' % title)


# =============================================================================
print('exact verification of "K_5 minus an edge is 2-conformally rigid but not')
print('1-conformally rigid" -- all arithmetic in Python integers and Fractions')
print('python %s' % sys.version.split()[0])

# -----------------------------------------------------------------------------
hr('Step 1: the object exhibited in the paper')
# -----------------------------------------------------------------------------
g6n, g6E = decode_graph6(GRAPH6)
check('graph6_decodes_to_five_vertices', g6n == 5, 'graph6 %r -> n = %d' % (GRAPH6, g6n))
check('graph6_has_nine_edges', len(g6E) == 9, '|E| = %d' % len(g6E))
allpairs5 = set(itertools.combinations(range(5), 2))
comp = allpairs5 - g6E
check('graph6_graph_is_K5_minus_exactly_one_edge', len(comp) == 1,
      'complement edge = %s (0-based)' % sorted(comp))

E = [tuple(sorted(e)) for e in EDGES]
check('printed_edge_list_has_nine_distinct_edges', len(set(E)) == 9 and len(E) == 9,
      '|E| = 9')
want = set(itertools.combinations(range(1, 6), 2)) - {MISSING_PAIR}
check('printed_edge_list_is_K5_on_1_to_5_minus_the_pair_4_5', set(E) == want,
      'missing pair = {4,5}')

mapped = set(tuple(sorted((G6_TO_PAPER[i], G6_TO_PAPER[j]))) for (i, j) in g6E)
check('printed_relabelling_carries_graph6_onto_the_paper_labelling', mapped == set(E),
      'sigma = %s' % {k: v for k, v in sorted(G6_TO_PAPER.items())})

deg = {v: sum(1 for e in E if v in e) for v in range(1, 6)}
check('degree_sequence_is_3_3_4_4_4', sorted(deg.values()) == [3, 3, 4, 4, 4],
      'degrees %s' % [deg[v] for v in range(1, 6)])

# -----------------------------------------------------------------------------
hr('Step 2: the unweighted Laplacian and the sums S_k(1), s_k(1)')
# -----------------------------------------------------------------------------
L = laplacian(E)
check('laplacian_is_symmetric', is_sym(L))
check('laplacian_has_zero_row_sums', all(sum(row, F(0)) == 0 for row in L))
check('laplacian_trace_is_2E_equals_18', trace(L) == 18, 'tr L = %s' % trace(L))

cp_L = ptrim(charpoly(L))
tgt_L, spec_L = roots_from_factored([(0, 1), (3, 1), (5, 3)])
check('laplacian_charpoly_is_x_times_x_minus_3_times_x_minus_5_cubed',
      cp_L == ptrim(tgt_L), 'coefficients (lowest first) %s' % [str(c) for c in cp_L])
check('laplacian_spectrum_is_0_3_5_5_5', spec_L == [F(0), F(3), F(5), F(5), F(5)],
      'spec L(1) = %s' % [str(x) for x in spec_L])

ones = [F(1)] * 5
d45 = [F(0), F(0), F(0), F(1), F(-1)]                    # e_4 - e_5
a_v = [F(1), F(-1), F(0), F(0), F(0)]
b_v = [F(0), F(1), F(-1), F(0), F(0)]
c_v = [F(2), F(2), F(2), F(-3), F(-3)]

check('allones_is_the_trivial_eigenvector', mvec(L, ones) == [F(0)] * 5)
check('e4_minus_e5_is_an_eigenvector_for_eigenvalue_3',
      mvec(L, d45) == [3 * x for x in d45], 'L(e4-e5) = 3(e4-e5)')
for nm, v in (('a', a_v), ('b', b_v), ('c', c_v)):
    check('vector_%s_is_an_eigenvector_for_eigenvalue_5' % nm,
          mvec(L, v) == [5 * x for x in v], '%s = %s' % (nm, [str(x) for x in v]))
B = [[ones[i], d45[i], a_v[i], b_v[i], c_v[i]] for i in range(5)]
check('the_five_exhibited_eigenvectors_are_linearly_independent', det(B) != 0,
      'det = %s' % det(B))


def S_k(spec, k):
    """Sum of the k largest NONTRIVIAL eigenvalues (the trivial 0 is dropped once)."""
    nz = sorted(spec)[1:]
    return sum(nz[len(nz) - k:], F(0))


def s_k(spec, k):
    """Sum of the k smallest nontrivial eigenvalues."""
    nz = sorted(spec)[1:]
    return sum(nz[:k], F(0))


for k, val in ((1, 5), (2, 10), (3, 15), (4, 18)):
    check('S_%d_allones_equals_%d' % (k, val), S_k(spec_L, k) == val,
          'S_%d(1) = %s' % (k, S_k(spec_L, k)))
check('s_1_allones_equals_lambda_2_equals_3', s_k(spec_L, 1) == 3,
      's_1(1) = %s' % s_k(spec_L, 1))
check('S_4_allones_equals_the_trace_so_k_equals_n_minus_1_is_trivial',
      S_k(spec_L, 4) == trace(L) == 18, 'S_4 = s_4 = tr L = 18')

# -----------------------------------------------------------------------------
hr('Step 3: half (i) -- the dual certificate X = M/45 proves upper-2 rigidity')
# -----------------------------------------------------------------------------
X = smul(F(1, X_DEN), mat(M_INT))
check('X_is_symmetric', is_sym(X))
check('trace_X_equals_2_equals_k', trace(X) == 2, 'tr X = %s' % trace(X))
check('X_annihilates_allones', mvec(X, ones) == [F(0)] * 5, 'X.1 = 0')

cp_X = ptrim(charpoly(X))
tgt_X, spec_X = roots_from_factored([(0, 2), (F(5, 9), 2), (F(8, 9), 1)])
check('X_charpoly_is_x2_times_x_minus_5_9_squared_times_x_minus_8_9',
      cp_X == ptrim(tgt_X), 'coefficients (lowest first) %s' % [str(c) for c in cp_X])
check('X_spectrum_is_0_0_5_9_5_9_8_9', spec_X == sorted([F(0), F(0), F(5, 9), F(5, 9), F(8, 9)]),
      'spec X = %s' % [str(x) for x in spec_X])
check('zero_le_X_le_I_holds_EXACTLY_not_to_a_tolerance',
      all(F(0) <= mu <= F(1) for mu in spec_X), 'every eigenvalue of X lies in [0,1]')

for nm, v, mu in (('a', a_v, F(5, 9)), ('b', b_v, F(5, 9)), ('c', c_v, F(8, 9)),
                  ('e4-e5', d45, F(0)), ('allones', ones, F(0))):
    check('X_acts_on_%s_by_%s' % (nm.replace('-', '_minus_'), str(mu).replace('/', '_over_')),
          mvec(X, v) == [mu * x for x in v], 'X %s = (%s) %s' % (nm, mu, nm))

check('L_times_X_equals_5X_so_the_range_of_X_lies_in_the_5_eigenspace',
      is_zero_mat(madd(mmul(L, X), X, F(-5))), 'L X - 5 X = 0')
check('trace_LX_equals_10_equals_S_2_allones',
      trace(mmul(L, X)) == 10 == S_k(spec_L, 2), 'tr(L X) = %s' % trace(mmul(L, X)))

ls = [Lstar(X, e) for e in E]
check('Lstar_X_is_the_constant_10_over_9_on_all_nine_edges',
      set(ls) == {F(10, 9)} and len(ls) == 9, 'L*(X)_e = %s for every e' % ls[0])
check('Lstar_X_edge_sum_is_10_equals_S_2_allones', sum(ls, F(0)) == 10,
      '9 * 10/9 = %s' % sum(ls, F(0)))

# A closed form for X that the paper does NOT print:  X = (5/9) P + (1/3) c c^T / 30, where
# P is the orthogonal projector onto the 5-eigenspace of L.  The paper prints only M and the
# eigen-facts of Proposition 2(2) (X.1 = Xd = 0, Xa = (5/9)a, Xb = (5/9)b, Xc = (8/9)c); this
# repackaging, and P itself, appear nowhere in the paper.  The next three checks therefore
# establish more than the paper claims; they are kept because they explain where M came from.
J5 = [[F(1)] * 5 for _ in range(5)]
dd = [[d45[i] * d45[j] for j in range(5)] for i in range(5)]
P = madd(madd(eye(5), J5, F(-1, 5)), dd, F(-1, 2))
cc = [[c_v[i] * c_v[j] for j in range(5)] for i in range(5)]
X_closed = madd(smul(F(5, 9), P), cc, F(1, 90))
check('closed_form_X_equals_five_ninths_P_plus_c_cT_over_90',
      is_zero_mat(madd(X, X_closed, F(-1))),
      'P = I - J/5 - (e4-e5)(e4-e5)^T/2,  |c|^2 = 30,  (1/3)/30 = 1/90; '
      'this closed form is NOT printed in the paper -- the paper gives M and the eigen-facts')
check('P_is_a_symmetric_idempotent_of_trace_3',
      is_sym(P) and is_zero_mat(madd(mmul(P, P), P, F(-1))) and trace(P) == 3,
      'tr P = %s' % trace(P))
check('L_times_P_equals_5P_so_P_projects_onto_the_5_eigenspace',
      is_zero_mat(madd(mmul(L, P), P, F(-5))))

# the pairing identity tr(L(w) X) = <w, L*(X)>, checked on exact rational weights
TESTW = [
    ('allones', lambda e: F(1)),
    ('w_star', w_star),
    ('ramp', lambda e: F(E.index(tuple(sorted(e))) + 1, 4)),
    ('lopsided', lambda e: F(7, 2) if tuple(sorted(e)) == (1, 2) else F(11, 16)),
    ('zero_on_a_triangle', lambda e: F(0) if set(e) <= {1, 2, 3} else F(3, 2)),
]
pair_ok = True
for nm, wf in TESTW:
    lhs = trace(mmul(laplacian(E, wf), X))
    rhs = sum((wf(e) * Lstar(X, e) for e in E), F(0))
    pair_ok = pair_ok and (lhs == rhs)
check('pairing_identity_trace_L_w_X_equals_w_dot_Lstar_X', pair_ok,
      'exact on %d rational weight vectors' % len(TESTW))
kyfan_ok = all(sum((wf(e) * Lstar(X, e) for e in E), F(0)) == F(10, 9) * sum((wf(e) for e in E), F(0))
               for _, wf in TESTW)
check('constancy_makes_the_bound_on_the_five_TESTW_vectors', kyfan_ok,
      'for each of the 5 TESTW vectors, <w, L*(X)> = (10/9) sum_e w_e')
note('Ky Fan: for X with 0 <= X <= I, tr X = k and X.1 = 0, tr(L(w)X) <= S_k(w).  '
     'The three checks above therefore give S_2(w) >= 10 = S_2(1) for every w in Delta_E, '
     'with equality at all-ones: K_5 - e is upper-2 rigid.')

# -----------------------------------------------------------------------------
hr('Step 4: half (ii) -- the weight w* kills lower-1 (equivalently upper-3).  The FACT that '
   'K_5 - e is not lower-1 rigid is not claimed as new by the paper: it is published '
   'numerical data (Niu, row id=450, lcr=False); only the exact rational w* is supplied here')
# -----------------------------------------------------------------------------
wvals = [w_star(e) for e in E]
check('w_star_is_nonnegative', all(v >= 0 for v in wvals),
      'values %s' % sorted(set(str(v) for v in wvals)))
check('w_star_sums_to_9_equals_E_so_it_lies_in_Delta_E', sum(wvals, F(0)) == 9,
      '6*(9/7) + 3*(3/7) = %s' % sum(wvals, F(0)))
check('w_star_is_9_over_7_on_the_six_edges_meeting_4_5_and_3_over_7_inside_1_2_3',
      sum(1 for e in E if w_star(e) == F(9, 7)) == 6
      and sum(1 for e in E if w_star(e) == F(3, 7)) == 3, '6 + 3 = 9 edges')

Lw = laplacian(E, w_star)
check('L_w_star_has_zero_row_sums', all(sum(row, F(0)) == 0 for row in Lw))
check('trace_L_w_star_equals_18_equals_2E', trace(Lw) == 18, 'tr L(w*) = %s' % trace(Lw))
cp_w = ptrim(charpoly(Lw))
tgt_w, spec_w = roots_from_factored([(0, 1), (F(27, 7), 3), (F(45, 7), 1)])
check('L_w_star_charpoly_is_x_times_x_minus_27_7_cubed_times_x_minus_45_7',
      cp_w == ptrim(tgt_w), 'coefficients (lowest first) %s' % [str(c) for c in cp_w])
check('L_w_star_spectrum_is_0_27_7_27_7_27_7_45_7',
      spec_w == sorted([F(0), F(27, 7), F(27, 7), F(27, 7), F(45, 7)]),
      'spec L(w*) = %s' % [str(x) for x in spec_w])
for nm, v, mu in (('a', a_v, F(27, 7)), ('b', b_v, F(27, 7)), ('e4-e5', d45, F(27, 7)),
                  ('c', c_v, F(45, 7))):
    check('w_star_eigenvector_%s_has_eigenvalue_%s'
          % (nm.replace('-', '_minus_'), str(mu).replace('/', '_over_')),
          mvec(Lw, v) == [mu * x for x in v], 'L(w*) %s = (%s) %s' % (nm, mu, nm))

check('s_1_w_star_is_27_over_7_strictly_greater_than_s_1_allones_equals_3',
      s_k(spec_w, 1) == F(27, 7) and s_k(spec_w, 1) > s_k(spec_L, 1),
      '27/7 = %s > 3' % s_k(spec_w, 1))
check('lower_1_rigidity_FAILS_so_K5_minus_e_is_NOT_1_conformally_rigid',
      s_k(spec_w, 1) > s_k(spec_L, 1),
      'all-ones does not maximise s_1 over Delta_E')
check('S_3_w_star_is_99_over_7_strictly_below_S_3_allones_equals_15',
      S_k(spec_w, 3) == F(99, 7) and S_k(spec_w, 3) < S_k(spec_L, 3),
      '99/7 = %s < 15' % S_k(spec_w, 3))
check('upper_3_rigidity_FAILS_by_the_same_weight', S_k(spec_w, 3) < S_k(spec_L, 3))
check('S_2_w_star_is_72_over_7_which_respects_the_certificate_bound_10',
      S_k(spec_w, 2) == F(72, 7) and S_k(spec_w, 2) >= 10,
      '72/7 = %s >= 10' % S_k(spec_w, 2))

dual_ok = True
for spec in (spec_L, spec_w):
    for k in range(1, N - 1):
        dual_ok = dual_ok and (S_k(spec, k) + s_k(spec, N - 1 - k) == 18)
check('duality_identity_S_k_plus_s_n_minus_1_minus_k_equals_2E', dual_ok,
      'checked for w = all-ones and w = w*, k = 1,2,3; hence upper-k iff lower-(n-1-k)')
check('at_n_5_and_k_2_the_duality_is_self_paired_so_upper_2_gives_2_conformal_rigidity',
      N - 1 - 2 == 2, 'n-1-k = 5-1-2 = 2 = k')

trivial_ok = all(S_k(spec, 4) == 18 for spec in (spec_L, spec_w))
check('S_4_is_the_constant_18_on_Delta_E_so_k_equals_4_is_trivially_rigid', trivial_ok,
      'S_4(w) = tr L(w) = 2|E| = 18 for every w in Delta_E')

# -----------------------------------------------------------------------------
hr('Step 5: material BEYOND THE PAPER -- the upper-1 certificate X/2, the values of L*(P), '
   'and the k-map at n = 5')
note('The five checks in this step correspond to nothing in the paper: the paper proves no '
     'upper-1 certificate (it cites Niu, ucr=True, for that fact), never mentions P, and '
     'states no k-map.  They are recomputations the program supplies on its own; they are '
     'consistent with the paper and establish more than it claims, and none of them is '
     'evidence for a claim the paper makes.')
# -----------------------------------------------------------------------------
Xh = smul(F(1, 2), X)
check('half_X_is_a_feasible_upper_1_certificate',
      trace(Xh) == 1 and mvec(Xh, ones) == [F(0)] * 5
      and all(F(0) <= mu / 2 <= F(1) for mu in spec_X),
      'tr(X/2) = 1, (X/2).1 = 0, spectrum {0,0,5/18,5/18,4/9} in [0,1]')
lsh = [Lstar(Xh, e) for e in E]
check('Lstar_half_X_is_the_constant_5_over_9_with_edge_sum_5_equals_S_1_allones',
      set(lsh) == {F(5, 9)} and sum(lsh, F(0)) == 5 == S_k(spec_L, 1),
      '9 * 5/9 = %s' % sum(lsh, F(0)))

lsP = {'inside_123': set(Lstar(P, e) for e in E if set(e) <= {1, 2, 3}),
       'meeting_45': set(Lstar(P, e) for e in E if set(e) & {4, 5})}
check('Lstar_P_is_2_inside_1_2_3_and_3_over_2_on_the_six_cross_edges_hence_NOT_constant',
      lsP['inside_123'] == {F(2)} and lsP['meeting_45'] == {F(3, 2)},
      '3*2 + 6*(3/2) = %s' % sum((Lstar(P, e) for e in E), F(0)))
check('Lstar_P_edge_sum_is_15_equals_S_3_allones',
      sum((Lstar(P, e) for e in E), F(0)) == 15 == S_k(spec_L, 3))

kmap = {1: False, 2: True, 3: False, 4: True}
derived = {
    1: (s_k(spec_w, 1) <= s_k(spec_L, 1)),                 # lower-1 fails -> False
    2: True,                                               # certificate X + self-duality
    3: (S_k(spec_w, 3) >= S_k(spec_L, 3)),                 # upper-3 fails -> False
    4: (S_k(spec_L, 4) == S_k(spec_w, 4) == 18),           # trivial
}
check('k_conformal_rigidity_at_n_5_holds_exactly_for_k_in_2_and_4', derived == kmap,
      'k-map %s ; the non-trivial content is {1:F, 2:T, 3:F}, non-monotone in both directions'
      ' -- computed here from the two witnesses; the paper states no k-map, asserting only '
      '2-conformal rigidity and the failure of lower-1 (equivalently upper-3)'
      % {k: ('T' if v else 'F') for k, v in sorted(derived.items())})

# -----------------------------------------------------------------------------
print()
note('SCOPE -- what this program does NOT cover.  It re-derives, in exact arithmetic, the '
     'quantities named in the PASS lines above for the single graph K_5 - e.  Its inputs -- '
     'the graph6 string, the edge list, the matrix M and the weight w* -- are TRANSCRIBED BY '
     'HAND from the paper and hard-coded here: this program never reads paper.tex, so neither '
     'the faithfulness of that transcription nor the completeness of the list above against '
     'what the paper states is machine-checked, and "every quantity the paper states" is not '
     'something this program verifies.  NOT RE-RUN here: (a) the Ky Fan maximum principle and '
     'the source\'s upper-k/lower-(n-1-k) duality, which are quoted theorems, not '
     'recomputations; (b) minimality and uniqueness -- the paper claims NEITHER, asserting '
     'neither that K_5 - e is the smallest graph with this behaviour nor that it is the only '
     '5-vertex witness; the census of the 994 connected graphs on 3..7 vertices that looked at '
     'those questions was float64 at a bound tolerance near 1e-7 and is NOT reproduced here, '
     'so nothing here bears on them either; (c) every order n >= 8, which was never swept; '
     '(d) any infinite family.  BEYOND THE PAPER, in the other direction: all of Step 5 (the '
     'upper-1 certificate X/2, the values of L*(P), the k-map) and the closed form for X in '
     'Step 3 are established here but are not stated or claimed by the paper.  This program '
     'certifies the two witnesses and nothing else.')

print()
if _FAILS:
    print('VERDICT: %d CHECK(S) FAILED: %s' % (len(_FAILS), ', '.join(_FAILS)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(_PASSES))
sys.exit(0)
