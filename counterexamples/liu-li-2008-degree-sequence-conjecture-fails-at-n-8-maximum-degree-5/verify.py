#!/usr/bin/env python3
"""verify.py -- computational checks for the note on (n, Delta) = (8, 5)

    The statement of the Liu-Li 2008 conjecture is NOT re-derived here and is taken from the
    paper, Section 2 (its predicted degree sequence at (8,5) is assumed to be (5,5,5,5,5,5,5,3)).

Python 3.9+, STANDARD LIBRARY ONLY (sys, itertools, fractions).  No third-party package, no external
data file, no network.  All decisions are made in exact integer or exact rational arithmetic; no
floating point value is ever compared.

The program reads the objects TRANSCRIBED INTO THIS FILE FROM THE PAPER -- the edge lists of H and
of C1, C2, C3, the graph6 strings, the integer Rayleigh vector -- and re-derives the quantities
transcribed here about them.  The paper itself is not parsed, so completeness of that list against
the paper is not verified here.  It then re-establishes, by its own complete enumeration, the two class statements the paper
makes: that there are exactly 21000 labeled connected graphs on 8 vertices with degree sequence
(5,5,5,5,5,5,5,3) and that every one of them satisfies lambda_1 < 97/20, and that there are exactly
2520 labeled graphs with degree sequence (5,5,5,5,5,5,5,1), all connected, all with lambda_1 > 49/10.

One line per check.  Exits 0 iff every check passed.
"""

import itertools
import sys
from fractions import Fraction

N = 8

# ----------------------------------------------------------------------------------------------
# THE OBJECTS, EXACTLY AS PRINTED IN THE PAPER
# ----------------------------------------------------------------------------------------------

# The witness H, labelling L1 (Section 3.1 of the paper), graph6 GI[z~g
H1_G6 = 'GI[z~g'
H1_EDGES = [(0, 7), (1, 2), (1, 3), (1, 4), (1, 6), (1, 7), (2, 4), (2, 5), (2, 6), (2, 7),
            (3, 4), (3, 5), (3, 6), (3, 7), (4, 5), (4, 6), (5, 6), (5, 7)]

# The same graph in the second labelling L2 used by the Rayleigh vector, graph6 GNz~s?
H2_G6 = 'GNz~s?'
H2_EDGES = [(0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
            (2, 3), (2, 4), (2, 5), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6)]

# The integer Rayleigh vector of the paper, in labelling L2
RAYLEIGH_X = [3251590, 3914936, 3914936, 3818824, 3818824, 3818824, 3818824, 663346]
RAYLEIGH_Q = Fraction(61272483059723, 12499997522989)

# The three isomorphism-class representatives of the PREDICTED class, degree sequence (5^7, 3)
C1_G6 = 'GFz~v?'
C1_EDGES = [(0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7),
            (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (3, 5), (3, 6), (4, 5), (4, 6)]
C2_G6 = 'GLv~v?'
C2_EDGES = [(0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (1, 2), (1, 4), (1, 5), (1, 6), (1, 7),
            (2, 3), (2, 5), (2, 6), (2, 7), (3, 4), (3, 5), (3, 6), (4, 5), (4, 6)]
C3_G6 = 'GZn]~?'
C3_EDGES = [(0, 2), (0, 4), (0, 5), (0, 6), (0, 7), (1, 2), (1, 3), (1, 5), (1, 6), (1, 7),
            (2, 3), (2, 4), (2, 7), (3, 4), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)]

# Characteristic polynomials as printed, coefficients of x^8 .. x^0
P_H = [1, 0, -18, -32, 5, 32, 12, 0, 0]
P_C1 = [1, 0, -19, -24, 12, 0, 0, 0, 0]
P_C2 = [1, 0, -19, -28, 24, 36, -8, 0, 0]
P_C3 = [1, 0, -19, -30, 28, 62, -3, -24, 1]

# Exact evaluations as printed
P_H_AT_49_10 = Fraction(-27065510199, 100000000)
P_H_AT_5 = 16800
P_H_AT_99_20 = Fraction(197988807248001, 25600000000)
P_C_AT_49_10 = [Fraction(846855031701, 100000000),
                Fraction(813058555701, 100000000),
                Fraction(784932497701, 100000000)]

# The eight leading principal minors of the integer matrix 97 I - 20 A, as printed
MINORS_97 = [
    [97, 9409, 912673, 77238481, 6396925057, 355826560529, 8841684881313, 28053247092561],
    [97, 9409, 873873, 77398481, 6308461057, 349728752529, 8491067905313, 23557099988561],
    [97, 9409, 835073, 72082881, 5633183857, 364595360529, 8202660545313, 18809429316561],
]

# Exact rational brackets for the four spectral radii, as printed in the paper
BRACKETS = [
    ('H', Fraction(49017, 10000), Fraction(49018, 10000)),
    ('C1', Fraction(48419, 10000), Fraction(48420, 10000)),
    ('C2', Fraction(48431, 10000), Fraction(48432, 10000)),
    ('C3', Fraction(48444, 10000), Fraction(48445, 10000)),
]

# ----------------------------------------------------------------------------------------------
# CHECK HARNESS
# ----------------------------------------------------------------------------------------------

_passes = []
_fails = []


def check(name, ok, detail=''):
    if ok:
        _passes.append(name)
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _fails.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(text):
    print('NOTE %s' % text)


# ----------------------------------------------------------------------------------------------
# GRAPH PRIMITIVES
# ----------------------------------------------------------------------------------------------

def adj(n, edges):
    A = [[0] * n for _ in range(n)]
    for u, v in edges:
        A[u][v] = 1
        A[v][u] = 1
    return A


def degrees(n, edges):
    d = [0] * n
    for u, v in edges:
        d[u] += 1
        d[v] += 1
    return d


def is_connected(n, edges):
    if n == 0:
        return True
    nbr = [[] for _ in range(n)]
    for u, v in edges:
        nbr[u].append(v)
        nbr[v].append(u)
    seen = {0}
    stack = [0]
    while stack:
        x = stack.pop()
        for y in nbr[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == n


def components(verts, edgeset):
    """Connected components of the induced subgraph on `verts`. Returns a list of frozensets."""
    verts = set(verts)
    nbr = {v: [] for v in verts}
    for u, v in edgeset:
        if u in verts and v in verts:
            nbr[u].append(v)
            nbr[v].append(u)
    out = []
    todo = set(verts)
    while todo:
        s = todo.pop()
        seen = {s}
        stack = [s]
        while stack:
            x = stack.pop()
            for y in nbr[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        todo -= seen
        out.append(frozenset(seen))
    return out


def decode_graph6(s):
    """graph6 -> (n, sorted edge list). Upper triangle, column major: (0,1),(0,2),(1,2),(0,3),..."""
    n = ord(s[0]) - 63
    bits = ''.join(bin(ord(c) - 63)[2:].zfill(6) for c in s[1:])
    need = n * (n - 1) // 2
    if len(bits) < need:
        raise ValueError('graph6 string too short')
    k = 0
    E = []
    for j in range(1, n):
        for i in range(j):
            if bits[k] == '1':
                E.append((i, j))
            k += 1
    return n, E


def encode_graph6(n, edges):
    es = set()
    for u, v in edges:
        es.add((min(u, v), max(u, v)))
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append('1' if (i, j) in es else '0')
    while len(bits) % 6:
        bits.append('0')
    out = chr(n + 63)
    for i in range(0, len(bits), 6):
        out += chr(int(''.join(bits[i:i + 6]), 2) + 63)
    return out


def complement_edges(n, edges):
    es = set((min(u, v), max(u, v)) for u, v in edges)
    return [(i, j) for i in range(n) for j in range(i + 1, n) if (i, j) not in es]


# ----------------------------------------------------------------------------------------------
# EXACT LINEAR ALGEBRA
# ----------------------------------------------------------------------------------------------

def matmul(X, Y):
    n = len(X)
    Yt = list(zip(*Y))
    return [[sum(a * b for a, b in zip(row, col)) for col in Yt] for row in X]


def charpoly(A):
    """Coefficients of det(x I - A), x^n .. x^0, by Newton's identities on integer traces.

    Deliberately a DIFFERENT algorithm from a determinant expansion: the traces of A^k are counts of
    closed walks, so every intermediate quantity is a nonnegative integer with a combinatorial
    meaning, and the Newton divisions are exact or the program raises.
    """
    n = len(A)
    P = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    tr = []
    for _ in range(n):
        P = matmul(P, A)
        tr.append(sum(P[i][i] for i in range(n)))
    e = [1] + [0] * n
    for k in range(1, n + 1):
        s = 0
        for i in range(1, k + 1):
            s += ((-1) ** (i - 1)) * e[k - i] * tr[i - 1]
        q, r = divmod(s, k)
        if r != 0:
            raise ArithmeticError('Newton identity produced a non-integer elementary symmetric function')
        e[k] = q
    return [((-1) ** k) * e[k] for k in range(n + 1)]


def poly_eval(coeffs, x):
    """Horner, on coefficients x^n .. x^0. Exact for int or Fraction x."""
    acc = 0
    for c in coeffs:
        acc = acc * x + c
    return acc


def leading_minors_positive(M):
    """(all_positive, minors_computed) for a SYMMETRIC integer matrix, by fraction-free (Bareiss)
    elimination WITHOUT pivoting -- so the k-th pivot IS the k-th leading principal minor.

    Aborts at the first minor that is not strictly positive, which is exactly Sylvester's criterion
    failing.  For a positive definite matrix no pivot can vanish, so no division by zero arises on
    the path that returns True.
    """
    n = len(M)
    A = [row[:] for row in M]
    minors = []
    prev = 1
    for k in range(n):
        p = A[k][k]
        minors.append(p)
        if p <= 0:
            return False, minors
        if k == n - 1:
            break
        for i in range(k + 1, n):
            aik = A[i][k]
            for j in range(k + 1, n):
                num = A[i][j] * p - aik * A[k][j]
                q, r = divmod(num, prev)
                if r != 0:
                    raise ArithmeticError('Bareiss elimination produced a non-exact division')
                A[i][j] = q
        prev = p
    return True, minors


def det_int(M):
    """Exact integer determinant, Bareiss with row pivoting, with an exact-rational fallback."""
    n = len(M)
    A = [row[:] for row in M]
    sign = 1
    prev = 1
    try:
        for k in range(n - 1):
            if A[k][k] == 0:
                for r in range(k + 1, n):
                    if A[r][k] != 0:
                        A[k], A[r] = A[r], A[k]
                        sign = -sign
                        break
                else:
                    return 0
            p = A[k][k]
            for i in range(k + 1, n):
                aik = A[i][k]
                for j in range(k + 1, n):
                    num = A[i][j] * p - aik * A[k][j]
                    q, r = divmod(num, prev)
                    if r != 0:
                        raise ArithmeticError('non-exact')
                    A[i][j] = q
            prev = p
        return sign * A[n - 1][n - 1]
    except ArithmeticError:
        return _det_fraction(M)


def _det_fraction(M):
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    det = Fraction(1)
    for k in range(n):
        if A[k][k] == 0:
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    A[k], A[r] = A[r], A[k]
                    det = -det
                    break
            else:
                return 0
        det *= A[k][k]
        inv = A[k][k]
        for i in range(k + 1, n):
            f = A[i][k] / inv
            if f:
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    assert det.denominator == 1
    return int(det)


def shifted_matrix(A, num, den):
    """The integer matrix num*I - den*A, whose positive definiteness is exactly lambda_1(A) < num/den."""
    n = len(A)
    return [[(num if i == j else 0) - den * A[i][j] for j in range(n)] for i in range(n)]


def lambda1_below(A, num, den):
    """Exact: True iff Sylvester's criterion certifies lambda_1(A) < num/den (num, den > 0)."""
    ok, _m = leading_minors_positive(shifted_matrix(A, num, den))
    return ok


def lambda1_above(A, num, den):
    """Exact, one-sided: det(num*I - den*A) < 0 forces an eigenvalue of A strictly above num/den.

    det(num*I - den*A) = den^n * p_A(num/den) and p_A is monic with all roots real, so a negative
    value at num/den puts the LARGEST root strictly above it.
    """
    n = len(A)
    d = det_int(shifted_matrix(A, num, den))
    if n % 2 == 0:
        return d < 0
    return d < 0  # den^n > 0, so the sign of the determinant is the sign of p_A(num/den)


# ----------------------------------------------------------------------------------------------
# ENUMERATION OF LABELED GRAPHS WITH A PRESCRIBED DEGREE FUNCTION
# ----------------------------------------------------------------------------------------------

def gen_by_degree(deg):
    """All labeled graphs on len(deg) vertices whose degree at vertex v is exactly deg[v].

    Vertices are processed in order; at vertex v every edge to a LOWER vertex has already been
    committed, so the residual demand rem[v] must be met exactly by neighbours above v.  Yields
    sorted edge tuples; each graph is produced once.
    """
    n = len(deg)
    rem = list(deg)
    edges = []

    def rec(v):
        if v == n:
            if all(r == 0 for r in rem):
                yield tuple(edges)
            return
        k = rem[v]
        if k == 0:
            yield from rec(v + 1)
            return
        cand = [u for u in range(v + 1, n) if rem[u] > 0]
        if len(cand) < k:
            return
        for S in itertools.combinations(cand, k):
            for u in S:
                rem[u] -= 1
            rem[v] = 0
            for u in S:
                edges.append((v, u))
            yield from rec(v + 1)
            del edges[len(edges) - k:]
            rem[v] = k
            for u in S:
                rem[u] += 1

    yield from rec(0)


def cycle_signature(n, edges):
    """For a graph all of whose degrees are 2 except one vertex v of even degree 2t: the multiset of
    the t cycle lengths glued at v, together with the multiset of the remaining component sizes.

    Removing v from its own component leaves t disjoint paths; a path on p vertices closes up with v
    into a cycle of length p+1.  Every other component has all degrees 2, hence is a single cycle,
    hence its size IS its length.  So the signature determines the graph up to isomorphism.
    """
    d = degrees(n, edges)
    hub = [v for v in range(n) if d[v] != 2]
    if len(hub) != 1:
        return None
    v = hub[0]
    es = [(a, b) for a, b in edges]
    comps = components(range(n), es)
    home = [c for c in comps if v in c][0]
    paths = components(home - {v}, [(a, b) for a, b in es if a != v and b != v])
    glued = tuple(sorted(len(p) + 1 for p in paths))
    others = tuple(sorted(len(c) for c in comps if c is not home))
    return (glued, others)


# ----------------------------------------------------------------------------------------------
print('verify.py -- computational checks for the note on (n, Delta) = (8, 5); the statement of the')
print('    Liu-Li 2008 conjecture is NOT re-derived here and is taken from the paper, Section 2')
print('    (its predicted degree sequence at (8,5) is assumed to be (5,5,5,5,5,5,5,3))')
print('exact integer / exact rational arithmetic only; standard library only')
print('python: %s' % sys.version.split()[0])
print('')

# ==============================================================================================
print('=== Step 1: the cell (n, Delta) = (8, 5) and what the conjecture predicts there')
# ==============================================================================================
n, Delta = 8, 5
check('ambient_clause_3_le_Delta_le_n_minus_2', 3 <= Delta <= n - 2, '3 <= 5 <= 6')
check('n_times_Delta_is_even', (n * Delta) % 2 == 0, 'n*Delta = 40')
predicted_delta = Delta - 2 if (n * Delta) % 2 == 0 else Delta - 1
check('predicted_minimum_degree_is_Delta_minus_2_equals_3', predicted_delta == 3,
      'delta = %d, predicted degree sequence (5,5,5,5,5,5,5,3)' % predicted_delta)
deg_sum = (n - 1) * Delta + predicted_delta
check('predicted_class_has_19_edges', deg_sum == 38 and deg_sum // 2 == 19,
      'degree sum %d, edges %d' % (deg_sum, deg_sum // 2))
comp_edges_count = n * (n - 1) // 2 - deg_sum // 2
comp_degs = tuple(sorted([n - 1 - Delta] * (n - 1) + [n - 1 - predicted_delta]))
check('predicted_class_complement_has_9_edges_and_degrees_4_and_seven_2s',
      comp_edges_count == 9 and comp_degs == (2, 2, 2, 2, 2, 2, 2, 4),
      'complement: %d edges, degrees %s' % (comp_edges_count, list(comp_degs)))
# parity: (n-1)*Delta = 35 is odd, so delta must be odd; non-regularity forces 1 <= delta <= 4
odd_deltas = [d for d in range(1, Delta) if ((n - 1) * Delta + d) % 2 == 0]
check('parity_leaves_only_delta_in_1_3_on_8_vertices', odd_deltas == [1, 3],
      'admissible delta values with an integral edge count: %s' % odd_deltas)

# ==============================================================================================
print('')
print('=== Step 2: the witness H, read from the edge list printed in the paper')
# ==============================================================================================
H1 = sorted(set((min(a, b), max(a, b)) for a, b in H1_EDGES))
check('H_edge_list_is_18_distinct_edges_no_loop',
      len(H1) == 18 and len(H1) == len(H1_EDGES) and all(a != b for a, b in H1),
      '%d distinct edges' % len(H1))
dH = degrees(N, H1)
check('H_degree_sequence_is_1_and_seven_5s', sorted(dH) == [1, 5, 5, 5, 5, 5, 5, 5], 'degrees %s' % dH)
check('H_maximum_degree_is_exactly_5_and_H_is_not_regular',
      max(dH) == 5 and min(dH) != max(dH), 'max %d, min %d' % (max(dH), min(dH)))
check('H_is_connected', is_connected(N, H1))
nH, EH = decode_graph6(H1_G6)
check('H_graph6_decodes_to_the_printed_edge_list',
      nH == 8 and sorted(EH) == H1, 'graph6 %s' % H1_G6)
check('H_graph6_re_encodes_to_the_same_string', encode_graph6(N, H1) == H1_G6,
      'encode -> %s' % encode_graph6(N, H1))
# structure: K_7 minus (P_3 + 2K_2), pendant on the P_3 centre
core = [v for v in range(N) if dH[v] == 5]
pend = [v for v in range(N) if dH[v] == 1]
missing = [(i, j) for i in core for j in core if i < j and (i, j) not in set(H1)]
check('H_core_is_the_seven_degree_5_vertices_and_the_pendant_is_unique',
      len(core) == 7 and len(pend) == 1, 'core %s, pendant %s' % (core, pend))
check('H_core_is_K7_minus_exactly_four_edges', len(missing) == 4, 'missing pairs %s' % missing)
mdeg = {}
for a, b in missing:
    mdeg[a] = mdeg.get(a, 0) + 1
    mdeg[b] = mdeg.get(b, 0) + 1
centre = [v for v, k in mdeg.items() if k == 2]
check('the_four_missing_edges_form_P3_plus_2K2',
      len(mdeg) == 7 and sorted(mdeg.values()) == [1, 1, 1, 1, 1, 1, 2] and len(centre) == 1
      and len(components(mdeg.keys(), missing)) == 3,
      'the four missing pairs span 7 vertices in 3 components with degrees %s; the P_3 centre is %s'
      % (sorted(mdeg.values()), centre))
nbr_of_pendant = [v for (a, b) in H1 for v in (a, b) if {a, b} & set(pend) and v not in pend]
check('the_pendant_hangs_on_the_P3_centre', nbr_of_pendant == centre,
      'pendant %d attached to %s; P_3 centre %s' % (pend[0], nbr_of_pendant, centre))

# second labelling
H2 = sorted(set((min(a, b), max(a, b)) for a, b in H2_EDGES))
d2 = degrees(N, H2)
check('H_second_labelling_is_18_edges_with_the_same_degree_sequence',
      len(H2) == 18 and sorted(d2) == [1, 5, 5, 5, 5, 5, 5, 5] and is_connected(N, H2),
      'degrees %s, connected' % d2)
n2, E2 = decode_graph6(H2_G6)
check('H_second_labelling_graph6_decodes_and_re_encodes',
      n2 == 8 and sorted(E2) == H2 and encode_graph6(N, H2) == H2_G6, 'graph6 %s' % H2_G6)

A_H1 = adj(N, H1)
A_H2 = adj(N, H2)
p_h1 = charpoly(A_H1)
p_h2 = charpoly(A_H2)
check('the_two_labellings_of_H_are_cospectral', p_h1 == p_h2, 'identical characteristic polynomials')
check('H_characteristic_polynomial_matches_the_paper', p_h1 == P_H, 'x^8..x^0 = %s' % p_h1)

# ==============================================================================================
print('')
print('=== Step 3: Fact 2 -- lambda_1(H) > 49/10, by two independent exact routes')
# ==============================================================================================
check('the_Rayleigh_vector_is_a_nonzero_integer_vector',
      all(isinstance(t, int) for t in RAYLEIGH_X) and any(t != 0 for t in RAYLEIGH_X),
      '8 integers, largest %d' % max(RAYLEIGH_X))
num = sum(RAYLEIGH_X[i] * A_H2[i][j] * RAYLEIGH_X[j] for i in range(N) for j in range(N))
den = sum(t * t for t in RAYLEIGH_X)
q = Fraction(num, den)
check('the_Rayleigh_quotient_equals_the_fraction_printed_in_the_paper', q == RAYLEIGH_Q,
      'x^T A x / x^T x = %d/%d' % (q.numerator, q.denominator))
check('the_Rayleigh_quotient_exceeds_49_over_10', q > Fraction(49, 10),
      '%d/%d > 49/10' % (q.numerator, q.denominator))
check('route_A_gives_lambda_1_of_H_above_49_over_10', q > Fraction(49, 10),
      "Rayleigh's principle: lambda_1 >= any quotient")
v_h = poly_eval(P_H, Fraction(49, 10))
check('p_H_at_49_over_10_equals_the_value_printed_in_the_paper', v_h == P_H_AT_49_10,
      'p_H(49/10) = %d/%d' % (v_h.numerator, v_h.denominator))
check('p_H_at_49_over_10_is_negative', v_h < 0, 'route B: a monic real-rooted polynomial negative '
                                                'at t has its largest root above t')
check('p_H_is_monic', P_H[0] == 1)
check('p_H_at_5_matches_the_paper', poly_eval(P_H, 5) == P_H_AT_5, 'p_H(5) = %d' % poly_eval(P_H, 5))
v99 = poly_eval(P_H, Fraction(99, 20))
check('p_H_at_99_over_20_matches_the_paper', v99 == P_H_AT_99_20,
      'p_H(99/20) = %d/%d > 0' % (v99.numerator, v99.denominator))

# ==============================================================================================
print('')
print('=== Step 4: the three representatives of the predicted class')
# ==============================================================================================
REPS = [('C1', C1_EDGES, C1_G6, P_C1), ('C2', C2_EDGES, C2_G6, P_C2), ('C3', C3_EDGES, C3_G6, P_C3)]
rep_polys = []
for idx, (nm, ed, g6, pp) in enumerate(REPS):
    E = sorted(set((min(a, b), max(a, b)) for a, b in ed))
    d = degrees(N, E)
    check('%s_is_19_edges_with_degree_sequence_5pow7_3_and_connected' % nm,
          len(E) == 19 and sorted(d) == [3, 5, 5, 5, 5, 5, 5, 5] and is_connected(N, E),
          'degrees %s' % d)
    nn, EE = decode_graph6(g6)
    check('%s_graph6_decodes_and_re_encodes' % nm,
          nn == 8 and sorted(EE) == E and encode_graph6(N, E) == g6, 'graph6 %s' % g6)
    A = adj(N, E)
    cp = charpoly(A)
    rep_polys.append(cp)
    check('%s_characteristic_polynomial_matches_the_paper' % nm, cp == pp, 'x^8..x^0 = %s' % cp)
    ok, minors = leading_minors_positive(shifted_matrix(A, 97, 20))
    check('%s_all_eight_leading_minors_of_97I_minus_20A_are_positive' % nm, ok and all(m > 0 for m in minors),
          'minors %s' % minors)
    check('%s_those_minors_are_exactly_the_ones_printed_in_the_paper' % nm, minors == MINORS_97[idx])
    val = poly_eval(cp, Fraction(49, 10))
    check('%s_value_at_49_over_10_matches_the_paper_and_is_positive' % nm,
          val == P_C_AT_49_10[idx] and val > 0, 'p(49/10) = %d/%d' % (val.numerator, val.denominator))
    # the reader's hand cross-check: the last minor is det(97I - 20A) = 20^8 * p(97/20)
    lhs = det_int(shifted_matrix(A, 97, 20))
    rhs = 20 ** 8 * poly_eval(cp, Fraction(97, 20))
    check('%s_det_97I_minus_20A_equals_20pow8_times_p_at_97_over_20' % nm,
          lhs == rhs == minors[-1], 'det = %d' % lhs)

# exact brackets, so the decimal values quoted in the paper are checked and not merely displayed
BR_A = {'H': A_H1, 'C1': adj(N, C1_EDGES), 'C2': adj(N, C2_EDGES), 'C3': adj(N, C3_EDGES)}
for nm, lo, hi in BRACKETS:
    A = BR_A[nm]
    above = lambda1_above(A, lo.numerator, lo.denominator)
    below = lambda1_below(A, hi.numerator, hi.denominator)
    check('lambda_1_of_%s_lies_strictly_between_the_two_rationals_printed' % nm, above and below,
          '%s < lambda_1(%s) < %s' % (lo, nm, hi))
check('the_three_representatives_are_strictly_ordered_C1_lt_C2_lt_C3',
      BRACKETS[1][2] <= BRACKETS[2][1] and BRACKETS[2][2] <= BRACKETS[3][1],
      'disjoint brackets, so lambda_1(C1) < lambda_1(C2) < lambda_1(C3) among these three '
      'representatives; that these three exhaust the conjecture-allowed class is established '
      'separately in Step 5, and the conjecture\'s own statement is not verified here')

# ==============================================================================================
print('')
print('=== Step 5: Fact 3 -- the COMPLETE predicted class, by our own enumeration')
# ==============================================================================================
print('NOTE enumerating every labeled graph on 8 vertices with degree sequence (5,5,5,5,5,5,5,3)')
print('NOTE     by generating its 9-edge complement, whose degrees are (4,2,2,2,2,2,2,2)')
pred_total = 0
pred_connected = 0
viol97 = 0
viol49 = 0
sig_counts = {}
rep_sigs = {}
for hub in range(N):
    dseq = [4 if v == hub else 2 for v in range(N)]
    for cedges in gen_by_degree(dseq):
        pred_total += 1
        E = complement_edges(N, cedges)
        d = degrees(N, E)
        if sorted(d) != [3, 5, 5, 5, 5, 5, 5, 5]:
            viol97 += 1  # cannot happen; counted so a silent mis-generation shows up as a violation
            continue
        if is_connected(N, E):
            pred_connected += 1
        else:
            continue
        A = adj(N, E)
        if not lambda1_below(A, 97, 20):
            viol97 += 1
        if not lambda1_below(A, 49, 10):
            viol49 += 1
        sg = cycle_signature(N, list(cedges))
        sig_counts[sg] = sig_counts.get(sg, 0) + 1

check('the_predicted_class_has_exactly_21000_labeled_members', pred_total == 21000,
      'labeled graphs with degree sequence (5^7,3): %d' % pred_total)
check('every_member_of_the_predicted_class_is_connected', pred_connected == pred_total,
      'connected: %d of %d' % (pred_connected, pred_total))
check('THE_DECISION_every_one_of_the_21000_satisfies_Sylvester_on_97I_minus_20A',
      viol97 == 0 and pred_total == 21000,
      'violations = %d, i.e. lambda_1(G) < 97/20 = 4.85 for EVERY G in the predicted class' % viol97)
check('the_weaker_matrix_49I_minus_10A_also_has_zero_violations', viol49 == 0,
      'violations = %d (this matrix certifies only lambda_1 < 49/10 = 4.90)' % viol49)

want_sigs = {((3, 3), (3,)): 840, ((3, 6), ()): 10080, ((4, 5), ()): 10080}
check('the_class_splits_into_exactly_three_isomorphism_types', len(sig_counts) == 3,
      'complement cycle types found: %s' % sorted(sig_counts))
check('the_three_orbit_sizes_are_840_10080_10080_and_sum_to_21000',
      sig_counts == want_sigs and sum(sig_counts.values()) == 21000,
      '%s' % {str(k): v for k, v in sorted(sig_counts.items())})
check('the_hand_orbit_formulas_reproduce_those_three_counts',
      8 * 35 * 3 == 840 and 8 * 21 * 60 == 10080 and 8 * 35 * 3 * 12 == 10080
      and 840 + 10080 + 10080 == 21000,
      '8*C(7,4)*3 = 840; 8*C(7,2)*(5!/2) = 10080; 8*C(7,3)*(3!/2)*(4!/2) = 10080')
for nm, ed, _g6, _pp in REPS:
    rep_sigs[nm] = cycle_signature(N, complement_edges(N, sorted(set((min(a, b), max(a, b)) for a, b in ed))))
check('C1_C2_C3_realise_the_three_distinct_types',
      len(set(rep_sigs.values())) == 3 and set(rep_sigs.values()) == set(want_sigs),
      '%s' % {k: str(v) for k, v in sorted(rep_sigs.items())})
check('C1_is_the_type_of_orbit_size_840', sig_counts[rep_sigs['C1']] == 840,
      'C1 complement type %s has orbit 840, consistent with |Aut| = 8!/840 = 48'
      % (rep_sigs['C1'],))

# ==============================================================================================
print('')
print('=== Step 6: the witness class (5,5,5,5,5,5,5,1) is a single isomorphism class')
# ==============================================================================================
wit_total = 0
wit_connected = 0
wit_sigs = {}
wit_above = 0
for hub in range(N):
    dseq = [6 if v == hub else 2 for v in range(N)]
    for cedges in gen_by_degree(dseq):
        wit_total += 1
        E = complement_edges(N, cedges)
        if sorted(degrees(N, E)) != [1, 5, 5, 5, 5, 5, 5, 5]:
            continue
        if is_connected(N, E):
            wit_connected += 1
        A = adj(N, E)
        if lambda1_above(A, 49, 10):
            wit_above += 1
        sg = cycle_signature(N, list(cedges))
        wit_sigs[sg] = wit_sigs.get(sg, 0) + 1

check('the_witness_class_has_exactly_2520_labeled_members', wit_total == 2520,
      'labeled graphs with degree sequence (5^7,1): %d' % wit_total)
check('every_member_of_the_witness_class_is_connected', wit_connected == wit_total,
      'connected: %d of %d' % (wit_connected, wit_total))
check('the_witness_class_is_a_SINGLE_isomorphism_type', len(wit_sigs) == 1 and list(wit_sigs.values()) == [2520],
      'complement cycle type %s, orbit 2520' % (list(wit_sigs)[0],))
check('2520_equals_8_factorial_over_16_so_the_automorphism_group_has_order_16',
      2520 * 16 == 40320, '8!/|Aut| = 40320/16 = 2520')
check('every_one_of_the_2520_has_lambda_1_above_49_over_10', wit_above == wit_total,
      '%d of %d have det(49I - 10A) < 0' % (wit_above, wit_total))
check('H_lies_in_that_class_hence_H_is_forced_by_its_degree_sequence',
      sorted(degrees(N, H1)) == [1, 5, 5, 5, 5, 5, 5, 5]
      and cycle_signature(N, complement_edges(N, H1)) in wit_sigs,
      'H is the unique (5^7,1) graph up to isomorphism, so it is not a discovery')

# ==============================================================================================
print('')
print('=== Step 7: controls -- the certifier must fire and must stay silent')
# ==============================================================================================
K8 = [(i, j) for i in range(8) for j in range(i + 1, 8)]
K6K2 = [(i, j) for i in range(6) for j in range(i + 1, 6)] + [(6, 7)]
P8 = [(i, i + 1) for i in range(7)]
C8 = [(i, (i + 1) % 8) for i in range(8)]
S8 = [(0, j) for j in range(1, 8)]
check('control_K8_is_REFUSED_by_Sylvester_at_97_over_20', not lambda1_below(adj(8, K8), 97, 20),
      'lambda_1(K_8) = 7 > 4.85, so the certifier must and does refuse it')
check('control_K6_plus_K2_is_disconnected_so_it_is_outside_the_class',
      not is_connected(8, K6K2) and max(degrees(8, K6K2)) == 5,
      'max degree 5, non-regular, lambda_1 = 5 exactly -- the anti-control a broken connectivity '
      'filter would report as a spectacular refutation')
check('control_K6_plus_K2_is_also_REFUSED_by_Sylvester_at_97_over_20',
      not lambda1_below(adj(8, K6K2), 97, 20), 'lambda_1 = 5 > 4.85')
check('control_P8_is_certified_at_97_over_20', lambda1_below(adj(8, P8), 97, 20), 'lambda_1(P_8) < 2')
check('control_C8_is_certified_at_97_over_20', lambda1_below(adj(8, C8), 97, 20), 'lambda_1(C_8) = 2')
check('control_K1_7_is_certified_at_97_over_20', lambda1_below(adj(8, S8), 97, 20),
      'lambda_1(K_{1,7}) = sqrt(7) = 2.6457...')
check('control_the_certifier_REFUSES_H_at_97_over_20', not lambda1_below(A_H1, 97, 20),
      'the certifier failing in the right direction on the right object')
check('control_the_lower_bound_test_stays_silent_on_C3_at_49_over_10',
      not lambda1_above(adj(N, C3_EDGES), 49, 10), 'C3 has p(49/10) > 0, so no root above 49/10 is claimed')

# ==============================================================================================
print('')
print('=== Step 8: the conclusion, assembled from the exact facts above')
# ==============================================================================================
check('97_over_20_is_strictly_below_49_over_10', Fraction(97, 20) < Fraction(49, 10),
      '4.85 < 4.90, so the two bounds separate strictly')
check('H_beats_every_graph_of_the_predicted_degree_sequence',
      q > Fraction(49, 10) > Fraction(97, 20) and viol97 == 0 and pred_total == 21000,
      'lambda_1(H) > 49/10 > 97/20 > lambda_1(G) for all 21000 labeled G with degseq (5^7,3)')
check('therefore_no_maximiser_at_8_5_can_carry_the_predicted_degree_sequence',
      pred_connected == 21000 and viol97 == 0 and q > Fraction(97, 20),
      'the family is finite and non-empty (H is in it), so a maximiser M exists and '
      'lambda_1(M) >= lambda_1(H) > 97/20 > lambda_1(G) for every predicted G')
check('reading_robustness_a_graph_of_maximum_degree_at_most_4_cannot_compete',
      Fraction(4) < Fraction(97, 20),
      'lambda_1(G) <= Delta(G) <= 4 < 4.85 < lambda_1(H), so "maximum degree at most Delta" '
      'changes nothing')

print('')
note('SCOPE. This program re-derives (i) the quantities transcribed into this file from the paper '
     'about H, C1, C2 and C3 -- their edge lists, graph6 strings, characteristic polynomials, the '
     'leading minors of 97I - 20A, the exact evaluations at 49/10, 5 and 99/20, the Rayleigh '
     'quotient, and the four spectral brackets (NOT RE-RUN: the paper itself is not parsed, so '
     'completeness of this list against the paper is not verified here), (ii) the completeness of the predicted class at (8,5) and the '
     'strict bound lambda_1 < 97/20 over ALL 21000 of its labeled members, and (iii) the '
     'single-isomorphism-class structure of the witness class. It establishes nothing about any '
     'cell other than (8,5) -- in particular nothing here concerns Delta = n-2, which belongs to '
     'Huang-Liu-Yang. Which graph ATTAINS lambda(8,5) is not certified here and is not needed: the '
     'refutation needs only that H beats every graph of the predicted degree sequence.')

print('')
if _fails:
    print('VERDICT: %d of %d CHECKS FAILED: %s' % (len(_fails), len(_passes) + len(_fails),
                                                   ', '.join(_fails)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(_passes))
sys.exit(0)
