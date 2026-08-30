#!/usr/bin/env python3
"""verify.py -- re-derives every quantity claimed in paper.tex, from the objects the paper prints.

Python 3.9+, STANDARD LIBRARY ONLY (itertools, fractions, sys).  No numpy, no sympy, no
networkx, no external data file, no eigenvalue routine, and no floating-point number takes
part in any decision: every comparison below is between exact integers or exact rationals.

What is read from the paper and not recomputed here:
  * the generator rule for r_k^eps (paper, Section 2);
  * the partition rule for V_{i,j} (paper, Section 2);
  * the eigenvector table v(d), d in Z_6 (paper, Theorem 1);
  * the matrix T and its blocks S, S' (paper, Section 4);
  * the separators t = 3, t = 47/10 and the table t_m (paper, Sections 3, 4, 5).
Everything else -- the graphs, the quotients, the characteristic polynomials, the eigenvalue
counts, the printed decimal digits -- is rebuilt here from those inputs.

Prints one PASS line per check and closes with a VERDICT line; exits 0 iff every check passed.
"""

import itertools
import sys
from fractions import Fraction

CHECKS = 0
FAILED = 0


def check(good, name, detail=""):
    global CHECKS, FAILED
    CHECKS += 1
    if not good:
        FAILED += 1
    print(("PASS  " if good else "FAIL  ") + name + (("  " + detail) if detail else ""))


# ----------------------------------------------------------------------------------------
# 1.  The graph P_m(n), exactly as the paper defines it.
# ----------------------------------------------------------------------------------------
# A vertex is a pair (w, c): w is a permutation word of [n] (position j carries letter w_j)
# and c in Z_m^n is the colour vector (position j carries colour c_j).  We store it as the
# tuple ((w_1, c_1), ..., (w_n, c_n)).
#
# r_k^eps reverses the prefix of length k and adds eps to the colour of every entry that sat
# in positions 1..k.  For m >= 3, R_m = {r_k^+, r_k^- : k in [n]}, so P_m(n) is 2n-regular.

def vertices(m, n):
    out = []
    for w in itertools.permutations(range(1, n + 1)):
        for c in itertools.product(range(m), repeat=n):
            out.append(tuple(zip(w, c)))
    return out


def r(sigma, k, eps, m):
    new = list(sigma)
    for j in range(k):
        letter, col = sigma[k - 1 - j]
        new[j] = (letter, (col + eps) % m)
    return tuple(new)


def pancake(m, n):
    V = vertices(m, n)
    index = {v: i for i, v in enumerate(V)}
    adj = []
    for v in V:
        adj.append([index[r(v, k, eps, m)] for k in range(1, n + 1) for eps in (1, -1)])
    return V, index, adj


def where(sigma, letter):
    """(position, colour) of `letter` in sigma; positions are 1-based."""
    for j, (lt, col) in enumerate(sigma):
        if lt == letter:
            return j + 1, col
    raise KeyError(letter)


def connected(adj):
    seen = {0}
    stack = [0]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == len(adj)


def basic_graph_checks(m, n, tag):
    V, index, adj = pancake(m, n)
    fact = 1
    for i in range(2, n + 1):
        fact *= i
    check(len(V) == m ** n * fact, "%s |V| = m^n n! = %d" % (tag, m ** n * fact),
          "built %d vertices" % len(V))
    check(all(len(set(a)) == 2 * n for a in adj) and all(i not in adj[i] for i in range(len(V))),
          "%s simple and %d-regular" % (tag, 2 * n),
          "every vertex has %d distinct neighbours, no loop" % (2 * n))
    check(all(i in adj[j] for i, a in enumerate(adj) for j in a),
          "%s adjacency is symmetric" % tag, "R_m is closed under inverses")
    check(connected(adj), "%s connected" % tag,
          "so lambda_1 = %d with multiplicity one" % (2 * n))
    return V, index, adj


# ----------------------------------------------------------------------------------------
# 2.  The partition P_{m,n} and its quotient, both rebuilt from the graph.
# ----------------------------------------------------------------------------------------
# V_{i,j} = { sigma : letter 1 sits at position j and carries colour i },  i in Z_m, j in [n].
# Cells are indexed i*n + (j-1).

def cell_labels(V, m, n):
    lab = []
    for sigma in V:
        j, col = where(sigma, 1)
        lab.append(col * n + (j - 1))
    return lab


def quotient(V, adj, m, n):
    """(Q, equitable, cell_sizes): Q[a][b] = # neighbours in cell b of a vertex in cell a."""
    lab = cell_labels(V, m, n)
    K = m * n
    rows = {}
    sizes = [0] * K
    equitable = True
    for i in range(len(V)):
        sizes[lab[i]] += 1
        cnt = [0] * K
        for j in adj[i]:
            cnt[lab[j]] += 1
        seen = rows.setdefault(lab[i], tuple(cnt))
        if seen != tuple(cnt):
            equitable = False
    return [list(rows[a]) for a in range(K)], equitable, sizes


def circ_form(m, n):
    """circ(2 D_n, E_n, O, ..., O, E_n) with D_n = diag(0..n-1), (E_n)_{ij} = [i+j <= n+1]."""
    twoD = [[2 * i if i == j else 0 for j in range(n)] for i in range(n)]
    E = [[1 if (i + 1) + (j + 1) <= n + 1 else 0 for j in range(n)] for i in range(n)]
    Z = [[0] * n for _ in range(n)]
    C = [twoD, E] + [Z] * (m - 3) + [E]
    K = m * n
    Q = [[0] * K for _ in range(K)]
    for a in range(m):
        for b in range(m):
            B = C[(b - a) % m]
            for j in range(n):
                for jp in range(n):
                    Q[a * n + j][b * n + jp] = B[j][jp]
    return Q


def kron_form(m, n):
    """I_m (x) 2 D_n  +  C_m (x) E_n, with C_m the adjacency matrix of the m-cycle."""
    twoD = [[2 * i if i == j else 0 for j in range(n)] for i in range(n)]
    E = [[1 if (i + 1) + (j + 1) <= n + 1 else 0 for j in range(n)] for i in range(n)]
    C = cycle(m)
    K = m * n
    Q = [[0] * K for _ in range(K)]
    for a in range(m):
        for j in range(n):
            for jp in range(n):
                Q[a * n + j][a * n + jp] += twoD[j][jp]
            for b in range(m):
                if C[a][b]:
                    for jp in range(n):
                        Q[a * n + j][b * n + jp] += C[a][b] * E[j][jp]
    return Q


def cycle(m):
    C = [[0] * m for _ in range(m)]
    for i in range(m):
        C[i][(i + 1) % m] += 1
        C[i][(i - 1) % m] += 1
    return C


# ----------------------------------------------------------------------------------------
# 3.  Exact linear algebra:  characteristic polynomials over Z, and eigenvalue COUNTS by
#     Sylvester's law of inertia.  No eigenvalue is ever approximated.
# ----------------------------------------------------------------------------------------

def charpoly(M):
    """Exact characteristic polynomial of an integer matrix (Faddeev-LeVerrier over Q;
    the result is integral).  Returns [a_0, ..., a_N] for a_0 x^N + ... + a_N."""
    N = len(M)
    A = [[Fraction(v) for v in row] for row in M]
    Mk = [[Fraction(int(i == j)) for j in range(N)] for i in range(N)]
    co = [Fraction(1)]
    for k in range(1, N + 1):
        AM = [[sum(A[i][t] * Mk[t][j] for t in range(N)) for j in range(N)] for i in range(N)]
        c = -Fraction(sum(AM[i][i] for i in range(N)), k)
        co.append(c)
        Mk = [[AM[i][j] + (c if i == j else 0) for j in range(N)] for i in range(N)]
    assert all(x.denominator == 1 for x in co)
    return [int(x) for x in co]


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def poly_prod(*fs):
    out = [1]
    for f in fs:
        out = poly_mul(out, f)
    return out


def poly_pow(f, e):
    out = [1]
    for _ in range(e):
        out = poly_mul(out, f)
    return out


def poly_at(coeffs, x):
    v = Fraction(0)
    for a in coeffs:
        v = v * x + a
    return v


def poly_str(c):
    N = len(c) - 1
    parts = []
    for i, a in enumerate(c):
        if a == 0:
            continue
        p = N - i
        if p == 0:
            parts.append(str(a))
        else:
            head = "" if a == 1 else ("-" if a == -1 else str(a))
            parts.append(head + ("x^%d" % p if p > 1 else "x"))
    s = " + ".join(parts).replace("+ -", "- ")
    return s or "0"


def inertia(M):
    """(#negative, #zero, #positive) eigenvalues of a symmetric matrix of Fractions.
    Exact symmetric LDL^T with pivoting; a 2x2 block with zero diagonal and nonzero
    off-diagonal b has eigenvalues +-b, contributing one of each sign."""
    S = [[Fraction(v) for v in row] for row in M]
    neg = pos = zero = 0
    while S:
        k = len(S)
        piv = next((i for i in range(k) if S[i][i] != 0), None)
        if piv is None:
            off = next(((i, j) for i in range(k) for j in range(i + 1, k) if S[i][j] != 0), None)
            if off is None:
                zero += k
                break
            i, j = off
            S[0], S[i] = S[i], S[0]
            for row in S:
                row[0], row[i] = row[i], row[0]
            S[1], S[j] = S[j], S[1]
            for row in S:
                row[1], row[j] = row[j], row[1]
            b = S[0][1]
            pos += 1
            neg += 1
            S = [[S[a][c] - (S[a][0] * S[1][c] + S[a][1] * S[0][c]) / b
                  for c in range(2, k)] for a in range(2, k)]
        else:
            S[0], S[piv] = S[piv], S[0]
            for row in S:
                row[0], row[piv] = row[piv], row[0]
            d = S[0][0]
            if d > 0:
                pos += 1
            else:
                neg += 1
            S = [[S[a][c] - S[a][0] * S[0][c] / d for c in range(1, k)] for a in range(1, k)]
    return neg, zero, pos


def above(M, t):
    """(number of eigenvalues of M strictly above the rational t, number equal to t)."""
    t = Fraction(t)
    N = [[(t if i == j else Fraction(0)) - M[i][j] for j in range(len(M))] for i in range(len(M))]
    neg, zero, _ = inertia(N)
    return neg, zero


# ----------------------------------------------------------------------------------------
# 4.  Z[omega], omega a primitive 6th root of unity:  omega^2 = omega - 1.
# ----------------------------------------------------------------------------------------
OMEGA = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]   # omega^k as a + b*omega


def zadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def zscale(k, x):
    return (k * x[0], k * x[1])


# ========================================================================================
print("=" * 88)
print("A  THE CELL (m,n) = (6,2):  P_6(2) ON 72 VERTICES")
print("=" * 88)

V62, IDX62, ADJ62 = basic_graph_checks(6, 2, "(6,2)")

# --- the partition and the quotient, rebuilt from the graph ------------------------------
Q62, EQ62, SZ62 = quotient(V62, ADJ62, 6, 2)
check(len(SZ62) == 12 and set(SZ62) == {6}, "(6,2) partition has mn = 12 cells of size 6",
      "6 * 12 = 72")
check(EQ62, "(6,2) partition is equitable",
      "neighbour-count vector constant on each cell, tested at every one of 72 vertices")
check(Q62 == circ_form(6, 2), "(6,2) quotient = circ(2 D_2, E_2, O, O, O, E_2)",
      "reproduces Greaves-Zhu lem:circ from the raw graph")
check(Q62 == kron_form(6, 2), "(6,2) quotient = I_6 (x) 2 D_2 + C_6 (x) E_2",
      "the identity the Fourier block-diagonalisation uses")
check(set(sum(row) for row in Q62) == {4}, "(6,2) Q has constant row sum 4",
      "so lambda_1(Q) = 4")

CP62 = charpoly(Q62)
WANT62 = poly_prod([1, 0], [1, -4], [1, 0, -8], poly_pow([1, -3, 1], 2), poly_pow([1, -1, -3], 2))
check(CP62 == WANT62, "(6,2) charpoly(Q) = x (x-4) (x^2-8) (x^2-3x+1)^2 (x^2-x-3)^2",
      "degree %d, computed exactly over Z" % (len(CP62) - 1))

# --- the eigenvector v ------------------------------------------------------------------
VAL = {0: 2, 1: 1, 2: -1, 3: -2, 4: -1, 5: 1}          # v(d) = 2 cos(pi d / 3), from the paper


def d_of(sigma):
    p, cp = where(sigma, 1)
    q, cq = where(sigma, 2)
    return (cp - cq) % 6


classes = [0] * 6
for s in V62:
    classes[d_of(s)] += 1
check(classes == [12] * 6, "(6,2) each of the 6 classes d = 0..5 holds 12 vertices",
      "so v is orthogonal to the all-ones vector")
vec = [VAL[d_of(s)] for s in V62]
check(sum(vec) == 0, "(6,2) sum of v over V(P_6(2)) is 0", "v = (2,1,-1,-2,-1,1) on d = 0..5")
check(any(x != 0 for x in vec), "(6,2) v is not the zero vector", "v(0) = 2")
bad = [i for i in range(len(V62)) if sum(vec[j] for j in ADJ62[i]) != 3 * vec[i]]
check(not bad, "(6,2) A v = 3 v at every one of the 72 vertices",
      "exact integer arithmetic, %d violations" % len(bad))

sigma0 = ((1, 0), (2, 0))
i0 = IDX62[sigma0]
nb = sorted(vec[j] for j in ADJ62[i0])
check(vec[i0] == 2 and sum(nb) == 6, "(6,2) the paper's hand row at w = 12, c = (0,0)",
      "neighbour values %s sum to 6 = 3 * 2" % ("+".join(str(x) for x in sorted(nb))))

# --- lambda_2(Q) = 2 sqrt 2, and the decision -------------------------------------------
n_above3, n_at3 = above(Q62, 3)
check((n_above3, n_at3) == (1, 0), "(6,2) exactly one eigenvalue of Q exceeds 3, none equals 3",
      "inertia of 3I - Q; hence lambda_2(Q) < 3 strictly")
LO2S2 = Fraction(2828427124746189, 10 ** 15)
HI2S2 = Fraction(2828427124746191, 10 ** 15)
check(LO2S2 ** 2 < 8 < HI2S2 ** 2, "(6,2) 2 sqrt 2 = 2.828427124746190 to 15 decimals",
      "exact rational squaring brackets sqrt 8")
check(above(Q62, LO2S2)[0] == 2 and above(Q62, HI2S2)[0] == 1,
      "(6,2) lambda_2(Q) lies in (2.828427124746189, 2.828427124746191]",
      "so lambda_2(Q) = 2 sqrt 2 and beta_{6,2} = 4 - 2 sqrt 2 = 1.171572875253810")
check(poly_at([1, 0, -8], Fraction(3)) == 1, "(6,2) x^2 - 8 at x = 3 equals 1 > 0",
      "the whole refutation is the integer inequality 9 > 8")
check(Fraction(4) - HI2S2 > 1, "(6,2) DECISION psi_{6,2} <= 1 < 4 - 2 sqrt 2 = beta_{6,2}",
      "strict deficit at least 3 - 2 sqrt 2 > 0.17157287525380")

print()
print("=" * 88)
print("B  THE CELL (m,n) = (6,3):  P_6(3) ON 1296 VERTICES")
print("=" * 88)

V63, IDX63, ADJ63 = basic_graph_checks(6, 3, "(6,3)")

Q63, EQ63, SZ63 = quotient(V63, ADJ63, 6, 3)
check(len(SZ63) == 18 and set(SZ63) == {72}, "(6,3) partition has mn = 18 cells of size 72",
      "18 * 72 = 1296")
check(EQ63, "(6,3) partition is equitable",
      "neighbour-count vector constant on each cell, tested at every one of 1296 vertices")
check(Q63 == circ_form(6, 3), "(6,3) quotient = circ(2 D_3, E_3, O, O, O, E_3)", "lem:circ again")
check(Q63 == kron_form(6, 3), "(6,3) quotient = I_6 (x) 2 D_3 + C_6 (x) E_3", "")
check(set(sum(row) for row in Q63) == {6}, "(6,3) Q has constant row sum 6", "so lambda_1(Q) = 6")

CP63 = charpoly(Q63)
WANT63 = poly_prod([1, 0], [1, -4], [1, -6],
                   poly_pow([1, -8, 17, -5], 2), poly_pow([1, -4, -3, 9], 2), [1, -2, -16, 16])
check(CP63 == WANT63,
      "(6,3) charpoly(Q) = x(x-4)(x-6) (x^3-8x^2+17x-5)^2 (x^3-4x^2-3x+9)^2 (x^3-2x^2-16x+16)",
      "degree %d, computed exactly over Z" % (len(CP63) - 1))

# --- the character-twisted invariant subspace, and T ------------------------------------
PAIRS = [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]
T_PAPER = [[1, 0, 2, 0, 0, 2],
           [0, 1, 0, 1, 2, 0],
           [2, 0, 1, 2, 0, 0],
           [0, 1, 2, 2, 0, 0],
           [0, 2, 0, 0, 1, 1],
           [2, 0, 0, 0, 1, 2]]

feat = []
for s in V63:
    p, cp = where(s, 1)
    q, cq = where(s, 2)
    feat.append(((p, q), (cp - cq) % 6))
F = [[OMEGA[d] if pk == pr else (0, 0) for (pk, d) in feat] for pr in PAIRS]

ortho = True
for row in range(6):
    tot = (0, 0)
    for x in F[row]:
        tot = zadd(tot, x)
    if tot != (0, 0):
        ortho = False
check(ortho, "(6,3) every f_{p,q} is orthogonal to the all-ones vector",
      "sums to 0 in Z[omega], omega^2 = omega - 1")

T_BUILT = [[0] * 6 for _ in range(6)]
invariant = True
for col in range(6):
    Af = []
    for i in range(len(V63)):
        tot = (0, 0)
        for j in ADJ63[i]:
            tot = zadd(tot, F[col][j])
        Af.append(tot)
    for row, pr2 in enumerate(PAIRS):
        i = next(i for i in range(len(V63)) if feat[i] == (pr2, 0))   # there omega^0 = 1
        val = Af[i]
        if val[1] != 0:
            invariant = False
        T_BUILT[row][col] = val[0]
    rebuilt = [(0, 0)] * len(V63)
    for row in range(6):
        c = T_BUILT[row][col]
        if c:
            for i in range(len(V63)):
                rebuilt[i] = zadd(rebuilt[i], zscale(c, F[row][i]))
    if rebuilt != Af:
        invariant = False
check(invariant, "(6,3) span{f_{p,q}} is A-invariant, verified at all 1296 vertices",
      "A f_{p,q} = sum_h T[h][(p,q)] f_h, exactly in Z[omega]; the step the review left asserted")
check(T_BUILT == T_PAPER, "(6,3) the matrix of A on that subspace is the T printed in the paper",
      "derived from the 1296-vertex graph, not transcribed")
check(all(T_PAPER[i][j] == T_PAPER[j][i] for i in range(6) for j in range(6)),
      "(6,3) T is symmetric", "so its eigenvalues are real")
check([sum(row) for row in T_PAPER] == [5, 4, 5, 5, 4, 5],
      "(6,3) T has row sums (5,4,5,5,4,5), NOT the degree 6",
      "T is NOT an equitable-partition quotient; it bounds lambda_2, not lambda_1")

PERM = {0: 2, 1: 4, 2: 0, 3: 5, 4: 1, 5: 3}                # (p,q) -> (q,p) = (1 3)(2 5)(4 6)
check(all(T_PAPER[i][j] == T_PAPER[PERM[i]][PERM[j]] for i in range(6) for j in range(6)),
      "(6,3) T commutes with the letter-swap involution (1 3)(2 5)(4 6)", "")

SYM = [[1, 0, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0], [0, 0, 0, 1, 0, 1]]
ASYM = [[1, 0, -1, 0, 0, 0], [0, 1, 0, 0, -1, 0], [0, 0, 0, 1, 0, -1]]


def restrict(basis):
    cols = []
    for k in range(3):
        img = [sum(T_PAPER[i][j] * basis[k][j] for j in range(6)) for i in range(6)]
        coef = []
        for l in range(3):
            piv = next(i for i in range(6) if basis[l][i] != 0)
            coef.append(Fraction(img[piv], basis[l][piv]))
        rebuilt = [sum(coef[l] * basis[l][i] for l in range(3)) for i in range(6)]
        assert rebuilt == [Fraction(x) for x in img]
        cols.append([int(x) for x in coef])
    return [[cols[k][l] for k in range(3)] for l in range(3)]


S = restrict(SYM)
SP = restrict(ASYM)
check(S == [[3, 0, 2], [0, 3, 1], [2, 1, 2]], "(6,3) symmetric block S = [[3,0,2],[0,3,1],[2,1,2]]",
      "restriction of T to the +1 eigenspace of the involution")
check(SP == [[-1, 0, -2], [0, -1, 1], [-2, 1, 2]],
      "(6,3) antisymmetric block S' = [[-1,0,-2],[0,-1,1],[-2,1,2]]", "")
CS, CSP = charpoly(S), charpoly(SP)
check(CS == poly_prod([1, -3], [1, -5, 1]), "(6,3) charpoly(S) = %s = (x-3)(x^2-5x+1)" % poly_str(CS), "")
check(CSP == poly_prod([1, 1], [1, -1, -7]), "(6,3) charpoly(S') = %s = (x+1)(x^2-x-7)" % poly_str(CSP), "")
check(charpoly(T_PAPER) == poly_mul(CS, CSP), "(6,3) charpoly(T) = charpoly(S) * charpoly(S')",
      "the involution splits T with no remainder")

SEP63 = Fraction(47, 10)
check(above(T_PAPER, SEP63)[0] >= 1, "(6,3) T has an eigenvalue above 47/10",
      "so lambda_2(A(P_6(3))) >= lambda_max(T) > 47/10")
check(poly_at([1, -5, 1], SEP63) == Fraction(-41, 100),
      "(6,3) x^2 - 5x + 1 at 47/10 is -41/100 < 0",
      "so (5 + sqrt 21)/2 > 47/10, i.e. 25 * 21 = 525 > 484 = 22^2")
LO21 = Fraction(4582575694955840, 10 ** 15)
HI21 = Fraction(4582575694955841, 10 ** 15)
check(LO21 ** 2 < 21 < HI21 ** 2, "(6,3) sqrt 21 = 4.582575694955840 to 15 decimals",
      "so lambda_max(T) = (5 + sqrt 21)/2 = 4.791287847477920")
check(above(T_PAPER, (5 + LO21) / 2)[0] == 1 and above(T_PAPER, (5 + HI21) / 2)[0] == 0,
      "(6,3) lambda_max(T) lies in ((5+4.582575694955840)/2, (5+4.582575694955841)/2]",
      "so psi_{6,3} <= 6 - lambda_max(T) = (7 - sqrt 21)/2 = 1.208712152522080")

n_ab, n_at = above(Q63, SEP63)
check((n_ab, n_at) == (1, 0), "(6,3) exactly one eigenvalue of Q exceeds 47/10, none equals it",
      "inertia of (47/10) I - Q; hence lambda_2(Q) < 47/10")
G3 = [1, -2, -16, 16]
LOL2 = Fraction(4685846165554339, 10 ** 15)
HIL2 = Fraction(4685846165554341, 10 ** 15)
check(poly_at(G3, LOL2) < 0 < poly_at(G3, HIL2),
      "(6,3) x^3 - 2x^2 - 16x + 16 changes sign in (4.685846165554339, 4.685846165554341)", "")
check(above(Q63, LOL2)[0] == 2 and above(Q63, HIL2)[0] == 1,
      "(6,3) lambda_2(Q) = 4.685846165554340 to 15 decimals",
      "beta_{6,3} = 6 - lambda_2(Q) = 1.314153834445660")
check(6 - (5 + HI21) / 2 < Fraction(13, 10) <= 6 - HIL2,
      "(6,3) DECISION psi_{6,3} <= 1.20871215252208 < 13/10 <= beta_{6,3} = 1.31415383444566",
      "the rational 47/10 separates the two spectra")

print()
print("=" * 88)
print("C  THE FAMILY (m,2), m >= 6:  THE ALGEBRAIC IDENTITIES OF THE PROOF")
print("=" * 88)
# Polynomials in one variable over Z, coefficient lists high degree first, as above.

def padd(*polys):
    """Sum of polynomials given as coefficient lists, high degree first."""
    d = max(len(p) for p in polys)
    out = [0] * d
    for p in polys:
        q = [0] * (d - len(p)) + list(p)
        out = [a + b for a, b in zip(out, q)]
    return out


def pneg(p):
    return [-a for a in p]


def is_zero(p):
    return all(a == 0 for a in p)


# C1.  Build 2 D_2 + 2c E_2 entrywise as a matrix over Z[c] and take its trace and
#      determinant symbolically:  the characteristic polynomial is x^2 - (2c+2) x + (4c - 4c^2).
TWO_C = [2, 0]                                                     # the polynomial 2c
MSYM = [[padd([0], TWO_C), padd([0], TWO_C)], [padd([0], TWO_C), [2]]]
TR = padd(MSYM[0][0], MSYM[1][1])                                  # trace
DT = padd(poly_mul(MSYM[0][0], MSYM[1][1]), pneg(poly_mul(MSYM[0][1], MSYM[1][0])))
check(TR == [2, 2] and DT == [-4, 4, 0],
      "family charpoly(2 D_2 + 2c E_2) = x^2 - (2c+2) x + (4c - 4c^2)",
      "trace and determinant of [[2c, 2c], [2c, 2]] taken symbolically over Z[c]")

# C2.  ((c+1) + s)^2 - (2c+2)((c+1) + s) + (4c - 4c^2) = 0  when s^2 = 5c^2 - 2c + 1,
#      i.e. F(c) = 4 - (c+1) - sqrt(5c^2-2c+1) really is 4 minus the top eigenvalue.
RAD = [5, -2, 1]                              # 5c^2 - 2c + 1
CP1 = [1, 1]                                  # c + 1
S_FREE = padd(poly_mul(CP1, CP1), RAD, pneg(poly_mul(TR, CP1)), DT)
S_LIN = padd(poly_mul([2], CP1), pneg(TR))    # 2(c+1) - (2c+2)
check(is_zero(S_FREE) and is_zero(S_LIN),
      "family (c+1) + sqrt(5c^2-2c+1) is a root of that polynomial",
      "s-free and s-linear parts both vanish identically in Z[c] after s^2 -> 5c^2-2c+1")

# C3.  (5c-1)^2 - (5c^2 - 2c + 1) = 20c^2 - 8c:  the shape lemma.  F'(c) = 0 forces
#      sqrt(5c^2-2c+1) = -(5c-1), hence 20c^2 - 8c = 0 with 5c - 1 < 0, so c = 0 only.
check(padd(poly_mul([5, -1], [5, -1]), pneg(RAD)) == [20, -8, 0],
      "family (5c-1)^2 - (5c^2-2c+1) = 20c^2 - 8c",
      "c = 0 is the only critical point of F; F increases on [-1,0], decreases on [0,1]")

# C4.  c = 1 - eps  =>  5c^2 - 2c + 1 = 4 - 8 eps + 5 eps^2.
ONE_MINUS = [-1, 1]                           # 1 - eps
RAD_EPS = padd(poly_mul([5], poly_mul(ONE_MINUS, ONE_MINUS)), poly_mul([-2], ONE_MINUS), [1])
check(RAD_EPS == [5, -8, 4], "family at c = 1 - eps the radicand is 4 - 8 eps + 5 eps^2",
      "not 4 - 6 eps + 5 eps^2, which is the slip an earlier draft made")

# C5.  (2 - eps)^2 - (4 - 8 eps + 5 eps^2) = 4 eps - 4 eps^2, so F(c_1) > 2 eps <=> eps < 1.
check(padd(poly_mul([-1, 2], [-1, 2]), pneg(RAD_EPS)) == [-4, 4, 0],
      "family (2-eps)^2 - (4 - 8 eps + 5 eps^2) = 4 eps - 4 eps^2",
      "branch (i): F(cos 2pi/m) > 2 eps  <=>  4 eps > 4 eps^2  <=>  eps < 1")

# C6.  F(-1) = 4 - sqrt 8, and 4 - 2 sqrt 2 > 2 eps  <=>  cos(2pi/m) > sqrt 2 - 1.
check(poly_at(RAD, Fraction(-1)) == 8, "family radicand at c = -1 is 8, so F(-1) = 4 - 2 sqrt 2",
      "branch (ii): F is increasing on [-1,0], so F(c) >= 4 - 2 sqrt 2 for c <= 0")
check(Fraction(3, 2) ** 2 > 2, "family 1/2 > sqrt 2 - 1, because (3/2)^2 = 9/4 > 2",
      "cos(2pi/6) = 1/2, so branch (ii) first holds at m = 6: 9 > 8 again")
check(all(Fraction(2, m) < Fraction(1, 2) for m in range(5, 200)),
      "family 2/m < 1/2 for every m in 5..199",
      "so 2pi/m < pi/2 and eps = 1 - cos(2pi/m) < 1: branch (i) holds from m = 5")
check(all(Fraction(2, m) <= Fraction(1, 3) for m in range(6, 200)),
      "family 2/m <= 1/3 for every m in 6..199",
      "so cos(2pi/m) >= cos(pi/3) = 1/2 > sqrt 2 - 1: branch (ii) holds from m = 6")

print()
print("=" * 88)
print("D  THE FAMILY (m,2):  AN INDEPENDENT EXACT CERTIFICATE AT EACH m = 6..20")
print("=" * 88)
# For each m the paper prints a rational separator t_m.  Two exact inertia counts decide the
# cell: Q_m has exactly one eigenvalue above t_m (namely lambda_1 = 4) and none equal to it,
# so beta_{m,2} = 4 - lambda_2(Q_m) > 4 - t_m; and C_m has at least two eigenvalues above
# t_m - 2, so 2 + 2cos(2pi/m) > t_m and psi_{m,2} <= 4 - (2 + 2cos(2pi/m)) < 4 - t_m.
SEPARATOR = {6: Fraction(29, 10), 7: Fraction(3), 8: Fraction(16, 5), 9: Fraction(7, 2),
             10: Fraction(7, 2), 11: Fraction(18, 5), 12: Fraction(37, 10), 13: Fraction(37, 10),
             14: Fraction(19, 5), 15: Fraction(19, 5), 16: Fraction(19, 5), 17: Fraction(19, 5),
             18: Fraction(387, 100), 19: Fraction(387, 100), 20: Fraction(39, 10)}

for m in range(6, 21):
    V, IDX, ADJ = pancake(m, 2)
    # (a) the neighbour law: the four neighbours of a vertex with colour difference d carry
    #     colour differences d+1, d-1, d, d.  Hence for ANY f on Z_m, (A f)(d) = f(d+1) +
    #     f(d-1) + 2 f(d); with f(d) = zeta^d + zeta^-d this is (2 + 2cos(2pi/m)) f.
    lawful = True
    for i, s in enumerate(V):
        p, cp = where(s, 1)
        q, cq = where(s, 2)
        d = (cp - cq) % m
        got = []
        for j in ADJ[i]:
            pp, cpp = where(V[j], 1)
            qq, cqq = where(V[j], 2)
            got.append((cpp - cqq) % m)
        if sorted(got) != sorted([(d + 1) % m, (d - 1) % m, d, d]):
            lawful = False
            break
    check(lawful, "family m = %2d neighbour law {d+1, d-1, d, d} at all %d vertices"
          % (m, len(V)), "so psi_{m,2} <= 2 - 2 cos(2pi/m)")
    Qm, eqm, szm = quotient(V, ADJ, m, 2)
    check(eqm and Qm == circ_form(m, 2) == kron_form(m, 2),
          "family m = %2d partition equitable and quotient = circ(2 D_2, E_2, O.., E_2)" % m,
          "cells of size %d" % szm[0])
    t = SEPARATOR[m]
    aQ = above(Qm, t)
    aC = above(cycle(m), t - 2)
    good = aQ == (1, 0) and aC[0] >= 2
    check(good, "family m = %2d REFUTED by the separator t = %s" % (m, t),
          "eigenvalues of Q above t: %d (=%d), of C_m above t-2: %d (>=2)"
          % (aQ[0], aQ[0], aC[0]))

print()
print("=" * 88)
print("E  NEGATIVE CONTROLS:  THE CELLS WHERE THIS TEST DOES NOT FIRE")
print("=" * 88)
# m = 5 is an exact tie and at m = 3, 4 the quotient bound is vacuous, so the eigenvector
# argument separates nothing at m <= 5.  That is a failure of the TEST, not a decision about
# the cell: the test bounds psi above and beta below, so it cannot show psi = beta either.
C5 = cycle(5)
check(charpoly(C5) == poly_prod([1, -2], poly_pow([1, 1, -1], 2)),
      "control charpoly(C_5) = (x-2)(x^2+x-1)^2",
      "so 2 cos(2pi/5) = (sqrt 5 - 1)/2 and 2 eps = 2 - 2cos(2pi/5) = (5 - sqrt 5)/2")
Q5 = circ_form(5, 2)
A5 = Fraction(2618033, 10 ** 6)
B5 = Fraction(2618035, 10 ** 6)
check(above(Q5, A5)[0] >= 2 and above(Q5, B5)[0] == 1,
      "control m = 5: lambda_2(Q_5) lies in (2.618033, 2.618035)", "")
check(above(C5, A5 - 2)[0] >= 2 and above(C5, B5 - 2)[0] <= 1,
      "control m = 5: 2 + 2cos(2pi/5) lies in (2.618033, 2.618035) TOO",
      "an exact tie at (3 + sqrt 5)/2, so no separator exists at m = 5 and this test decides "
      "the cell in neither direction")
Q3, Q4 = circ_form(3, 2), circ_form(4, 2)
check(above(Q3, Fraction(2))[0] >= 2 and above(cycle(3), Fraction(0))[0] <= 1,
      "control m = 3: lambda_2(Q_3) > 2 > 1 = 2 + 2cos(2pi/3)",
      "beta_{3,2} = (7 - sqrt 13)/2; the eigenvector gives no refutation here")
check(above(Q4, Fraction(5, 2))[0] >= 2 and above(cycle(4), Fraction(1, 2))[0] <= 1,
      "control m = 4: lambda_2(Q_4) > 5/2 > 2 = 2 + 2cos(2pi/4)",
      "no refutation at m = 4 either")
check(above(circ_form(6, 2), Fraction(29, 10))[0] == 1
      and above(cycle(6), Fraction(9, 10))[0] >= 2,
      "control the same test applied at m = 6 does fire",
      "so the m = 3,4,5 negatives are not an artefact of the test")

print()
print("NOT RE-RUN: the exact value psi_{6,2} = 1.  This program proves only psi_{6,2} <= 1,")
print("NOT RE-RUN:   which is all the refutation needs; the matching lower bound is Blanco's")
print("NOT RE-RUN:   published theorem (arXiv:2608.15398v1) and is quoted, not verified here.")
print("NOT RE-RUN: lambda_2(A(P_6(2))) and lambda_2(A(P_6(3))) are never computed -- only")
print("NOT RE-RUN:   bounded below, by the exhibited eigenvector and by lambda_max(T).")
print("NOT RE-RUN: the cells m = 3, 4, 5 at n = 2.  Section E shows only that the separator")
print("NOT RE-RUN:   test does not fire there.  The test bounds psi_{m,2} above and beta_{m,2}")
print("NOT RE-RUN:   below, so a negative is NOT evidence that psi_{m,2} = beta_{m,2}; those")
print("NOT RE-RUN:   cells are left undecided in both directions.")
print("NOT RE-RUN: MINIMALITY.  Nothing here shows (6,2) is the smallest counterexample. The")
print("NOT RE-RUN:   cell (2,3), on 48 vertices, is not decided by this program at all, and")
print("NOT RE-RUN:   m = 2 uses a different connection set (q_2 = 1) that is not built here.")
print("NOT RE-RUN: cells with n >= 4, and m > 20 at n = 2, are not examined; the family claim")
print("NOT RE-RUN:   for all m >= 6 rests on section C's identities, not on section D's table.")
print("NOT RE-RUN: the asymptotic expansion beta_{m,2} = 3 eps - eps^2/4 + O(eps^3).")

print()
if FAILED:
    print("VERDICT: %d OF %d CHECKS FAILED" % (FAILED, CHECKS))
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % CHECKS)
sys.exit(0)
