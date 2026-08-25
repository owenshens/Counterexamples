#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- verification program for

    "Counterexamples to Wirdane's Cyclotomic Catalan Conjecture"
    (disproof of Conjecture 6.6 of arXiv:2605.14682v1)

Run:   python3 verify.py          (Python 3.9+, standard library only)
Exit:  0 if every check passes, 1 otherwise.

The paper claims that the Fuerlinger--Hofbauer q-Catalan value C_n(zeta_mu)
does NOT equal (1/mu)*binom(n+1,(n+1)/mu) whenever mu | n+1, mu >= 2, apart
from the single pair (n,mu) = (1,2); and it exhibits three infinite families
of counterexamples, (2m+1,2), (3r+2,3) and (p-1,p).

--------------------------------------------------------------------------
TAKEN FROM THE PAPER (definitions and exhibited data -- inputs, not results)
--------------------------------------------------------------------------
  P1. Equation (1) of the paper, the definition of the q-Catalan numbers:
        C_0(q) = 1,  C_n(q) = sum_{k=0}^{n-1} q^k C_k(q) C_{n-1-k}(q).
  P2. The right-hand side of the conjecture, the paper's equation (4):
        (1/mu) * binom(n+1, (n+1)/mu),   for mu | n+1.
  P3. The combinatorial definitions quoted from [Wirdane]:
        S_n(312); S_{n,k}(312) = {pi in S_n(312) : pi_{k+1} = n};
        S'_{n,k}(312) = disjoint union of S_{n,j}(312) for j = 0..k;
        inv_i(pi) = #{(a,b): a<b, pi_a>pi_b, pi_a = i mod mu}   (Def 6.1);
        C_{n,k}(q_1..q_mu) = sum_{pi in S'_{n,k}} prod_i q_i^{inv_i(pi)}
                                                            (Def 6.2);
        C^{(mu)}_{n,k} = C_{n,k}(zeta_mu, ..., zeta_mu).
  P4. Wirdane's own q-Catalan-triangle recurrence (used to re-derive the
      bridge a second, independent way):
        C_{n,0}(q) = q^binom(n,2),
        C_{n,k}(q) = C_{n,k-1}(q) + q^{n-k-1} C_{n-1,k}(q)   (1 <= k <= n-1).
  P5. The single tabulated polynomial quoted from [Wirdane]:
        C_{3,3}(q) = 1 + 2q + q^2 + q^3.
  P6. The Gaussian q-Catalan number, for the last Remark:
        Cat~_n(q) = qbinom(2n,n) / [n+1]_q.
  P7. Remark 6.7 of [Wirdane], as quoted:
        C^{(2)}_{n,k} = (-1)^{floor(n/2)} C_{n,k}.
  P8. Numerical figures held here ONLY as the targets that the derived
      values are compared against (never as seeds).  All of them are printed
      in the paper EXCEPT the first, which is flagged because provenance
      matters:
        C_n(-1) for n=0..13 equals 1,1,0,-1,0,2,0,-5,0,14,0,-42,0,132.
          *** NOT printed in the paper. ***  The paper prints only the closed
          form, equation (6); this 14-term table is a comparison target
          assembled here from that closed form.  Check B1 is therefore a
          comparison against a target fixed here, not against a figure the
          paper prints; the paper's own statement is what B2/B3 decide.
        81 admissible pairs with n <= 29; the conjecture holds at exactly
        1 of them, namely (1,2), and fails at 80; the three families
        account for 32 of the 81; non-reality alone disposes of 66;
        the odd-prime congruence reaches 9; Cat~_n(-1) at odd n <= 11
        takes the values 1, 3, 10, 35, 126, 462; C_2(zeta_3) = 1 + zeta_3;
        Cat~_2(zeta_3) = 1 + zeta_3^2 = 1/2 - i*sqrt(3)/2.

--------------------------------------------------------------------------
DERIVED HERE (everything the checks actually decide)
--------------------------------------------------------------------------
  D1.  C_n(q) for n <= 30 as exact integer coefficient lists, from P1 alone.
  D2.  Cat_n independently as binom(2n,n)/(n+1); the identity C_n(1)=Cat_n.
  D3.  Equation (5), F_q(z) = 1 + z F_q(z) F_q(qz), verified coefficientwise
       by actually multiplying the two truncated series.
  D4.  C_n(-1) for n <= 30, and equation (6) (vanishing at even n, and
       (-1)^m Cat_m at n = 2m+1); the identity (1/2)binom(2m+2,m+1) =
       (2m+1)Cat_m; hence the family (2m+1,2).
  D5.  C_n(omega) in Z[omega] (omega^2 = -1-omega) by exact reduction modulo
       Phi_3, and the closed form 1 + omega*sum_{j<=(n-2)/3} Cat_j for ALL
       n <= 30 -- stronger than the paper's n = 3r+2 statement; the
       omega-coefficient is then shown nonzero, which is exactly
       "Im > 0 while the right-hand side is real".
  D6.  Cat_{p-1} = -1 mod p, binom(2p-1,p-1) = 1 mod p, and the exact
       integer identity (2p-1)Cat_{p-1} = binom(2p-1,p-1); the image of
       C_{p-1}(zeta_p) - 1 in Z[zeta_p]/(1-zeta_p) = F_p.
  D7.  Phi_mu(q) for every mu used, by exact integer polynomial division of
       q^mu - 1 by the product of the proper divisor cyclotomics.
  D8.  The full sweep: the set of admissible pairs, its cardinality, the
       exact reduction of C_n(q) - (1/mu)binom(n+1,(n+1)/mu) modulo Phi_mu
       for every pair, and the counts of vanishing / non-vanishing pairs.
  D9.  Reality of C_n(zeta_mu), decided exactly by comparing the residue of
       C_n(q) with the residue of C_n(q^{mu-1}) modulo Phi_mu (complex
       conjugation), and the count of non-real pairs.
  D10. The three families as explicit pair sets, their pairwise
       intersections, and the cardinality of their union.
  D11. Cat~_n(q) by exact polynomial division, its integrality, and its
       residue modulo Phi_mu at every admissible pair.
  D12. S_n(312) by brute-force enumeration for n <= 8: its cardinality,
       S_{n,n}(312) = empty and equation (2), sum_i inv_i = inv, the [FH]
       interpretation C_n(q) = sum q^inv, the multivariate C_{n,k} and the
       bridge, equation (3), C^{(mu)}_{n,n} = C_n(zeta_mu) for mu = 2,3,4,5.
  D13. Wirdane's triangle iterated to n <= 11, C_{n,n}(q) = C_n(q)
       (independently supplying the source's unproved Proposition 2.4), and
       agreement of the triangle with the brute-force C_{n,k}(q), n <= 8.
  D14. The Remark 6.7 conflict, and the three colliding values at n = k = 3.
  D15. The degenerate order mu = 1, and minimality of n = 2 among mu >= 2.
  D16. The paper's CONCLUSION in one check (block J): Conjecture 6.6 is false,
       shown at (n,mu) = (2,3) with its left side built from Definitions
       6.1/6.2 by enumeration and never from equation (1).

A section headed GAPS is printed before the verdict.  It names the steps
between these facts and the paper's claim that no check covers: the bridge to
C^{(mu)}_{n,n} is machine-verified only for n <= 8 (mu <= 5) and n <= 11
(univariately); the three families are infinite but only finitely tested; the
transcription of the source's definitions is a document fact; five checks
assert more than the paper does; and four checks are structural, i.e. cannot
fail and are evidence about this program only, not about the paper.

Arithmetic is exact throughout: Python ints and fractions.Fraction only.
No floating point occurs anywhere in any certificate.  The one place the
paper speaks analytically ("its imaginary part is positive") is discharged
algebraically: 1 and omega are a Z-basis of Z[omega], so a + b*omega with
b != 0 is not rational, and Im(a + b*omega) = b*sqrt(3)/2 has the sign of b.
"""

import sys
from fractions import Fraction
from itertools import permutations

# ---------------------------------------------------------------- parameters
N_MAX = 30      # C_n(q) built for 0 <= n <= N_MAX
SWEEP_N = 29    # sweep over admissible pairs with n <= SWEEP_N  (paper: 81)
BRUTE_N = 8     # brute-force enumeration of S_n(312) for n <= BRUTE_N
TRI_N = 11      # Wirdane triangle iterated for n <= TRI_N
PRIME_MAX = 31  # odd primes p with 3 <= p <= PRIME_MAX

# ------------------------------------------------------- comparison targets
# NOTE: this 14-term table is a comparison target assembled here, not a figure
# the paper prints -- the paper prints only the closed form, equation (6), from
# which these values follow.  See P8 above.
TARGET_MINUS1 = [1, 1, 0, -1, 0, 2, 0, -5, 0, 14, 0, -42, 0, 132]  # n = 0..13
PAPER_C33 = [1, 2, 1, 1]                      # C_{3,3}(q) = 1+2q+q^2+q^3
PAPER_N_PAIRS = 81
PAPER_HOLDS_AT = (1, 2)
PAPER_N_HOLDS = 1
PAPER_N_FAILS = 80
PAPER_FAMILY_UNION = 32
PAPER_NONREAL = 66
PAPER_PRIME_PAIRS = 9
PAPER_GAUSS_MU2 = {1: 1, 3: 3, 5: 10, 7: 35, 9: 126, 11: 462}

# The statements of [Wirdane] this paper cites by number (a document fact, not
# a computation; printed for the reader, deliberately NOT counted as a check).
SOURCE_STATEMENTS = [
    ("the cyclotomic Catalan conjecture", "Conjecture 6.6"),
    ("the definition of inv_i", "Definition 6.1"),
    ("the definition of C_{n,k}(q_1,...,q_mu)", "Definition 6.2"),
    ("C_{n,k}(q) = sum over S'_{n,k}(312) of q^inv", "Theorem 3.3"),
    ("the [FH] interpretation of C_n(q)", "Remark 3.4"),
    ("C_{n,n}(q) = C_n(q), its proof deferred there", "Proposition 2.4"),
    ("specialization of all q_i to a common q", "Corollary 6.5"),
    ("the mu = 2 sign rule contradicted here", "Remark 6.7"),
]

RESULTS = []    # list of (bool, text)


def check(ok, text):
    """Record one check and print its PASS/FAIL line immediately."""
    RESULTS.append((bool(ok), text))
    print(("PASS " if ok else "FAIL ") + text)
    return bool(ok)


def info(text):
    print("      " + text)


# =========================================================================
# 1.  Exact polynomial arithmetic.  A polynomial is a list of coefficients,
#     index = degree, over Z (ints) or Q (Fraction).  [] is the zero poly.
# =========================================================================

def p_trim(a):
    i = len(a)
    while i > 0 and a[i - 1] == 0:
        i -= 1
    return a[:i]


def p_add(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i, c in enumerate(a):
        out[i] += c
    for i, c in enumerate(b):
        out[i] += c
    return p_trim(out)


def p_sub(a, b):
    n = max(len(a), len(b))
    out = [0] * n
    for i, c in enumerate(a):
        out[i] += c
    for i, c in enumerate(b):
        out[i] -= c
    return p_trim(out)


def p_mul(a, b):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        if ca:
            for j, cb in enumerate(b):
                if cb:
                    out[i + j] += ca * cb
    return p_trim(out)


def p_shift(a, k):
    """multiply by q^k"""
    return ([0] * k + list(a)) if a else []


def p_eval(a, x):
    """Horner evaluation at a ring element x supporting * and +; exact."""
    acc = 0
    for c in reversed(a):
        acc = acc * x + c
    return acc


def p_str(a, var="q"):
    if not a:
        return "0"
    parts = []
    for i, c in enumerate(a):
        if c == 0:
            continue
        if i == 0:
            parts.append(str(c))
        elif i == 1:
            parts.append(("%s*%s" % (c, var)) if c not in (1,) else var)
        else:
            parts.append(("%s*%s^%d" % (c, var, i)) if c not in (1,)
                         else "%s^%d" % (var, i))
    return " + ".join(parts)


def p_divmod(a, b):
    """Exact division with remainder over Q.  b must be nonzero.
    Returns (quotient, remainder) with Fraction coefficients."""
    a = [Fraction(c) for c in p_trim(list(a))]
    b = [Fraction(c) for c in p_trim(list(b))]
    if not b:
        raise ZeroDivisionError("polynomial division by zero")
    q = [Fraction(0)] * max(0, len(a) - len(b) + 1)
    lead = b[-1]
    while len(a) >= len(b) and a:
        shift = len(a) - len(b)
        coef = a[-1] / lead
        q[shift] = coef
        for i, cb in enumerate(b):
            a[shift + i] -= coef * cb
        a = p_trim(a)
    return p_trim(q), a


def p_exact_div(a, b):
    """Divide, insisting the remainder is exactly zero."""
    q, r = p_divmod(a, b)
    if r:
        raise ValueError("inexact polynomial division, remainder " + p_str(r))
    return q


def p_mod(a, m):
    """Reduce a modulo m (over Q).  Canonical representative, degree < deg m."""
    return p_divmod(a, m)[1]


def to_int_poly(a):
    """Convert a Fraction-coefficient poly to ints, insisting it is integral."""
    out = []
    for c in a:
        f = Fraction(c)
        if f.denominator != 1:
            raise ValueError("non-integral coefficient " + str(f))
        out.append(int(f))
    return p_trim(out)


_PHI_CACHE = {}


def cyclotomic(mu):
    """Phi_mu(q), derived (not tabulated) as (q^mu - 1) divided by the product
    of Phi_d(q) over the proper divisors d of mu."""
    if mu in _PHI_CACHE:
        return _PHI_CACHE[mu]
    num = [-1] + [0] * (mu - 1) + [1]          # q^mu - 1
    den = [1]
    for d in range(1, mu):
        if mu % d == 0:
            den = p_mul(den, cyclotomic(d))
    phi = to_int_poly(p_exact_div(num, den))
    _PHI_CACHE[mu] = phi
    return phi


# =========================================================================
# 2.  The objects of the paper.
# =========================================================================

def binom(n, k):
    """Exact integer binomial coefficient, built multiplicatively."""
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    num, den = 1, 1
    for j in range(1, k + 1):
        num *= n - k + j
        den *= j
    return num // den


def catalan(m):
    """Cat_m = binom(2m,m)/(m+1), derived independently of the recurrence."""
    b = binom(2 * m, m)
    assert b % (m + 1) == 0
    return b // (m + 1)


def q_catalan(nmax):
    """C_n(q) for 0 <= n <= nmax, from the paper's equation (1),
         C_0 = 1,  C_n = sum_{k=0}^{n-1} q^k C_k C_{n-1-k}.
    Returns a list of integer coefficient lists."""
    C = [[1]]
    for n in range(1, nmax + 1):
        acc = []
        for k in range(n):
            acc = p_add(acc, p_shift(p_mul(C[k], C[n - 1 - k]), k))
        C.append(acc)
    return C


def conj_rhs(n, mu):
    """The right-hand side of equation (4), (1/mu)binom(n+1,(n+1)/mu),
    as an exact rational.  Requires mu | n+1."""
    assert (n + 1) % mu == 0
    return Fraction(binom(n + 1, (n + 1) // mu), mu)


def gaussian_pascal(nmax, kmax):
    """Table T[m][j] = qbinom(m,j) for 0 <= j <= min(m,kmax), 0 <= m <= nmax,
    from the q-Pascal recurrence qbinom(m,j) = qbinom(m-1,j-1)
    + q^j qbinom(m-1,j).  One table, shared by every n."""
    T = [[[1]] + [[] for _ in range(kmax)]]
    for m in range(1, nmax + 1):
        prev = T[-1]
        row = [[1]] + [[] for _ in range(kmax)]
        for j in range(1, min(kmax, m) + 1):
            row[j] = p_add(prev[j - 1], p_shift(prev[j], j))
        T.append(row)
    return T


def gaussian_pascal_alt(nmax, kmax):
    """The SAME table by the OTHER q-Pascal recurrence,
         qbinom(m,j) = q^{m-j} qbinom(m-1,j-1) + qbinom(m-1,j).
    Used so that check F1 has something independent to fail against: the
    quotient Cat~_n(q) is produced by dividing the first table, so
    re-multiplying it can only reproduce that same table by construction."""
    T = [[[1]] + [[] for _ in range(kmax)]]
    for m in range(1, nmax + 1):
        prev = T[-1]
        row = [[1]] + [[] for _ in range(kmax)]
        for j in range(1, min(kmax, m) + 1):
            row[j] = p_add(p_shift(prev[j - 1], m - j), prev[j])
        T.append(row)
    return T


def gaussian_q_catalans(nmax):
    """Cat~_n(q) = qbinom(2n,n) / [n+1]_q for 0 <= n <= nmax, by exact
    polynomial division (which also proves integrality).
    Returns (list of Cat~_n, the q-Pascal table it was divided out of)."""
    T = gaussian_pascal(2 * nmax, nmax)
    out = []
    for n in range(nmax + 1):
        bracket = [1] * (n + 1)                # [n+1]_q = 1 + q + ... + q^n
        out.append(to_int_poly(p_exact_div(T[2 * n][n], bracket)))
    return out, T


# =========================================================================
# 3.  Working inside Z[zeta_mu] = Z[q]/(Phi_mu(q)).  Every element has a
#     unique representative of degree < deg Phi_mu, so equality of residues
#     is equality of algebraic numbers, and a residue of degree > 0 is a
#     provably IRRATIONAL value -- which is all that is needed against the
#     rational right-hand side of the conjecture.  Irrational does not mean
#     non-real (zeta_5 + zeta_5^4 = 2cos(2pi/5) is real, residue q + q^4), so
#     reality is decided separately and exactly, by the conjugation test
#     is_real_at below.
# =========================================================================

def residue(poly, mu):
    """The canonical representative of poly mod Phi_mu, integer coefficients."""
    return to_int_poly(p_mod(poly, cyclotomic(mu)))


def subst_power(a, e):
    """Substitute q -> q^e."""
    if not a:
        return []
    out = [0] * ((len(a) - 1) * e + 1)
    for i, c in enumerate(a):
        out[i * e] += c
    return p_trim(out)


def conjugate_residue(poly, mu):
    """Complex conjugation on Z[zeta_mu] is zeta_mu -> zeta_mu^{-1} =
    zeta_mu^{mu-1}, i.e. q -> q^{mu-1} followed by reduction."""
    return residue(subst_power(residue(poly, mu), mu - 1), mu)


def is_real_at(poly, mu):
    """Exact decision: is the value of poly at zeta_mu fixed by conjugation?"""
    return residue(poly, mu) == conjugate_residue(poly, mu)


def omega_coords(poly):
    """Write the value at omega = zeta_3 as a + b*omega (omega^2 = -1-omega)."""
    r = residue(poly, 3)
    r = list(r) + [0, 0]
    return r[0], r[1]


def diff_mod_phi(poly, rhs, mu):
    """Reduce poly(q) - rhs modulo Phi_mu(q) over Q; returns the residue.
    Empty list means the conjecture's two sides agree at (n,mu)."""
    lhs = [Fraction(c) for c in poly]
    return p_mod(p_sub(lhs, [Fraction(rhs)]), cyclotomic(mu))


# =========================================================================
# 4.  The combinatorial side: S_n(312), inv, inv_i.
# =========================================================================

def avoids_312(pi):
    """True iff pi contains no occurrence of the pattern 312, i.e. no
    a < b < c with pi_b < pi_c < pi_a."""
    n = len(pi)
    pref = 0                    # max of pi_1..pi_{b-1}
    for b in range(n):
        if pref > 0:
            for c in range(b + 1, n):
                if pi[b] < pi[c] < pref:
                    return False
        if pi[b] > pref:
            pref = pi[b]
    return True


def gen_312(n):
    """S_n(312) by brute force over all n! permutations of {1,...,n}."""
    if n == 0:
        return [()]
    return [pi for pi in permutations(range(1, n + 1)) if avoids_312(pi)]


def inv_count(pi):
    """The inversion number inv(pi)."""
    n = len(pi)
    return sum(1 for a in range(n) for b in range(a + 1, n) if pi[a] > pi[b])


def inv_residues(pi, mu):
    """(inv_1(pi), ..., inv_mu(pi)) of Definition 6.1: inversion (a,b) is
    credited to the class i in {1,...,mu} with pi_a = i (mod mu)."""
    out = [0] * mu
    n = len(pi)
    for a in range(n):
        va = pi[a]
        i = va % mu
        if i == 0:
            i = mu                      # residue system is 1,...,mu
        for b in range(a + 1, n):
            if va > pi[b]:
                out[i - 1] += 1
    return tuple(out)


def max_position(pi):
    """k with pi_{k+1} = n, i.e. the index of the largest entry (0-based)."""
    return pi.index(max(pi))


# =========================================================================
# 5.  Wirdane's q-Catalan triangle, iterated from his own recurrence (P4).
#     The recurrence step is applied for 1 <= k <= n-1; the last entry uses
#     S_{n,n}(312) = empty (equation (2)), i.e. C_{n,n} = C_{n,n-1}.  Note that
#     a literal k = n step would need q^{n-k-1} = q^{-1}, which is not a
#     polynomial -- so the empty-set reading is the only consistent one, and
#     it is confirmed independently against brute force below.
# =========================================================================

def wirdane_triangle(nmax):
    """tri[n][k] = C_{n,k}(q) for 0 <= k <= n <= nmax, integer coefficients."""
    tri = [[[1]]]                                  # C_{0,0}(q) = q^0 = 1
    for n in range(1, nmax + 1):
        row = [p_shift([1], binom(n, 2))]          # C_{n,0}(q) = q^binom(n,2)
        for k in range(1, n):
            row.append(p_add(row[k - 1],
                             p_shift(tri[n - 1][k], n - k - 1)))
        row.append(list(row[n - 1]))               # C_{n,n} = C_{n,n-1}
        tri.append(row)
    return tri


# =========================================================================
# 6.  Definition 6.2 evaluated by enumeration, and the bridge, equation (3).
# =========================================================================

def brute_triangle_row(n, s312):
    """C_{n,k}(q) = sum over S'_{n,k}(312) of q^inv, for k = 0..n, built by
    bucketing S_n(312) on the position of its largest entry and taking
    prefix sums.  Also returns the buckets' sizes."""
    buckets = [[] for _ in range(n + 1)]
    sizes = [0] * (n + 1)
    for pi in s312:
        k = max_position(pi) if n else 0
        sizes[k] += 1
        buckets[k] = p_add(buckets[k], p_shift([1], inv_count(pi)))
    row, acc = [], []
    for k in range(n + 1):
        acc = p_add(acc, buckets[k])
        row.append(acc)
    return row, sizes


def res_mul(a, b, mu):
    return residue(p_mul(a, b), mu)


def zeta_pow(e, mu):
    """zeta_mu^e as a canonical residue in Z[q]/Phi_mu."""
    return residue(p_shift([1], e), mu)


def multivar_diagonal(n, mu, s312):
    """C^{(mu)}_{n,n} = sum over S'_{n,n}(312) of prod_{i=1}^{mu}
    zeta_mu^{inv_i(pi)}, formed one monomial at a time from the exponent
    vector (inv_1,...,inv_mu) -- never collapsed to inv beforehand."""
    total = []
    for pi in s312:
        term = [1]
        for e in inv_residues(pi, mu):
            term = res_mul(term, zeta_pow(e, mu), mu)
        total = residue(p_add(total, term), mu)
    return total


# =========================================================================
# 7.  CHECK BLOCK A -- the q-Catalan polynomials themselves.
# =========================================================================

def check_block_a(C):
    print("--- A.  eq. (1), C_n(1) = Cat_n, eq. (5) -------------------------")
    info("C_1(q) = %s" % p_str(C[1]))
    info("C_2(q) = %s" % p_str(C[2]))
    info("C_3(q) = %s" % p_str(C[3]))
    info("C_4(q) = %s" % p_str(C[4]))
    info("deg C_n = binom(n,2) for n <= %d : %s"
         % (N_MAX, all(len(C[n]) - 1 == binom(n, 2) for n in range(1, N_MAX + 1))))

    # A1: C_n(1) = Cat_n, with Cat_n derived from binom(2n,n)/(n+1).
    ones = [p_eval(C[n], 1) for n in range(N_MAX + 1)]
    cats = [catalan(n) for n in range(N_MAX + 1)]
    info("C_n(1) for n=0..10 : %s" % ones[:11])
    info("Cat_n  for n=0..10 : %s" % cats[:11])
    info("C_%d(1) = %d ; Cat_%d = %d" % (N_MAX, ones[N_MAX], N_MAX, cats[N_MAX]))
    check(ones == cats,
          "A1  C_n(1) = Cat_n for all 0 <= n <= %d (%d values)"
          % (N_MAX, N_MAX + 1))

    # A2: equation (5), F(z) = 1 + z F(z) F(qz), coefficientwise: the z^n
    #     coefficient of the product is sum_{l<n} C_l(q) q^{n-1-l} C_{n-1-l}(q).
    #     HONESTY NOTE (do not read this check as evidence about the paper):
    #     substituting k = n-1-l turns that sum into sum_k q^k C_k C_{n-1-k},
    #     which is verbatim the recurrence C was built from, so A2 CANNOT FAIL
    #     for any input.  The paper's assertion "equation (1) is equivalent to
    #     equation (5)" is itself this reindexing.  A2 and A3 below therefore
    #     validate the implementation (commutativity of p_mul, the shift
    #     bookkeeping), not the mathematics.
    bad = []
    for n in range(N_MAX + 1):
        if n == 0:
            rhs = [1]
        else:
            rhs = []
            for l in range(n):
                rhs = p_add(rhs, p_shift(p_mul(C[l], C[n - 1 - l]), n - 1 - l))
        if rhs != C[n]:
            bad.append(n)
    check(not bad,
          "A2  eq. (5), F_q(z) = 1 + z F_q(z) F_q(qz), holds "
          "coefficientwise through z^%d  [TAUTOLOGICAL: this is eq. (1) "
          "reindexed by k = n-1-l, so it cannot fail; implementation check "
          "only]" % N_MAX)
    if bad:
        info("first failure at n = %d" % bad[0])

    # A3: the same recurrence read the other way round.  This is equation (1)
    #     read with l = n-1-k; it agrees with the build order by commutativity
    #     of Z[q], so it is a consistency check on the implementation only.
    bad3 = [n for n in range(1, N_MAX + 1)
            if C[n] != _reindexed(C, n)]
    check(not bad3,
          "A3  sum_{l<n} q^l C_{n-1-l}(q) C_l(q) = C_n(q) for 1 <= n <= %d"
          "  [TAUTOLOGICAL, as A2: implementation check only]" % N_MAX)


def _reindexed(C, n):
    acc = []
    for l in range(n):
        acc = p_add(acc, p_shift(p_mul(C[n - 1 - l], C[l]), l))
    return acc


# =========================================================================
# 8.  CHECK BLOCK B -- order two:  equation (6) and the family (2m+1,2).
# =========================================================================

def check_block_b(C):
    print("--- B.  order two: C_n(-1) and the family (2m+1,2) ---------------")
    vals = [p_eval(C[n], -1) for n in range(N_MAX + 1)]
    info("derived C_n(-1), n = 0..13 : %s" % vals[:14])
    info("comparison target (P8)     : %s" % TARGET_MINUS1)
    info("(the PAPER prints no such table -- only the closed form, "
         "equation (6), which is what B2/B3 decide)")
    check(vals[:14] == TARGET_MINUS1,
          "B1  derived C_n(-1) for n = 0..13 matches the 14-term comparison "
          "target fixed here (NOT a figure printed in the paper)")

    # B2/B3: equation (6), derived for the whole range.
    ev = [(m, vals[2 * m]) for m in range(1, N_MAX // 2 + 1)]
    check(all(v == 0 for _, v in ev),
          "B2  C_{2m}(-1) = 0 for all 1 <= m <= %d" % (N_MAX // 2))
    od = [(m, vals[2 * m + 1], (-1) ** m * catalan(m))
          for m in range((N_MAX - 1) // 2 + 1)]
    check(all(v == w for _, v, w in od),
          "B3  C_{2m+1}(-1) = (-1)^m Cat_m for all 0 <= m <= %d"
          % ((N_MAX - 1) // 2))

    # B4: consistency with the mu = 2 residue used by the sweep: Phi_2 = q+1.
    info("Phi_2(q) = %s" % p_str(cyclotomic(2)))
    check(all(residue(C[n], 2) == p_trim([vals[n]]) for n in range(N_MAX + 1)),
          "B4  reduction of C_n(q) mod Phi_2 agrees with evaluation at q = -1")

    # B5: the arithmetic identity the paper uses for the conjectured value.
    ids = [(m, conj_rhs(2 * m + 1, 2), Fraction((2 * m + 1) * catalan(m)))
           for m in range(0, (N_MAX - 1) // 2 + 1)]
    info("(1/2)binom(2m+2,m+1) for m=0..6 : %s" % [str(a) for _, a, _ in ids[:7]])
    info("(2m+1)Cat_m          for m=0..6 : %s" % [str(b) for _, _, b in ids[:7]])
    check(all(a == b for _, a, b in ids),
          "B5  (1/2)binom(2m+2,m+1) = (2m+1)Cat_m for 0 <= m <= %d"
          % ((N_MAX - 1) // 2))

    # the family itself.  Every m >= 1 is a counterexample; m = 0 is not.
    mmax = (SWEEP_N - 1) // 2
    fam = [(m, vals[2 * m + 1], conj_rhs(2 * m + 1, 2)) for m in range(mmax + 1)]
    for m, lhs, rhs in fam[:6]:
        info("  n = %2d : C_n(-1) = %-6s   conjectured %-6s  %s"
             % (2 * m + 1, lhs, rhs, "differ" if lhs != rhs else "AGREE"))
    check(all(lhs != rhs for m, lhs, rhs in fam if m >= 1),
          "B6  C_{2m+1}(-1) != (1/2)binom(2m+2,m+1) for every 1 <= m <= %d "
          "(family (2m+1,2))" % mmax)
    check(fam[0][1] == fam[0][2],
          "B7  the two sides do agree at m = 0, i.e. at (n,mu) = (1,2)")


# =========================================================================
# 9.  CHECK BLOCK C -- order three:  equations (7) and (8), family (3r+2,3).
# =========================================================================

def check_block_c(C):
    print("--- C.  order three: C_n(omega) in Z[omega], family (3r+2,3) ------")
    info("Phi_3(q) = %s   (omega^2 = -1 - omega)" % p_str(cyclotomic(3)))
    coords = [omega_coords(C[n]) for n in range(N_MAX + 1)]
    for n in range(13):
        a, b = coords[n]
        info("  C_%-2d(omega) = %d + %d*omega" % (n, a, b))

    # C1: the closed form implied by equation (7), for EVERY n <= 30.  (The
    #     paper states equation (8) only along n = 3r+2; the coefficient
    #     extraction from equation (7) gives it for all n, so that is checked.)
    predicted = []
    for n in range(N_MAX + 1):
        top = (n - 2) // 3 if n >= 2 else -1
        predicted.append((1, sum(catalan(j) for j in range(top + 1))))
    check(coords == predicted,
          "C1  C_n(omega) = 1 + omega*sum_{j <= (n-2)/3} Cat_j for ALL "
          "0 <= n <= %d (from eq. (7), stronger than eq. (8))" % N_MAX)

    # C2: non-reality.  1 and omega are a Z-basis of Z[omega], so b != 0 means
    #     the value is not rational; concretely Im(a+b*omega) = b*sqrt(3)/2.
    check(all(coords[n][1] > 0 for n in range(2, N_MAX + 1)),
          "C2  the omega-coefficient of C_n(omega) is strictly positive for "
          "every 2 <= n <= %d, i.e. Im C_n(omega) = b*sqrt(3)/2 > 0" % N_MAX)
    check(all(is_real_at(C[n], 3) == (coords[n][1] == 0)
              for n in range(N_MAX + 1)),
          "C3  the conjugation test (q -> q^2 mod Phi_3) calls C_n(omega) real "
          "exactly when its omega-coefficient vanishes")

    # C4: the family (3r+2,3) itself, against the rational right-hand side.
    rmax = (SWEEP_N - 2) // 3
    ok, shown, seen = True, 0, 0
    for r in range(rmax + 1):
        n = 3 * r + 2
        a, b = coords[n]
        rhs = conj_rhs(n, 3)
        d = diff_mod_phi(C[n], rhs, 3)
        seen += 1
        if b == 0 or not d:
            ok = False
        if shown < 5:
            info("  n = %2d : C_n(zeta_3) = %d + %d*omega   conjectured %s   "
                 "residue of difference = %s"
                 % (n, a, b, rhs, p_str(d)))
            shown += 1
    # guard against a vacuous pass: an empty family would satisfy "ok".
    expect_members = len([n for n in range(SWEEP_N + 1) if n % 3 == 2])
    info("family members tested = %d (independently counted as the n <= %d "
         "with n = 2 mod 3: %d)" % (seen, SWEEP_N, expect_members))
    check(ok and seen == expect_members and seen > 0,
          "C4  for every r with 0 <= r <= %d -- all %d of them, count checked "
          "so the quantifier is not vacuous -- the value C_{3r+2}(zeta_3) is "
          "non-real and differs from (1/3)binom(3r+3,r+1) (family (3r+2,3))"
          % (rmax, expect_members))
    check(coords[2] == (1, 1) and conj_rhs(2, 3) == 1,
          "C5  the smallest instance: C_2(zeta_3) = 1 + zeta_3 against the "
          "conjectured (1/3)binom(3,1) = 1")


def is_prime(m):
    """Integer-only trial division (no float square root anywhere)."""
    if m < 2:
        return False
    d = 2
    while d * d <= m:
        if m % d == 0:
            return False
        d += 1
    return True


def odd_primes(upto):
    return [p for p in range(3, upto + 1) if is_prime(p)]


# =========================================================================
# 10.  CHECK BLOCK D -- odd prime orders:  the family (p-1,p).
# =========================================================================

def check_block_d(C):
    print("--- D.  odd prime orders: the family (p-1,p) ----------------------")
    ps = odd_primes(min(PRIME_MAX, N_MAX + 1))
    info("primes used: %s" % ps)

    # D1/D2: the two congruences, derived from the integers themselves.
    c1 = [(p, catalan(p - 1) % p) for p in ps]
    c2 = [(p, binom(2 * p - 1, p - 1) % p) for p in ps]
    info("Cat_{p-1} mod p       : %s" % [v for _, v in c1])
    info("p-1                   : %s" % [p - 1 for p in ps])
    info("binom(2p-1,p-1) mod p : %s" % [v for _, v in c2])
    check(all(v == p - 1 for p, v in c1),
          "D1  Cat_{p-1} = p-1 = -1 (mod p) for every odd prime p <= %d"
          % ps[-1])
    check(all(v == 1 for p, v in c2),
          "D2  binom(2p-1,p-1) = 1 (mod p) for every odd prime p <= %d"
          % ps[-1])
    check(all((2 * p - 1) * catalan(p - 1) == binom(2 * p - 1, p - 1)
              for p in ps),
          "D3  the exact identity (2p-1)Cat_{p-1} = binom(2p-1,p-1)")

    # D4: the ideal-theoretic step, made exact.  Z[zeta_p]/(1-zeta_p) = F_p
    # sends zeta_p to 1, so the image of C_{p-1}(zeta_p) - 1 is
    # C_{p-1}(1) - 1 = Cat_{p-1} - 1 = -2 (mod p), which is nonzero for p odd.
    ok4 = True
    for p in ps:
        R = residue(p_sub(C[p - 1], [1]), p)     # C_{p-1}(q) - 1 mod Phi_p
        img_reduced = p_eval(R, 1) % p
        img_direct = (catalan(p - 1) - 1) % p
        if img_reduced != img_direct or img_reduced != (-2) % p or img_reduced == 0:
            ok4 = False
        if p <= 7:
            info("  p = %2d : C_{p-1}(q)-1 mod Phi_p = %s ; image in F_p = %d "
                 "( = -2 mod %d )" % (p, p_str(R), img_reduced, p))
    check(ok4,
          "D4  the image of C_{p-1}(zeta_p) - 1 in Z[zeta_p]/(1-zeta_p) = F_p "
          "equals -2 != 0 for every odd prime p <= %d" % ps[-1])

    # D5: the sharper, purely algebraic statement -- the residue itself is a
    # nonzero element of Z[q]/Phi_p, so the two sides are unequal outright.
    nz = [p for p in ps if residue(p_sub(C[p - 1], [1]), p)]
    check(len(nz) == len(ps),
          "D5  C_{p-1}(q) - 1 has nonzero residue mod Phi_p for all %d odd "
          "primes p <= %d, so C_{p-1}(zeta_p) != 1 = (1/p)binom(p,1)"
          % (len(ps), ps[-1]))


# =========================================================================
# 11.  CHECK BLOCK E -- the exhaustive sweep (first Remark of the paper).
# =========================================================================

def admissible_pairs(nmax):
    """Every (n,mu) with 0 <= n <= nmax, mu >= 2 and mu | n+1."""
    out = []
    for n in range(nmax + 1):
        for mu in range(2, n + 2):
            if (n + 1) % mu == 0:
                out.append((n, mu))
    return out


def check_block_e(C):
    print("--- E.  exhaustive sweep over the admissible pairs, n <= %d -------"
          % SWEEP_N)
    pairs = admissible_pairs(SWEEP_N)
    info("derived number of admissible pairs = %d   (paper: %d)"
         % (len(pairs), PAPER_N_PAIRS))
    info("pairs with n = 0 : %d   (no mu >= 2 divides 1)"
         % sum(1 for n, _ in pairs if n == 0))
    check(len(pairs) == PAPER_N_PAIRS,
          "E1  there are exactly %d pairs (n,mu) with mu >= 2, mu | n+1, "
          "n <= %d" % (PAPER_N_PAIRS, SWEEP_N))

    holds, fails = [], []
    for (n, mu) in pairs:
        d = diff_mod_phi(C[n], conj_rhs(n, mu), mu)
        (holds if not d else fails).append((n, mu))
    info("pairs where the conjecture HOLDS : %s" % holds)
    info("number of failures = %d   (paper: %d)" % (len(fails), PAPER_N_FAILS))
    for (n, mu) in fails[:6]:
        r = residue(C[n], mu)
        info("  (n,mu) = (%2d,%2d) : C_n(zeta_mu) = [%s] ,  conjectured %s"
             % (n, mu, p_str(r), conj_rhs(n, mu)))
    check(holds == [PAPER_HOLDS_AT] and len(holds) == PAPER_N_HOLDS,
          "E2  C_n(q) - (1/mu)binom(n+1,(n+1)/mu) vanishes mod Phi_mu at "
          "exactly one pair, and that pair is (n,mu) = (1,2)")
    check(len(fails) == PAPER_N_FAILS,
          "E3  the conjecture fails at the remaining %d pairs"
          % PAPER_N_FAILS)

    # E4: the three families, as explicit sets, and their union.
    famA = set((2 * m + 1, 2) for m in range(1, (SWEEP_N - 1) // 2 + 1))
    famB = set((3 * r + 2, 3) for r in range((SWEEP_N - 2) // 3 + 1))
    famC = set((p - 1, p) for p in odd_primes(SWEEP_N + 1))
    union = famA | famB | famC
    info("|(2m+1,2)| = %d, |(3r+2,3)| = %d, |(p-1,p)| = %d"
         % (len(famA), len(famB), len(famC)))
    info("pairwise intersections: A&B = %s, A&C = %s, B&C = %s"
         % (sorted(famA & famB), sorted(famA & famC), sorted(famB & famC)))
    info("derived |union| = %d   (paper: %d)" % (len(union), PAPER_FAMILY_UNION))
    check(len(union) == PAPER_FAMILY_UNION,
          "E4  the three families account for exactly %d of the %d pairs"
          % (PAPER_FAMILY_UNION, PAPER_N_PAIRS))
    check(union <= set(fails),
          "E5  every pair of the three families is one of the %d failures"
          % len(fails))
    check(len(famC) == PAPER_PRIME_PAIRS,
          "E6  the odd-prime congruence reaches exactly %d pairs"
          % PAPER_PRIME_PAIRS)

    # E7: the non-reality observation, decided by exact conjugation.
    nonreal = [(n, mu) for (n, mu) in pairs if not is_real_at(C[n], mu)]
    real = [(n, mu) for (n, mu) in pairs if is_real_at(C[n], mu)]
    info("derived non-real pairs = %d   (paper: %d)"
         % (len(nonreal), PAPER_NONREAL))
    info("real pairs = %d, all with mu = 2 : %s"
         % (len(real), all(mu == 2 for _, mu in real)))
    check(len(nonreal) == PAPER_NONREAL,
          "E7  C_n(zeta_mu) is non-real at exactly %d of the %d pairs, so "
          "non-reality alone disposes of %d" % (PAPER_NONREAL, len(pairs),
                                               PAPER_NONREAL))
    check(all((n, mu) in set(fails) for (n, mu) in nonreal),
          "E8  every non-real pair is a failure (the right-hand side is "
          "rational, hence real)")
    return pairs, fails


# =========================================================================
# 12.  CHECK BLOCK F -- the Gaussian q-Catalan numbers (last Remark).
# =========================================================================

def check_block_f(pairs):
    print("--- F.  the Gaussian q-Catalan number Cat~_n(q) -------------------")
    G, T = gaussian_q_catalans(SWEEP_N)          # exact division => integral
    # A second, independent q-Pascal recurrence.  Re-multiplying G against the
    # table it was divided out of would be vacuous (and a non-integral quotient
    # raises inside to_int_poly rather than reaching a check), so F1 compares
    # against the ALTERNATIVE table: this can and will fail if either
    # recurrence, the division, or [n+1]_q is wrong.
    T2 = gaussian_pascal_alt(2 * SWEEP_N, SWEEP_N)
    same_tables = all(T[2 * n][n] == T2[2 * n][n] for n in range(SWEEP_N + 1))
    remult = all(p_mul(G[n], [1] * (n + 1)) == T2[2 * n][n]
                 for n in range(SWEEP_N + 1))
    intg = all(all(isinstance(c, int) for c in G[n]) for n in range(SWEEP_N + 1))
    info("qbinom(2n,n) from the two q-Pascal recurrences agree: %s"
         % same_tables)
    info("qbinom(4,2) = %s ; Cat~_2(q)*[3]_q = %s"
         % (p_str(T2[4][2]), p_str(p_mul(G[2], [1, 1, 1]))))
    check(remult and intg and same_tables,
          "F1  Cat~_n(q)*[n+1]_q = qbinom(2n,n) with Cat~_n(q) in Z[q], for "
          "every 0 <= n <= %d, the Gaussian binomial being recomputed by the "
          "independent recurrence qbinom(m,j) = q^{m-j}qbinom(m-1,j-1) "
          "+ qbinom(m-1,j)" % SWEEP_N)
    info("Cat~_2(q) = %s   (paper: 1 + q^2)" % p_str(G[2]))
    info("Cat~_3(q) = %s" % p_str(G[3]))
    check(G[2] == [1, 0, 1],
          "F2  Cat~_2(q) = 1 + q^2, as printed in the paper")

    # F3: at mu = 2 the Gaussian normalization reproduces the conjectured value.
    mu2 = [(n, mu) for (n, mu) in pairs if mu == 2]
    vals = [(n, p_eval(G[n], -1), conj_rhs(n, 2)) for (n, mu) in mu2]
    info("Cat~_n(-1) for odd n = 1..11 : %s" % [v for n, v, _ in vals[:6]])
    info("paper's printed values        : %s"
         % [PAPER_GAUSS_MU2[n] for n in sorted(PAPER_GAUSS_MU2)])
    check([v for n, v, _ in vals[:6]]
          == [PAPER_GAUSS_MU2[n] for n in sorted(PAPER_GAUSS_MU2)],
          "F3  Cat~_n(-1) at odd n <= 11 takes the six values the paper "
          "prints (1, 3, 10, 35, 126, 462)")
    check(all(v == r for _, v, r in vals),
          "F4  Cat~_n(-1) = (1/2)binom(n+1,(n+1)/2) at every one of the %d "
          "admissible pairs with mu = 2 and n <= %d" % (len(mu2), SWEEP_N))

    # F5: at mu >= 3 it never reproduces it.
    hi = [(n, mu) for (n, mu) in pairs if mu >= 3]
    agree = [(n, mu) for (n, mu) in hi
             if not diff_mod_phi(G[n], conj_rhs(n, mu), mu)]
    info("mu >= 3 pairs = %d ; of these the Gaussian value agrees at %d"
         % (len(hi), len(agree)))
    check((not agree) and len(hi) == PAPER_NONREAL,
          "F5  Cat~_n(zeta_mu) differs from (1/mu)binom(n+1,(n+1)/mu) at all "
          "%d admissible pairs with mu >= 3 (the universe is size-checked "
          "against the %d of the paper's Remark, so the quantifier is not "
          "vacuous).  NOTE: the paper claims this only at (2,3); the general "
          "statement is an extension made here, so a FAIL here need not be a "
          "defect in the paper" % (len(hi), PAPER_NONREAL))

    # F6: the paper's explicit order-three value.
    a, b = omega_coords(G[2])
    info("Cat~_2(zeta_3) = %d + %d*omega  = 1 + zeta_3^2 = 1/2 - sqrt(3)/2*i"
         % (a, b))
    alt = residue(p_add([1], zeta_pow(2, 3)), 3)   # 1 + zeta_3^2, independently
    info("independent route: 1 + zeta_3^2 reduces to %s" % p_str(alt))
    check((a, b) == (0, -1) and alt == [0, -1] and conj_rhs(2, 3) == 1,
          "F6  Cat~_2(zeta_3) = 1 + zeta_3^2 = -omega, whose omega-coefficient "
          "is -1 != 0, so it is non-real and differs from the conjectured 1")
    return G


# =========================================================================
# 13.  CHECK BLOCK G -- the bridge, equation (3), way 1: brute force.
#      This is the step the paper does not prove outright, so it is checked
#      here directly from Definitions 6.1 and 6.2.
# =========================================================================

def check_block_g(C):
    print("--- G.  brute force S_n(312), eq. (2), [FH], eq. (3) --------------")
    sets = [gen_312(n) for n in range(BRUTE_N + 1)]
    sizes = [len(s) for s in sets]
    info("|S_n(312)| for n = 0..%d : %s" % (BRUTE_N, sizes))
    info("Cat_n       for n = 0..%d : %s"
         % (BRUTE_N, [catalan(n) for n in range(BRUTE_N + 1)]))
    check(sizes == [catalan(n) for n in range(BRUTE_N + 1)],
          "G1  |S_n(312)| = Cat_n for 0 <= n <= %d (brute force over n! "
          "permutations)" % BRUTE_N)

    rows, bucket_sizes = [], []
    for n in range(BRUTE_N + 1):
        r, b = brute_triangle_row(n, sets[n])
        rows.append(r)
        bucket_sizes.append(b)
    info("|S_{8,k}(312)| for k = 0..8 : %s" % bucket_sizes[BRUTE_N])
    check(all(bucket_sizes[n][n] == 0 for n in range(1, BRUTE_N + 1)),
          "G2  S_{n,n}(312) is empty for 1 <= n <= %d  [DEFINITIONAL: pi_{n+1} "
          "does not exist for pi in S_n, and the buckets here are indexed by a "
          "position in 0..n-1, so this cannot fail]" % BRUTE_N)
    check(all(sum(bucket_sizes[n][:n]) == sizes[n]
              for n in range(1, BRUTE_N + 1)),
          "G3  S_{n,0}(312),...,S_{n,n-1}(312) partition S_n(312), so "
          "eq. (2), S'_{n,n}(312) = S_n(312), holds (1 <= n <= %d)" % BRUTE_N)
    check(all(bucket_sizes[n][0] == 1 for n in range(1, BRUTE_N + 1)),
          "G4  |S_{n,0}(312)| = 1, the decreasing permutation, matching the "
          "base case C_{n,0}(q) = q^binom(n,2)")

    # G5: the residue refinement of Definition 6.1 really refines inv.
    ok = True
    for n in range(BRUTE_N + 1):
        for pi in sets[n]:
            iv = inv_count(pi)
            for mu in range(2, 6):
                if sum(inv_residues(pi, mu)) != iv:
                    ok = False
    check(ok,
          "G5  sum_{i=1}^{mu} inv_i(pi) = inv(pi) for every pi in S_n(312), "
          "n <= %d, mu = 2..5  [STRUCTURAL: inv_residues walks the same (a,b) "
          "loop as inv_count and credits each inversion to exactly one class, "
          "so this cannot fail; it is the paper's own one-line argument]"
          % BRUTE_N)

    # G6: the [FH] interpretation C_n(q) = sum over S_n(312) of q^inv.
    check(all(rows[n][n] == C[n] for n in range(BRUTE_N + 1)),
          "G6  sum_{pi in S_n(312)} q^inv(pi) = C_n(q) coefficientwise for "
          "0 <= n <= %d (the [FH] interpretation, verified not assumed)"
          % BRUTE_N)

    # G7: equation (3) itself, from Definition 6.2 with all q_i = zeta_mu.
    bad = []
    for mu in range(2, 6):
        for n in range(BRUTE_N + 1):
            lhs = multivar_diagonal(n, mu, sets[n])
            rhs = residue(C[n], mu)
            if lhs != rhs:
                bad.append((n, mu))
    info("eq. (3) tested at %d (n,mu) combinations; mismatches: %s"
         % (4 * (BRUTE_N + 1), bad if bad else "none"))
    d3 = multivar_diagonal(3, 3, sets[3])
    info("C^{(3)}_{3,3} = %s  and  C_3(omega) = %s"
         % (p_str(d3), p_str(residue(C[3], 3))))
    check(not bad,
          "G7  C^{(mu)}_{n,n} = C_n(zeta_mu) for all 0 <= n <= %d and "
          "mu = 2,3,4,5 (eq. (3), the bridge the paper does not prove "
          "outright).  Its substance is G6: given G5, which is structural, "
          "eq. (3) reduces to the [FH] interpretation" % BRUTE_N)
    return rows


# =========================================================================
# 14.  CHECK BLOCK H -- the bridge, way 2: Wirdane's own triangle recurrence.
#      Independently supplies the source's unproved Proposition 2.4.
# =========================================================================

def check_block_h(C, rows):
    print("--- H.  Wirdane's triangle recurrence, C_{n,n}(q) = C_n(q) --------")
    tri = wirdane_triangle(TRI_N)
    for n in range(5):
        info("row %d : %s" % (n, [p_str(t) for t in tri[n]]))

    # H1: the recurrence's own diagonal reproduces C_n(q) -- Proposition 2.4.
    bad = [n for n in range(TRI_N + 1) if tri[n][n] != C[n]]
    check(not bad,
          "H1  C_{n,n}(q) = C_n(q) for 0 <= n <= %d, iterated from "
          "C_{n,0}=q^binom(n,2), C_{n,k}=C_{n,k-1}+q^{n-k-1}C_{n-1,k} "
          "(supplies [Wirdane, Prop. 2.4])" % TRI_N)
    if bad:
        info("first failure at n = %d" % bad[0])

    # H2: the recurrence agrees with Definition 6.2 everywhere brute force
    #     reaches, entry by entry -- not only on the diagonal.
    mism = [(n, k) for n in range(BRUTE_N + 1) for k in range(n + 1)
            if tri[n][k] != rows[n][k]]
    info("triangle entries compared: %d ; mismatches: %s"
         % (sum(n + 1 for n in range(BRUTE_N + 1)), mism if mism else "none"))
    check(not mism,
          "H2  every entry C_{n,k}(q) of the recurrence equals the "
          "enumeration sum_{pi in S'_{n,k}(312)} q^inv(pi), 0 <= k <= n <= %d"
          % BRUTE_N)

    # H3: the one polynomial the paper quotes from [Wirdane]'s table.
    info("derived C_{3,3}(q) = %s   (tabulated in [Wirdane]: %s)"
         % (p_str(tri[3][3]), p_str(PAPER_C33)))
    check(tri[3][3] == PAPER_C33 and rows[3][3] == PAPER_C33,
          "H3  C_{3,3}(q) = 1 + 2q + q^2 + q^3, the value tabulated in "
          "[Wirdane], from both the recurrence and the enumeration")
    return tri


# =========================================================================
# 15.  CHECK BLOCK I -- the remaining Remarks: the internal conflict with
#      [Wirdane, Rem. 6.7], the triple collision, mu = 1, minimality.
# =========================================================================

def check_block_i(C, tri):
    print("--- I.  Remark 6.7, the n=k=3 collision, mu=1, minimality ---------")

    # I1: setting every q_i = 1 in Definition 6.2 gives the classical Catalan
    #     triangle, whose diagonal is Cat_n.
    diag1 = [p_eval(tri[n][n], 1) for n in range(TRI_N + 1)]
    info("C_{n,n}(1) for n=0..%d : %s" % (TRI_N, diag1))
    check(diag1 == [catalan(n) for n in range(TRI_N + 1)],
          "I1  C_{n,n} = C_n(1) = Cat_n on the diagonal of the classical "
          "Catalan triangle, 0 <= n <= %d" % TRI_N)

    # I2: Remark 6.7 of [Wirdane] asserts C_n(-1) = (-1)^{floor(n/2)}Cat_n;
    #     false for every even n >= 2 since C_n(-1) = 0 but Cat_n != 0.
    ev = [(n, p_eval(C[n], -1), (-1) ** (n // 2) * catalan(n))
          for n in range(2, N_MAX + 1, 2)]
    info("even n : C_n(-1) = %s" % [v for _, v, _ in ev][:8])
    info("         Rem 6.7 predicts %s" % [w for _, _, w in ev][:8])
    check(all(v == 0 and w != 0 for _, v, w in ev),
          "I2  [Wirdane, Rem. 6.7] fails for every even n with 2 <= n <= %d: "
          "C_n(-1) = 0 while (-1)^floor(n/2) Cat_n != 0" % N_MAX)
    check(ev[0] == (2, 0, -2),
          "I3  the smallest instance is n = 2, where Rem. 6.7 predicts "
          "-Cat_2 = -2 against C_2(-1) = 0")

    # I4: three values collide at n = k = 3.
    s3 = gen_312(3)
    true_val = multivar_diagonal(3, 2, s3)          # C^{(2)}_{3,3}
    true_int = p_eval(C[3], -1)
    rem67 = (-1) ** (3 // 2) * catalan(3)
    conj66 = conj_rhs(3, 2)
    info("C^{(2)}_{3,3} (from Def 6.2) = %s ; C_3(-1) = %d ; "
         "Rem 6.7 -> %d ; Conj 6.6 -> %s"
         % (p_str(true_val), true_int, rem67, conj66))
    check(p_trim([true_int]) == true_val and true_int == -1
          and rem67 == -5 and conj66 == 3
          and len({true_int, rem67, Fraction(conj66)}) == 3,
          "I4  at n = k = 3 the three values are -1 (true), -5 (Rem. 6.7) and "
          "3 (Conj. 6.6), pairwise distinct")

    # I5: the degenerate order mu = 1 (zeta_1 = 1, right-hand side = 1).
    mu1 = [(n, catalan(n), Fraction(binom(n + 1, n + 1), 1))
           for n in range(N_MAX + 1)]
    firstbad = min(n for n, l, r in mu1 if l != r)
    info("mu = 1 : C_n(1) vs binom(n+1,n+1) agrees at n = %s, first failure "
         "at n = %d" % ([n for n, l, r in mu1 if l == r], firstbad))
    check(firstbad == 2,
          "I5  admitting mu = 1, the formula holds at n = 0,1 and already "
          "fails at n = 2")

    # I6: minimality of n = 2 among mu >= 2.
    low = admissible_pairs(1)
    info("admissible pairs with n <= 1 : %s" % low)
    check(low == [(1, 2)] and not diff_mod_phi(C[1], conj_rhs(1, 2), 2)
          and bool(diff_mod_phi(C[2], conj_rhs(2, 3), 3)),
          "I6  no mu >= 2 divides n+1 at n = 0; the only pair at n = 1 is "
          "(1,2), where both sides are 1; and (2,3) fails -- so n = 2 is the "
          "least n in a counterexample")

    # I7: the overlap the last Remark points out.  Decided on the two families
    #     as DERIVED sets -- the prime family from the primality test, the
    #     order-three family from the r-formula -- and not on constant-true
    #     comparisons against the literal pair (2,3), which no computed
    #     quantity could have falsified.
    famB = set((3 * r + 2, 3) for r in range((SWEEP_N - 2) // 3 + 1))
    famC = set((p - 1, p) for p in odd_primes(SWEEP_N + 1))
    overlap = famB & famC
    hits = [r for r in range((SWEEP_N - 2) // 3 + 1)
            if (3 * r + 2, 3) in famC]
    info("order-three family & prime family = %s ; the r values that land in "
         "the prime family = %s" % (sorted(overlap), hits))
    check(overlap == set([(2, 3)]) and hits == [0]
          and bool(diff_mod_phi(C[2], conj_rhs(2, 3), 3)),
          "I7  the prime family and the order-three family meet in exactly one "
          "pair: (n,mu) = (2,3), which is p = 3 of (p-1,p) and r = 0 of "
          "(3r+2,3), and the conjecture does fail there")


# =========================================================================
# 15b. CHECK BLOCK J -- the paper's CONCLUSION, in one place.
#      Blocks A..I verify many true facts one at a time, but none of them on
#      its own says "Conjecture 6.6 is false", which is the paper's result.
#      This one does, end to end and at the one pair where every link of the
#      chain lies inside the machine-verified range: the left side is built
#      from Definitions 6.1/6.2 by enumeration (not from equation (1)), the
#      right side from the conjecture's own formula.
# =========================================================================

def check_block_j():
    print("--- J.  the paper's CONCLUSION, end to end at the headline pair ----")
    s2 = gen_312(2)
    lhs = multivar_diagonal(2, 3, s2)          # C^{(3)}_{2,2} from Def 6.2
    rhs = conj_rhs(2, 3)                       # (1/3)binom(3,1)
    phi3 = cyclotomic(3)
    info("S'_{2,2}(312) = S_2(312) = %s" % [list(pi) for pi in s2])
    info("(inv_1,inv_2,inv_3) per permutation = %s"
         % [inv_residues(pi, 3) for pi in s2])
    info("C^{(3)}_{2,2} = %s   (that is 1 + zeta_3)" % p_str(lhs))
    info("Conjecture 6.6 predicts (1/3)binom(3,1) = %s" % rhs)
    d = p_mod(p_sub([Fraction(c) for c in lhs], [Fraction(rhs)]), phi3)
    info("difference in Z[q]/Phi_3 = %s ; nonzero: %s" % (p_str(d), bool(d)))
    info("conjugate of the left side = %s , so it is not real: %s"
         % (p_str(conjugate_residue(lhs, 3)), not is_real_at(lhs, 3)))
    check(lhs == [1, 1] and rhs == 1 and bool(d) and not is_real_at(lhs, 3),
          "J1  CONCLUSION: Conjecture 6.6 of [Wirdane] is FALSE.  Its left "
          "side at (n,mu) = (2,3), computed from Definitions 6.1 and 6.2 by "
          "enumerating S'_{2,2}(312) and never from equation (1), is "
          "1 + zeta_3; that is not real, hence cannot equal the rational "
          "(1/3)binom(3,1) = 1 the conjecture predicts")


# =========================================================================
# 16.  What these checks do NOT establish.  Printed, not hidden.
# =========================================================================

def print_gaps():
    print("--- GAPS: what passing these checks does NOT establish -----------")
    print("  G-A. Conjecture 6.6's left-hand side is C^{(mu)}_{n,n}, not")
    print("       C_n(zeta_mu).  The paper bridges them with its eq. (3),")
    print("       whose last step is the [FH] interpretation cited from a")
    print("       paywalled 1985 paper.  This program verifies that bridge")
    print("       only for 0 <= n <= %d with mu in {2,3,4,5} (check G7), and"
          % BRUTE_N)
    print("       univariately for n <= %d (check H1 with [Wirdane, Thm 3.3])."
          % TRI_N)
    print("       The paper's headline counterexample (n,mu) = (2,3) lies")
    print("       INSIDE that range, so the headline result is fully certified.")
    print("       But the 80-failure sweep and the family members with n > %d"
          % TRI_N)
    print("       or mu > 5 are certified here only as statements about")
    print("       C_n(zeta_mu); their transfer to Conjecture 6.6 rests on a")
    print("       citation this program cannot check.")
    print("  G-B. The three families are INFINITE.  Only m <= %d, r <= %d and"
          % ((SWEEP_N - 1) // 2, (SWEEP_N - 2) // 3))
    print("       p <= %d are tested.  The paper's own proofs of the general" % PRIME_MAX)
    print("       statements -- F(z)+F(-z)=2 at order two, the G(z) identity")
    print("       behind eq. (7), and the ideal argument mod (1-zeta_p) --")
    print("       are NOT machine-checked; only their consequences are, and")
    print("       only in the finite ranges above.")
    print("  G-C. Source fidelity is a document fact.  Beyond the statement")
    print("       -> number map printed below, the TRANSCRIPTIONS of Def 6.1,")
    print("       Def 6.2, Rem 6.7 and of the statement of Conjecture 6.6")
    print("       itself are unverifiable without the e-print.  One mitigation")
    print("       is recorded: eq. (3) is insensitive to which endpoint of")
    print("       an inversion assigns its residue class, since any assignment")
    print("       of each inversion to exactly one class gives sum_i inv_i =")
    print("       inv, so a slip in Def 6.1 of that kind would not affect it.")
    print("       Exactly ONE datum in this program cross-validates the")
    print("       transcription at all: check H3, where the C_{3,3}(q)")
    print("       tabulated in [Wirdane] is reproduced independently, both")
    print("       from the triangle recurrence and from the enumeration of")
    print("       S'_{3,3}(312).  Every other quoted string -- Conj 6.6,")
    print("       Rem 6.7, Def 6.1, Def 6.2, and the attributions of Prop 2.4")
    print("       and Rem 6.4 -- is uncorroborated here, and the recipe below")
    print("       is the only route to checking it.")
    print("  G-D. Some checks assert MORE than the paper does: B1 (a table")
    print("       fixed here, not printed by the paper), C1 and C2 (all")
    print("       n <= %d, the paper says n = 3r+2), F4 (n <= %d, the paper"
          % (N_MAX, SWEEP_N))
    print("       says n <= 11), F5 (all mu >= 3 pairs, the paper says (2,3))")
    print("       and E7.  A FAIL in those is a defect in a target fixed here")
    print("       or in an unproved extrapolation, not necessarily in the paper.")
    print("  G-E. Four checks are structural and carry no evidence about the")
    print("       paper: A2 and A3 (eq. (1) reindexed), G2 (definitional")
    print("       emptiness), G5 (a partition of the same loop).  G7 is forced")
    print("       by G5 together with G6.  Discount them: of the %d checks,"
          % len(RESULTS))
    print("       %d are implementation checks." % 4)
    print("")


# =========================================================================
# 17.  main
# =========================================================================

def main():
    print("=" * 74)
    print("verify.py -- Counterexamples to Wirdane's Cyclotomic Catalan")
    print("                Conjecture (Conjecture 6.6 of arXiv:2605.14682v1)")
    print("=" * 74)
    print("Exact arithmetic only: Python ints and fractions.Fraction.")
    print("No floating point appears in any certificate.  Values in")
    print("Z[zeta_mu] are canonical residues in Z[q]/(Phi_mu(q)), so equality")
    print("of residues is equality of algebraic numbers.")
    print("")

    print("Building C_n(q) for 0 <= n <= %d from the recurrence, eq. (1) ..."
          % N_MAX)
    C = q_catalan(N_MAX)
    print("done; deg C_%d = %d, %d nonzero coefficients."
          % (N_MAX, len(C[N_MAX]) - 1, sum(1 for c in C[N_MAX] if c)))
    print("")

    check_block_a(C)
    print("")
    check_block_b(C)
    print("")
    check_block_c(C)
    print("")
    check_block_d(C)
    print("")
    pairs, _fails = check_block_e(C)
    print("")
    check_block_f(pairs)
    print("")
    rows = check_block_g(C)
    print("")
    tri = check_block_h(C, rows)
    print("")
    check_block_i(C, tri)
    print("")
    check_block_j()
    print("")

    print_gaps()

    print("--- source fidelity (a document fact, NOT one of the checks) ------")
    print("      The statement -> number map the paper relies on cannot be")
    print("      decided by computation; a reader confirms it against")
    print("      the numbered statements of arXiv:2605.14682v1, either by")
    print("      reading them off the abs page / PDF at")
    print("        https://arxiv.org/abs/2605.14682v1")
    print("      or from the LaTeX source, whose main file name is not known")
    print("      to this program and must be located rather than assumed:")
    print("        curl -sL https://arxiv.org/e-print/2605.14682v1 -o w.tar.gz")
    print("        mkdir w && tar xzf w.tar.gz -C w")
    print("        grep -l documentclass w/*.tex        # the main file, MAIN")
    print("        <any LaTeX engine> MAIN.tex          # then read off the")
    print("                                             # numbering below")
    for what, num in SOURCE_STATEMENTS:
        print("        %-44s expected: %s" % (what, num))
    print("")

    nfail = sum(1 for ok, _ in RESULTS if not ok)
    n = len(RESULTS)
    if nfail:
        print("failed checks:")
        for ok, t in RESULTS:
            if not ok:
                print("   " + t)
        print("VERDICT: %d OF %d CHECKS FAILED" % (nfail, n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
