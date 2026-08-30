#!/usr/bin/env python3
"""verify.py -- checks every computational claim of

    "A Counterexample to the Singh--Singh Width-2 ROABP Modular-Stability Proposition"

against the objects PRINTED IN THAT PAPER.  Python 3.9+, standard library only, no network, no
external data file, exact integer arithmetic throughout; there is no floating-point decision
anywhere below (the two decimal numbers 80.45... and 518.74... that the paper prints for
131^(9/10) and 1039^(9/10) are checked in the equivalent exact integer form (r-1)^10 > r^9).

WHAT IS TRANSCRIBED, AND FROM WHERE.  Everything this program consumes is a literal in this file
and every one of those literals is printed in the paper:

  * the substitution of the source's Definition 4.7.1, quoted in Section 1 of the paper:
        Gamma_g : x_i |--> lambda^(g^i)   in   R_r = K[lambda]/<lambda^r - 1>,
    implemented in gamma_vector() as a FULL length-r coefficient vector over the basis
    lambda^0, ..., lambda^(r-1).  Nothing here uses a derived criterion for C_g = 0; the derived
    criterion (exponent congruence) is implemented SEPARATELY, in bad_set_exp(), and the two are
    compared against each other.
  * the objects: C = x_1 - x_r  (Theorem 1 / Corollary 3 of the paper) and
                 C = x_1 - x_(1+k), k = (r-1)/2  (Remark, the n < r witness);
  * the primes 31, 131, 1039, 4099, 10007; the class parameters w = 2, d = 1; the source's own
    illustrative threshold R(w,d) = (w*d)^7 = 128; epsilon = 1/10; and the six primes
    1009, 1013, 1019, 1021, 1031, 1033 at which the halved witness fails.

The paper's Theorem 1 is a two-line argument in the cyclic group Z_r^*; this program is
corroboration, not a premise.  It prints one `PASS <name>` line per check, closes with
`VERDICT: ALL <n> CHECKS PASS`, and exits 0 if and only if every check passed.
"""

import sys
from math import gcd

# --------------------------------------------------------------------------------------------
# 0. THE LITERALS, ALL PRINTED IN THE PAPER
# --------------------------------------------------------------------------------------------
W, D = 2, 1                      # width and individual degree of the witness ROABP
THRESHOLD = (W * D) ** 7         # the source's own illustrative R(w,d); = 128 at w=2, d=1
EPS_DEN = 10                     # epsilon = 1/EPS_DEN
PRIMES_RING = (31, 131, 1039)    # primes at which the full ring engine is run
PRIMES_EXP = (31, 131, 1039, 4099, 10007)   # primes at which the congruence route is run
CHARS = (0, 2, 3)                # characteristic 0 (Q / Z) and the two smallest primes
QR_PRIMES = (1039, 4099)         # the n < r witness of the paper's Remark
QR_FAILING_PRIMES = (1009, 1013, 1019, 1021, 1031, 1033)
LEAST_R_FOR_EPS = 7             # paper: (r-1)^10 > r^9 holds for every integer r >= 7
LEAST_PRIME_ABOVE_THRESHOLD = 131
NOT_RE_RUN = []

# --------------------------------------------------------------------------------------------
# 1. THE CHECK HARNESS
# --------------------------------------------------------------------------------------------
_passed = 0
_failed = 0


def check(name, ok, detail=''):
    """One check.  A check is a claim of the paper; `detail` prints the numbers behind it."""
    global _passed, _failed
    if ok:
        _passed += 1
        print('PASS %s  %s' % (name, detail))
    else:
        _failed += 1
        print('FAIL %s  %s' % (name, detail))


# --------------------------------------------------------------------------------------------
# 2. POLYNOMIALS
# --------------------------------------------------------------------------------------------
# A polynomial is a list of (coefficient, monomial) pairs; a monomial is a dict {var: exponent}
# with 1-based variable indices.  Coefficients are exact Python ints.

def poly_diff_of_vars(i, j):
    """The paper's object x_i - x_j."""
    return [(1, {i: 1}), (-1, {j: 1})]


def poly_sum_of_vars(idxs):
    """x_{i1} + x_{i2} + ... -- used for the controls."""
    return [(1, {i: 1}) for i in idxs]


def monomial_count(poly):
    return len(poly)


def is_zero_polynomial(poly, char):
    """Is `poly` the zero polynomial in K[x_1..x_n]?  Collect equal monomials, reduce, compare."""
    acc = {}
    for c, mon in poly:
        key = tuple(sorted(mon.items()))
        acc[key] = acc.get(key, 0) + c
    for v in acc.values():
        if (v % char if char else v) != 0:
            return False
    return True


# --------------------------------------------------------------------------------------------
# 3. THE RING ENGINE -- Definition 4.7.1, implemented literally
# --------------------------------------------------------------------------------------------
def gamma_vector(poly, r, g, char, nested=True):
    """C_g as the FULL length-r coefficient vector over the basis lambda^0, ..., lambda^(r-1).

    nested=True  -> x_i |--> lambda^(g^i)   (the source's Definition 4.7.1, printed p.26)
    nested=False -> x_i |--> lambda^(g*i)   (the alternative reading; the paper's reading control)
    """
    vec = [0] * r
    for c, mon in poly:
        e = 0
        for v, a in mon.items():
            step = pow(g, v, r) if nested else (g * v) % r
            e = (e + a * step) % r
        vec[e] += c
    if char:
        vec = [x % char for x in vec]
    return vec


def is_zero_in_ring(poly, r, g, char, nested=True):
    return not any(gamma_vector(poly, r, g, char, nested))


def bad_set_ring(poly, r, char, nested=True):
    """B_r(C) computed by the ring engine: brute force over all of Z_r^*."""
    return set(g for g in range(1, r) if is_zero_in_ring(poly, r, g, char, nested))


def bad_set_exp(i, j, r):
    """B_r(x_i - x_j) by the derived congruence criterion g^i = g^j mod r.  Independent route."""
    return set(g for g in range(1, r) if pow(g, i, r) == pow(g, j, r))


# --------------------------------------------------------------------------------------------
# 4. SYMBOLIC 2x2 DIAGONAL ROABP -- Proposition 2 of the paper
# --------------------------------------------------------------------------------------------
def sp_mul(a, b):
    """Multiply two symbolic polynomials held as {monomial_key: coeff}, monomial_key sorted tuple."""
    out = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            m = dict(ka)
            for v, e in kb:
                m[v] = m.get(v, 0) + e
            k = tuple(sorted(m.items()))
            out[k] = out.get(k, 0) + ca * cb
    return {k: c for k, c in out.items() if c != 0}


def sp_add_scaled(a, b, s):
    out = dict(a)
    for k, c in b.items():
        out[k] = out.get(k, 0) + s * c
    return {k: c for k, c in out.items() if c != 0}


def roabp_width2_diagonal(i, j, n):
    """Run the paper's explicit construction and return (computed polynomial, max entry degree).

    M_t = diag(alpha_t, beta_t) with alpha_t = x_t if t == i else 1, beta_t = x_t if t == j else 1,
    u = (1,1), v = (1,-1).  The product of diagonal matrices is entrywise, so we carry the two
    diagonal entries as symbolic polynomials and combine at the end.
    """
    one = {(): 1}
    prod1, prod2 = dict(one), dict(one)
    maxdeg = 0
    for t in range(1, n + 1):
        alpha = {((t, 1),): 1} if t == i else dict(one)
        beta = {((t, 1),): 1} if t == j else dict(one)
        for entry in (alpha, beta):
            for k in entry:
                maxdeg = max(maxdeg, max([e for _v, e in k] or [0]))
        prod1 = sp_mul(prod1, alpha)
        prod2 = sp_mul(prod2, beta)
    return sp_add_scaled(prod1, prod2, -1), maxdeg


def sp_of_poly(poly):
    out = {}
    for c, mon in poly:
        k = tuple(sorted(mon.items()))
        out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c != 0}


# --------------------------------------------------------------------------------------------
# 5. SMALL NUMBER THEORY (exact, standard library only)
# --------------------------------------------------------------------------------------------
def is_prime(m):
    if m < 2:
        return False
    if m % 2 == 0:
        return m == 2
    f = 3
    while f * f <= m:
        if m % f == 0:
            return False
        f += 2
    return True


def primes_up_to(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b'\x00\x00'
    p = 2
    while p * p <= n:
        if sieve[p]:
            sieve[p * p::p] = bytearray(len(sieve[p * p::p]))
        p += 1
    return [i for i in range(n + 1) if sieve[i]]


def beats_bound(size, r, eps_den):
    """size > r^(1 - 1/eps_den), in EXACT integers: size^eps_den > r^(eps_den - 1)."""
    return size ** eps_den > r ** (eps_den - 1)


def quadratic_residues(r):
    return set(pow(x, 2, r) for x in range(1, r))


def is_subgroup(S, r):
    if not S or 1 not in S:
        return False
    for a in S:
        for b in S:
            if (a * b) % r not in S:
                return False
    return True


# --------------------------------------------------------------------------------------------
# 6. THE CHECKS
# --------------------------------------------------------------------------------------------
def main():
    print('# objects, primes and constants are the literals of block 0 below, all printed in the paper')
    print('# w = %d, d = %d, R(w,d) = (w*d)^7 = %d, epsilon = 1/%d'
          % (W, D, THRESHOLD, EPS_DEN))
    print()

    # --- 6.1  the ring R_r and its basis ---------------------------------------------------
    r0 = 31
    ok = True
    npairs = 0
    for char in CHARS:
        for a in range(r0):
            for b in range(a + 1, r0):
                v = [0] * r0
                v[a] += 1
                v[b] -= 1
                if char:
                    v = [x % char for x in v]
                if not any(v):
                    ok = False
                npairs += 1
        # and the zero vector must read as zero, in the same code path
        if any([0] * r0):
            ok = False
    check('ring-basis-independence',
          ok and npairs == 3 * (r0 * (r0 - 1) // 2),
          'lambda^a - lambda^b != 0 in R_31 at all %d = 465 x 3 combinations of a pair '
          '0 <= a < b < 31 with a characteristic in {0, 2, 3}, so for a two-term C the condition '
          'C_g = 0 is exactly an exponent congruence mod r, in every characteristic' % npairs)

    # --- 6.2  the reading of Definition 4.7.1 is load-bearing -------------------------------
    prod_reading_sizes = {}
    for r in PRIMES_EXP:
        C = poly_diff_of_vars(1, r)
        prod_reading_sizes[r] = sum(
            1 for g in range(1, r) if ((g * 1) % r) == ((g * r) % r))
    check('reading-is-nested-not-product',
          all(v == 0 for v in prod_reading_sizes.values()),
          'under the alternative reading x_i -> lambda^(g*i) the object x_1 - x_r has an EMPTY bad '
          'set at every r in %s (g*r = 0 mod r never equals g), so the nested exponent of '
          'Definition 4.7.1 is load-bearing and this note stands or falls with it'
          % (list(PRIMES_EXP),))
    small = 131
    C = poly_diff_of_vars(1, small)
    check('reading-control-ring-agrees',
          bad_set_ring(C, small, 2, nested=False) == set()
          and bad_set_ring(C, small, 2, nested=True) == set(range(1, small)),
          'ring engine at r = %d over F_2: bad set is empty under the product reading and all of '
          'Z_r^* (%d elements) under the nested reading' % (small, small - 1))

    # --- 6.3  the witness is in the class Proposition 4.7.3 speaks about --------------------
    cases = [(1, 31, 31), (1, 131, 131), (1, 520, 520)]
    ok, details = True, []
    for i, j, n in cases:
        got, maxdeg = roabp_width2_diagonal(i, j, n)
        want = sp_of_poly(poly_diff_of_vars(i, j))
        if got != want or maxdeg != D:
            ok = False
        details.append('n=%d j=%d' % (n, j))
    check('roabp-width-2-diagonal-realisation', ok,
          'the explicit product u^T diag-M_1 ... diag-M_n v of Proposition 2 evaluates symbolically '
          'to exactly x_i - x_j for %s, every matrix entry of degree <= %d in its own variable, '
          'width %d, diagonal' % (', '.join(details), D, W))
    C131 = poly_diff_of_vars(1, 131)
    check('witness-nonzero-and-two-term',
          all(not is_zero_polynomial(C131, c) for c in CHARS) and monomial_count(C131) == 2,
          'x_1 - x_131 is a nonzero polynomial over char 0, 2 and 3 and has exactly 2 monomials')
    check('witness-support-density',
          monomial_count(C131) == 2 and 2 < 2 ** 131,
          '2 monomials out of the 2^n multilinear monomials in n = 131 variables: the object lies '
          'in the constant-width, extremely-low-density class that the source\'s own p.27 Remark '
          'excludes from the conjecture, which is why the correction is to Proposition 4.7.3 (p.28) '
          'and to the printed quantifier form, not to the conjecture as described')

    # --- 6.4  the main theorem, ring engine, three characteristics --------------------------
    ring_sets = {}
    for r in PRIMES_RING:
        C = poly_diff_of_vars(1, r)
        for char in CHARS:
            assert char == 0 or r % char != 0, 'Definition 4.7.1 needs gcd(r, char K) = 1'
            B = bad_set_ring(C, r, char)
            ring_sets[(r, char)] = B
            check('ring-bad-set-r%d-char%d' % (r, char),
                  B == set(range(1, r)) and len(B) == r - 1,
                  'C = x_1 - x_%d, n = %d: B_r(C) = Z_r^* exactly, |B| = %d = r - 1, the maximum '
                  'possible' % (r, r, len(B)))

    # --- 6.5  the main theorem, independent congruence route --------------------------------
    exp_sets = {}
    for r in PRIMES_EXP:
        B = bad_set_exp(1, r, r)
        exp_sets[r] = B
        check('congruence-bad-set-r%d' % r,
              B == set(range(1, r)) and len(B) == r - 1,
              'g^1 = g^r mod %d for all %d elements of Z_r^*, so |B_r(x_1 - x_r)| = r - 1 = %d'
              % (r, r - 1, r - 1))
    check('two-routes-agree',
          all(ring_sets[(r, char)] == exp_sets[r] for r in PRIMES_RING for char in CHARS),
          'the full length-r coefficient-vector engine and the exponent-congruence criterion give '
          'the SAME bad set at every (r, char) in %s x %s'
          % (list(PRIMES_RING), list(CHARS)))
    check('fermat-mechanism',
          all(pow(g, r, r) == g % r for r in PRIMES_EXP for g in range(1, r)),
          'g^r = g mod r for every g in Z_r^* at r in %s -- Fermat\'s little theorem, the whole '
          'mechanism of Theorem 1' % (list(PRIMES_EXP),))

    # --- 6.6  the subgroup lemma, over every admissible index difference -------------------
    for r in PRIMES_RING:
        ok = all(len(bad_set_exp(1, 1 + delta, r)) == gcd(delta, r - 1)
                 for delta in range(1, r))
        check('subgroup-order-lemma-r%d' % r, ok,
              '|B_%d(x_i - x_j)| = gcd(j - i, r - 1) for every one of the %d index differences '
              'j - i = 1..%d' % (r, r - 1, r - 1))
    ok = True
    for i in range(1, 32):
        for j in range(i + 1, 32):
            if len(bad_set_ring(poly_diff_of_vars(i, j), 31, 2)) != gcd(j - i, 30):
                ok = False
    check('subgroup-order-lemma-ring-r31', ok,
          'the ring engine confirms |B_31(x_i - x_j)| = gcd(j - i, 30) for all %d pairs '
          '1 <= i < j <= 31 over F_2' % (31 * 30 // 2))
    check('bad-set-is-a-subgroup',
          all(is_subgroup(bad_set_exp(1, 1 + delta, 131), 131) for delta in range(1, 131)),
          'B_131(x_i - x_j) is closed under multiplication and contains 1 for every index '
          'difference, as the kernel of g -> g^(j-i) must be')

    # --- 6.7  the inequalities the paper prints --------------------------------------------
    check('inequality-at-r131-exact',
          beats_bound(130, 131, EPS_DEN)
          and 130 ** 10 == 1378584918490000000000
          and 131 ** 9 == 11361656654439817571,
          '130^10 = %d > 131^9 = %d, i.e. |B| = 130 > 131^(9/10), in exact integers'
          % (130 ** 10, 131 ** 9))
    check('threshold-and-least-admissible-prime',
          THRESHOLD == 128
          and is_prime(LEAST_PRIME_ABOVE_THRESHOLD)
          and LEAST_PRIME_ABOVE_THRESHOLD >= THRESHOLD
          and not any(is_prime(m) for m in range(THRESHOLD, LEAST_PRIME_ABOVE_THRESHOLD)),
          '(w*d)^7 = %d, and %d is prime and is the LEAST prime >= %d (128, 129, 130 are '
          'composite), so the witness clears the source\'s own illustrative threshold'
          % (THRESHOLD, LEAST_PRIME_ABOVE_THRESHOLD, THRESHOLD))
    fails_below = [r for r in range(2, LEAST_R_FOR_EPS) if beats_bound(r - 1, r, EPS_DEN)]
    holds_above = all(beats_bound(r - 1, r, EPS_DEN) for r in range(LEAST_R_FOR_EPS, 20001))
    check('inequality-holds-from-r7',
          not fails_below and holds_above,
          '(r-1)^10 > r^9 fails for every integer 2 <= r <= 6 and holds for every integer '
          '7 <= r <= 20000, so the least r is exactly %d -- far below the sufficient condition '
          '2^(1/eps) = 1024 the paper quotes' % LEAST_R_FOR_EPS)
    sieve = primes_up_to(25000)
    pset = set(sieve)
    ok, tested = True, 0
    for m in range(1, 15):
        lo = 2 ** m
        for r in range(lo + 1, min(lo + 3001, 25001)):
            if r in pset:
                tested += 1
                if not (r - 1) ** m > r ** (m - 1):
                    ok = False
    check('sufficient-condition-r-above-2-to-the-1-over-eps', ok and tested > 0,
          'for eps = 1/m with m = 1..14, (r-1)^m > r^(m-1) at every one of the %d primes r with '
          '2^m < r <= 2^m + 3000 -- the proof step r - 1 >= r/2 > r^(1-eps) of Corollary 3'
          % tested)

    # --- 6.8  the n < r witness of the paper's Remark --------------------------------------
    for r in QR_PRIMES:
        k = (r - 1) // 2
        B = bad_set_exp(1, 1 + k, r)
        QR = quadratic_residues(r)
        check('qr-witness-r%d' % r,
              len(B) == k == gcd(k, r - 1) and B == QR and beats_bound(len(B), r, EPS_DEN),
              'C = x_1 - x_%d on n = %d variables: |B| = %d = (r-1)/2 = gcd(k, r-1), B is exactly '
              'the set of quadratic residues mod %d (Euler\'s criterion, computed independently), '
              'and %d^10 > %d^9' % (1 + k, k + 1, len(B), r, len(B), r))
    r = 1039
    k = (r - 1) // 2
    check('qr-witness-ring-r1039',
          bad_set_ring(poly_diff_of_vars(1, 1 + k), r, 2) == quadratic_residues(r),
          'the full ring engine at r = 1039 over F_2 returns exactly the 519 quadratic residues '
          'for C = x_1 - x_520, matching the congruence route and Euler\'s criterion')
    check('qr-witness-least-prime-is-1039',
          all(not beats_bound((p - 1) // 2, p, EPS_DEN) for p in QR_FAILING_PRIMES)
          and beats_bound(519, 1039, EPS_DEN)
          and all(not beats_bound((p - 1) // 2, p, EPS_DEN)
                  for p in sieve if 3 <= p < 1039),
          'for the HALVED witness ((r-1)/2)^10 > r^9 fails at every odd prime below 1039, '
          'including %s -- two of them above 2^10 = 1024 -- and first holds at r = 1039, where '
          '519 > 1039^(9/10). The naive threshold 2^(1/eps) is NOT correct for this object.'
          % (list(QR_FAILING_PRIMES),))
    check('qr-witness-leaves-algorithm-1-intact',
          all((p - 1) // 2 <= p - 2 for p in QR_PRIMES)
          and all(p - 1 > p - 2 for p in QR_PRIMES),
          'a good parameter still exists for the halved witness, since (r-1)/2 <= r-2 at r in %s, '
          'while the n >= r object of Theorem 1 has |B| = r-1 > r-2 and therefore removes EVERY '
          'good parameter; only the latter breaks the single step of Algorithm 1 (printed p.40)'
          % (list(QR_PRIMES),))

    # --- 6.9  the controls -----------------------------------------------------------------
    B = bad_set_ring(poly_sum_of_vars([1, 2]), 1039, 2)
    check('control-forced-positive-x1-plus-x2',
          B == {1} and len(B) == 1,
          'C = x_1 + x_2 at r = 1039 over F_2 MUST have B = {1}, since g^1 = g^2 mod r iff g = 1; '
          'the engine returned B = {1}, |B| = 1')
    ok = all(bad_set_ring([(1, {1: 1})], 1039, c) == set() for c in CHARS)
    check('control-proved-silent-x1', ok,
          'C = x_1 can never vanish under Gamma_g (it maps to a single basis vector), and the '
          'engine returned |B_1039(x_1)| = 0 over char 0, 2 and 3')
    C3 = poly_sum_of_vars([1, 2, 3])
    b2, b3, b0 = (bad_set_ring(C3, 31, 2), bad_set_ring(C3, 31, 3), bad_set_ring(C3, 31, 0))
    check('control-characteristic-sensitivity',
          b3 == {1} and b2 == set() and b0 == set(),
          'C = x_1 + x_2 + x_3 at r = 31 sends all three monomials to lambda^1 at g = 1 with '
          'coefficient sum 3, so g = 1 must be bad over F_3 and good over F_2 and char 0: the '
          'engine returned |B| = 1 with B = {1} over F_3, and 0 over both others -- the engine is '
          'not blind to characteristic')
    B = bad_set_ring(poly_diff_of_vars(1, 2), 31, 2)
    check('control-forced-negative-difference-1',
          B == {1} and gcd(1, 30) == 1,
          'C = x_1 - x_2 has index difference 1, so the lemma forces |B| = gcd(1, 30) = 1; the '
          'engine returned B = {1}. A decider tuned to manufacture large bad sets would not '
          'return a singleton here.')

    # --- 6.10  scope ----------------------------------------------------------------------
    NOT_RE_RUN.append(
        'NOT RE-RUN: the full length-r ring engine was run only at r = 31, 131 and 1039. The two '
        'largest primes the paper prints, r = 4099 and r = 10007, are checked by the exponent-'
        'congruence route alone; the two routes are shown to agree at the three smaller primes, '
        'but no coefficient vector of length 4099 or 10007 is built here.')
    NOT_RE_RUN.append(
        'NOT RE-RUN: nothing here tests any ROABP of width >= 3, and nothing here tests any object '
        'of non-degenerate coefficient support density. The regime the source DESCRIBES in its '
        'p.27 Remark is therefore untouched by this program as it is by the paper, and remains '
        'open. This program exhibits no counterexample there and refutes nothing there.')
    NOT_RE_RUN.append(
        'NOT RE-RUN: no search over a family of ROABPs is performed, and in particular the small-'
        'parameter census at r = 31, n = 5 that appears in this result\'s internal record is NOT '
        'reproduced here and is NOT offered as evidence anywhere. A sweep at n << r is '
        'structurally blind to the family of Theorem 1, so it could not bear on the surviving '
        'regime in either direction; Section 3 of the paper says so.')
    NOT_RE_RUN.append(
        'NOT RE-RUN: the source PDF (arXiv:2602.13449v1, 1,086,923 bytes, PDF-only) is not fetched '
        'or parsed. Conjecture 4.7.2, Proposition 4.7.3, the p.27 Remark, Example 4.7.3, '
        'Definition 4.7.1, Algorithm 1 and the Appendix B table are TRANSCRIPTIONS from the '
        'printed pages named in the paper; this program cannot re-check a transcription, and a '
        'reader who doubts one should compare it against the PDF by eye.')
    NOT_RE_RUN.append(
        'NOT RE-RUN: no prior-art search is performed here. The paper states that the witness '
        'OBJECT is the source\'s own Example 4.7.3 at S = {1}, S\' = {r}, and it claims novelty '
        'for neither the object nor the count |B_r(C)| = r - 1. This program checks the count, '
        'not the literature.')

    for line in NOT_RE_RUN:
        print(line)

    if _failed:
        print('VERDICT: %d of %d CHECKS FAILED' % (_failed, _passed + _failed))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _passed)
    return 0


if __name__ == '__main__':
    sys.exit(main())
