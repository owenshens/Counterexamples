#!/usr/bin/env python3
"""
verify.py -- checks every computational claim of the accompanying paper

    "Polynomial Certificates for the First Two Partial Quotients of ( (3)-S_N)/2, and the Third for N 2artial Quotients of
     sqrt((zeta(3) - S_N)/2)"

Python 3.9+, STANDARD LIBRARY ONLY (fractions, math.isqrt).  No third-party
package, no external data file, no floating point in any decision: every
comparison below is between exact integers or exact Fractions.

WHAT IT READS.  The block `PAPER` immediately below is transcribed verbatim
from the paper -- the two antidifference polynomials A and B, the four
certificate coefficient vectors, the residue table of Section 4, the four
small-index brackets of Section 5, and the control values of Section 6.  The
program then re-derives each of those objects from the definitions alone and
compares.  Nothing is read from disk.

WHAT IT PRINTS.  One `PASS <name> <detail>` line per check, a block of
`NOT RE-RUN:` lines saying what is *not* covered, and a closing
`VERDICT: ALL <n> CHECKS PASS`.  Exit status is 0 if and only if every check
passed.

Runtime: about one second on one core.
"""

import sys
from fractions import Fraction as F
from math import isqrt

# ===================================================================== #
#  THE OBJECT, AS PRINTED IN THE PAPER                                  #
# ===================================================================== #

PAPER = {
    # Section 2: the two antidifferences, as coefficient lists in ascending
    # powers of x.   psi(x) = A(x)/(12 x^6),  psi3(x) = B(x)/(12 x^8).
    "A": [-1, 0, 3, 6, 6],                       # 6x^4 + 6x^3 + 3x^2 - 1
    "B": [1, 0, -1, 0, 3, 6, 6],                 # 6x^6 + 6x^5 + 3x^4 - x^2 + 1

    # Section 2: the two sign polynomials.
    #   num(12 D)  over j^6 (j+1)^6,  D(j)  = psi(j)  - psi(j+1)  - 1/j^3
    #   num(12 D3) over j^8 (j+1)^8,  D3(j) = psi3(j) - psi3(j+1) - 1/j^3
    "num12D":  [-1, -6, -12, -8],                # = -(2j+1)^3, all <= 0
    "num12D3": [1, 8, 27, 48, 45, 18],           # all > 0 -- the opposite sign

    # Section 2: the tail bound on D3 and its summed form.
    "D3_pointwise_bound": (49, 4, 11),           # D3(j) <= 49/(4 j^11)
    "D3_tail_bound": (49, 40, 10),               # sum_{j>N} D3(j) <= 49/(40 N^10)

    # Section 3: the two decisive certificates, as coefficient lists in t
    # where N = 1 + t.
    "I2_at_1_plus_t":   [659, 1380, 1055, 342, 38],
    "I1_3_at_1_plus_t": [492, 3448, 7664, 8328, 5032, 1728, 316, 24],
    "I2_values_1_to_5":   [659, 3474, 10983, 26606, 54675],
    "I1_3_values_1_to_5": [492, 27032, 263772, 1415016, 5415116],

    # Section 4: the 19 residue classes.  (rho, floor((12rho+6)/19), s0)
    # where the class N = 19 s + rho is certified for every s >= s0.
    "residue_table": [
        (0, 0, 1), (1, 0, 1), (2, 1, 1), (3, 2, 1), (4, 2, 1),
        (5, 3, 0), (6, 4, 0), (7, 4, 0), (8, 5, 0), (9, 6, 0),
        (10, 6, 0), (11, 7, 0), (12, 7, 0), (13, 8, 0), (14, 9, 0),
        (15, 9, 0), (16, 10, 0), (17, 11, 0), (18, 11, 0),
    ],
    "certificate_covers_from": 5,                # every N >= 5
    "uncovered": [1, 2, 3, 4],
    # The delicate class rho = 9 (where 19 | 12N+6), coefficients in s.
    "rho9_J1": [12491451, 142067432, 671815989, 1690280508,
                2385550453, 1789910922, 557513238],
    "rho9_J2_len": 18,
    "rho9_J2_c0": 37139362880053056000,

    # Section 5: the four small-index brackets, (N, J, lo, hi, cf) with
    # lo < R_N < hi certified and cf the determined initial quotients.
    "small": [
        (1, 1, (155, 768),         (207, 1024),            [0, 3, 6, 1]),
        (2, 2, (337, 4374),        (6067, 78732),          [0, 5, 10, 1]),
        (3, 3, (1967, 49152),      (10491, 262144),        [0, 7, 14, 2]),
        (4, 5, (1707247, 69984000), (61461017, 2519424000), [0, 9, 18, 2]),
    ],
    "a3_small_actual":    {1: 1, 2: 1, 3: 2, 4: 2},
    "a3_small_predicted": {1: 0, 2: 1, 3: 2, 4: 2},

    # Section 1 and 5: the Pell data and the target cell.
    "target_cell": {"k": 1, "N": 3, "M": 5, "Q": 7, "cf": [0, 7, 14, 2]},
    "N_k_first_four": [3, 20, 119, 696],

    # Section 6: the controls.
    "sqrt2_cf": [1, 2, 2, 2, 2, 2, 2, 2],
    "thm_main_k": [3, 4, 5],
    "phidef": {2: (1, 2), 3: (1, 2), 4: (1, 4), 6: (-1, 12), 8: (1, 12),
               10: (-3, 20), 12: (5, 12), 14: (-691, 420), 16: (35, 4),
               18: (-3617, 60)},
    "eqNum": [25319, 455742, 3870132, 20594154, 76914556, 213946830,
              458881964, 774825162, 1041773850, 1119835402, 959399980,
              647615118, 336386804, 128519066, 32900250, 4386700],
    "zeta3_decimal": "1.20205690315959428539973816151144999076498629234",
    "zeta3_bracket_N": [1, 2, 3, 5, 20, 119],
    "sweep_max": 60,
    "boundary": {9: 6, 28: 18, 47: 30},

    # Section 7: the superseded four-term enclosure, recorded as a negative
    # measurement.  psi(N+1) < R_N < psi(N+1) + 9/(32 N^8).
    "four_term_I1_values_1_to_4": [-3337600, -219443968, -1361791872, 71258831168],
    "four_term_rho9_deg": 12,
    "four_term_rho9_c0": -1616647295583072,
    "four_term_rho9_clead": -3454978457596019520,
}

# ===================================================================== #
#  CHECK HARNESS                                                        #
# ===================================================================== #

_N_OK = 0
_BAD = []


def check(name, ok, detail=""):
    global _N_OK
    if ok:
        _N_OK += 1
        print("PASS %s %s" % (name, detail))
    else:
        _BAD.append(name)
        print("CHECK-FAILED %s %s" % (name, detail))


def head(t):
    print("")
    print("=== %s ===" % t)


# ===================================================================== #
#  EXACT POLYNOMIAL ARITHMETIC OVER Fraction                            #
#  A polynomial is a list of Fractions, ascending powers.               #
#  A rational function is a pair (numerator, denominator).              #
# ===================================================================== #

def P(cs):
    return [F(c) for c in cs]


def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(*ps):
    n = max(len(p) for p in ps)
    r = [F(0)] * n
    for p in ps:
        for i, c in enumerate(p):
            r[i] += c
    return ptrim(r)


def psc(s, p):
    return ptrim([F(s) * c for c in p])


def pmul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return ptrim(r)


def ppow(p, n):
    r = [F(1)]
    for _ in range(n):
        r = pmul(r, p)
    return r


def pev(p, x):
    s = F(0)
    for c in reversed(p):
        s = s * x + c
    return s


def pcomp(p, q):
    """p(q(x))."""
    r = [F(0)]
    for i, c in enumerate(p):
        if c:
            r = padd(r, psc(c, ppow(q, i)))
    return r


def pdiv_exact(p, q):
    """p/q, asserting an exact division."""
    p = list(p)
    out = [F(0)] * (len(p) - len(q) + 1)
    for i in range(len(p) - len(q), -1, -1):
        c = p[i + len(q) - 1] / q[-1]
        out[i] = c
        for j, y in enumerate(q):
            p[i + j] -= c * y
    assert all(x == 0 for x in p), "inexact polynomial division"
    return ptrim(out)


def allpos(p):
    return all(c >= 0 for c in p) and any(c > 0 for c in p)


def radd(A, B):
    return (padd(pmul(A[0], B[1]), pmul(B[0], A[1])), pmul(A[1], B[1]))


def rsub(A, B):
    return radd(A, (psc(-1, B[0]), B[1]))


def rmul(A, B):
    return (pmul(A[0], B[0]), pmul(A[1], B[1]))


def rinv(A):
    return (A[1], A[0])


def req(A, B):
    """equality of rational functions, by cross-multiplication."""
    return pmul(A[0], B[1]) == pmul(B[0], A[1])


def rc(c):
    return ([F(c)], [F(1)])


X = P([0, 1])
X1 = P([1, 1])


def shift1(p):
    """p(x+1)."""
    return pcomp(p, X1)


def cert_from(poly, lo, hi):
    """least s0 in [lo,hi] with poly(s0 + t) having all coefficients >= 0."""
    for s0 in range(lo, hi + 1):
        sh = pcomp(poly, P([s0, 1]))
        if allpos(sh):
            return s0, sh
    return None, None


# ===================================================================== #
#  EXACT NUMERIC ROUTE (integers only; used as an INDEPENDENT check)    #
# ===================================================================== #

def psi_at(x):
    return F(6 * x ** 4 + 6 * x ** 3 + 3 * x ** 2 - 1, 12 * x ** 6)


def psi3_at(x):
    return F(6 * x ** 6 + 6 * x ** 5 + 3 * x ** 4 - x ** 2 + 1, 12 * x ** 8)


def iroot(n, k):
    if n <= 0:
        return 0
    x = 1 << ((n.bit_length() + k - 1) // k + 1)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            return x
        x = y


def Renc(N, digits):
    """(lo, hi) with lo < R_N < hi and hi - lo < 10^-digits.

    R_N = sum_{j=N+1}^{J} j^-3 + R_J and psi(J+1) < R_J < psi3(J+1).  The
    partial sum is accumulated in fixed point with integer floor division, so
    the running value is a certified under-estimate and adding one unit per
    term gives a certified over-estimate.  No floating point anywhere.
    """
    J = max(N, iroot(10 ** digits // 12, 8) + 2)
    prec = (digits * 10) // 3 + 80 + (J - N).bit_length()
    S = 1 << prec
    acc = 0
    for j in range(N + 1, J + 1):
        acc += S // j ** 3
    return F(acc, S) + psi_at(J + 1), F(acc + (J - N), S) + psi3_at(J + 1), prec


def sqrt_lo(x, prec):
    return F(isqrt((x.numerator << (2 * prec)) // x.denominator), 1 << prec)


def sqrt_hi(x, prec):
    return F(isqrt((x.numerator << (2 * prec)) // x.denominator) + 1, 1 << prec)


def cf_bracket(lo, hi, nmax):
    """The initial regular-continued-fraction digits COMMON to every real in
    [lo, hi]; stops as soon as a digit is not determined by the bracket."""
    out = []
    for _ in range(nmax):
        a = lo.numerator // lo.denominator
        b = hi.numerator // hi.denominator
        if a != b:
            break
        out.append(a)
        fl, fh = lo - a, hi - a
        if fl <= 0 or fh <= 0:
            break
        lo, hi = 1 / fh, 1 / fl
    return out


def cf_of_sqrt_R_over(N, d, digits, nmax):
    lo, hi, prec = Renc(N, digits)
    return cf_bracket(sqrt_lo(lo / d, prec), sqrt_hi(hi / d, prec), nmax)


def cf_from_bracket(lo, hi, nmax, prec=400):
    return cf_bracket(sqrt_lo(lo / 2, prec), sqrt_hi(hi / 2, prec), nmax)


# ===================================================================== #
#  1.  THE ANTIDIFFERENCES psi AND psi3                                 #
# ===================================================================== #

head("1. the two antidifferences and their signs (Section 2 of the paper)")

A = P(PAPER["A"])
B = P(PAPER["B"])

# psi(x) = 1/(2x^2) + 1/(2x^3) + 1/(4x^4) - 1/(12x^6); check A = 12 x^6 psi.
A_derived = padd(psc(6, ppow(X, 4)), psc(6, ppow(X, 3)), psc(3, ppow(X, 2)), P([-1]))
check("A-equals-12x6-psi", A == A_derived,
      "A(x) = 12 x^6 (1/(2x^2) + 1/(2x^3) + 1/(4x^4) - 1/(12x^6)) has "
      "coefficient vector %s" % [int(c) for c in A_derived])

B_derived = padd(pmul(ppow(X, 2), A), P([1]))
check("B-equals-x2A-plus-1", B == B_derived,
      "B(x) = x^2 A(x) + 1 = %s, i.e. psi3 = psi + 1/(12 x^8)"
      % [int(c) for c in B_derived])

PSI = (A, psc(12, ppow(X, 6)))
PSI3 = (B, psc(12, ppow(X, 8)))
PSI_p1 = (shift1(A), psc(12, ppow(X1, 6)))
PSI3_p1 = (shift1(B), psc(12, ppow(X1, 8)))
CUBE = ([F(1)], ppow(X, 3))

D = rsub(rsub(PSI, PSI_p1), CUBE)
D3 = rsub(rsub(PSI3, PSI3_p1), CUBE)

num12D = pdiv_exact(pmul(D[0], psc(12, pmul(ppow(X, 6), ppow(X1, 6)))), D[1])
num12D3 = pdiv_exact(pmul(D3[0], psc(12, pmul(ppow(X, 8), ppow(X1, 8)))), D3[1])

check("D-numerator-as-printed", num12D == P(PAPER["num12D"]),
      "12 j^6 (j+1)^6 (psi(j)-psi(j+1)-1/j^3) = %s" % [int(c) for c in num12D])
check("D-numerator-is-minus-2j-plus-1-cubed",
      num12D == psc(-1, ppow(P([1, 2]), 3)),
      "the vector equals -(2j+1)^3 identically")
check("D-all-coefficients-nonpositive", all(c <= 0 for c in num12D),
      "so D(j) < 0 for every j >= 1, hence R_N > psi(N+1)")

check("D3-numerator-as-printed", num12D3 == P(PAPER["num12D3"]),
      "12 j^8 (j+1)^8 (psi3(j)-psi3(j+1)-1/j^3) = %s" % [int(c) for c in num12D3])
check("D3-all-coefficients-positive", all(c > 0 for c in num12D3),
      "sign OPPOSITE to D, so D3(j) > 0 and R_N < psi3(N+1) with NO tail term")

# The two hand checks quoted in the paper.
check("D-hand-check-j-1", psi_at(1) - psi_at(2) - F(1, 1) == F(-27, 768),
      "psi(1) - psi(2) - 1 = 7/6 - 155/768 - 1 = -27/768 = -(2*1+1)^3/(12*1^6*2^6)")
check("D3-hand-check-j-1", psi3_at(1) - psi3_at(2) - F(1, 1) == F(49, 1024),
      "psi3(1) - psi3(2) - 1 = 49/1024 = 147/(12*1^8*2^8) > 0")

# Pointwise tail bound  D3(j) <= 49/(4 j^11),  as a coefficientwise statement:
#   588 (j+1)^8 - 4 j^3 num12D3(j)  has all coefficients >= 0.
cn, cd, ce = PAPER["D3_pointwise_bound"]
gap = padd(psc(cd * 3 * cn, ppow(X1, 8)), psc(-cd, pmul(ppow(X, 3), num12D3)))
check("D3-pointwise-tail-bound", allpos(gap),
      "%d (j+1)^8 - %d j^3 num(12 D3) has all coefficients >= 0, so "
      "D3(j) <= %d/(%d j^%d) for every j >= 1"
      % (cd * 3 * cn, cd, cn, cd, ce))

# Summed tail bound  sum_{j>N} j^-11 <= 1/(10 N^10), via the telescoping
# inequality  1/(10 j^10) - 1/(10 (j+1)^10) >= 1/(j+1)^11, i.e.
#   (j+1)^11 - j^11 - j^10 - 10 j^10  has all coefficients >= 0.
tele = padd(ppow(X1, 11), psc(-1, ppow(X, 11)), psc(-11, ppow(X, 10)))
check("tail-sum-telescopes", allpos(tele),
      "(j+1)^11 - j^11 - 11 j^10 has all coefficients >= 0, so "
      "1/(10 j^10) - 1/(10 (j+1)^10) >= 1/(j+1)^11 and, summing, "
      "sum_{j>N} j^-11 <= 1/(10 N^10)")
tn, td, te = PAPER["D3_tail_bound"]
check("D3-summed-tail-bound",
      F(cn, cd) * F(1, 10) == F(tn, td) and te == ce - 1,
      "%d/(%d j^%d) summed over j > N is at most %d/%d * 1/(10 N^%d) = "
      "%d/(%d N^%d), so R_N > psi3(N+1) - %d/(%d N^%d) for every N >= 1"
      % (cn, cd, ce, cn, cd, te, tn, td, te, tn, td, te))

# The enclosures, checked numerically against an independent exact route.
lo_ok = hi_ok = tight_ok = True
for N in range(1, 13):
    lo, hi, _ = Renc(N, 24)
    lo_ok = lo_ok and psi_at(N + 1) < lo
    hi_ok = hi_ok and hi < psi3_at(N + 1)
    tight_ok = tight_ok and (psi3_at(N + 1) - F(tn, td * N ** te)) < lo
check("enclosure-lower-holds", lo_ok,
      "psi(N+1) < R_N for N = 1..12, against an independent exact partial-sum route")
check("enclosure-upper-holds", hi_ok, "R_N < psi3(N+1) for N = 1..12")
check("enclosure-tight-lower-holds", tight_ok,
      "psi3(N+1) - 49/(40 N^10) < R_N for N = 1..12")

# ===================================================================== #
#  2.  THE LEVEL-2 CERTIFICATE:  a_1 = 2N+1 AND a_2 = 4N+2              #
# ===================================================================== #

head("2. the certificate for a_1 = 2N+1 and a_2 = 4N+2 (Section 3)")

A1 = (P([1, 2]), [F(1)])          # 2N + 1
A2 = (P([2, 4]), [F(1)])          # 4N + 2
loV = radd(A1, rinv(radd(A2, rc(1))))    # (2N+1) + 1/(4N+3)
hiV = radd(A1, rinv(A2))                 # (2N+1) + 1/(4N+2)

check("loV-closed-form", req(loV, (P([4, 10, 8]), P([3, 4]))),
      "loV = (2N+1) + 1/(4N+3) = (8N^2+10N+4)/(4N+3)")
check("hiV-closed-form", req(hiV, (P([3, 8, 8]), P([2, 4]))),
      "hiV = (2N+1) + 1/(4N+2) = (8N^2+8N+3)/(4N+2)")

I2 = rsub(rmul(rmul(hiV, hiV), PSI_p1), rc(2))     # hiV^2 psi(N+1) - 2
I1_3 = rsub(rc(2), rmul(rmul(loV, loV), PSI3_p1))  # 2 - loV^2 psi3(N+1)

I2_sh = pcomp(I2[0], X1)
I1_3_sh = pcomp(I1_3[0], X1)

check("I2-shifted-vector-as-printed", I2_sh == P(PAPER["I2_at_1_plus_t"]),
      "numerator of I2 at N = 1+t is %s (degree %d)"
      % ([int(c) for c in I2_sh], len(I2_sh) - 1))
check("I2-all-coefficients-positive", allpos(I2_sh),
      "so I2 > 0 for every N >= 1")
check("I2-denominator-positive", allpos(pcomp(I2[1], X1)),
      "denominator 12 (N+1)^6 (4N+2)^2 has all coefficients >= 0 after N -> 1+t")
check("I2-denominator-as-printed",
      I2[1] == psc(12, pmul(ppow(X1, 6), ppow(P([2, 4]), 2))),
      "denominator is exactly 12 (N+1)^6 (4N+2)^2")
check("I2-values-1-to-5",
      [pev(I2[0], F(n)) for n in range(1, 6)] == P(PAPER["I2_values_1_to_5"]),
      "numerator of I2 at N = 1..5 is %s" % PAPER["I2_values_1_to_5"])

check("I1-3-shifted-vector-as-printed", I1_3_sh == P(PAPER["I1_3_at_1_plus_t"]),
      "numerator of I1_3 at N = 1+t is %s (degree %d)"
      % ([int(c) for c in I1_3_sh], len(I1_3_sh) - 1))
check("I1-3-all-coefficients-positive", allpos(I1_3_sh),
      "so I1_3 > 0 for every N >= 1")
check("I1-3-denominator-positive", allpos(pcomp(I1_3[1], X1)),
      "denominator 12 (N+1)^8 (4N+3)^2 has all coefficients >= 0 after N -> 1+t")
check("I1-3-denominator-as-printed",
      I1_3[1] == psc(12, pmul(ppow(X1, 8), ppow(P([3, 4]), 2))),
      "denominator is exactly 12 (N+1)^8 (4N+3)^2")
check("I1-3-values-1-to-5",
      [pev(I1_3[0], F(n)) for n in range(1, 6)] == P(PAPER["I1_3_values_1_to_5"]),
      "numerator of I1_3 at N = 1..5 is %s" % PAPER["I1_3_values_1_to_5"])
check("I1-3-hand-check-N-1", 2 * 12 * 256 * 49 - 484 * 621 == 492,
      "2*12*256*49 - 484*621 = 301056 - 300564 = 492, the N=1 value by hand")
check("I2-numerator-primitive",
      all(659 % p for p in range(2, 26)) and __import__("math").gcd(659, 38) == 1,
      "659 is prime and gcd(659, 38) = 1, so the I2 vector is primitive")

# The two inequalities that turn positivity into partial quotients.
check("loV-above-2N-plus-1",
      allpos(pcomp(psc(-1, padd(pmul(P([1, 2]), loV[1]), psc(-1, loV[0]))), X1)),
      "loV - (2N+1) = 1/(4N+3) > 0")
check("hiV-below-2N-plus-2",
      allpos(pcomp(padd(pmul(P([2, 2]), hiV[1]), psc(-1, hiV[0])), X1)),
      "(2N+2) - hiV = 1 - 1/(4N+2) > 0, so floor(V) = 2N+1")
check("loV-below-hiV", allpos(pcomp(rsub(hiV, loV)[0], X1))
      and allpos(pcomp(rsub(hiV, loV)[1], X1)),
      "hiV - loV = 1/(4N+2) - 1/(4N+3) > 0")
check("second-digit-interval-is-4N-plus-2",
      req(rinv(rsub(hiV, A1)), A2) and req(rinv(rsub(loV, A1)), radd(A2, rc(1))),
      "1/(hiV-(2N+1)) = 4N+2 and 1/(loV-(2N+1)) = 4N+3 identically, so for any "
      "V in (loV, hiV) the reciprocal of its fractional part lies strictly in "
      "(4N+2, 4N+3) and a_2 = 4N+2")

# ===================================================================== #
#  3.  THE LEVEL-3 CERTIFICATE:  19 RESIDUE CLASSES                     #
# ===================================================================== #

head("3. the 19 residue-class certificates for a_3 (Section 4)")

LO3 = rsub(PSI3_p1, ([F(tn)], psc(td, ppow(X, te))))   # psi3(N+1) - 49/(40 N^10)
HI3 = PSI3_p1                                          # psi3(N+1)

rows = []
rho9 = None
for rho, base_p, s0_p in PAPER["residue_table"]:
    base = (12 * rho + 6) // 19
    sub = P([rho, 19])                                  # N = 19 s + rho
    a1r = (padd(psc(2, sub), [F(1)]), [F(1)])
    a2r = (padd(psc(4, sub), [F(2)]), [F(1)])
    a3r = (P([base, 12]), [F(1)])                       # A_3 = 12 s + base
    a3p = (P([base + 1, 12]), [F(1)])
    loV3 = radd(a1r, rinv(radd(a2r, rinv(a3r))))
    hiV3 = radd(a1r, rinv(radd(a2r, rinv(a3p))))
    HIs = (pcomp(HI3[0], sub), pcomp(HI3[1], sub))
    LOs = (pcomp(LO3[0], sub), pcomp(LO3[1], sub))
    J1 = rsub(rc(2), rmul(rmul(loV3, loV3), HIs))       # > 0  =>  V > loV3
    J2 = rsub(rmul(rmul(hiV3, hiV3), LOs), rc(2))       # > 0  =>  V < hiV3
    s1, _ = cert_from(J1[0], 0, 40)
    s2, _ = cert_from(J2[0], 0, 40)
    d1, _ = cert_from(J1[1], 0, 40)
    d2, _ = cert_from(J2[1], 0, 40)
    got = None if (s1 is None or s2 is None) else max(s1, s2, d1 or 0, d2 or 0)
    rows.append((rho, base, got))
    ok = (base == base_p) and (got == s0_p)
    check("class-rho-%02d" % rho, ok,
          "N = 19s+%d, A_3 = 12s+%d: numerator and denominator of both J_1 and "
          "J_2 have all coefficients >= 0 after s -> %s + t (paper says %s)"
          % (rho, base, got, s0_p))
    if rho == 9:
        rho9 = (J1, J2)

check("all-19-classes-certified",
      all(g is not None for _, _, g in rows),
      "19 of 19 residue classes carry a positive-coefficient certificate")

J1_9, J2_9 = rho9
check("rho-9-is-the-integer-class",
      all((12 * (19 * s + 9) + 6) % 19 == 0 for s in range(6))
      and [r for r in range(19) if (12 * r + 6) % 19 == 0] == [9],
      "19 | 12N+6 exactly when N = 9 mod 19, so rho = 9 is the unique class in "
      "which the floor sits ON an integer")
check("rho-9-J1-vector-as-printed",
      J1_9[0] == P(PAPER["rho9_J1"]),
      "J_1 numerator in class rho = 9 is %s (degree %d, all >= 0)"
      % ([int(c) for c in J1_9[0]], len(J1_9[0]) - 1))
check("rho-9-J2-shape-as-printed",
      len(J2_9[0]) == PAPER["rho9_J2_len"] and int(J2_9[0][0]) == PAPER["rho9_J2_c0"]
      and allpos(J2_9[0]),
      "J_2 numerator in class rho = 9 has %d coefficients, all positive, the "
      "first being %d" % (len(J2_9[0]), int(J2_9[0][0])))
check("rho-1-is-the-18-over-19-class",
      [r for r in range(19) if (12 * r + 6) % 19 == 18] == [1],
      "the fractional part of (12N+6)/19 equals 18/19 exactly when N = 1 mod 19, "
      "which is why N = 1 and not N = 20 is the exceptional index")

# What the 19 certificates cover.
covered = set()
for rho, base, s0 in rows:
    for s in range(s0, s0 + 600):
        covered.add(19 * s + rho)
uncovered = [n for n in range(1, 1 + 19 * 500) if n not in covered]
check("certificate-covers-every-N-at-least-5",
      all(n >= PAPER["certificate_covers_from"] or n in PAPER["uncovered"]
          for n in range(1, 200))
      and all(n in covered for n in range(PAPER["certificate_covers_from"], 200)),
      "every N with 5 <= N <= 199 lies in a certified class")
check("uncovered-set-is-exactly-1-2-3-4",
      uncovered == PAPER["uncovered"],
      "the only indices no class certifies are N in %s" % PAPER["uncovered"])

# ===================================================================== #
#  4.  THE FOUR UNCOVERED INDICES, DECIDED EXACTLY                      #
# ===================================================================== #

head("4. the four uncovered indices, decided from printed rational brackets "
     "(Section 5)")

for N, J, lo_t, hi_t, cf_t in PAPER["small"]:
    part = sum(F(1, j ** 3) for j in range(N + 1, J + 1))
    lo = part + psi_at(J + 1)
    hi = part + psi3_at(J + 1)
    check("bracket-N-%d" % N,
          (lo.numerator, lo.denominator) == lo_t
          and (hi.numerator, hi.denominator) == hi_t,
          "R_%d in (%d/%d, %d/%d), from the exact partial sum to j = %d and "
          "psi(%d) < R_%d < psi3(%d)"
          % (N, lo_t[0], lo_t[1], hi_t[0], hi_t[1], J, J + 1, J, J + 1))
    got = cf_from_bracket(lo, hi, 4)
    check("cf-N-%d" % N, got[:4] == cf_t,
          "that bracket determines g(%d) = [0; %s]"
          % (N, ", ".join(str(x) for x in cf_t[1:])))
    a3 = got[3]
    check("a3-N-%d" % N, a3 == PAPER["a3_small_actual"][N],
          "a_3(g(%d)) = %d, against the displayed prediction "
          "floor((12*%d+6)/19) = %d -- %s"
          % (N, a3, N, (12 * N + 6) // 19,
             "AGREE" if a3 == (12 * N + 6) // 19 else "DISAGREE"))

check("N-1-is-the-only-small-disagreement",
      [N for N in (1, 2, 3, 4)
       if PAPER["a3_small_actual"][N] != (12 * N + 6) // 19] == [1],
      "of the four uncovered indices only N = 1 disagrees with the display: "
      "a_3(g(1)) = 1 while floor(18/19) = 0")
check("floor-18-over-19-is-zero", (12 * 1 + 6) // 19 == 0,
      "floor((12*1+6)/19) = floor(18/19) = 0, which is not a legal partial "
      "quotient at a position >= 1")

# ===================================================================== #
#  5.  THE TARGET CELL AND THE PELL SPECIALISATION                      #
# ===================================================================== #

head("5. the target cell k = 1 and eq:gcf-pell (Sections 1 and 5)")

Pell = [0, 1]
while len(Pell) < 24:
    Pell.append(2 * Pell[-1] + Pell[-2])
Qn = [1] + [Pell[i] + Pell[i - 1] for i in range(1, 24)]

check("pell-recursion",
      Pell[:8] == [0, 1, 2, 5, 12, 29, 70, 169] and Qn[:8] == [1, 1, 3, 7, 17, 41, 99, 239],
      "P = 0,1,2,5,12,29,70,169,... and Q_k = P_k + P_{k-1} = 1,1,3,7,17,41,99,239,...")
check("pell-negative-unit",
      all(Qn[2 * k + 1] ** 2 - 2 * Pell[2 * k + 1] ** 2 == -1 for k in range(1, 10)),
      "Q_{2k+1}^2 - 2 M_k^2 = -1 with M_k = P_{2k+1}, for k = 1..9")
Nk = [(Qn[2 * k + 1] - 1) // 2 for k in range(1, 5)]
check("N-k-values", Nk == PAPER["N_k_first_four"],
      "N_k = (Q_{2k+1}-1)/2 = %s" % Nk)

tc = PAPER["target_cell"]
check("target-cell-parameters",
      tc["N"] == (Qn[3] - 1) // 2 and tc["M"] == Pell[3] and tc["Q"] == Qn[3],
      "k = 1 gives M_1 = P_3 = %d, Q_3 = %d, N_1 = %d -- the index the source "
      "names as the site of its gap" % (Pell[3], Qn[3], (Qn[3] - 1) // 2))
N, J, lo_t, hi_t, _ = PAPER["small"][2]
part = sum(F(1, j ** 3) for j in range(N + 1, J + 1))
cf3 = cf_from_bracket(part + psi_at(J + 1), part + psi3_at(J + 1), 4)
check("target-cell-continued-fraction", cf3 == tc["cf"],
      "g(3) = sqrt((zeta(3) - 251/216)/2) = [0; 7, 14, 2], matching the three "
      "quotients the source displays at N = 3")
check("eq-gcf-pell-at-every-k",
      all(2 * ((Qn[2 * k + 1] - 1) // 2) + 1 == Qn[2 * k + 1]
          and 4 * ((Qn[2 * k + 1] - 1) // 2) + 2 == 2 * Qn[2 * k + 1]
          for k in range(1, 10)),
      "2N_k+1 = Q_{2k+1} and 4N_k+2 = 2 Q_{2k+1} identically, so the certificate "
      "of Section 3 gives eq:gcf-pell for every k >= 1")
check("S-3-exact", sum(F(1, j ** 3) for j in range(1, 4)) == F(251, 216),
      "S_3 = 1 + 1/8 + 1/27 = 251/216")

# ===================================================================== #
#  6.  CONTROLS                                                         #
# ===================================================================== #

head("6. controls -- forced positives, forced negatives, and the source's own "
     "printed integers (Section 6)")

check("control-decider-sqrt-2",
      cf_bracket(sqrt_lo(F(2), 300), sqrt_hi(F(2), 300), 8) == PAPER["sqrt2_cf"],
      "the bracket-continued-fraction routine returns %s for sqrt(2)"
      % PAPER["sqrt2_cf"])

for k in PAPER["thm_main_k"]:
    M = Pell[2 * k + 1]
    r = (10 * M) % 261
    pred = [0, M - 1, 1, 6 * M ** 3 + 12 * M - 2, 1, (10 * M - 261) // 261,
            261 // r + (-1 if 261 % r == 0 else 0)]
    Nk_ = (Qn[2 * k + 1] - 1) // 2
    got = []
    for dig in (24, 32, 40, 48, 56):
        got = cf_of_sqrt_R_over(Nk_, 1, dig, 7)
        if len(got) >= 7:
            break
    check("control-published-thm-main-k-%d" % k, got[:len(pred)] == pred,
          "our routine returns the source's own six printed quotients for "
          "sqrt(R_{N_%d}) at N = %d, M = %d: %s" % (k, Nk_, M, pred))

neg = []
for d in (3, 4):
    for Nn in (3, 20, 119):
        c = cf_of_sqrt_R_over(Nn, d, 24, 3)
        neg.append((d, Nn, c[1:3], [2 * Nn + 1, 4 * Nn + 2]))
check("control-negative-polarity", all(a != b for _, _, a, b in neg),
      "the same pipeline on sqrt(R_N/3) and sqrt(R_N/4) reproduces the /2 "
      "prediction in 0 of 6 cells: %s"
      % "; ".join("d=%d N=%d gives %s not %s" % t for t in neg))

# The source's own antidifference and defect numerator, reproduced.
BER = {2: F(1, 6), 4: F(-1, 30), 6: F(1, 42), 8: F(-1, 30), 10: F(5, 66),
       12: F(-691, 2730), 14: F(7, 6), 16: F(-3617, 510), 18: F(43867, 798)}
phi = {2: F(1, 2), 3: F(1, 2)}
for k in range(1, 9):
    phi[2 * k + 2] = phi.get(2 * k + 2, F(0)) + F(2 * k + 1, 2) * BER[2 * k]
check("control-source-phidef",
      phi == {p: F(a, b) for p, (a, b) in PAPER["phidef"].items()},
      "the Bernoulli-derived antidifference reproduces all 10 coefficients the "
      "source prints in its eq:phidef")
Pdeg = max(max(phi), 3)
Num = [F(0)]
for p, cp in phi.items():
    Num = padd(Num, psc(cp, padd(pmul(ppow(X1, Pdeg), ppow(X, Pdeg - p)),
                                 psc(-1, pmul(ppow(X, Pdeg), ppow(X1, Pdeg - p))))))
Num = padd(Num, psc(-1, pmul(ppow(X1, Pdeg), ppow(X, Pdeg - 3))))
check("control-source-eqNum",
      [-420 * c for c in Num] == P(PAPER["eqNum"]),
      "and all 16 integers of its eq:Num, under the source's own normalisation "
      "D = -Num/(420 j^18 (j+1)^18); degree 15, 2P-d = %d matching its stated "
      "j^-21 majorant" % (2 * Pdeg - (len(Num) - 1)))
check("control-source-defect-one-sign",
      len(set(1 if c > 0 else -1 for c in Num if c)) == 1,
      "that defect has one sign throughout, as the source asserts")

Z = PAPER["zeta3_decimal"]
z3lo = F(int(Z.replace(".", "")), 10 ** (len(Z) - 2))
z3hi = z3lo + F(1, 10 ** (len(Z) - 2))
zok = True
for Nn in PAPER["zeta3_bracket_N"]:
    SN = sum(F(1, j ** 3) for j in range(1, Nn + 1))
    lo, hi, _ = Renc(Nn, 24)
    zok = zok and lo < z3lo - SN and z3hi - SN < hi
check("control-zeta3-decimal-bracket", zok,
      "our enclosure brackets the published decimal of zeta(3) minus the exact "
      "rational S_N at N = %s, 6 of 6 -- the only place a decimal of zeta(3) "
      "appears anywhere" % PAPER["zeta3_bracket_N"])

f1 = f2 = 0
f3 = []
for Nn in range(1, PAPER["sweep_max"] + 1):
    c = cf_of_sqrt_R_over(Nn, 2, 20, 5)
    if len(c) < 4:
        f3.append((Nn, "undetermined"))
        continue
    if c[1] != 2 * Nn + 1:
        f1 += 1
    if c[2] != 4 * Nn + 2:
        f2 += 1
    if c[3] != (12 * Nn + 6) // 19:
        f3.append((Nn, c[3], (12 * Nn + 6) // 19))
check("sweep-a1-a2-agree", f1 == 0 and f2 == 0,
      "an independent exact sweep over N = 1..%d finds 0 disagreements for a_1 "
      "and 0 for a_2" % PAPER["sweep_max"])
check("sweep-a3-single-failure", f3 == [(1, 1, 0)],
      "and exactly one for a_3, at N = 1, where the true value is 1 and the "
      "display predicts 0")
bok = all(cf_of_sqrt_R_over(Nn, 2, 24, 5)[3] == v
          for Nn, v in PAPER["boundary"].items())
check("sweep-boundary-indices", bok,
      "the three class boundaries N = 9, 28, 47 give a_3 = 6, 18, 30 as the "
      "quasi-polynomial predicts")

# ===================================================================== #
#  7.  THE SUPERSEDED FOUR-TERM ENCLOSURE (a NEGATIVE measurement)      #
# ===================================================================== #

head("7. why the shorter enclosure does not do it (Section 7)")

HI4 = radd(PSI_p1, ([F(9)], psc(32, ppow(X, 8))))     # psi(N+1) + 9/(32 N^8)
I1_4 = rsub(rc(2), rmul(rmul(loV, loV), HI4))
vals = [pev(I1_4[0], F(n)) for n in range(1, 5)]
check("four-term-I1-negative-at-1-2-3",
      vals == P(PAPER["four_term_I1_values_1_to_4"]) and all(v < 0 for v in vals[:3]),
      "the four-term upper bound gives an I_1 numerator %s -- negative at "
      "exactly the three smallest indices, so it cannot pin a_1 there"
      % [int(v) for v in vals])
sub9 = P([9, 19])
a1r = (padd(psc(2, sub9), [F(1)]), [F(1)])
a2r = (padd(psc(4, sub9), [F(2)]), [F(1)])
a3r = (P([6, 12]), [F(1)])
loV3_9 = radd(a1r, rinv(radd(a2r, rinv(a3r))))
HI4s = (pcomp(HI4[0], sub9), pcomp(HI4[1], sub9))
J1_4 = rsub(rc(2), rmul(rmul(loV3_9, loV3_9), HI4s))
check("four-term-rho-9-all-negative",
      len(J1_4[0]) - 1 == PAPER["four_term_rho9_deg"]
      and all(c < 0 for c in J1_4[0])
      and int(J1_4[0][0]) == PAPER["four_term_rho9_c0"]
      and int(J1_4[0][-1]) == PAPER["four_term_rho9_clead"],
      "and in the class N = 9 mod 19 its J_1 numerator is degree %d with EVERY "
      "coefficient negative (constant %d, leading %d), so that enclosure never "
      "pins a_3 there for any s >= 0"
      % (PAPER["four_term_rho9_deg"], PAPER["four_term_rho9_c0"],
         PAPER["four_term_rho9_clead"]))
check("five-term-succeeds-where-four-term-fails",
      allpos(J1_9[0]) and all(c < 0 for c in J1_4[0]),
      "the same class, same target value, same reduction: five-term all "
      "non-negative, four-term all negative")

# ===================================================================== #
#  SCOPE AND VERDICT                                                    #
# ===================================================================== #

print("")
print("=== scope: what this program does NOT cover ===")
print("NOT RE-RUN: the trailing ellipsis of the source's display. Nothing here")
print("  addresses a_4, a_5 or their quasi-periods, so the source's Conjecture")
print("  prop:gcf as a whole is NOT verified and NOT claimed; only the three")
print("  quotients its equation eq:gcf displays are settled.")
print("NOT RE-RUN: any Lean formalisation. No Lean was written and none was")
print("  inspected; the source's development is not public, so its single")
print("  `sorry` could not be diffed against any artifact. What is discharged is")
print("  the mathematical statement the source describes, not a machine proof.")
print("NOT RE-RUN: the search that found the antidifferences psi and psi3.")
print("  Nothing depends on it -- a positive-coefficient certificate is")
print("  self-checking once the polynomial is exhibited.")
print("NOT COVERED: irrationality, transcendence or normality of zeta(3), and")
print("  the growth of the partial quotients. A decimal of zeta(3) enters this")
print("  program in exactly one place, the bracket control above, and no proved")
print("  statement rests on it.")
print("NOT COVERED: the source's Theorem thm:main. It concerns sqrt(R_{N_k}),")
print("  a different real number, and is used here only as a forced-positive")
print("  control on the continued-fraction routine; it is not reproved.")
print("NOT COVERED: minimality of the certificate. No claim is made that the")
print("  five-term enclosure is the shortest that works, only that the")
print("  four-term one measurably does not.")
print("")

if _BAD:
    print("VERDICT: %d CHECK(S) DID NOT PASS: %s" % (len(_BAD), ", ".join(_BAD)))
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % _N_OK)
sys.exit(0)
