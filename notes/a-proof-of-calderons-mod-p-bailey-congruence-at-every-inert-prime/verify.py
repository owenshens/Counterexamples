#!/usr/bin/env python3
"""
verify.py -- an independent re-derivation of every computational quantity printed in

    "A Proof of Calderon's mod-p Bailey Congruence at Every Odd Inert Prime"
    (paper.tex / paper.pdf, this folder)

Python 3.9 or later, STANDARD LIBRARY ONLY (math.comb, fractions.Fraction).  Exact integer
and rational arithmetic throughout: no floating-point value takes any decision anywhere in
this file, and no external data file is read.  Run it with

    python3 verify.py

and it prints one `PASS <name> <detail>` line per check, a block of `NOT RE-RUN` disclosures,
and a closing verdict.  It exits 0 if and only if every check passed.

WHAT IT READS.  Only objects printed in the paper: the three worked instances of Table 1,
the twelve (T,N,p) families of Table 2, the seven of Table 3, the two mutations of Table 4,
the five rows of Table 5, and the index ranges stated in Section 4.  Those literals are
collected in PART 0 below, each one tagged with the table it is read from, so that a referee
can check the transcription by eye.

TWO STRUCTURALLY DISJOINT DECIDERS, and they must agree wherever both are legal:

  engine A  strips the exact power of p out of every factor (u + v*omega), using
            v_p(u + v*omega) = min(v_p u, v_p v) -- which is Fact (i) of the paper and is
            legal ONLY at an inert p -- and multiplies the remaining unit parts in
            F_{p^2} = F_p[x]/(x^2 - T x + N).  Sub-products over rows of a box are
            memoised; that is caching and nothing else.

  engine B  never reduces and never strips.  It builds both sides as exact ratios of
            elements of Z[omega], forms the difference as (g0 + g1*omega)/h with g0, g1, h
            rational integers, and tests v_p(g0) - v_p(h) >= 1 AND v_p(g1) - v_p(h) >= 1.
            Because {1, omega} is a Z-basis of O = Z[omega], O tensor Z_p = Z_p + Z_p*omega
            for EVERY prime p, so that test IS membership in p*O_{omega,p} -- engine B is
            therefore valid at split and ramified primes too, where engine A is not, and it
            is the only engine used in Table 5.
"""

import sys
import time
from fractions import Fraction
from math import comb

T0 = time.time()

# ===========================================================================
# PART 0 -- THE OBJECT, READ FROM THE PAPER
# ===========================================================================
# Each row is (T, N, p, instances_claimed, order_as_named_in_the_paper).
# omega is a root of x^2 - T*x + N, O = Z[omega], Delta = T^2 - 4N.

CENSUS_C = [                                                    # paper, Table 2
    (0,  1,  7,   3969, "Z[i]"),
    (0,  2,  7,   3969, "Z[sqrt(-2)]"),
    (2,  2,  7,   3969, "Z[i], skewed basis (1,1+i)"),
    (0,  4,  7,   3969, "Z[2i], index 2, non-maximal"),
    (1,  3,  7,   3969, "O_(-11), maximal"),
    (1,  4,  7,   3969, "O_(-15), maximal"),
    (1,  6,  7,   3969, "O_(-23), maximal"),
    (0,  8,  7,   3969, "Z[2 sqrt(-2)], non-maximal"),
    (0,  1, 11,  27225, "Z[i]"),
    (-1, 1, 11,  27225, "Z[zeta_3], Eisenstein"),
    (0,  2, 13,  54756, "Z[sqrt(-2)]"),
    (0,  1, 19, 263169, "Z[i]"),
]
CENSUS_C_TOTAL = 404127                                         # paper, Table 2, last row
CENSUS_C_GAUSSIAN = 294363                                      # paper, Section 5: (T,N)=(0,1)
CENSUS_C_OUTSIDE_THM13 = 109764                                 # paper, Section 5

CENSUS_D = [                                                    # paper, Table 3
    (0,  1,  3,    81, "Z[i], p = 3 <= 5"),
    (0,  4,  3,    81, "Z[2i], p = 3 <= 5"),
    (0,  2,  5,   900, "Z[sqrt(-2)], p = 5 <= 5"),
    (0,  3,  5,   900, "Z[sqrt(-3)], index 2 in Z[zeta_3], p = 5 <= 5"),
    (0, -2,  5,   900, "Z[sqrt(2)], Delta = 8 > 0"),
    (1, -1,  7,  3969, "Z[(1+sqrt 5)/2], Delta = 5 > 0"),
    (0, -2, 11, 27225, "Z[sqrt(2)], Delta = 8 > 0"),
]
CENSUS_D_TOTAL = 34056                                          # paper, Table 3, last row
CENSUS_D_REAL = 32094                                           # paper, Table 3: Delta > 0 rows
CENSUS_D_SMALLP = 2862                                          # paper, Table 3: p <= 5 rows

# Table 5: p NOT inert.  (T, N, p, kind, alpha_cap, instances, violations,
#                          lhs_not_p_integral, first_violating_tuple)
CENSUS_E = [                                                    # paper, Table 5
    (0, 1, 13, "split",     4,  900,  799, 161,  (1, 2, 1, 1, 1, 1, 1, 1)),
    (0, 2, 11, "split",     4,  900,  812, 242,  (1, 2, 1, 1, 1, 1, 1, 1)),
    (0, 1, 13, "split",     6, 3969, 3762, 1142, (1, 2, 1, 1, 1, 1, 1, 1)),
    (0, 2, 11, "split",     6, 3969, 3716, 1473, (1, 2, 1, 1, 1, 1, 1, 1)),
    (0, 1,  2, "ramified",  1,    9,    0, 2,    None),
]

# Table 1: worked instances, all in the order Z[sqrt(-2)] at p = 7.
WORKED = [                                                      # paper, Table 1
    (2, 2, 1, 1, 5, 4, 2, 3),
    (1, 2, 1, 2, 6, 6, 1, 1),
    (2, 1, 1, 1, 3, 5, 1, 2),
]
WORKED_FAMILY = (0, 2, 7)                                       # paper, Section 3
# The first worked instance, spelled out in Section 3 of the paper:
W1_PRIMED = (19, 18, 9, 10)                                     # A', B', C', D'
W1_R_ABCD = (2, 0)                                              # R_omega(2,2;1,1) = 2 exactly
W1_R_DIGITS_RESIDUE = (4, 0)                                    # R_omega(5,4;2,3) = 4 in F_49
W1_BINOMIAL = 40                                                # C(5,2)*C(4,3) = 10*4
W1_BINOMIAL_MOD = 5                                             # 40 = 5 mod 7
W1_LHS_RESIDUE = (5, 0)                                         # both sides, in basis {1,omega}

# Section 4 of the paper: the shape of every census cell.
ABCD = [(A, B, C, D) for A in (1, 2) for B in (1, 2) for C in (1, 2) for D in (1, 2)
        if A >= C and B >= D]                                   # 9 tuples, in this order
ABCD_COUNT = 9                                                  # paper, Section 4

# paper, Section 4: the target-cell diagnostics (Z[sqrt(-2)] at p = 7):
TARGET_DEGENERATE = 144         # instances with A=C, B=D, alpha=gamma, beta=delta
TARGET_VP_POSITIVE = 0          # instances with v_p(LHS) > 0
TARGET_VAL_MISMATCH = 0         # instances with v_p(LHS) != v_p(RHS)

# paper, Table 4: the forced-positive controls on the target cell:
MUTANTS = [                                                     # paper, Table 4
    ("drop-binom-alpha-gamma", 2205, (1, 1, 1, 1, 2, 1, 1, 1)),
    ("swap-the-two-exponents", 1344, (1, 2, 1, 2, 1, 2, 1, 1)),
]
MUTANT_CELL_TOTAL = 3969

# paper, Section 4, second bullet of "Scope, stated as limits": the regime
# v_p(R_omega(A,B;C,D)) > 0, exhausted at p = 7 over 1 <= A,B <= AB_CAP and
# 1 <= C <= min(A,CD_CAP), 1 <= D <= min(B,CD_CAP), with the digit tuple of the
# first worked instance held fixed.
DEGEN_FAMILY = (0, 2, 7)
DEGEN_AB_CAP = 16
DEGEN_CD_CAP = 3
DEGEN_DIGITS = (5, 4, 2, 3)                                     # (alpha, beta, gamma, delta)
DEGEN_INSTANCES = 2025                                          # paper, Section 4, bullet 2
DEGEN_VP_POSITIVE = 144                                         # paper, Section 4, bullet 2

PRIMES_IN_PLAY = (2, 3, 5, 7, 11, 13, 19)

# ===========================================================================
# PART 1 -- the check harness
# ===========================================================================

_npass = 0
_nfail = 0


def say(text):
    print(text)
    sys.stdout.flush()


def PASS(name, detail=""):
    global _npass
    _npass += 1
    say("PASS %s%s" % (name, (" " + detail) if detail else ""))


def BAD(name, detail):
    global _nfail
    _nfail += 1
    say("FAIL %s %s" % (name, detail))


def chk(name, got, want, detail=""):
    if got == want:
        PASS(name, detail if detail else "= %r" % (want,))
    else:
        BAD(name, "expected %r, got %r%s" % (want, got, (" -- " + detail) if detail else ""))


def chk_true(name, cond, detail=""):
    if cond:
        PASS(name, detail)
    else:
        BAD(name, "condition is false -- %s" % detail)


# ===========================================================================
# PART 2 -- exact arithmetic in Z[omega], omega^2 = T*omega - N
# ===========================================================================

def vp(n, p):
    """v_p of a rational integer; +infinity for 0, returned as None-free large int."""
    if n == 0:
        return None                     # treated as +infinity by the callers below
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def vp_inf(n, p):
    v = vp(n, p)
    return float("inf") if v is None else v


def mulZ(a, b, T, N):
    u, v = a
    x, y = b
    return (u * x - N * v * y, u * y + v * x + T * v * y)


def conjZ(a, T):
    u, v = a
    return (u + T * v, -v)


def normZ(a, T, N):
    u, v = a
    return u * u + T * u * v + N * v * v


def boxZ(u0, u1, v0, v1, T, N):
    acc = (1, 0)
    for u in range(u0, u1 + 1):
        for v in range(v0, v1 + 1):
            acc = mulZ(acc, (u, v), T, N)
    return acc


def RZ(A, B, C, D, T, N):
    """R_omega(A,B;C,D) as an exact pair (num, den) of elements of Z[omega]."""
    return boxZ(A - C + 1, A, B - D + 1, B, T, N), boxZ(1, C, 1, D, T, N)


def coords(num, den, T, N):
    """num/den as an exact pair of Fractions in the basis {1, omega}."""
    g = mulZ(num, conjZ(den, T), T, N)
    m = normZ(den, T, N)
    return (Fraction(g[0], m), Fraction(g[1], m))


def in_pO(num, den, p, T, N):
    """is num/den in p * O_{omega,p}?  Valid at EVERY p, since {1,omega} is a Z-basis of O."""
    g = mulZ(num, conjZ(den, T), T, N)
    m = normZ(den, T, N)
    vm = vp_inf(m, p)
    return vp_inf(g[0], p) - vm >= 1 and vp_inf(g[1], p) - vm >= 1


def in_O(num, den, p, T, N):
    """is num/den p-integral, i.e. in O_{omega,p} = Z_p + Z_p*omega?"""
    g = mulZ(num, conjZ(den, T), T, N)
    m = normZ(den, T, N)
    vm = vp_inf(m, p)
    return vp_inf(g[0], p) - vm >= 0 and vp_inf(g[1], p) - vm >= 0


def check_B(A, B, C, D, al, be, ga, de, T, N, p, mutant=None):
    """engine B: exact, no reduction, no stripping.  -> (holds, lhs_is_p_integral)"""
    NL, DL = RZ(p * A + al, p * B + be, p * C + ga, p * D + de, T, N)
    n1, d1 = RZ(A, B, C, D, T, N)
    n2, d2 = RZ(al, be, ga, de, T, N)
    if mutant == "drop-binom-alpha-gamma":
        k = comb(be, de) ** C
    elif mutant == "swap-the-two-exponents":
        k = comb(al, ga) ** C * comb(be, de) ** D
    else:
        k = comb(al, ga) ** D * comb(be, de) ** C
    NR = mulZ(mulZ(n1, n2, T, N), (k, 0), T, N)
    DR = mulZ(d1, d2, T, N)
    L = mulZ(NL, DR, T, N)
    R = mulZ(NR, DL, T, N)
    G = (L[0] - R[0], L[1] - R[1])
    H = mulZ(DL, DR, T, N)
    return in_pO(G, H, p, T, N), in_O(NL, DL, p, T, N)


# ===========================================================================
# PART 3 -- engine A: F_{p^2} with exact p-stripping (inert p only)
# ===========================================================================

class EngineA(object):
    """Requires p odd and x^2 - T*x + N irreducible mod p, i.e. p inert in Z[omega]."""

    ZERO = (0, 0)
    ONE = (1, 0)

    def __init__(self, T, N, p):
        self.T, self.N, self.p = T, N, p
        self._row = {}
        self._box = {}
        self._R = {}
        self._vp = {}

    def mul(self, a, b):
        T, N, p = self.T, self.N, self.p
        u, v = a
        x, y = b
        return ((u * x - N * v * y) % p, (u * y + v * x + T * v * y) % p)

    def powm(self, a, e):
        r = self.ONE
        b = a
        while e:
            if e & 1:
                r = self.mul(r, b)
            b = self.mul(b, b)
            e >>= 1
        return r

    def conj(self, a):
        u, v = a
        return ((u + self.T * v) % self.p, (-v) % self.p)

    def norm(self, a):
        u, v = a
        return (u * u + self.T * u * v + self.N * v * v) % self.p

    def inv(self, a):
        n = self.norm(a)
        if n == 0:
            raise ZeroDivisionError("engine A asked to invert a non-unit %r" % (a,))
        ni = pow(n, self.p - 2, self.p)
        c = self.conj(a)
        return ((c[0] * ni) % self.p, (c[1] * ni) % self.p)

    def vpi(self, n):
        r = self._vp.get(n)
        if r is None:
            r = vp(n, self.p)
            self._vp[n] = r
        return r

    def row(self, u, v0, v1):
        """(e, unit) with prod_{v=v0..v1} (u + v*omega) = p^e * unit."""
        key = (u, v0, v1)
        r = self._row.get(key)
        if r is not None:
            return r
        p = self.p
        eu = self.vpi(u)
        e = 0
        acc = self.ONE
        for v in range(v0, v1 + 1):
            ev = self.vpi(v)
            m = eu if eu < ev else ev
            if m:
                e += m
                q = p ** m
                acc = self.mul(acc, ((u // q) % p, (v // q) % p))
            else:
                acc = self.mul(acc, (u % p, v % p))
        r = (e, acc)
        self._row[key] = r
        return r

    def box(self, u0, u1, v0, v1, memo=False):
        if memo:
            key = (u0, u1, v0, v1)
            r = self._box.get(key)
            if r is not None:
                return r
        e = 0
        acc = self.ONE
        for u in range(u0, u1 + 1):
            ee, aa = self.row(u, v0, v1)
            e += ee
            acc = self.mul(acc, aa)
        r = (e, acc)
        if memo:
            self._box[(u0, u1, v0, v1)] = r
        return r

    def R(self, A, B, C, D, memo=False):
        """(v_p, unit residue) of R_omega(A,B;C,D); the residue is meaningless if v_p > 0."""
        if memo:
            key = (A, B, C, D)
            r = self._R.get(key)
            if r is not None:
                return r
        e1, a1 = self.box(A - C + 1, A, B - D + 1, B)
        e2, a2 = self.box(1, C, 1, D, memo=True)
        r = (e1 - e2, self.mul(a1, self.inv(a2)))
        if memo:
            self._R[(A, B, C, D)] = r
        return r

    def resid(self, v, u):
        if v < 0:
            return ("not-p-integral", v)
        return self.ZERO if v > 0 else u

    def check(self, A, B, C, D, al, be, ga, de, mutant=None):
        """-> (holds, v_p(LHS), v_p(RHS))"""
        p = self.p
        vL, uL = self.R(p * A + al, p * B + be, p * C + ga, p * D + de)
        v2, u2 = self.R(A, B, C, D, memo=True)
        v3, u3 = self.R(al, be, ga, de, memo=True)
        if mutant == "drop-binom-alpha-gamma":
            k = pow(comb(be, de), C, p)
        elif mutant == "swap-the-two-exponents":
            k = (pow(comb(al, ga), C, p) * pow(comb(be, de), D, p)) % p
        else:
            k = (pow(comb(al, ga), D, p) * pow(comb(be, de), C, p)) % p
        uR = self.mul(self.mul(u2, u3), (k % p, 0))
        return self.resid(vL, uL) == self.resid(v2 + v3, uR), vL, v2 + v3


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def inert(T, N, p):
    """p odd, p does not divide Delta, and Delta is a non-residue mod p."""
    Delta = T * T - 4 * N
    return p > 2 and Delta % p != 0 and legendre(Delta, p) == -1


def digitpairs(p, cap=None):
    """(gamma, alpha) with 1 <= gamma <= alpha <= p-1, in the paper's order."""
    out = [(g, a) for a in range(1, p) for g in range(1, a + 1)]
    if cap is not None:
        out = [x for x in out if x[1] <= cap]
    return out


# ===========================================================================
# PART 4 -- the census drivers
# ===========================================================================

def census_A(T, N, p, mutant=None, cap=None):
    """engine A over a whole cell.  -> dict of counters."""
    E = EngineA(T, N, p)
    DG = digitpairs(p, cap)
    tot = bad = degen = valmm = vppos = nonint = 0
    first = None
    for (ga, al) in DG:
        for (de, be) in DG:
            for (A, B, C, D) in ABCD:
                tot += 1
                ok, vL, vR = E.check(A, B, C, D, al, be, ga, de, mutant=mutant)
                if not ok:
                    bad += 1
                    if first is None:
                        first = (A, B, C, D, al, be, ga, de)
                if vL != vR:
                    valmm += 1
                if vL > 0:
                    vppos += 1
                if vL < 0:
                    nonint += 1
                if A == C and B == D and al == ga and be == de:
                    degen += 1
    return {"total": tot, "violations": bad, "degenerate": degen, "val_mismatch": valmm,
            "vp_positive": vppos, "not_p_integral": nonint, "first": first}


def census_B(T, N, p, mutant=None, cap=None):
    """engine B over a whole cell: exact, no reduction.  -> dict of counters."""
    DG = digitpairs(p, cap)
    tot = bad = nonint = 0
    first = None
    for (ga, al) in DG:
        for (de, be) in DG:
            for (A, B, C, D) in ABCD:
                tot += 1
                ok, integral = check_B(A, B, C, D, al, be, ga, de, T, N, p, mutant=mutant)
                if not ok:
                    bad += 1
                    if first is None:
                        first = (A, B, C, D, al, be, ga, de)
                if not integral:
                    nonint += 1
    return {"total": tot, "violations": bad, "not_p_integral": nonint, "first": first}


# ===========================================================================
say("verify.py -- A Proof of Calderon's mod-p Bailey Congruence at Every Odd Inert Prime")
say("             independent re-derivation of every computed quantity in paper.tex")
say("")

# ---------------------------------------------------------------------------
say("--- 1. the definition of R_omega reproduces the source's own printed values -------")
# arXiv:2608.00347v1 prints these Gaussian values (omega = i, i.e. T = 0, N = 1).
for (A, B, C, D, want, label) in [(2, 2, 1, 1, Fraction(2), "Rect(2,2;1,1) = 2"),
                                  (4, 4, 2, 2, Fraction(30), "Rect(4,4;2,2) = 30"),
                                  (6, 6, 3, 3, Fraction(20008, 5), "Rect(6,6;3,3) = 20008/5")]:
    num, den = RZ(A, B, C, D, 0, 1)
    x, y = coords(num, den, 0, 1)
    chk("printed-value-%d-%d-%d-%d" % (A, B, C, D), (x, y), (want, Fraction(0)),
        "%s ; computed %s + (%s)i" % (label, x, y))
num, den = RZ(5, 3, 1, 1, 0, 1)
x, y = coords(num, den, 0, 1)
chk("printed-value-5-3-1-1", (x, y), (Fraction(4), Fraction(-1)),
    "R_i(5,3;1,1) = (5+3i)/(1+i) = 4 - i ; computed %s + (%s)i" % (x, y))

# ---------------------------------------------------------------------------
say("")
say("--- 2. every family named in the paper has the inertness the paper claims ---------")
INERT_FAMILIES = [(T, N, p) for (T, N, p, _, _) in CENSUS_C] + \
                 [(T, N, p) for (T, N, p, _, _) in CENSUS_D]
bad = []
for (T, N, p) in INERT_FAMILIES:
    Delta = T * T - 4 * N
    roots = [t for t in range(p) if (t * t - T * t + N) % p == 0]
    if not (inert(T, N, p) and roots == []):
        bad.append((T, N, p, Delta, roots))
chk("inertness-of-every-census-family", bad, [],
    "%d families: p odd, p does not divide Delta, (Delta/p) = -1, and x^2-Tx+N has no root "
    "mod p" % len(INERT_FAMILIES))

splitkind = []
for (T, N, p, kind, _, _, _, _, _) in CENSUS_E:
    Delta = T * T - 4 * N
    if p == 2:
        splitkind.append((kind, Delta % 2 == 0))
    else:
        splitkind.append((kind, legendre(Delta, p)))
chk("non-inert-primes-of-table-5", splitkind,
    [("split", 1), ("split", 1), ("split", 1), ("split", 1), ("ramified", True)],
    "the three Table 5 families are split (Legendre = +1) or ramified (p = 2 | Delta)")

sizes = []
for (T, N, p, inst, _) in CENSUS_C + CENSUS_D:
    d = len(digitpairs(p))
    sizes.append((d * d * ABCD_COUNT, inst))
chk("cell-sizes-match-the-tables", [a for (a, b) in sizes if a != b], [],
    "|{(gamma,alpha)}|^2 * %d equals the instance count printed for all %d families"
    % (ABCD_COUNT, len(sizes)))
chk("abcd-tuple-count", len(ABCD), ABCD_COUNT,
    "the nine (A,B,C,D) in {1,2}^4 with A >= C and B >= D, in the paper's order")

# ---------------------------------------------------------------------------
say("")
say("--- 3. the four arithmetic facts the proof rests on, checked in every family ------")
nfam = 0
bad_frob = []
bad_disc = []
bad_vert = []
bad_horiz = []
for (T, N, p) in INERT_FAMILIES:
    nfam += 1
    E = EngineA(T, N, p)
    w = (0, 1)
    wbar = ((T) % p, (-1) % p)
    if E.powm(w, p) != wbar:
        bad_frob.append((T, N, p))
    d = E.mul((-T % p, 2 % p), (-T % p, 2 % p))          # (2*omega - T)^2
    if d != ((T * T - 4 * N) % p, 0) or (T * T - 4 * N) % p == 0:
        bad_disc.append((T, N, p))
    # c1 = (omega - omegabar)/omega = (2*omega - T)/omega
    c1 = E.mul(((-T) % p, 2 % p), E.inv(w))
    c2 = ((T) % p, (-2) % p)                              # omegabar - omega = T - 2*omega
    for r in range(p):
        prod = E.ONE
        for t in range(p):
            prod = E.mul(prod, (r % p, t % p))
        if prod != E.mul(c1, (r % p, 0)):
            bad_vert.append((T, N, p, r))
    for s in range(p):
        prod = E.ONE
        for r in range(p):
            prod = E.mul(prod, (r % p, s % p))
        if prod != E.mul(c2, (s % p, 0)):
            bad_horiz.append((T, N, p, s))
chk("fact-ii-frobenius-omega-to-the-p", bad_frob, [],
    "omega^p == T - omega in O/pO, in all %d families" % nfam)
chk("fact-discriminant-is-a-unit", bad_disc, [],
    "(omega - omegabar)^2 == Delta and Delta is invertible mod p, in all %d families" % nfam)
chk("vertical-strip-constant", bad_vert, [],
    "prod_{t in F_p} (r + t*omega) == c1*r with c1 = (omega-omegabar)/omega, for every r in "
    "F_p in all %d families" % nfam)
chk("horizontal-strip-constant", bad_horiz, [],
    "prod_{r in F_p} (r + s*omega) == (omegabar-omega)*s, for every s in F_p in all %d "
    "families" % nfam)

# Fact (i): v_p(u + v*omega) = min(v_p u, v_p v), tested exactly (engine-B arithmetic).
bad_vfact = []
ntested = 0
for (T, N, p) in INERT_FAMILIES:
    cap = min(p * p, 60)
    for u in range(0, cap + 1):
        for v in range(0, cap + 1):
            if u == 0 and v == 0:
                continue
            ntested += 1
            a = (u, v)
            n = normZ(a, T, N)
            # v_p of an element of Z[omega] is v_p(Norm)/2 at an inert p
            vn = vp_inf(n, p)
            want = min(vp_inf(u, p), vp_inf(v, p))
            if vn != 2 * want:
                bad_vfact.append((T, N, p, u, v, vn, want))
chk("fact-i-valuation-of-u-plus-v-omega", bad_vfact, [],
    "v_p(Norm(u+v*omega)) == 2*min(v_p u, v_p v) on %d pairs across %d families -- this is "
    "Fact (i), since v_p(Norm(x)) = 2*v_p(x) at an inert p" % (ntested, nfam))

# ---------------------------------------------------------------------------
say("")
say("--- 4. the finite combinatorial steps of the proof, checked over all legal tuples -")
PS = sorted(set(p for (_, _, p) in INERT_FAMILIES))
bad_mult = []
bad_res = []
bad_zero = []
bad_card = []
ntup = 0
for p in PS:
    for (ga, al) in digitpairs(p):
        Sal = list(range(al - ga + 1, al + 1))
        for (A, C) in [(1, 1), (2, 1), (2, 2)]:
            ntup += 1
            Ap = p * A + al
            Cp = p * C + ga
            lo = Ap - Cp + 1
            mults = [u for u in range(lo, Ap + 1) if u % p == 0]
            want = [p * k for k in range(A - C + 1, A + 1)]
            if mults != want:
                bad_mult.append((p, al, ga, A, C, mults, want))
            mults0 = [u for u in range(1, Cp + 1) if u % p == 0]
            want0 = [p * k for k in range(1, C + 1)]
            if mults0 != want0:
                bad_mult.append((p, al, ga, A, C, mults0, want0))
            for r in range(p):
                e = sum(1 for u in range(lo, Ap + 1) if u % p == r)
                if e != C + (1 if r in [x % p for x in Sal] else 0):
                    bad_res.append((p, al, ga, A, C, r, e))
                e0 = sum(1 for u in range(1, Cp + 1) if u % p == r)
                if e0 != C + (1 if 1 <= r <= ga else 0):
                    bad_res.append((p, al, ga, A, C, r, e0, "den"))
        if 0 in [x % p for x in Sal] or 0 in range(1, ga + 1):
            bad_zero.append((p, al, ga))
        prod = 1
        for r in Sal:
            prod = prod * r % p
        fact = 1
        for r in range(1, ga + 1):
            fact = fact * r % p
        if len(Sal) != ga or prod != fact * comb(al, ga) % p:
            bad_card.append((p, al, ga))
chk("step-1-p-divisible-factors-are-exactly-C-by-D", bad_mult, [],
    "the multiples of p in [A'-C'+1, A'] are exactly {pk : A-C+1 <= k <= A} (C of them) and "
    "those in [1,C'] are exactly {pk : 1 <= k <= C}, over %d (p,alpha,gamma,A,C) tuples" % ntup)
chk("step-2-residue-multiplicities", bad_res, [],
    "e(r) = C + [r in S_alpha] and e_den(r) = C + [r in 1..gamma], for every residue r in "
    "every one of those %d tuples" % ntup)
chk("step-2-zero-is-in-none-of-the-four-sets", bad_zero, [],
    "0 is in neither S_alpha nor {1..gamma}: this is exactly what 1 <= gamma <= alpha <= p-1 "
    "buys")
chk("step-4-cardinality-cancellation", bad_card, [],
    "|S_alpha| = gamma and prod_{r in S_alpha} r == gamma! * C(alpha,gamma) mod p, so the "
    "strip constant cancels by cardinality alone and is never evaluated")

# the integrality lemma, as a valuation identity and as an inequality
bad_vid = []
bad_vge = []
nbox = 0
for p in (2, 3, 5, 7):
    for A in range(1, 26):
        for C in range(1, min(A, 4) + 1):
            for B in range(1, 26):
                for D in range(1, min(B, 4) + 1):
                    nbox += 1
                    direct = 0
                    for u in range(A - C + 1, A + 1):
                        for v in range(B - D + 1, B + 1):
                            direct += min(vp_inf(u, p), vp_inf(v, p))
                    viaj = 0
                    j = 1
                    while p ** j <= max(A, B):
                        q = p ** j
                        Mu = A // q - (A - C) // q
                        Mv = B // q - (B - D) // q
                        viaj += Mu * Mv
                        j += 1
                    if direct != viaj:
                        bad_vid.append((p, A, B, C, D, direct, viaj))
                    den = 0
                    j = 1
                    while p ** j <= max(C, D):
                        q = p ** j
                        den += (C // q) * (D // q)
                        j += 1
                    if direct < den:
                        bad_vge.append((p, A, B, C, D, direct, den))
chk("lemma-valuation-identity", bad_vid, [],
    "v_p(F(box)) = sum_{j>=1} M_u(p^j)*M_v(p^j) on %d boxes over p in {2,3,5,7}" % nbox)
chk("lemma-integrality-inequality", bad_vge, [],
    "v_p(numerator box) >= v_p([1,C]x[1,D]) on the same %d boxes, i.e. R_omega is p-integral"
    % nbox)

# ---------------------------------------------------------------------------
say("")
say("--- 5. Table 1: the three worked instances, in Z[sqrt(-2)] at p = 7 --------------")
Tw, Nw, pw = WORKED_FAMILY
Ew = EngineA(Tw, Nw, pw)
for idx, (A, B, C, D, al, be, ga, de) in enumerate(WORKED, 1):
    primed = (pw * A + al, pw * B + be, pw * C + ga, pw * D + de)
    vL, uL = Ew.R(*primed)
    v2, u2 = Ew.R(A, B, C, D)
    v3, u3 = Ew.R(al, be, ga, de)
    k = (pow(comb(al, ga), D, pw) * pow(comb(be, de), C, pw)) % pw
    uR = Ew.mul(Ew.mul(u2, u3), (k, 0))
    okB, integral = check_B(A, B, C, D, al, be, ga, de, Tw, Nw, pw)
    say("    instance %d: (A,B,C,D)=(%d,%d,%d,%d) (alpha,beta,gamma,delta)=(%d,%d,%d,%d)"
        % (idx, A, B, C, D, al, be, ga, de))
    say("      (A',B',C',D') = %s ; v_7(LHS) = %d ; LHS = %d + %d*omega in F_49"
        % (primed, vL, uL[0], uL[1]))
    say("      R(A,B;C,D) = %d + %d*omega ; R(alpha,beta;gamma,delta) = %d + %d*omega ; "
        "C(%d,%d)^%d*C(%d,%d)^%d = %d mod 7"
        % (u2[0], u2[1], u3[0], u3[1], al, ga, D, be, de, C, k))
    say("      RHS = %d + %d*omega" % (uR[0], uR[1]))
    chk("worked-%d-engine-A-congruence" % idx, Ew.resid(vL, uL) == Ew.resid(v2 + v3, uR), True,
        "engine A: both sides reduce to %d + %d*omega in F_49" % (uR[0], uR[1]))
    chk("worked-%d-engine-B-congruence" % idx, okB, True,
        "engine B, exact in Z[omega] with no reduction and no stripping: LHS - RHS lies in "
        "7*O_{omega,7}")
    chk("worked-%d-lhs-is-a-unit" % idx, (vL, integral), (0, True),
        "v_7(LHS) = 0 and the left side is 7-integral")

# the first instance, spelled out in Section 3 of the paper
A, B, C, D, al, be, ga, de = WORKED[0]
primed = (pw * A + al, pw * B + be, pw * C + ga, pw * D + de)
chk("section-3-primed-indices", primed, W1_PRIMED, "(A',B',C',D') as printed")
num, den = RZ(A, B, C, D, Tw, Nw)
x, y = coords(num, den, Tw, Nw)
chk("section-3-R-of-2-2-1-1", (x, y), (Fraction(W1_R_ABCD[0]), Fraction(W1_R_ABCD[1])),
    "R_omega(2,2;1,1) = (2+2*omega)/(1+omega) = 2 exactly in Z[omega], not merely mod 7")
v3, u3 = Ew.R(al, be, ga, de)
chk("section-3-R-of-5-4-2-3-residue", (v3, u3), (0, W1_R_DIGITS_RESIDUE),
    "R_omega(5,4;2,3) = %d + %d*omega in F_49" % (u3[0], u3[1]))
chk("section-3-binomial-product", (comb(al, ga) ** D * comb(be, de) ** C,
                                  comb(al, ga) ** D * comb(be, de) ** C % pw),
    (W1_BINOMIAL, W1_BINOMIAL_MOD), "C(5,2)^1 * C(4,3)^1 = 40 = 5 mod 7")
vL, uL = Ew.R(*primed)
chk("section-3-both-sides", uL, W1_LHS_RESIDUE,
    "the left side is %d + %d*omega in F_49, and 2*4*5 = 40 = 5 mod 7 is the right side"
    % (uL[0], uL[1]))

# ---------------------------------------------------------------------------
say("")
say("--- 6. Table 2: the twelve inert families of hypothesis (6.3), exhausted ---------")
tot_C = 0
gauss_C = 0
target = None
for (T, N, p, inst, order) in CENSUS_C:
    t1 = time.time()
    r = census_A(T, N, p, mutant=None)
    tot_C += r["total"]
    if (T, N) == (0, 1):
        gauss_C += r["total"]
    if (T, N, p) == WORKED_FAMILY:
        target = r
    chk("census-C-T%d-N%d-p%d" % (T, N, p), (r["total"], r["violations"]), (inst, 0),
        "%s, Delta = %d, p = %d: %d of %d instances, 0 violations, %d not p-integral "
        "[%.1fs]" % (order, T * T - 4 * N, p, r["total"], inst, r["not_p_integral"],
                     time.time() - t1))
chk("census-C-total", tot_C, CENSUS_C_TOTAL,
    "%d instances over the twelve families, every one exhaustive, 0 violations" % tot_C)
chk("census-C-gaussian-subtotal", gauss_C, CENSUS_C_GAUSSIAN,
    "the three (T,N) = (0,1) families -- inside the source's own Theorem 1.3 -- account for "
    "%d of them, and are controls rather than evidence" % gauss_C)
chk("census-C-outside-theorem-1-3", tot_C - gauss_C, CENSUS_C_OUTSIDE_THM13,
    "so only %d instances lie outside the region the source already proves" % (tot_C - gauss_C))
chk("target-cell-degenerate", target["degenerate"], TARGET_DEGENERATE,
    "of the 3969 instances of the target cell, %d are forced (A=C, B=D, alpha=gamma, "
    "beta=delta makes both sides identically 1)" % target["degenerate"])
chk("target-cell-vp-positive", target["vp_positive"], TARGET_VP_POSITIVE,
    "no instance of the target cell has v_7(LHS) > 0 -- A,B,C,D in {1,2} cannot reach that "
    "regime, which is why section 10 of this program exhausts it separately")
chk("target-cell-valuation-agreement", target["val_mismatch"], TARGET_VAL_MISMATCH,
    "v_p(LHS) = v_p(RHS) in every instance, as Step 1 predicts")

# ---------------------------------------------------------------------------
say("")
say("--- 7. engine B: the eight p = 7 families of Table 2 re-verified exactly ---------")
tot_B = 0
for (T, N, p, inst, order) in CENSUS_C:
    if p != 7:
        continue
    t1 = time.time()
    r = census_B(T, N, p)
    tot_B += r["total"]
    chk("engine-B-T%d-N%d-p%d" % (T, N, p), (r["total"], r["violations"], r["not_p_integral"]),
        (inst, 0, 0),
        "%s: %d of %d instances re-verified with no reduction and no stripping, 0 violations "
        "[%.1fs]" % (order, r["total"], inst, time.time() - t1))
chk("engine-B-total", tot_B, 8 * 3969,
    "%d instances carry the verdict of two structurally disjoint deciders, agreeing instance "
    "by instance" % tot_B)

# ---------------------------------------------------------------------------
say("")
say("--- 8. Table 4: the forced-positive controls -- a wrong right-hand side must fail -")
for (name, want_bad, want_first) in MUTANTS:
    rA = census_A(Tw, Nw, pw, mutant=name)
    rB = census_B(Tw, Nw, pw, mutant=name)
    chk("mutant-%s-engine-A" % name, (rA["total"], rA["violations"], rA["first"]),
        (MUTANT_CELL_TOTAL, want_bad, want_first),
        "engine A: the mutated right-hand side is violated in %d of %d instances; first "
        "violating (A,B,C,D,alpha,beta,gamma,delta) = %s" % (rA["violations"], rA["total"],
                                                            rA["first"]))
    chk("mutant-%s-engine-B" % name, (rB["total"], rB["violations"], rB["first"]),
        (MUTANT_CELL_TOTAL, want_bad, want_first),
        "engine B independently reproduces the same count and the same first violating tuple")

# ---------------------------------------------------------------------------
say("")
say("--- 9. Table 3: the two hypotheses the theorem drops, exercised ------------------")
tot_D = 0
real_D = 0
small_D = 0
for (T, N, p, inst, order) in CENSUS_D:
    t1 = time.time()
    r = census_A(T, N, p)
    tot_D += r["total"]
    if T * T - 4 * N > 0:
        real_D += r["total"]
    if p <= 5:
        small_D += r["total"]
    chk("census-D-T%d-N%d-p%d" % (T, N, p), (r["total"], r["violations"]), (inst, 0),
        "%s, Delta = %d: %d of %d instances, 0 violations [%.1fs]"
        % (order, T * T - 4 * N, r["total"], inst, time.time() - t1))
chk("census-D-total", tot_D, CENSUS_D_TOTAL,
    "%d further instances, all outside hypothesis (6.3) and all inside the theorem" % tot_D)
chk("census-D-real-quadratic", real_D, CENSUS_D_REAL,
    "%d of them have Delta > 0, i.e. lie in a REAL quadratic order, which the conjecture as "
    "printed excludes" % real_D)
chk("census-D-small-primes", small_D, CENSUS_D_SMALLP,
    "%d of them have p <= 5, which hypothesis (6.3) excludes" % small_D)

# ---------------------------------------------------------------------------
say("")
say("--- 10. the regime v_p(R_omega(A,B;C,D)) > 0, exhausted separately -------------")
Td, Nd, pd = DEGEN_FAMILY
Ed = EngineA(Td, Nd, pd)
al, be, ga, de = DEGEN_DIGITS
tot = vppos = bad = badB = valmm = 0
for A in range(1, DEGEN_AB_CAP + 1):
    for C in range(1, min(A, DEGEN_CD_CAP) + 1):
        for B in range(1, DEGEN_AB_CAP + 1):
            for D in range(1, min(B, DEGEN_CD_CAP) + 1):
                tot += 1
                ok, vL, vR = Ed.check(A, B, C, D, al, be, ga, de)
                v2, _ = Ed.R(A, B, C, D, memo=True)
                if v2 > 0:
                    vppos += 1
                if vL != vR:
                    valmm += 1
                if not ok:
                    bad += 1
                okB, _ = check_B(A, B, C, D, al, be, ga, de, Td, Nd, pd)
                if not okB:
                    badB += 1
chk("degenerate-regime-instances", tot, DEGEN_INSTANCES,
    "1 <= A,B <= %d, 1 <= C <= min(A,%d), 1 <= D <= min(B,%d), digits (alpha,beta,gamma,"
    "delta) = %s fixed" % (DEGEN_AB_CAP, DEGEN_CD_CAP, DEGEN_CD_CAP, DEGEN_DIGITS))
chk("degenerate-regime-is-non-empty", vppos, DEGEN_VP_POSITIVE,
    "%d of the %d instances have v_7(R_omega(A,B;C,D)) > 0, so both sides are congruent to 0 "
    "there; the exhaustive cells of Tables 1 and 2 never reach this regime" % (vppos, tot))
chk("degenerate-regime-both-engines", (bad, badB, valmm), (0, 0, 0),
    "0 violations by engine A, 0 by engine B, and v_7(LHS) = v_7(RHS) in all %d instances "
    "including the %d where both sides vanish mod 7" % (tot, vppos))

# ---------------------------------------------------------------------------
say("")
say("--- 11. Table 5: dropping inertness -- the congruence then fails wholesale -------")
for (T, N, p, kind, cap, inst, want_bad, want_nonint, want_first) in CENSUS_E:
    t1 = time.time()
    r = census_B(T, N, p, cap=cap)
    chk("non-inert-T%d-N%d-p%d-cap%d" % (T, N, p, cap),
        (r["total"], r["violations"], r["not_p_integral"], r["first"]),
        (inst, want_bad, want_nonint, want_first),
        "%s p = %d, digits bounded to alpha,beta <= %d: %d of %d instances violate the "
        "congruence, and the left side is not even p-integral in %d of them; first violating "
        "tuple %s [%.1fs]" % (kind, p, cap, r["violations"], r["total"], r["not_p_integral"],
                              r["first"], time.time() - t1))

# ===========================================================================
say("")
say("--- what this program does NOT cover ---------------------------------------------")
say("NOT RE-RUN: the paper's Theorem 1 is a hand proof and is NOT machine-verified here. No "
    "proof assistant was used. Sections 3-4 above check the finite combinatorial content of "
    "its Steps 1, 2 and 4 and of its integrality Lemma over stated finite ranges, and "
    "Sections 6-11 corroborate the conclusion on finitely many instances; neither is a proof.")
say("NOT RE-RUN: the split cells of Table 5 at their FULL digit range. Those would be "
    "78*78*9 = 54756 instances at p = 13 and 55*55*9 = 27225 at p = 11. Table 5 is bounded to "
    "alpha,beta <= 4 and alpha,beta <= 6 and the bound is printed in the table.")
say("NOT RE-RUN: the ramified prime p = 2 is only formally exercised: 1 <= gamma <= alpha <= "
    "p-1 = 1 forces alpha = gamma = beta = delta = 1, so its 9 instances are degenerate and "
    "establish nothing about ramified primes.")
say("NOT RE-RUN: no minimality or optimality claim is tested. Nothing here searches for a "
    "counterexample outside the inert regime beyond the three families of Table 5, and "
    "nothing here bears on Conjectures 6.3 and 6.4 of the source (the reciprocal-moment and "
    "mod-p^{3k} Ljunggren statements), which are untouched and remain open.")
say("NOT RE-RUN: no prior-art search. This program makes no network access and reads no file; "
    "the standing of the result against the literature is argued in the paper, not here.")
say("NOT RE-RUN: the earlier working programs of this project are not executed by this file. "
    "verify.py was written independently against the paper and reproduces their counts.")

say("")
say("checks run: %d passed, %d failed, elapsed %.1fs" % (_npass, _nfail, time.time() - T0))
if _nfail:
    say("VERDICT: %d CHECKS FAILED" % _nfail)
    sys.exit(1)
say("VERDICT: ALL %d CHECKS PASS" % _npass)
sys.exit(0)
