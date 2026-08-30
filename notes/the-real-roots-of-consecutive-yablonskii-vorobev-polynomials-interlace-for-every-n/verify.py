#!/usr/bin/env python3
"""verify.py -- re-derives the computational content of paper.tex of

    "Gap-Only Interlacing of the Real Roots of Consecutive Yablonskii-Vorob'ev
     Polynomials, Given the Published Real-Root Count"

from the objects the paper exhibits, and from nothing else.

WHAT IT READS FROM THE PAPER.  One thing, printed in the paper and hard-coded below
so a referee can diff it against the typeset text:

  * SEED + RECURRENCE   Y_0 = 1, Y_1 = z and
                        Y_{n+1} Y_{n-1} = z Y_n^2 - 4 ( Y_n Y_n'' - (Y_n')^2 ),
                        which is the displayed recurrence of section 1 of the paper.
                        Every polynomial used below is generated from these three
                        lines by exact division in Z[z]; no coefficient table is
                        read in.

WHAT IT ALSO CARRIES, AS INTERNAL REGRESSION LITERALS.  Three tables of literals live
in this file and are compared against the generated objects.  They are NOT displayed in
the paper and no check below claims that they are; they pin this program's own output so
that an accidental change to the arithmetic above cannot pass unnoticed:

  * the coefficient lists of Y_2, Y_3, Y_4, Y_5;
  * isolating rational intervals for the real roots of Y_1..Y_5;
  * the merged root words of (Y_n, Y_{n+1}) for n = 1..12.

WHAT IT DOES NOT DO.  It does not prove the theorem.  The theorem of the paper holds
for every n and is proved by hand there; a referee needs no machine.  This program is
a CONTROL on that proof.  It checks the published inputs (F1)-(F4) as exact algebra,
it checks the INTERIOR of the induction -- the occupancy of the gaps of Z_n by the
roots of Y_{n-1} and of Y_{n+1}, and the resulting merged word -- as an exact
statement about isolated real roots for 1 <= n <= 40, it checks the purely
arithmetical "zero slack" core for
1 <= n <= 200000, and it runs an ANTI-CONTROL showing that the parity lemma ALONE
proves nothing, so the load-bearing role of the published root count is visible
rather than asserted.  Its own gaps are printed as `NOT RE-RUN:` lines at the end.

ARITHMETIC.  Exact integers and exact rationals only (fractions.Fraction).  There is
no floating-point number anywhere in this file and no tolerance of any kind: every
decision is an integer sign or an integer equality.  Standard library only -- no
numpy, no sympy, no external data file.

Python 3.9+.  Run as:  python3 verify.py       (about a minute on a laptop)
"""

import sys
from fractions import Fraction as Fr

# ---------------------------------------------------------------------------
# 0.  bookkeeping: one PASS line per check, and the verdict counts them
# ---------------------------------------------------------------------------
_PASSES = []
_FAILS = []


def check(name, cond, detail=''):
    if cond:
        _PASSES.append(name)
        print('PASS %s%s' % (name, (' -- ' + detail) if detail else ''), flush=True)
    else:
        _FAILS.append(name)
        print('FAILED CHECK %s%s' % (name, (' -- ' + detail) if detail else ''), flush=True)
    return bool(cond)


# ---------------------------------------------------------------------------
# 1.  exact integer polynomial arithmetic.  coefficient lists, LOW -> HIGH.
# ---------------------------------------------------------------------------
def padd(a, b):
    r = [0] * max(len(a), len(b))
    for i, c in enumerate(a):
        r[i] += c
    for i, c in enumerate(b):
        r[i] += c
    while r and r[-1] == 0:
        r.pop()
    return r


def pscal(a, k):
    r = [c * k for c in a]
    while r and r[-1] == 0:
        r.pop()
    return r


def pmul(a, b):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    while r and r[-1] == 0:
        r.pop()
    return r


def pder(a):
    return [i * a[i] for i in range(1, len(a))]


def pshift(a, k):
    return ([0] * k + a) if a else []


def pdivexact(a, b):
    """a / b in Z[z].  Raises ArithmeticError unless the division is exact."""
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        num = a[i + len(b) - 1]
        if num % b[-1] != 0:
            raise ArithmeticError('non-exact leading division')
        c = num // b[-1]
        q[i] = c
        if c:
            for j, y in enumerate(b):
                a[i + j] -= c * y
    if any(a):
        raise ArithmeticError('non-exact division: nonzero remainder')
    while q and q[-1] == 0:
        q.pop()
    return q


def deg(a):
    return len(a) - 1


def sgn(x):
    return 0 if x == 0 else (1 if x > 0 else -1)


# ---------------------------------------------------------------------------
# 2.  THE OBJECT, generated from the seed and recurrence PRINTED IN THE PAPER
# ---------------------------------------------------------------------------
NMAX = 40                 # the induction interior is verified for n = 1..NMAX
NTOP = NMAX + 1           # so Y_{NMAX+1} and its roots are needed
NALG = NMAX + 2           # algebraic identities are checked up to here

SEED0 = [1]               # Y_0 = 1        (paper, section 1)
SEED1 = [0, 1]            # Y_1 = z        (paper, section 1)

# INTERNAL REGRESSION LITERALS, not displayed in the paper: coefficient lists of
# Y_2..Y_5, LOW -> HIGH.
REGRESSION_COEFFS = {
    2: [4, 0, 0, 1],
    3: [-80, 0, 0, 20, 0, 0, 1],
    4: [0, 11200, 0, 0, 0, 0, 0, 60, 0, 0, 1],
    5: [-6272000, 0, 0, -3136000, 0, 0, 78400, 0, 0, 2800, 0, 0, 140, 0, 0, 1],
}

# INTERNAL REGRESSION LITERALS, not displayed in the paper: isolating rational
# intervals for the real roots of Y_1..Y_5, in increasing order.  `EXACT0` denotes the
# exact root z = 0.
EXACT0 = 'EXACT0'
REGRESSION_INTERVALS = {
    1: [EXACT0],
    2: [(Fr(-13, 8), Fr(-3, 2))],
    3: [(Fr(-23, 8), Fr(-11, 4)), (Fr(3, 2), Fr(13, 8))],
    4: [(Fr(-4), Fr(-31, 8)), EXACT0],
    5: [(Fr(-5), Fr(-39, 8)), (Fr(-5, 4), Fr(-9, 8)), (Fr(21, 8), Fr(11, 4))],
}

# INTERNAL REGRESSION LITERALS, not displayed in the paper: the merged root word of the
# pair (Y_n, Y_{n+1}); '1' marks a root of Y_{n+1} and '0' a root of Y_n, read left to
# right along the real line.
REGRESSION_WORDS = {
    1: '10', 2: '101', 3: '1010', 4: '10101', 5: '101010', 6: '1010101',
    7: '10101010', 8: '101010101', 9: '1010101010', 10: '10101010101',
    11: '101010101010', 12: '1010101010101',
}


def build(N):
    """Y_0 .. Y_N from the printed seed and the printed recurrence, in Z[z]."""
    Y = [SEED0[:], SEED1[:]]
    for n in range(1, N):
        d1 = pder(Y[n])
        d2 = pder(d1)
        num = padd(pshift(pmul(Y[n], Y[n]), 1),
                   pscal(padd(pmul(Y[n], d2), pscal(pmul(d1, d1), -1)), -4))
        Y.append(pdivexact(num, Y[n - 1]))     # raises unless exact
    return Y


print('=== 1. the object, generated from the recurrence printed in the paper ===')
Y = build(NALG + 1)

check('recurrence-exact-division', len(Y) == NALG + 2,
      'Y_2..Y_%d were produced by Y_{n+1} = (z Y_n^2 - 4(Y_n Y_n\'\' - (Y_n\')^2)) / Y_{n-1} '
      'with zero remainder in Z[z] at every one of the %d steps'
      % (NALG + 1, NALG))

bad = [n for n, c in REGRESSION_COEFFS.items() if Y[n] != c]
check('regression-coefficients', not bad,
      'the coefficient lists of Y_2, Y_3, Y_4, Y_5 held as literals inside this program '
      'agree with the generated ones -- an internal regression check on the generator; '
      'the paper displays no coefficient list and none is claimed here'
      if not bad else 'MISMATCH at n=%s' % bad)

bad = [n for n in range(NALG + 2) if not (deg(Y[n]) == n * (n + 1) // 2 and Y[n][-1] == 1)]
check('monic-and-degree', not bad,
      '(F4) Y_n is monic of degree n(n+1)/2 for 0 <= n <= %d; deg Y_%d = %d'
      % (NALG + 1, NALG + 1, deg(Y[NALG + 1])) if not bad else 'MISMATCH at n=%s' % bad)

# The sign lemma of the paper (Lemma 3), as an exact identity in Z[z]:
#     Y_{n+1} Y_{n-1} - 4 (Y_n')^2  =  Y_n * ( z Y_n - 4 Y_n'' ).
# At a root a of Y_n the right side vanishes, so Y_{n+1}(a) Y_{n-1}(a) = 4 Y_n'(a)^2.
bad = []
for n in range(1, NALG + 1):
    lhs = padd(pmul(Y[n + 1], Y[n - 1]), pscal(pmul(pder(Y[n]), pder(Y[n])), -4))
    rhs = pmul(Y[n], padd(pshift(Y[n], 1), pscal(pder(pder(Y[n])), -4)))
    if lhs != rhs:
        bad.append(n)
check('lemma-3-sign-identity', not bad,
      'Y_{n+1}Y_{n-1} - 4(Y_n\')^2 = Y_n (z Y_n - 4 Y_n\'\') exactly in Z[z] for '
      'n = 1..%d, so Y_{n+1}(a)Y_{n-1}(a) = 4 Y_n\'(a)^2 at every root a of Y_n'
      % NALG if not bad else 'MISMATCH at n=%s' % bad)

# the Fukutani-Okamoto-Umemura Wronskian identity
bad = []
for n in range(1, NALG + 1):
    lhs = padd(pmul(pder(Y[n + 1]), Y[n - 1]), pscal(pmul(Y[n + 1], pder(Y[n - 1])), -1))
    rhs = pscal(pmul(Y[n], Y[n]), 2 * n + 1)
    if lhs != rhs:
        bad.append(n)
check('wronskian-identity', not bad,
      "Y'_{n+1}Y_{n-1} - Y_{n+1}Y'_{n-1} = (2n+1) Y_n^2 exactly in Z[z] for n = 1..%d"
      % NALG if not bad else 'MISMATCH at n=%s' % bad)


def reduce_cube(q):
    """q(z) = z^e * F(z^3) -> (e, F).  Raises if q is not of that shape."""
    e = 0
    while q[e] == 0:
        e += 1
    for k in range(len(q)):
        if (k - e) % 3 != 0 and q[k] != 0:
            raise ArithmeticError('not of the form z^e * F(z^3)')
    return e, [q[k] for k in range(e, len(q), 3)]


RED = [reduce_cube(Y[n]) for n in range(NALG + 2)]
bad = [n for n in range(NALG + 2) if RED[n][0] != (1 if n % 3 == 1 else 0)]
check('cube-structure', not bad,
      'Y_n = z^{e_n} F_n(z^3) with F_n in Z[u], F_n(0) != 0 and e_n = [n = 1 mod 3], '
      'for 0 <= n <= %d' % (NALG + 1) if not bad else 'MISMATCH at n=%s' % bad)

# --- (F1): simple roots, verified as squarefreeness over a finite field -----
# gcd_{F_p}(F, F') = 1 IMPLIES gcd_{Q}(F, F') = 1: F is monic, so it does not drop
# degree mod p, and a nontrivial rational common factor can be taken primitive in
# Z[u] with unit leading coefficient, hence stays nontrivial mod p.  The implication
# therefore runs in the direction we need, and the check is a proof rather than a
# sample.  Two different primes are used so no single unlucky prime carries it.
P1 = (1 << 31) - 1                     # 2147483647
P2 = (1 << 61) - 1                     # 2305843009213693951


def pgcd_modp(a, b, p):
    a = [c % p for c in a]
    b = [c % p for c in b]
    while a and a[-1] == 0:
        a.pop()
    while b and b[-1] == 0:
        b.pop()
    while b:
        inv = pow(b[-1], p - 2, p)
        lb = len(b)
        while len(a) >= lb:
            c = (a[-1] * inv) % p
            sh = len(a) - lb
            if c:
                a[sh:sh + lb] = [(a[sh + j] - c * y) % p for j, y in enumerate(b)]
            while a and a[-1] == 0:
                a.pop()
            if not a:
                break
        a, b = b, a
    return a


def coprime_modp(f, g, p):
    return deg(pgcd_modp(f, g, p)) == 0


for tag, p in (('F_%d' % P1, P1), ('F_%d' % P2, P2)):
    bad = [n for n in range(1, NTOP + 1)
           if deg(RED[n][1]) >= 1 and not coprime_modp(RED[n][1], pder(RED[n][1]), p)]
    check('simple-roots-F1-over-' + tag, not bad,
          '(F1) gcd(F_n, F_n\') = 1 over %s for n = 1..%d, hence over Q, hence every '
          'root of Y_n is simple' % (tag, NTOP) if not bad else 'MISMATCH at n=%s' % bad)

bad = [n for n in range(1, NTOP + 1) if not coprime_modp(RED[n][1], RED[n + 1][1], P1)]
check('coprime-consecutive-F2', not bad,
      '(F2) gcd(F_n, F_{n+1}) = 1 over F_%d for n = 1..%d; since e_n and e_{n+1} are '
      'never both 1 this gives gcd(Y_n, Y_{n+1}) = 1, so no root of Y_{n+1} is a root '
      'of Y_n' % (P1, NTOP) if not bad else 'MISMATCH at n=%s' % bad)

bad = [n for n in range(1, NTOP + 1) if not coprime_modp(RED[n - 1][1], RED[n + 1][1], P1)]
check('coprime-index-2', not bad,
      'gcd(F_{n-1}, F_{n+1}) = 1 over F_%d for n = 1..%d, so the interval endpoints used '
      'in the parity-transfer lemma (Lemma 4) are non-roots of Y_{n-1} as well as of '
      'Y_{n+1}' % (P1, NTOP)
      if not bad else 'MISMATCH at n=%s' % bad)

# --- signs at infinity, and the odd degree gap ------------------------------
bad = [n for n in range(1, NALG + 1) if Y[n][-1] != 1]
gap = [n for n in range(1, NALG + 1) if deg(Y[n + 1]) - deg(Y[n - 1]) != 2 * n + 1]
check('signs-at-infinity-and-odd-gap', (not bad) and (not gap),
      'Y_n monic gives sign Y_n(+inf) = +1 and sign Y_n(-inf) = (-1)^{n(n+1)/2}; and '
      'deg Y_{n+1} - deg Y_{n-1} = 2n+1, an ODD number, for n = 1..%d -- the single '
      'source of the left-ray asymmetry in the parity-transfer lemma (Lemma 4)' % NALG
      if (not bad) and (not gap) else 'MISMATCH bad=%s gap=%s' % (bad, gap))

# ---------------------------------------------------------------------------
# 3.  EXACT REAL-ROOT ISOLATION.  no floating point, no tolerance.
# ---------------------------------------------------------------------------
def ev_num(c, p, q):
    """q^deg(c) * c(p/q) as an exact integer, for integers p and q >= 1."""
    D = len(c) - 1
    qp = [1] * (D + 1)
    for k in range(1, D + 1):
        qp[k] = qp[k - 1] * q
    v = c[D]
    for i in range(D - 1, -1, -1):
        v = v * p + c[i] * qp[D - i]
    return v


def sgn_F_at_z(F, z):
    """sign of F(z^3) at a rational z.  den(z^3) > 0, so no sign flip is possible."""
    u = Fr(z) ** 3
    return sgn(ev_num(F, u.numerator, u.denominator))


def fujiwara_z_bound(F):
    """An integer Z with |z| < Z for every real z satisfying F(z^3) = 0.

    Fujiwara's bound: every complex root u of a monic F of degree D obeys
      |u| <= 2 max( |c_{D-1}|, |c_{D-2}|^{1/2}, ..., |c_1|^{1/(D-1)}, |c_0/2|^{1/D} ).
    The maximum is bracketed by INTEGER comparisons |c_{D-j}| <= B^j -- no roots and
    no logarithms are taken -- and then |z| = |u|^{1/3} <= (2B)^{1/3}.
    """
    D = len(F) - 1
    lc = abs(F[D])
    B = 1
    while True:
        good, bj = True, B
        for j in range(1, D + 1):
            v = abs(F[D - j])
            if j == D:
                v = (v + 1) // 2
            if v > lc * bj:
                good = False
                break
            bj *= B
        if good:
            break
        B *= 2
    ub = 2 * B
    z = 1
    while z ** 3 <= ub:
        z += 1
    return z


def isolate(F, want, tag, rounds=40):
    """Isolating intervals for the real roots of z |-> F(z^3), sorted increasingly.

    Each entry is a pair (a, b) of rationals with a < b holding exactly one root, or
    (r, r) for an exact rational root.  `want` is the number of real roots predicted
    by the paper's cited input (F3); the routine REFUSES to return unless it has
    exactly that many.  That is what makes the isolation complete: a cell showing a
    sign change holds an odd number of roots and a cell showing none holds an even
    number, so once the number of sign-change cells reaches `want` every one of them
    holds exactly one root and every other cell holds none.
    """
    if deg(F) == 0:
        if want != 0:
            raise ArithmeticError('%s: F is constant but %d real roots wanted' % (tag, want))
        return []
    Z = fujiwara_z_bound(F)
    grid = [Fr(j) for j in range(-Z, Z + 1)]
    sig = [sgn_F_at_z(F, g) for g in grid]
    exact = [g for g, s in zip(grid, sig) if s == 0]
    cells = [[grid[i], grid[i + 1], sig[i], sig[i + 1]] for i in range(len(grid) - 1)]
    for _ in range(rounds):
        change = [c for c in cells if c[2] * c[3] < 0]
        if len(change) + len(exact) == want:
            out = [(c[0], c[1]) for c in change] + [(r, r) for r in exact]
            out.sort()
            return out
        nxt = []
        for c in cells:
            if c[2] * c[3] < 0:
                nxt.append(c)
                continue
            m = (c[0] + c[1]) / 2
            s = sgn_F_at_z(F, m)
            if s == 0:
                exact.append(m)
            nxt.append([c[0], m, c[2], s])
            nxt.append([m, c[1], s, c[3]])
        cells = nxt
    raise ArithmeticError('%s: isolation did not reach the %d roots (F3) predicts'
                          % (tag, want))


def refine(F, a, b):
    """One bisection of an isolating interval, keeping the root strictly inside."""
    if a == b:
        return a, b
    m = (a + b) / 2
    s = sgn_F_at_z(F, m)
    if s == 0:
        return m, m
    return (a, m) if s * sgn_F_at_z(F, a) < 0 else (m, b)


print()
print('=== 2. exact real-root isolation and the published count (F3) ===')

M = {n: (n + 1) // 2 for n in range(NALG + 2)}          # (F3): m_n = floor((n+1)/2)

ROOTS = {0: []}
for n in range(1, NTOP + 1):
    e, F = RED[n]
    iv = isolate(F, M[n] - e, 'Y_%d' % n)
    if e:
        iv = sorted(iv + [(Fr(0), Fr(0))])
    ROOTS[n] = iv

bad = [n for n in range(1, NTOP + 1) if len(ROOTS[n]) != M[n]]
check('root-count-F3', not bad,
      '(F3) Y_n has exactly floor((n+1)/2) real roots for n = 1..%d, each confined to '
      'its own isolating rational interval; the count is what certifies the isolation '
      'COMPLETE and not merely nonempty' % NTOP if not bad else 'MISMATCH at n=%s' % bad)


def sturm_real_root_count(F):
    """Number of distinct real roots of F over Q, by Sturm's theorem."""
    def rem(a, b):
        a = a[:]
        while len(a) >= len(b) and any(a):
            c = a[-1] / b[-1]
            sh = len(a) - len(b)
            for j, y in enumerate(b):
                a[sh + j] -= c * y
            while a and a[-1] == 0:
                a.pop()
        return a

    chain = [[Fr(x) for x in F], [Fr(i) * Fr(F[i]) for i in range(1, len(F))]]
    while chain[-1] and deg(chain[-1]) > 0:
        r = rem(chain[-2], chain[-1])
        if not r:
            break
        chain.append([-x for x in r])
    chain = [c for c in chain if c]

    def variations(at_plus):
        s = []
        for c in chain:
            g = sgn(c[-1])
            if not at_plus and deg(c) % 2 == 1:
                g = -g
            s.append(g)
        s = [x for x in s if x]
        return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])

    return variations(False) - variations(True)


NSTURM = 12
bad = []
for n in range(1, NSTURM + 1):
    e, F = RED[n]
    got = (sturm_real_root_count(F) if deg(F) >= 1 else 0) + e
    if got != M[n]:
        bad.append((n, got, M[n]))
check('sturm-cross-check', not bad,
      'Sturm\'s theorem, computed with no reference to the grid and with no use of '
      '(F3), returns the same real-root count floor((n+1)/2) for n = 1..%d' % NSTURM
      if not bad else 'MISMATCH %s' % bad)

bad = []
for n, want in REGRESSION_INTERVALS.items():
    e, F = RED[n]
    got = ROOTS[n]
    if len(got) != len(want):
        bad.append((n, 'count %d vs %d' % (len(want), len(got))))
        continue
    for w, (ga, gb) in zip(want, got):
        if w is EXACT0:
            if not (ga == 0 and gb == 0):
                bad.append((n, 'expected the exact root 0'))
            continue
        wa, wb = w
        # An exact sign change of Y_n across the literal endpoints -- so the literal
        # interval really does contain a root -- plus containment of the literal
        # interval in the (coarser) interval this program isolated independently, so
        # that the root it contains is the SAME root.  Here
        # sign Y_n(z) = sign(z)^{e_n} * sign F_n(z^3).
        sa = sgn_F_at_z(F, wa) * (sgn(wa) ** e)
        sb = sgn_F_at_z(F, wb) * (sgn(wb) ** e)
        if not (sa * sb < 0 and ga <= wa and wb <= gb):
            bad.append((n, str(w)))
check('regression-intervals', not bad,
      'every rational interval held as a literal inside this program carries an exact '
      'sign change of Y_n and sits inside the interval this program isolated '
      'independently, so it brackets the same root; checked for n = 1..5 -- an internal '
      'regression check on the isolator, the paper displays no such interval and none '
      'is claimed here' if not bad else 'MISMATCH %s' % bad)

# ---------------------------------------------------------------------------
# 4.  THE INTERIOR OF THE INDUCTION, exactly, for n = 1..NMAX
# ---------------------------------------------------------------------------
print()
print('=== 3. the interior of the induction: gap occupancy and the merged word, '
      'n = 1..%d ===' % NMAX)


def merged_order(n):
    """Tags of Z_{n-1} u Z_n u Z_{n+1} in increasing order of the root.

    Intervals belonging to DIFFERENT polynomials are bisected until pairwise
    disjoint, so the order is decided by exact rational comparison and never by a
    numerical estimate.  The coprimality checked above is what makes the roots of two
    different Y's distinct, hence what makes the separation terminate.
    """
    items = [[t, a, b] for t in (n - 1, n, n + 1) for (a, b) in ROOTS[t]]
    for _ in range(400):
        items.sort(key=lambda r: (r[1], r[2]))
        clash = None
        for i in range(len(items) - 1):
            A, B = items[i], items[i + 1]
            if A[0] != B[0] and A[2] > B[1]:
                clash = (i, i + 1)
                break
        if clash is None:
            return [r[0] for r in items]
        for i in clash:
            r = items[i]
            if r[1] != r[2]:
                r[1], r[2] = refine(RED[r[0]][1], r[1], r[2])
    raise ArithmeticError('n=%d: the roots could not be separated' % n)


rows, bad1, bad3, bad4, badw = [], [], [], [], []
for n in range(1, NMAX + 1):
    order = merged_order(n)
    idx = [i for i, t in enumerate(order) if t == n]
    if not idx:
        raise ArithmeticError('Z_n is empty at n=%d' % n)
    even = (n % 2 == 0)
    occ = {}
    for t in (n - 1, n + 1):
        occ['L', t] = sum(1 for i in range(0, idx[0]) if order[i] == t)
        occ['R', t] = sum(1 for i in range(idx[-1] + 1, len(order)) if order[i] == t)
        occ['G', t] = [sum(1 for i in range(idx[k] + 1, idx[k + 1]) if order[i] == t)
                       for k in range(len(idx) - 1)]
    s1 = (occ['L', n - 1] == 0 and all(c == 1 for c in occ['G', n - 1])
          and occ['R', n - 1] == (1 if even else 0))
    s3 = (occ['L', n + 1] == 1 and all(c == 1 for c in occ['G', n + 1])
          and occ['R', n + 1] == (1 if even else 0))
    word = ''.join('1' if t == n + 1 else '0' for t in order if t in (n, n + 1))
    s4 = (all(word[i] != word[i + 1] for i in range(len(word) - 1))
          and word.startswith('1') and len(word) == n + 1)
    if not s1:
        bad1.append(n)
    if not s3:
        bad3.append(n)
    if not s4:
        bad4.append(n)
    if n in REGRESSION_WORDS and REGRESSION_WORDS[n] != word:
        badw.append(n)
    print('   n=%2d  (|Z_{n-1}|,|Z_n|,|Z_{n+1}|) = (%2d,%2d,%2d)   '
          'N_{n-1}: L=%d G=%s R+=%d   N_{n+1}: L=%d G=%s R+=%d   word=%s'
          % (n, len(ROOTS[n - 1]), len(ROOTS[n]), len(ROOTS[n + 1]),
             occ['L', n - 1], ''.join(map(str, occ['G', n - 1])) or '-', occ['R', n - 1],
             occ['L', n + 1], ''.join(map(str, occ['G', n + 1])) or '-', occ['R', n + 1],
             word), flush=True)
    rows.append(word)

check('induction-hypothesis-occupancy', not bad1,
      'N_{n-1}(L) = 0, N_{n-1}(G_i) = 1 in every gap and N_{n-1}(R+) = [n even], for '
      'n = 1..%d -- this is S(n-1) read off against the roots of Y_n' % NMAX
      if not bad1 else 'MISMATCH at n=%s' % bad1)
check('induction-forced-occupancy', not bad3,
      'N_{n+1}(L) = 1, N_{n+1}(G_i) = 1 in every gap and N_{n+1}(R+) = [n even], for '
      'n = 1..%d -- exactly the occupancy the proof forces, with nothing to spare'
      % NMAX if not bad3 else 'MISMATCH at n=%s' % bad3)
check('induction-merged-word', not bad4,
      'the merged word of (Y_n, Y_{n+1}) strictly alternates, begins with a root of '
      'Y_{n+1} and has length n+1, for n = 1..%d' % NMAX
      if not bad4 else 'MISMATCH at n=%s' % bad4)
check('regression-words', not badw,
      'the merged words held as literals inside this program are reproduced exactly for '
      'n = 1..%d -- an internal regression check on the merge; the paper displays only '
      'the general parity-dependent form of the word, not this list, and no agreement '
      'with a list in the paper is claimed here' % max(REGRESSION_WORDS)
      if not badw else 'MISMATCH at n=%s' % badw)

# ---------------------------------------------------------------------------
# 5.  THE ARITHMETICAL CORE, over a range no census could reach
# ---------------------------------------------------------------------------
print()
print('=== 4. the counting core of the induction, n = 1..200000 ===')
NBIG = 200000


def m_(k):
    return (k + 1) // 2


bad = [n for n in range(1, NBIG + 1) if m_(n + 1) != m_(n) + (1 if n % 2 == 0 else 0)]
check('count-identity', not bad,
      'floor((n+2)/2) = floor((n+1)/2) + [n even] for n = 1..%d' % NBIG
      if not bad else 'MISMATCH at n=%s' % bad[:5])

bad = [n for n in range(1, NBIG + 1)
       if 1 + (m_(n) - 1) + (1 if n % 2 == 0 else 0) != m_(n + 1)]
check('zero-slack', not bad,
      'the forced minimum 1 + (m_n - 1) + [n even] equals m_{n+1} EXACTLY, in both '
      'parities, for n = 1..%d: the sum of the lower bounds already exhausts the '
      'published count, so every interval count is pinned at its minimum' % NBIG
      if not bad else 'MISMATCH at n=%s' % bad[:5])

bad = [n for n in range(1, NBIG + 1)
       if (abs(m_(n + 1) - m_(n)) == 1) != (n % 2 == 0)]
check('stronger-definition-variant', not bad,
      '|m_{n+1} - m_n| = 1 holds for EVEN n and fails for every ODD n (n = 1..%d): '
      'under the stronger variant of the definition of interlacing described in section '
      '1 of the paper -- the variant that also demands that the two root counts differ '
      'by one -- the statement is true for every even n and FALSE for every odd n' % NBIG
      if not bad else 'MISMATCH at n=%s' % bad[:5])

# ---------------------------------------------------------------------------
# 6.  CONTROLS IN BOTH POLARITIES
# ---------------------------------------------------------------------------
print()
print('=== 5. controls ===')


def ev_fr(c, x):
    s = Fr(0)
    for i, k in enumerate(c):
        s += Fr(k) * Fr(x) ** i
    return s


def alternates(seq):
    return all(seq[i] != seq[i + 1] for i in range(len(seq) - 1))


check('decider-forced-yes-and-no',
      alternates('abab') and not alternates('aabb'),
      'the alternation decider accepts the forced-YES word abab and rejects the '
      'forced-NO word aabb')

# THE ANTI-CONTROL.  M = (x-1)(x-2) plays the role of Y_n.  A = 2x-3 has ONE root in
# the gap (1,2); B = (5x-6)(5x-7)(5x-8) has THREE roots, all of them BELOW 1.  A and B
# take the same sign at x=1 and the same sign at x=2, so the hypothesis AND the
# conclusion of the parity-transfer lemma (Lemma 4 of the paper) hold for both -- yet
# (M,A) interlace and (M,B) do not.  That lemma alone therefore proves nothing; the
# published count (F3) is load-bearing.
Mc = pmul([-1, 1], [-2, 1])
A = [-3, 2]
B = pmul(pmul([-6, 5], [-7, 5]), [-8, 5])
same = (sgn(ev_fr(A, 1)) == sgn(ev_fr(B, 1))) and (sgn(ev_fr(A, 2)) == sgn(ev_fr(B, 2)))
zeroM = [Fr(1), Fr(2)]
zeroA = [Fr(3, 2)]
zeroB = [Fr(6, 5), Fr(7, 5), Fr(8, 5)]


def interlaces(p_roots, q_roots):
    """between any two roots of p lies a root of q, and conversely (the gap-only
    definition quoted in section 1 of the paper)."""
    for lo, hi in ((p_roots, q_roots), (q_roots, p_roots)):
        for i in range(len(lo) - 1):
            if not any(lo[i] < r < lo[i + 1] for r in hi):
                return False
    return True


check('anti-control-lemma-4-alone-is-insufficient',
      same and interlaces(zeroM, zeroA) and not interlaces(zeroM, zeroB)
      and all(sgn(ev_fr(Mc, r)) != 0 for r in zeroA + zeroB),
      'A = 2x-3 and B = (5x-6)(5x-7)(5x-8) take the same signs as each other at both '
      'roots of M = (x-1)(x-2), so both satisfy the hypothesis and the conclusion of '
      'the parity-transfer lemma (Lemma 4 of the paper); but (M,A) interlace and (M,B) '
      'do NOT.  That lemma alone is strictly insufficient and the published count (F3) '
      'is load-bearing.')

# A FORCED-NEGATIVE control on the isolation routine itself: perturb the constant term
# of Y_3 from -80 to +200.  The perturbed polynomial has NO real root, so the count
# floor((n+1)/2) = 2 becomes unattainable, and the routine must report that rather
# than return a pass.
tampered = Y[3][:]
tampered[0] = 200
te, tF = reduce_cube(tampered)
reported = False
try:
    isolate(tF, M[3] - te, 'tampered-Y_3', rounds=4)
except ArithmeticError:
    reported = True
check('forced-negative-control',
      reported and sturm_real_root_count(tF) == 0 and te == 0,
      'replacing the constant term of Y_3 by +200 leaves a polynomial with no real '
      'root at all (Sturm: 0), so the count (F3) predicts is unattainable and '
      'isolate() raises instead of returning; the same routine that reports the passes '
      'above does report a failure when there is one')

# ---------------------------------------------------------------------------
# 7.  WHAT THIS PROGRAM DOES NOT COVER
# ---------------------------------------------------------------------------
print()
print('NOT RE-RUN: the published inputs are used as inputs.  The program re-derives '
      '(F1), (F2) and the COUNT asserted by (F3) for n <= %d by exact algebra, but it '
      'does not verify the PROOFS of Fukutani-Okamoto-Umemura (Nagoya Math. J. 159 '
      '(2000) 179-200), of Clarkson (Semin. Congr. 14 (2006) 21-52) or of Roffelsen '
      '(SIGMA 8 (2012) 099).  The theorem of the paper is conditional on those, exactly '
      'as it says.' % NTOP)
print('NOT RE-RUN: the theorem holds for EVERY n and is proved by hand in the paper.  '
      'The root-order data above is a control on the interior of that induction for '
      '1 <= n <= %d only.  n > %d is covered by the proof and by nothing in this file.'
      % (NMAX, NMAX))
print('NOT RE-RUN: interlacing means throughout the GAP-ONLY definition quoted in '
      'section 1 of the paper.  Nothing in this file establishes the statement under the '
      'stronger variant that also demands |m_{n+1} - m_n| = 1; the '
      'stronger-definition-variant check above records instead that under that variant '
      'the statement FAILS for every odd n.')
print('NOT RE-RUN: no claim of minimality, extremality or uniqueness is tested; no '
      'polynomial family other than the Yablonskii-Vorob\'ev family is examined; and '
      'the generalised Okamoto polynomials that are the subject of the source paper are '
      'not touched.')
print('NOT RE-RUN: the bibliographic material of section 1 -- the quoted sentence, the '
      'definition of interlacing and the numbered results attributed to '
      'arXiv:2402.15887v1 and to the other sources -- is not machine-checked here.  It '
      'was read off those sources by hand, and the journal version of record was not '
      'reachable, so no preprint-to-journal diff of the quoted sentence was performed.')

print()
print('=== VERDICT ===')
if _FAILS:
    print('VERDICT: %d OF %d CHECKS FAILED: %s'
          % (len(_FAILS), len(_FAILS) + len(_PASSES), ', '.join(_FAILS)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(_PASSES))
sys.exit(0)
