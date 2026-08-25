#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- referee verification program for

    "The Least Counterexample to Conjecture 4.3 of Dolfi, Hafezieh, and Spiga"

Conjecture 4.3 of [DHS], as quoted in the paper:

    gcd(2^n - 1, n) = 1  and  omega(n) = l >= 3   ==>   omega(2^n - 1) >= 2^l .

The paper refutes it and claims n = 170 is the LEAST counterexample.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs; never used as evidence for themselves)
--------------------------------------------------------------------------
  * PAPER_N            = 170, with claimed factorisation 2 * 5 * 17.
  * PAPER_FACTORS      = the seven claimed prime factors of 2^170 - 1 in
                         equation (1) of the paper.
  * PAPER_SMALL_TRIAL  = the five factors the paper calls "prime by trial
                         division"; PAPER_BASE_PRIMES = the five auxiliary
                         base primes 709, 2437, 2857, 6367, 27953 of
                         Proposition 5.
  * PAPER_POCK_P/_Q    = the six rows (N, prime factors of F, witness a) of
                         the certificate table in the proof of Prop. 5.
  * PAPER_SQFREE_LIST  = the 13 squarefree three-prime products below 170
                         listed in the minimality proof.
  * PAPER_CENSUS       = {70, 130, 154, 165}, the survivors the paper claims.
  * PAPER_EXCL_ORDERS  = the paper's exclusion data ord_7(2) = 3 for n = 105
                         and ord_11(2) = 10 for n = 110.
  * PAPER_TABLE_FACS   = the cyclotomic factorisations Phi_35(2), Phi_130(2),
                         Phi_11(2) displayed in the minimality table.
  * PAPER_OMEGA_N = 3, PAPER_OMEGA_MERSENNE = 7, PAPER_TAU_N = 8.

Nothing else is read from the paper; in particular no factorisation, order,
census or primality verdict below is copied from it.

--------------------------------------------------------------------------
WHAT IS DERIVED HERE (the checks; all exact integer arithmetic)
--------------------------------------------------------------------------
  * 2^170 - 1 is computed as an integer and compared with the product of the
    seven claimed factors; the small factors are proved prime by bare trial
    division and the two 19/20-digit factors twice: once by executing the
    paper's Pocklington chains, and once independently by a Lucas n-1
    certificate built from a factorisation of N-1 computed here.
  * Hence omega(2^170 - 1) is DERIVED to be exactly 7, and the whole
    factorisation is re-derived from scratch out of the cyclotomic values
    Phi_d(2), d | 170, and compared with the paper's list.  Each Phi_d(2) is
    computed twice by structurally different routes -- exact big-integer
    division of 2^d - 1, and long division of x^d - 1 in Z[x] followed by
    evaluation at 2 -- and the two must agree, with deg Phi_d = phi(d) and
    Phi_d monic.  (On their own, "prod_{d|n} Phi_d(2) = 2^n - 1" and
    "Phi_d(2) | 2^n - 1" are identities that hold by construction and can
    never fail; the cross-check is what makes this route falsifiable.)
  * gcd(2^170 - 1, 170) is computed with big-integer gcd; the paper's
    mechanism (ord_5(2) = 4, ord_17(2) = 8, neither divides 170) is computed.
  * The refutation 7 < 8 = 2^omega(170) is computed, not asserted.
  * MINIMALITY: every n in [1, 169] is enumerated, omega(n) and
    gcd(2^n - 1, n) computed, and for each survivor 2^n - 1 is factored
    completely (cyclotomic decomposition + Pollard-Brent + Lucas
    certificates) so that omega(2^n - 1) >= 2^omega(n) is verified by
    computation. No probabilistic primality result is load-bearing:
    Miller-Rabin is used only to steer the factoriser, and every prime that
    enters a verdict carries a trial-division or Lucas certificate.
  * The paper's auxiliary numeric claims (the 13-element squarefree list, the
    exclusion reasons for 6 | n, 105, 110, the table factorisations, the
    non-squarefree tau bound) are recomputed.
  * The whole theorem is recomputed independently in section F: n is scanned
    upward from 1 and the first counterexample found must be 170.  The scan
    FAILS CLOSED: any exponent whose factorisation could not be completed is
    reported as undecided and makes F1 fail, so an exhausted factoriser can
    never be read as "no counterexample here".
  * A mutation self-test corrupts the exhibited object (a wrong factor, a
    product-PRESERVING corruption that is not a prime factorisation, a wrong
    P AND a wrong Q -- each rejected by its Pocklington row and by the
    independent Lucas certifier -- a non-counterexample exponent, known
    composites) and confirms that the
    load-bearing predicates then report failure.  The self-test calls the same
    predicate the C1 verdict uses, not a weaker stand-in.
  * Section G is SUPPLEMENTARY and goes beyond the paper: Lemma 2 is tested
    on every exponent n <= 100 satisfying the gcd hypothesis.

Standard library only; no external data files; no floating point in any
decision.  Runtime: about half a second.  Output: one PASS/FAIL line per
check and a final VERDICT line.  Exit status 0 iff every check passes.
"""

import math
import random
import sys
import time

random.seed(20260822)

_CHECKS = []


def check(name, ok, detail=""):
    """Record one check.  ok must be a real bool computed from data."""
    _CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return bool(ok)


def verdict():
    n = len(_CHECKS)
    bad = [nm for nm, ok in _CHECKS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


# ----------------------------------------------------------------------
# Arithmetic toolbox (exact integers only)
# ----------------------------------------------------------------------

def is_prime_trial(m):
    """Bare trial division.  Rigorous; used for the small numbers the paper
    itself calls 'prime by trial division'."""
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


def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i != n // i:
                ds.append(n // i)
        i += 1
    return sorted(ds)


def tau(n):
    return len(divisors(n))


def prime_factors_trial(n):
    """Complete prime factorisation by trial division; returns dict p->e.
    Only called on n < 10**7 in this program."""
    out = {}
    m = n
    d = 2
    while d * d <= m:
        while m % d == 0:
            out[d] = out.get(d, 0) + 1
            m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


def omega_small(n):
    return len(prime_factors_trial(n))


def is_squarefree_small(n):
    return all(e == 1 for e in prime_factors_trial(n).values())


def multiplicative_order(a, m):
    """Order of a modulo m by direct iteration (m small: m < 10**6 here)."""
    if math.gcd(a, m) != 1:
        return 0
    x = a % m
    k = 1
    while x != 1:
        x = (x * a) % m
        k += 1
        if k > m:
            return 0
    return k


_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)


def mr_probable_prime(n):
    """Miller-Rabin on fixed bases.  NOT load-bearing: it is only used to
    decide when the factoriser should stop splitting.  Every number that
    enters a verdict is afterwards certified by is_prime_trial or
    lucas_certified_prime, either of which would reject a pseudoprime."""
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def pollard_brent(n):
    """Return a nontrivial factor of composite n."""
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        y = random.randrange(1, n)
        c = random.randrange(1, n)
        m = 128
        g = r = q = 1
        x = ys = y
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2
        if g == n:
            g = 1
            y = ys
            while g == 1:
                y = (y * y + c) % n
                g = math.gcd(abs(x - y), n)
        if g != n:
            return g


def factorize(n):
    """Complete factorisation, returned as a dict p->e.  The primality of the
    returned parts is verified separately by lucas_certified_prime; the
    product identity prod p^e == n is checked by the callers."""
    out = {}
    m = n
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        while m % p == 0:
            out[p] = out.get(p, 0) + 1
            m //= p
    d = 49
    while d * d <= m and d < 200000:
        if m % d == 0:
            while m % d == 0:
                out[d] = out.get(d, 0) + 1
                m //= d
        d += 2
    if m == 1:
        return out
    stack = [m]
    while stack:
        cur = stack.pop()
        if cur == 1:
            continue
        if mr_probable_prime(cur):
            out[cur] = out.get(cur, 0) + 1
            continue
        f = pollard_brent(cur)
        stack.append(f)
        stack.append(cur // f)
    return out


_LUCAS_MEMO = {}


def lucas_certified_prime(n):
    """Deterministic primality proof (converse of Fermat / Lucas 1876):
    n is prime iff there is a with a^(n-1) = 1 mod n and a^((n-1)/q) != 1
    mod n for every prime q | n-1.  The factorisation of n-1 is computed
    here and each q is certified recursively, so a True answer is a proof
    and no probabilistic test is trusted."""
    if n in _LUCAS_MEMO:
        return _LUCAS_MEMO[n]
    if n < 2:
        return False
    if n < 10 ** 7:
        res = is_prime_trial(n)
        _LUCAS_MEMO[n] = res
        return res
    if pow(2, n - 1, n) != 1:
        # Fermat witness: n is composite.  Checked before factoring n-1 so
        # that a corrupted (composite) input is rejected immediately instead
        # of sending the factoriser after a large hard number.
        _LUCAS_MEMO[n] = False
        return False
    fac = factorize(n - 1)
    prod = 1
    for q, e in fac.items():
        prod *= q ** e
    if prod != n - 1:
        _LUCAS_MEMO[n] = False
        return False
    for q in fac:
        if not lucas_certified_prime(q):
            _LUCAS_MEMO[n] = False
            return False
    res = False
    for a in range(2, 500):
        if pow(a, n - 1, n) != 1:
            res = False
            break
        if all(pow(a, (n - 1) // q, n) != 1 for q in fac):
            res = True
            break
    _LUCAS_MEMO[n] = res
    return res


_PHI2 = {}


def phi_at_2(d):
    """Phi_d(2), the d-th cyclotomic polynomial evaluated at 2, computed
    exactly from  2^d - 1 = prod_{e | d} Phi_e(2)  by exact division.
    Raises if any division is inexact (which would signal a coding error)."""
    if d in _PHI2:
        return _PHI2[d]
    v = (1 << d) - 1
    for e in divisors(d):
        if e == d:
            continue
        w = phi_at_2(e)
        if v % w != 0:
            raise ArithmeticError("inexact cyclotomic division at d=%d" % d)
        v //= w
    _PHI2[d] = v
    return v


# --- a SECOND, structurally different implementation of Phi_d(2) -------------
# phi_at_2 divides big integers; the routine below builds the integer
# polynomial Phi_d(x) in Z[x] by long division of x^d - 1 and only then
# evaluates at x = 2 by Horner.  The two routes share no code path, so a
# coding error in either one is caught by their disagreement.  (The gates
# "prod_{d|n} Phi_d(2) == 2^n - 1" and "Phi_d(2) | 2^n - 1" are, on their own,
# algebraic identities that hold by construction of phi_at_2 and can never
# fire; this cross-check is what actually makes the cyclotomic route
# falsifiable.)

def _poly_exact_div(a, b):
    """Exact division in Z[x], coefficients low-degree first.  Raises if the
    quotient is not integral or the remainder is nonzero."""
    a = list(a)
    if not b or b[-1] == 0:
        raise ArithmeticError("zero divisor polynomial")
    out = [0] * max(0, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        if a[-1] == 0:
            a.pop()
            continue
        q, r = divmod(a[-1], b[-1])
        if r != 0:
            raise ArithmeticError("non-integral polynomial quotient")
        sh = len(a) - len(b)
        out[sh] = q
        for i, c in enumerate(b):
            a[sh + i] -= q * c
        while a and a[-1] == 0:
            a.pop()
    if any(a):
        raise ArithmeticError("nonzero polynomial remainder")
    return out


_CYCLO_POLY = {}


def cyclotomic_poly(d):
    """Phi_d(x) as a coefficient list, low-degree first, over Z."""
    if d in _CYCLO_POLY:
        return _CYCLO_POLY[d]
    num = [-1] + [0] * (d - 1) + [1]            # x^d - 1
    for e in divisors(d):
        if e == d:
            continue
        num = _poly_exact_div(num, cyclotomic_poly(e))
    _CYCLO_POLY[d] = num
    return num


def phi_at_2_by_polynomial(d):
    v = 0
    for c in reversed(cyclotomic_poly(d)):
        v = v * 2 + c
    return v


def euler_phi(n):
    r = n
    for p in prime_factors_trial(n):
        r = r // p * (p - 1)
    return r


def pocklington_row_ok(N, f_primes, a):
    """Execute the paper's Pocklington criterion (Lemma 4) for one row.
    f_primes is the list of
    primes whose product is F (taken from the paper's table).  Returns
    (ok, reason).  Every clause is computed."""
    if N <= 1 or N % 2 == 0:
        return False, "N not odd > 1"
    if len(set(f_primes)) != len(f_primes):
        return False, "F not squarefree"
    F = 1
    for q in f_primes:
        F *= q
    if (N - 1) % F != 0:
        return False, "F does not divide N-1"
    if F * F <= N:
        return False, "F^2 <= N"
    if pow(a, N - 1, N) != 1:
        return False, "a^(N-1) != 1 mod N"
    for q in f_primes:
        if math.gcd(pow(a, (N - 1) // q, N) - 1, N) != 1:
            return False, "gcd(a^((N-1)/q)-1, N) != 1 for q=%d" % q
    return True, "F=%d, F^2-N=%d" % (F, F * F - N)


# ----------------------------------------------------------------------
# INPUTS: everything below is transcribed from the paper
# ----------------------------------------------------------------------

PAPER_N = 170
PAPER_N_PRIMES = [2, 5, 17]                      # "170 = 2 . 5 . 17"
PAPER_OMEGA_N = 3                                # l = omega(170) = 3
PAPER_TAU_N = 8                                  # tau(170) = 8
PAPER_OMEGA_MERSENNE = 7                         # omega(2^170 - 1) = 7

PAPER_FACTORS = [3, 11, 31, 43691, 131071,
                 9520972806333758431, 26831423036065352611]   # equation (1)
PAPER_SMALL_TRIAL = [3, 11, 31, 43691, 131071]   # "prime by trial division"
PAPER_P = 9520972806333758431
PAPER_Q = 26831423036065352611
PAPER_BASE_PRIMES = [709, 2437, 2857, 6367, 27953]

# Certificate table of Proposition 5: (N, [primes of F], witness a)
PAPER_POCK_P = [(2598660833, [709, 2437], 3),
                (PAPER_P, [27953, 2598660833], 15)]
PAPER_POCK_Q = [(114281, [2857], 6),
                (10285291, [114281], 2),
                (204710635813423, [6367, 10285291], 3),
                (PAPER_Q, [204710635813423], 3)]

# Minimality proof
PAPER_SQFREE_LIST = [30, 42, 66, 70, 78, 102, 105,
                     110, 114, 130, 138, 154, 165]
PAPER_CENSUS = [70, 130, 154, 165]
PAPER_EXCL_ORDERS = {105: (7, 3), 110: (11, 10)}   # n -> (prime p, ord_p(2))
# (n, d, claimed factorisation of Phi_d(2)) from the displayed table
PAPER_TABLE_FACS = [(70, 35, [71, 122921]),
                    (130, 130, [131, 409891, 7623851]),
                    (154, 11, [23, 89]),
                    (165, 11, [23, 89])]


# ----------------------------------------------------------------------
# DERIVED: helpers used by the checks
# ----------------------------------------------------------------------

def certified_prime(m):
    """Rigorous primality: trial division below 10^7, Lucas n-1 above."""
    if m < 10 ** 7:
        return is_prime_trial(m)
    return lucas_certified_prime(m)


def order_mod_prime(a, r):
    """Multiplicative order of a modulo the prime r, by reducing the exponent
    r-1 through the prime factors of r-1.  Exact; no iteration bound.
    Returns 0 if the factorisation of r-1 cannot be verified, so that a
    caller can never be handed an order resting on a probabilistic test."""
    fac = factorize(r - 1)
    prod = 1
    for q, e in fac.items():
        prod *= q ** e
    if prod != r - 1:
        return 0
    for q in fac:
        if not certified_prime(q):
            return 0
    e = r - 1
    for q in fac:
        while e % q == 0 and pow(a, e // q, r) == 1:
            e //= q
    return e


def complete_prime_factorisation_ok(factors, M):
    """True iff `factors` is a COMPLETE prime factorisation of M: the entries
    are pairwise distinct, each is > 1 and carries a rigorous primality
    certificate, and their product is exactly M.  Used both by the C1 verdict
    and by the F2 self-test, so the self-test exercises the real predicate."""
    if len(set(factors)) != len(factors):
        return False
    if not all(f > 1 for f in factors):
        return False
    prod = 1
    for f in factors:
        prod *= f
    if prod != M:
        return False
    return all(certified_prime(f) for f in factors)


def mersenne_full_factorization(n):
    """Factor 2^n - 1 completely through its cyclotomic decomposition.
    Returns (ok, reason, prime_to_d) where prime_to_d maps each prime
    divisor found to the unique d | n with the prime dividing Phi_d(2)."""
    M = (1 << n) - 1
    prod = 1
    prime_to_d = {}
    for d in divisors(n):
        v = phi_at_2(d)
        # falsifiable cross-check: the big-integer route and the Z[x] route
        # must agree, and deg Phi_d must be Euler's phi(d).
        if v != phi_at_2_by_polynomial(d):
            return False, ("Phi_%d(2): integer route %d disagrees with the "
                           "polynomial route %d"
                           % (d, v, phi_at_2_by_polynomial(d))), {}
        cp = cyclotomic_poly(d)
        if len(cp) - 1 != euler_phi(d) or cp[-1] != 1:
            return False, "Phi_%d is not monic of degree phi(%d)" % (d, d), {}
        prod *= v
        if d == 1:
            continue
        fac = factorize(v)
        pv = 1
        for p, e in fac.items():
            pv *= p ** e
        if pv != v:
            return False, "factorisation of Phi_%d(2) is not exact" % d, {}
        for p in fac:
            if not certified_prime(p):
                return False, "uncertified factor %d of Phi_%d(2)" % (p, d), {}
            if M % p != 0:
                return False, "%d does not divide 2^%d-1" % (p, n), {}
            if p in prime_to_d:
                return False, "prime %d shared by d=%d and d=%d" % (
                    p, prime_to_d[p], d), {}
            prime_to_d[p] = d
    if prod != M:
        return False, "prod_{d|n} Phi_d(2) != 2^%d - 1" % n, {}
    return True, "omega=%d" % len(prime_to_d), prime_to_d


# ----------------------------------------------------------------------
# A. the exhibited object is well formed and is what the paper says
# ----------------------------------------------------------------------

def section_A():
    n = PAPER_N
    M = (1 << n) - 1
    print("# n = %d, 2^n - 1 = %d (%d decimal digits)"
          % (n, M, len(str(M))))

    prod = 1
    for p in PAPER_N_PRIMES:
        prod *= p
    check("A1.n_equals_product_of_its_primes",
          prod == n and all(is_prime_trial(p) for p in PAPER_N_PRIMES)
          and len(set(PAPER_N_PRIMES)) == len(PAPER_N_PRIMES),
          "%s -> %d" % (" * ".join(map(str, PAPER_N_PRIMES)), prod))

    computed = sorted(prime_factors_trial(n))
    check("A2.independent_factorisation_of_n_matches_paper",
          computed == sorted(PAPER_N_PRIMES)
          and omega_small(n) == PAPER_OMEGA_N,
          "trial division gives %s, omega=%d" % (computed, omega_small(n)))

    prod = 1
    for f in PAPER_FACTORS:
        prod *= f
    check("A3.product_of_paper_factors_equals_2^n-1", prod == M,
          "product has %d digits, 2^%d-1 has %d"
          % (len(str(prod)), n, len(str(M))))

    check("A4.paper_factors_distinct_and_gt_1",
          len(set(PAPER_FACTORS)) == len(PAPER_FACTORS)
          and all(f > 1 for f in PAPER_FACTORS)
          and list(PAPER_FACTORS) == sorted(PAPER_FACTORS),
          "%d factors, ascending" % len(PAPER_FACTORS))

    bad = [f for f in PAPER_SMALL_TRIAL if not is_prime_trial(f)]
    # also tie the paper's phrase "the first five factors" to the actual list
    # in (1), so the five certified numbers cannot drift from the five factors.
    aligned = list(PAPER_SMALL_TRIAL) == list(PAPER_FACTORS[:5])
    check("A5.five_small_factors_prime_by_trial_division",
          not bad and aligned,
          ("non-prime: %s" % bad if bad else
           "3,11,31,43691,131071 all prime")
          + ("" if aligned else "; NOT the first five entries of (1)"))

    bad = [p for p in PAPER_BASE_PRIMES if not is_prime_trial(p)]
    check("A6.five_certificate_base_primes_prime_by_trial_division", not bad,
          "non-prime: %s" % bad if bad else "709,2437,2857,6367,27953")


def run_chain(rows):
    """Execute a Pocklington chain.  A row is accepted only if every prime of
    its F is already proved prime (by trial division on the paper's base
    primes, or by an earlier row of the chain)."""
    proved = set(p for p in PAPER_BASE_PRIMES if is_prime_trial(p))
    details = []
    for (N, f_primes, a) in rows:
        unproved = [q for q in f_primes if q not in proved]
        if unproved:
            return False, "F contains unproved prime(s) %s for N=%d" % (
                unproved, N), proved
        ok, why = pocklington_row_ok(N, f_primes, a)
        if not ok:
            return False, "row N=%d: %s" % (N, why), proved
        details.append("%d ok (%s)" % (N, why))
        proved.add(N)
    return True, "; ".join(details), proved


# ----------------------------------------------------------------------
# B. primality of the two large factors, proved twice
# ----------------------------------------------------------------------

def section_B():
    okP, whyP, provedP = run_chain(PAPER_POCK_P)
    check("B1.pocklington_chain_proves_P", okP and PAPER_P in provedP, whyP)

    okQ, whyQ, provedQ = run_chain(PAPER_POCK_Q)
    check("B2.pocklington_chain_proves_Q", okQ and PAPER_Q in provedQ, whyQ)

    check("B3.chain_endpoints_are_the_two_large_factors_of_2^n-1",
          PAPER_POCK_P[-1][0] == PAPER_FACTORS[-2]
          and PAPER_POCK_Q[-1][0] == PAPER_FACTORS[-1],
          "P=%d, Q=%d" % (PAPER_P, PAPER_Q))

    # Independent proof, using no datum from the paper's table.
    lp = lucas_certified_prime(PAPER_P)
    lq = lucas_certified_prime(PAPER_Q)
    check("B4.independent_lucas_certificates_for_P_and_Q", lp and lq,
          "P prime=%s, Q prime=%s (n-1 factorisations computed here)"
          % (lp, lq))


# ----------------------------------------------------------------------
# C. omega(2^170 - 1) = 7, derived twice
# ----------------------------------------------------------------------

def section_C():
    n = PAPER_N
    M = (1 << n) - 1

    complete = complete_prime_factorisation_ok(PAPER_FACTORS, M)
    omega_M = len(set(PAPER_FACTORS)) if complete else -1
    check("C1.omega_2^n-1_is_exactly_7",
          complete and omega_M == PAPER_OMEGA_MERSENNE,
          "complete factorisation into %d certified distinct primes"
          % omega_M)

    ok, why, prime_to_d = mersenne_full_factorization(n)
    check("C2.independent_refactorisation_reproduces_the_seven_primes",
          ok and sorted(prime_to_d) == sorted(PAPER_FACTORS),
          "cyclotomic sieve gives %s (%s)" % (sorted(prime_to_d), why))

    vals = {}
    for d in divisors(n):
        vals[d] = phi_at_2(d)
    print("# Phi_d(2) for d | %d: %s"
          % (n, ", ".join("%d:%d" % (d, vals[d]) for d in sorted(vals))))
    nontrivial = sorted(vals[d] for d in vals if d > 1)
    check("C3.nontrivial_cyclotomic_values_are_the_paper_factors",
          nontrivial == sorted(PAPER_FACTORS) and vals[1] == 1,
          "each Phi_d(2), d>1, is prime here so the 7 factors are the 7 "
          "cyclotomic values")

    # "Phi_d(2) | 2^n-1" alone is an algebraic identity that holds by
    # construction of phi_at_2 and can never fail; the falsifiable content is
    # the agreement of the two independent implementations, plus the
    # structural facts deg Phi_d = phi(d) and Phi_d monic.
    bad = [d for d in vals if M % vals[d] != 0]
    disagree = [d for d in vals if vals[d] != phi_at_2_by_polynomial(d)]
    misshapen = [d for d in vals
                 if len(cyclotomic_poly(d)) - 1 != euler_phi(d)
                 or cyclotomic_poly(d)[-1] != 1]
    check("C4.cyclotomic_values_reconstructed_in_Z[x]_and_divide_2^n-1",
          not bad and not disagree and not misshapen,
          "non-divisors: %s; route disagreements: %s; wrong shape: %s"
          % (bad, disagree, misshapen) if (bad or disagree or misshapen)
          else "all %d values divide, and agree with Phi_d(x) built by long "
               "division in Z[x] (monic, deg = phi(d))" % len(vals))


# ----------------------------------------------------------------------
# D. hypotheses of Conjecture 4.3 hold at n = 170, its conclusion fails
# ----------------------------------------------------------------------

def section_D():
    n = PAPER_N
    M = (1 << n) - 1
    g = math.gcd(M, n)
    check("D1.gcd_2^n-1_n_equals_1", g == 1, "gcd = %d" % g)

    o5 = multiplicative_order(2, 5)
    o17 = multiplicative_order(2, 17)
    mech = (M % 2 == 1 and o5 == 4 and o17 == 8
            and n % o5 != 0 and n % o17 != 0
            and M % 5 != 0 and M % 17 != 0)
    check("D2.gcd_mechanism_orders_of_2_mod_5_and_17", mech,
          "2^n-1 odd; ord_5(2)=%d, ord_17(2)=%d, neither divides %d; "
          "5 and 17 do not divide 2^n-1" % (o5, o17, n))

    ell = omega_small(n)
    check("D3.hypothesis_omega_n_at_least_3", ell >= 3 and ell == PAPER_OMEGA_N,
          "omega(%d) = %d" % (n, ell))

    # omega is taken from the factorisation DERIVED here, not from the
    # paper's factor list, so this check stands on its own.
    fok, fwhy, fprimes = mersenne_full_factorization(n)
    omega_M = len(fprimes) if fok else -1
    bound = 2 ** ell
    # The detail string must not assert the refutation when the underlying
    # factorisation did not complete: a FAIL line reading "Conjecture 4.3 is
    # false" would be worse than no line at all.
    check("D4.CONCLUSION_FAILS_omega_2^n-1_lt_2^omega_n",
          fok and 0 < omega_M < bound,
          ("omega(2^%d-1) = %d < %d = 2^%d  ==> Conjecture 4.3 is false"
           % (n, omega_M, bound, ell)) if fok and 0 < omega_M < bound else
          ("NO CONCLUSION: factorisation of 2^%d-1 not established (%s)"
           % (n, fwhy)) if not fok else
          ("omega(2^%d-1) = %d is NOT < %d = 2^%d" % (n, omega_M, bound, ell)))

    t = tau(n)
    tight = fok and t == PAPER_TAU_N and omega_M == t - 1
    check("D5.tau_n_and_lemma_2_bound_are_tight", tight,
          ("tau(%d) = %d, omega(2^n-1) = %d = tau(n) - 1" % (n, t, omega_M))
          if tight else
          ("tau(%d) = %d (paper says %d), omega(2^n-1) = %s"
           % (n, t, PAPER_TAU_N, omega_M if fok else "not established")))

    ok, why, prime_to_d = mersenne_full_factorization(n)
    orders = {}
    good = ok
    for p, d in sorted(prime_to_d.items()):
        o = order_mod_prime(2, p)
        orders[p] = o
        if o != d or n % p == 0:
            good = False
    check("D6.lemma_mechanism_every_prime_of_Phi_d(2)_has_order_d", good,
          "d -> ord_r(2) for the prime r | Phi_d(2): "
          + ", ".join("%d->%d" % (prime_to_d[p], orders[p])
                      for p in sorted(orders, key=lambda q: prime_to_d[q])))


# ----------------------------------------------------------------------
# E. minimality: the exhaustive census of n < 170
# ----------------------------------------------------------------------

def census_below(limit):
    """All n < limit with omega(n) >= 3 and gcd(2^n - 1, n) = 1, computed by
    exhaustive enumeration and exact big-integer gcd."""
    out = []
    for n in range(1, limit):
        if omega_small(n) < 3:
            continue
        if math.gcd((1 << n) - 1, n) == 1:
            out.append(n)
    return out


def section_E1():
    n0 = PAPER_N
    cand = [n for n in range(1, n0) if omega_small(n) >= 3]
    cen = census_below(n0)
    print("# %d integers below %d have omega >= 3; %d of them satisfy the gcd "
          "hypothesis" % (len(cand), n0, len(cen)))
    check("E1.exhaustive_census_of_hypothesis_satisfiers_below_n",
          cen == sorted(PAPER_CENSUS),
          "computed %s, paper claims %s" % (cen, sorted(PAPER_CENSUS)))

    sq3 = [n for n in range(1, n0)
           if omega_small(n) == 3 and is_squarefree_small(n)]
    check("E2.squarefree_three_prime_products_below_n_match_paper_list",
          sq3 == sorted(PAPER_SQFREE_LIST),
          "%d values, computed == paper: %s" % (len(sq3),
                                                sq3 == sorted(PAPER_SQFREE_LIST)))

    big = [n for n in range(1, n0) if omega_small(n) >= 4]
    check("E3.no_integer_below_n_has_four_distinct_prime_factors",
          not big and 2 * 3 * 5 * 7 > n0,
          "found %s; 2*3*5*7 = %d > %d so l = 3 is forced"
          % (big, 2 * 3 * 5 * 7, n0))

    nonsq = [m for m in range(1, n0)
             if omega_small(m) >= 3 and not is_squarefree_small(m)]
    bad = [m for m in nonsq
           if not (tau(m) >= 3 * 2 ** (omega_small(m) - 1)
                   and tau(m) - 1 >= 2 ** omega_small(m))]
    ineq = all(3 * 2 ** (l - 1) - 1 >= 2 ** l for l in range(1, 65))
    check("E4.non_squarefree_case_is_killed_by_the_tau_bound",
          not bad and ineq and len(nonsq) > 0,
          "%d non-squarefree candidates, all with tau(n)-1 >= 2^omega(n); "
          "3*2^(l-1)-1 >= 2^l for l <= 64" % len(nonsq))


def section_E2():
    # The paper's stated exclusion reasons, recomputed one value at a time.
    # The list swept here is the one COMPUTED in E2, not the paper's transcript,
    # so E5/E7 are not paper-against-paper bookkeeping.
    sq3 = [m for m in range(1, PAPER_N)
           if omega_small(m) == 3 and is_squarefree_small(m)]
    six = [m for m in sq3 if m % 6 == 0]
    bad = [m for m in six
           if not (((1 << m) - 1) % 3 == 0 and m % 3 == 0
                   and math.gcd((1 << m) - 1, m) % 3 == 0)]
    check("E5.every_multiple_of_6_fails_the_gcd_hypothesis_via_3",
          not bad and len(six) == 7 and len(sq3) == 13,
          "%d of the %d computed values are divisible by 6; 3 | gcd for each: "
          "%s" % (len(six), len(sq3), six))

    reasons = []
    ok = True
    for m, (p, o) in sorted(PAPER_EXCL_ORDERS.items()):
        oc = multiplicative_order(2, p)
        good = (oc == o and m % oc == 0 and m % p == 0
                and ((1 << m) - 1) % p == 0
                and math.gcd((1 << m) - 1, m) % p == 0)
        ok = ok and good
        reasons.append("n=%d: ord_%d(2)=%d divides n, so %d | gcd"
                       % (m, p, oc, p))
    check("E6.the_values_105_and_110_fail_via_orders_of_2_mod_7_and_11",
          ok, "; ".join(reasons))

    rest = sorted(set(sq3) - set(six) - set(PAPER_EXCL_ORDERS))
    check("E7.the_stated_exclusions_leave_exactly_the_paper_census",
          rest == sorted(PAPER_CENSUS) and rest == census_below(PAPER_N),
          "%d - %d (mult. of 6) - 2 (105,110) = %s, and this equals the "
          "independently computed census"
          % (len(sq3), len(six), rest))


def section_E3():
    # the load-bearing minimality computation
    cen = census_below(PAPER_N)
    all_ok = True
    lemma_ok = True
    for m in cen:
        ok, why, prime_to_d = mersenne_full_factorization(m)
        om = len(prime_to_d) if ok else -1
        need = 2 ** omega_small(m)
        good = ok and om >= need
        all_ok = all_ok and good
        lemma_ok = lemma_ok and ok and om >= tau(m) - 1
        print("#   n=%d: omega(n)=%d, tau(n)=%d, omega(2^n-1)=%d, need >= %d"
              " -> %s" % (m, omega_small(m), tau(m), om, need,
                          "no counterexample" if good else "COUNTEREXAMPLE?"))
    check("E8.MINIMALITY_no_n_below_170_violates_the_conjectured_bound",
          all_ok and len(cen) == 4,
          "all %d hypothesis satisfiers below %d have omega(2^n-1) >= 2^3 = 8"
          % (len(cen), PAPER_N))
    check("E9.lemma_2_bound_holds_on_every_census_member",
          lemma_ok, "omega(2^n-1) >= tau(n) - 1 for n in %s" % cen)

    ok = True
    rows = []
    for (m, d, facs) in PAPER_TABLE_FACS:
        v = phi_at_2(d)
        prod = 1
        for f in facs:
            prod *= f
        good = (m % d == 0 and prod == v and len(set(facs)) == len(facs)
                and all(is_prime_trial(f) for f in facs) and len(facs) >= 2)
        ok = ok and good
        rows.append("n=%d: Phi_%d(2)=%d=%s" % (m, d, v,
                                               "*".join(map(str, facs))))
    check("E10.displayed_cyclotomic_factorisations_of_the_table_are_correct",
          ok, "; ".join(rows))


# ----------------------------------------------------------------------
# F. the theorem, recomputed from scratch, and a mutation self-test
# ----------------------------------------------------------------------

def is_counterexample(m):
    """Decide whether m satisfies the hypotheses of Conjecture 4.3 and violates
    its conclusion.  Everything is computed: omega(m), the gcd, and the
    complete factorisation of 2^m - 1.

    Returns (hit, why, decided).  `decided` is False exactly when the verdict
    could NOT be established -- i.e. the factoriser failed to produce a
    complete, certified factorisation of 2^m - 1.  An undecided m must never
    be silently read as "not a counterexample": that would let the minimality
    scan return 170 while knowing nothing about a smaller exponent.  Callers
    are required to propagate `decided`."""
    if omega_small(m) < 3:
        return False, "omega(n) < 3", True
    if math.gcd((1 << m) - 1, m) != 1:
        return False, "gcd != 1", True
    ok, why, prime_to_d = mersenne_full_factorization(m)
    if not ok:
        return False, "UNDECIDED, factorisation failed: " + why, False
    om = len(prime_to_d)
    need = 2 ** omega_small(m)
    return (om < need,
            "omega(2^n-1)=%d vs 2^omega(n)=%d" % (om, need), True)


def section_F():
    least = None
    undecided = []
    for m in range(1, PAPER_N + 1):
        hit, why, decided = is_counterexample(m)
        if not decided:
            undecided.append(m)
            continue
        if hit:
            least = m
            break
    # Fail closed: an exponent the program could not decide invalidates the
    # scan, even if the first DECIDED counterexample happens to be 170.
    check("F1.THEOREM_least_counterexample_recomputed_is_170",
          least == PAPER_N and not undecided,
          "smallest counterexample found by scanning n = 1..%d is %s; "
          "undecided exponents: %s"
          % (PAPER_N, least, undecided if undecided else "none"))

    # --- mutation self-tests: the load-bearing predicates must reject junk ---
    M = (1 << PAPER_N) - 1
    perturbed = list(PAPER_FACTORS)
    perturbed[3] = perturbed[3] + 2       # 43691 -> 43693
    # A product-PRESERVING corruption: merge 3 * 11 into the composite 33.
    # The naive "does the product still equal 2^170-1" test cannot see this
    # one, so it is put through the same completeness predicate that C1 uses.
    merged = sorted([33] + [f for f in PAPER_FACTORS if f not in (3, 11)])
    merged_prod = 1
    for f in merged:
        merged_prod *= f
    check("F2.selftest_corrupted_factorisations_are_rejected",
          (not complete_prime_factorisation_ok(perturbed, M))
          and merged_prod == M
          and (not complete_prime_factorisation_ok(merged, M))
          and complete_prime_factorisation_ok(PAPER_FACTORS, M),
          "43691 -> 43693 rejected; the product-preserving merge 3*11 -> 33 "
          "is also rejected (same product, not a prime factorisation); the "
          "true list is accepted")

    # Both large factors are perturbed, not just P: each of P+2 and Q+2 is put
    # through the final row of its own Pocklington chain and through the
    # independent Lucas certifier, and both must reject it, while the true P
    # and Q must still be accepted (so the mutation test is not vacuously
    # rejecting everything).  Q+2 is composite: the digit sum of Q is 67, so
    # Q = 1 mod 3 and 3 | Q+2; and F = 204710635813423 divides Q-1, hence
    # cannot divide (Q+2)-1 = Q+1, since that would force F | 2.
    badP = PAPER_P + 2
    badQ = PAPER_Q + 2
    rowP_ok, whyP = pocklington_row_ok(badP, PAPER_POCK_P[-1][1],
                                       PAPER_POCK_P[-1][2])
    rowQ_ok, whyQ = pocklington_row_ok(badQ, PAPER_POCK_Q[-1][1],
                                       PAPER_POCK_Q[-1][2])
    check("F3.selftest_corrupted_P_and_Q_rejected_by_both_primality_proofs",
          (not rowP_ok) and (not lucas_certified_prime(badP))
          and (not rowQ_ok) and (not lucas_certified_prime(badQ))
          and lucas_certified_prime(PAPER_P)
          and lucas_certified_prime(PAPER_Q),
          "P+2 rejected (%s) and Q+2 rejected (%s) by their Pocklington rows, "
          "and both by the Lucas certifier; the true P and Q still certify "
          "prime" % (whyP, whyQ))

    hit70, why70, dec70 = is_counterexample(70)
    hit30, why30, dec30 = is_counterexample(30)
    # `decided` is required: "not a counterexample" must be a verdict, not an
    # admission that the factoriser gave up.
    check("F4.selftest_violation_predicate_is_false_on_non_counterexamples",
          dec70 and dec30 and (not hit70) and (not hit30),
          "n=70 (%s), n=30 (%s)" % (why70, why30))

    composites = [(1 << PAPER_N) - 1, 409368176241571, 3215031751,
                  PAPER_P * 3]
    check("F5.selftest_certifier_rejects_known_composites",
          all(not certified_prime(c) for c in composites),
          "rejected 2^170-1, Phi_130(2), the strong pseudoprime 3215031751, "
          "and 3P")


LEMMA_SWEEP_LIMIT = 100


def section_G():
    """Supplementary (beyond the paper's claim): Lemma 2 is a proved
    statement, not a computation, but the minimality argument leans on it, so
    it is tested on every exponent up to LEMMA_SWEEP_LIMIT whose gcd
    hypothesis holds."""
    bad = []
    tested = 0
    for m in range(1, LEMMA_SWEEP_LIMIT + 1):
        if math.gcd((1 << m) - 1, m) != 1:
            continue
        ok, why, prime_to_d = mersenne_full_factorization(m)
        tested += 1
        if not ok or len(prime_to_d) < tau(m) - 1:
            bad.append((m, ok, len(prime_to_d), tau(m) - 1))
    check("G1.lemma_2_holds_for_every_exponent_up_to_100_with_gcd_1",
          tested > 0 and not bad,
          "omega(2^n-1) >= tau(n)-1 verified for all %d such n <= %d; "
          "violations: %s" % (tested, LEMMA_SWEEP_LIMIT, bad))

    hyp = [m for m in range(PAPER_N + 1, 1001)
           if omega_small(m) >= 3 and math.gcd((1 << m) - 1, m) == 1]
    print("# supplementary: %d integers in (170, 1000] also satisfy the two "
          "hypotheses (%s ...); the paper makes no claim about them and they "
          "are not factored here" % (len(hyp), hyp[:6]))


def main():
    t0 = time.time()
    print("# verify.py -- Conjecture 4.3 of Dolfi-Hafezieh-Spiga, "
          "least counterexample n = 170")
    print("# python %s, standard library only, exact integer arithmetic"
          % sys.version.split()[0])
    print("# ---- A. the exhibited object ----")
    section_A()
    print("# ---- B. primality of the two large factors ----")
    section_B()
    print("# ---- C. omega(2^170-1) = 7 ----")
    section_C()
    print("# ---- D. hypotheses hold, conclusion fails ----")
    section_D()
    print("# ---- E. minimality census over 1 <= n < 170 ----")
    section_E1()
    section_E2()
    section_E3()
    print("# ---- F. theorem recomputed, and mutation self-tests ----")
    section_F()
    print("# ---- G. supplementary, beyond the paper's claim ----")
    section_G()
    print("# SCOPE: the census in E and F is the paper's COMPLETE minimality")
    print("#   claim: all 169 integers below 170 are enumerated, the 23 with")
    print("#   omega >= 3 are tested against the gcd hypothesis, and each of")
    print("#   the 4 survivors has 2^n-1 factored completely.  Nothing in the")
    print("#   paper's claim is left unverified or narrowed.")
    print("# NOT re-run: (i) the supplementary sweep G1 stops at n = 100,")
    print("#   because factoring 2^n-1 for prime exponents near 170 (e.g.")
    print("#   n = 167) is far outside the budget -- G1 is an extra, not part")
    print("#   of the paper's claim, whose Lemma 2 is proved, not computed;")
    print("#   (ii) the bibliographic attributions in Remark 3 ([LLRS]")
    print("#   Example 3.8, OEIS A046800) are not checkable offline.")
    print("# elapsed %.2f s" % (time.time() - t0))
    return verdict()


if __name__ == "__main__":
    sys.exit(main())
