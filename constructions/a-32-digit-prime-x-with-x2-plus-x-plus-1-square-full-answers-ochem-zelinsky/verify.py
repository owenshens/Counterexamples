#!/usr/bin/env python3
# verify.py -- referee check of the note
#
#   A 32-digit prime x with x^2+x+1 square-full: Question sqfull of
#   Ochem and Zelinsky has the answer yes
#
# Everything below is read from the DECIMAL STRINGS printed in the paper and
# re-derived from them. Python 3.9+, standard library only (no numpy, no sympy,
# no gmpy2, no external data file, no network). Exact integer arithmetic only;
# no floating-point value ever decides a check. No randomness anywhere: the two
# primality routines are trial division and a Pratt (n-1) certificate whose
# factorisations and bases are transcribed from the paper, so the program has no
# search of its own to get lucky in.
#
# Run:   python3 verify.py
# It prints one `PASS <name> [detail]` line per check, some `NOTE` lines, and
# closes with `VERDICT: ALL <n> CHECKS PASS`, exiting 0 iff every check passed.

import sys
from array import array
from math import isqrt

# ---------------------------------------------------------------------------
# 0. THE OBJECT, exactly as printed in the paper
# ---------------------------------------------------------------------------
X = 47116128896261596331524418282573                    # the witness (32 digits)
N_PRINTED = 2219929602169136991765718912586970986546248112210754339293782903
A = 2544031832644333879196802002161                     # the square part of n/7^3
N_FACTORS = [(7, 3), (43, 2), (3416683, 2), (12353161, 2), (1401752363863729, 2)]
A_FACTORS = [43, 3416683, 12353161, 1401752363863729]

# x-1, factored, as printed in the paper
X_MINUS_1 = [(2, 2), (6299, 1), (9701563, 1), (192750846182877882539, 1)]

# The residues the paper prints for the Lucas certificate of x, base 2.
X_RESIDUES_PRINTED = {
    2:                     47116128896261596331524418282572,
    6299:                  42436530233755792060895430447505,
    9701563:               21277589209484664623743014586336,
    192750846182877882539: 19695209395513058587964373225481,
}

# The two subsidiary certificates printed in the paper, so that NO primality
# claim anywhere rests on a probable-prime test or on a published bound.
PRATT = {
    X:                     (2, X_MINUS_1),
    192750846182877882539: (2, [(2, 1), (17, 2), (27779, 1), (12004714807399, 1)]),
    1401752363863729:      (7, [(2, 4), (3, 1), (37, 1), (107, 1), (7376401679, 1)]),
}

# The Pell data printed in the paper
U_PRINTED = 94232257792523192663048836565147               # u = 2x+1
D = 1372                                                   # 4 * 7^3
T, S = 34100354867927167, 920623046934552                   # T^2 - 1372 S^2 = 1
FUND = (8, 3)                                              # the unit 8 + 3 sqrt 7
FUND_POWER = 14                                            # lcm(2,7)

# The two composite near-misses and the census target. NEITHER IS AN EXAMPLE IN
# THE SOURCE PAPER: x=18 is a(2) of OEIS A296376 and x=88916 occurs in no cited
# source at all -- it is a hit of the census in Step 7 below.
NEAR_MISSES = [18, 88916]
CENSUS_XMAX = 10 ** 8
ANTIPRUNE_XMAX = 10 ** 6
BRUTE_XMAX = 3000
# A118896: the count of powerful numbers <= 10^k, k = 1..8
A118896 = [4, 14, 54, 185, 619, 2027, 6553, 21044]

# ---------------------------------------------------------------------------
# 1. HARNESS
# ---------------------------------------------------------------------------
_state = {'pass': 0, 'fail': 0}


def chk(name, cond, detail=''):
    if cond:
        _state['pass'] += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _state['fail'] += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(s):
    print('NOTE %s' % s)


def head(s):
    print()
    print('=== %s' % s)


# ---------------------------------------------------------------------------
# 2. PRIMALITY: trial division, then Pratt (n-1) certificates
# ---------------------------------------------------------------------------
SIEVE_LIMIT = 4 * 10 ** 6          # >= sqrt(12004714807399), the largest trial case


def small_primes(limit):
    sieve = bytearray(b'\x01') * (limit + 1)
    sieve[0:2] = b'\x00\x00'
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p::p] = bytearray(len(sieve[p * p::p]))
    return [i for i in range(limit + 1) if sieve[i]]


PRIMES = small_primes(SIEVE_LIMIT)
TRIAL_MAX = SIEVE_LIMIT * SIEVE_LIMIT


def prime_by_trial(n):
    """Exhaustive trial division. Complete only for n <= TRIAL_MAX = 1.6e13."""
    assert 0 < n <= TRIAL_MAX, n
    if n < 2:
        return False
    for p in PRIMES:
        if p * p > n:
            return True
        if n % p == 0:
            return n == p
    return True


def prove_prime(n, trace=None):
    """A PROOF, never a probable-prime test.

    n <= TRIAL_MAX  -> exhaustive trial division.
    otherwise       -> the Pratt certificate transcribed from the paper: n-1 is
                       re-multiplied from its printed factorisation, every
                       printed factor is itself proved prime by this same
                       routine, and the base a from the paper is checked to
                       satisfy a^(n-1) = 1 and a^((n-1)/q) != 1 (mod n) for each
                       prime q | n-1. That is the converse of Fermat's theorem
                       (Lucas), so it certifies primality outright.
    """
    if trace is None:
        trace = []
    if n <= TRIAL_MAX:
        r = prime_by_trial(n)
        trace.append(('trial', n, r))
        return r
    if n not in PRATT:
        trace.append(('no-certificate', n, False))
        return False
    a, fac = PRATT[n]
    prod = 1
    for q, e in fac:
        prod *= q ** e
    if prod != n - 1:
        trace.append(('bad-factorisation', n, False))
        return False
    if pow(a, n - 1, n) != 1:
        trace.append(('fermat-fails', n, False))
        return False
    for q, _ in fac:
        if pow(a, (n - 1) // q, n) == 1:
            trace.append(('order-too-small', (n, q), False))
            return False
    for q, _ in fac:
        if not prove_prime(q, trace):
            trace.append(('factor-not-prime', (n, q), False))
            return False
    trace.append(('pratt', n, True))
    return True


# ---------------------------------------------------------------------------
# 3. SQUARE-FULLNESS, from first principles
# ---------------------------------------------------------------------------
def factor_small(n):
    """Complete factorisation by trial division; used only for n <= 10^7 or so."""
    out = []
    m = n
    for p in PRIMES:
        if p * p > m:
            break
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            out.append((p, e))
    if m > 1:
        out.append((m, 1))
    return out


def square_full_small(n):
    """(rad n)^2 | n, decided from a complete factorisation. n >= 1; 1 counts."""
    if n < 1:
        return False
    return all(e >= 2 for _, e in factor_small(n))


def perfect_power(n):
    """True iff n = r^t for some integers r >= 2, t >= 2. Exact, no floats."""
    if n < 4:
        return False
    for t in range(2, n.bit_length() + 1):
        lo, hi = 1, 1 << ((n.bit_length() + t) // t + 1)
        while lo < hi:
            mid = (lo + hi) // 2
            if mid ** t < n:
                lo = mid + 1
            else:
                hi = mid
        if lo ** t == n and lo >= 2:
            return True
    return False


# ---------------------------------------------------------------------------
# 4. THE CENSUS: every x <= XMAX with x^2+x+1 square-full
# ---------------------------------------------------------------------------
def icbrt(n):
    """The integer cube root of n >= 0, by bisection. No floats."""
    lo, hi = 0, 1
    while hi ** 3 <= n:
        hi *= 2
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** 3 <= n:
            lo = mid
        else:
            hi = mid - 1
    return lo


def sieve_flags(limit):
    """(admissible, squarefree) bytearrays over 0..limit.

    admissible[m] = 1 iff every prime factor of m is 1 mod 3 (so m is odd and 3
    does not divide m); m = 1 qualifies vacuously. Both are bytearrays, so the
    census uses about 3 bytes per integer of the range and stays well inside a
    hundred megabytes at the bounds used here.
    """
    isp = bytearray(b'\x01') * (limit + 1)
    isp[0:2] = b'\x00\x00'
    for p in range(2, isqrt(limit) + 1):
        if isp[p]:
            isp[p * p::p] = bytearray(len(range(p * p, limit + 1, p)))
    adm = bytearray(b'\x01') * (limit + 1)
    adm[0] = 0
    for p in range(2, limit + 1):
        if isp[p] and p % 3 != 1:
            adm[p::p] = bytearray(len(range(p, limit + 1, p)))
    sqfree = bytearray(b'\x01') * (limit + 1)
    sqfree[0] = 0
    for p in range(2, isqrt(limit) + 1):
        if isp[p]:
            pp = p * p
            sqfree[pp::pp] = bytearray(len(range(pp, limit + 1, pp)))
    return adm, sqfree


def census(xmax, restrict=True):
    """All x in [1, xmax] with x^2+x+1 square-full, as (x, a, b) with n = a^2 b^3.

    Every powerful n has a unique representation n = a^2 b^3 with b squarefree
    (Golomb), and n = x^2+x+1 iff 4n-3 = (2x+1)^2, so the search is over pairs
    (a, b) with a^2 b^3 <= nmax and 4 a^2 b^3 - 3 a perfect square.

    restrict=True applies the sound prune: if x^2+x+1 is square-full then it is
    odd, 3 does not divide it (v_3 <= 1 always), and every prime factor is
    1 mod 3 -- so a and b are both admissible in the sense above, and the least
    usable b is 7.
    restrict=False drops the prune entirely and is the anti-prune control; there
    the least usable b is 2, so the a-range is larger.
    """
    nmax = xmax * xmax + xmax + 1
    bmin = 7 if restrict else 2
    bmax = icbrt(nmax)
    amax_all = isqrt(nmax // bmin ** 3) + 1
    limit = max(amax_all, bmax) + 1
    adm, sqfree = sieve_flags(limit)
    if restrict:
        avals = array('q', (m for m in range(1, amax_all + 1) if adm[m]))
        bvals = [b for b in range(bmin, bmax + 1) if adm[b] and sqfree[b]]
    else:
        avals = array('q', range(1, amax_all + 1))
        bvals = [b for b in range(bmin, bmax + 1) if sqfree[b]]
    hits = []
    pairs = 0
    for b in bvals:
        b3 = b ** 3
        amax = isqrt(nmax // b3)
        for a in avals:
            if a > amax:
                break
            pairs += 1
            n = a * a * b3
            v = 4 * n - 3
            u = isqrt(v)
            if u * u == v and u % 2 == 1:
                x = (u - 1) // 2
                if 1 <= x <= xmax:
                    hits.append((x, a, b))
    return sorted(hits), pairs, len(avals), len(bvals)


def pell_fundamental(d):
    """The fundamental solution of t^2 - d s^2 = 1, by the continued-fraction
    expansion of sqrt(d). Exact integer arithmetic; d not a perfect square."""
    a0 = isqrt(d)
    assert a0 * a0 != d
    m, q, a = 0, 1, a0
    num_prev, num = 1, a0
    den_prev, den = 0, 1
    while num * num - d * den * den != 1:
        m = q * a - m
        q = (d - m * m) // q
        a = (a0 + m) // q
        num_prev, num = num, a * num + num_prev
        den_prev, den = den, a * den + den_prev
    return num, den


def powerful_count(limit):
    """#{ n <= limit : n powerful }, by the same a^2 b^3 enumeration. Includes 1."""
    amax = isqrt(limit) + 1
    squarefree = bytearray(b'\x01') * (amax + 1)
    for p in PRIMES:
        if p * p > amax:
            break
        for j in range(p * p, amax + 1, p * p):
            squarefree[j] = 0
    total = 0
    b = 1
    while b ** 3 <= limit:
        if squarefree[b]:
            total += isqrt(limit // b ** 3)
        b += 1
    return total


# ===========================================================================
#                                   CHECKS
# ===========================================================================
print('verify.py -- a 32-digit prime x with x^2+x+1 square-full;')
print('Question sqfull of Ochem and Zelinsky (arXiv:2607.19746) has the answer YES.')
print('exact integer arithmetic only, standard library only, no randomness')

# --- 1 -------------------------------------------------------------------
head('Step 1: the object read from the paper')
N = X * X + X + 1
note('x = %d' % X)
note('n = x^2+x+1 = %d' % N)
chk('x_has_32_decimal_digits', len(str(X)) == 32, '%d digits' % len(str(X)))
chk('n_equals_the_decimal_string_printed_in_the_paper', N == N_PRINTED)
chk('n_has_64_decimal_digits', len(str(N)) == 64, '%d digits' % len(str(N)))
chk('A_equals_the_product_of_the_four_factors_printed_for_it',
    A == 43 * 3416683 * 12353161 * 1401752363863729, 'A = %d' % A)
chk('n_equals_seven_cubed_times_A_squared', N == 7 ** 3 * A * A)
_prod = 1
for _p, _e in N_FACTORS:
    _prod *= _p ** _e
chk('n_equals_7^3_43^2_3416683^2_12353161^2_1401752363863729^2', N == _prod)
chk('seven_does_not_divide_A_so_v7_of_n_is_exactly_three', A % 7 == 3,
    'A mod 7 = %d' % (A % 7))

# --- 2 -------------------------------------------------------------------
head('Step 2: n is square-full, by the source\'s own definition (rad n)^2 | n')
chk('every_prime_dividing_n_other_than_7_divides_A',
    all(A % p == 0 for p, _ in N_FACTORS if p != 7),
    'so v_p(n) = 2 v_p(A) >= 2 for p != 7, and v_7(n) = 3')
for _p, _e in N_FACTORS:
    chk('listed_factor_%d_is_prime_with_exponent_%d_at_least_2' % (_p, _e),
        prove_prime(_p) and _e >= 2)
_rad = 1
for _p, _e in N_FACTORS:
    _rad *= _p
chk('rad_n_squared_divides_n', N % (_rad * _rad) == 0, 'rad n = %d' % _rad)
chk('n_is_square_full', all(e >= 2 for _, e in N_FACTORS))
chk('omega_of_n_is_5', len(N_FACTORS) == 5)
chk('n_is_odd_and_3_does_not_divide_n', N % 2 == 1 and N % 3 != 0)
chk('every_prime_factor_of_n_is_1_mod_3', all(p % 3 == 1 for p, _ in N_FACTORS),
    'residues ' + ','.join(str(p % 3) for p, _ in N_FACTORS))
chk('n_is_NOT_a_perfect_power', not perfect_power(N),
    'so the witness is not excluded by the source\'s Lemma rtt (Bennett-Levin), '
    'which closes the case x^2+x+1 = r^t with r prime')

# --- 3 -------------------------------------------------------------------
head('Step 3: x is prime -- Lucas (n-1) certificate, base 2, no probable-prime step')
_pm = 1
for _p, _e in X_MINUS_1:
    _pm *= _p ** _e
chk('x_minus_1_equals_2^2_6299_9701563_192750846182877882539', _pm == X - 1,
    'x-1 = %d' % (X - 1))
chk('two_to_the_x_minus_1_is_1_mod_x', pow(2, X - 1, X) == 1)
for _q, _ in X_MINUS_1:
    _r = pow(2, (X - 1) // _q, X)
    chk('two_to_the_x_minus_1_over_%d_is_not_1_mod_x' % _q, _r != 1, '= %d' % _r)
    chk('that_residue_equals_the_one_printed_in_the_paper_for_q_%d' % _q,
        _r == X_RESIDUES_PRINTED[_q])
for _q in (6299, 9701563):
    chk('certificate_prime_%d_is_prime_by_exhaustive_trial_division' % _q,
        prime_by_trial(_q), 'sqrt = %d' % isqrt(_q))
# the one certificate prime that itself needs a certificate
_a2, _f2 = PRATT[192750846182877882539]
_pm2 = 1
for _p, _e in _f2:
    _pm2 *= _p ** _e
chk('192750846182877882539_minus_1_equals_2_17^2_27779_12004714807399',
    _pm2 == 192750846182877882539 - 1)
chk('base_2_certifies_192750846182877882539',
    pow(_a2, 192750846182877882539 - 1, 192750846182877882539) == 1
    and all(pow(_a2, (192750846182877882539 - 1) // q, 192750846182877882539) != 1
            for q, _ in _f2))
for _q in (17, 27779, 12004714807399):
    chk('subsidiary_prime_%d_is_prime_by_exhaustive_trial_division' % _q,
        prime_by_trial(_q), 'sqrt = %d' % isqrt(_q))
_trace = []
chk('prove_prime_returns_TRUE_for_x_with_a_fully_grounded_certificate',
    prove_prime(X, _trace), '%d certificate steps, every leaf trial-divided'
    % len(_trace))
chk('no_prime_below_10^6_divides_x_independent_sanity_check',
    all(X % p != 0 for p in PRIMES if p < 10 ** 6))
note('x IS PRIME unconditionally: x-1 is completely factored into proven primes '
     'and base 2 has full order modulo x.')

# --- 4 -------------------------------------------------------------------
head('Step 4: the auxiliary full factorisation of n (not load-bearing)')
_a3, _f3 = PRATT[1401752363863729]
_pm3 = 1
for _p, _e in _f3:
    _pm3 *= _p ** _e
chk('1401752363863729_minus_1_equals_2^4_3_37_107_7376401679',
    _pm3 == 1401752363863729 - 1)
chk('base_7_certifies_1401752363863729_base_2_does_NOT',
    pow(7, 1401752363863729 - 1, 1401752363863729) == 1
    and all(pow(7, (1401752363863729 - 1) // q, 1401752363863729) != 1 for q, _ in _f3)
    and pow(2, (1401752363863729 - 1) // 2, 1401752363863729) == 1,
    'base 2 is a quadratic residue there, which is why the paper prints base 7')
chk('7376401679_is_prime_by_exhaustive_trial_division', prime_by_trial(7376401679),
    'sqrt = %d' % isqrt(7376401679))
chk('all_four_factors_of_A_are_prime', all(prove_prime(p) for p in A_FACTORS))
note('the square-fullness argument in Step 2 needs NONE of Step 4: it uses only '
     'n = 7^3 A^2 and 7 not dividing A.')

# --- 5 -------------------------------------------------------------------
head('Step 5: the norm-form family the witness comes from')
U = 2 * X + 1
chk('u_equals_2x_plus_1_as_printed', U == U_PRINTED, 'u = %d' % U)
chk('four_n_minus_3_equals_u_squared', 4 * N - 3 == U * U)
chk('u_squared_minus_1372_A_squared_equals_minus_3', U * U - D * A * A == -3)
chk('u_squared_minus_7_times_14A_squared_equals_minus_3',
    U * U - 7 * (14 * A) ** 2 == -3, 'the same equation over Z[sqrt 7]')
chk('T_squared_minus_1372_S_squared_equals_1', T * T - D * S * S == 1)
chk('T_S_is_the_FUNDAMENTAL_solution_of_t^2_minus_1372_s^2_equals_1',
    pell_fundamental(D) == (T, S),
    'computed independently by the continued-fraction expansion of sqrt(1372)')
_p, _q = 1, 0
for _ in range(FUND_POWER):
    _p, _q = FUND[0] * _p + 7 * FUND[1] * _q, FUND[1] * _p + FUND[0] * _q
chk('the_automorph_is_the_14th_power_of_the_unit_8_plus_3_sqrt_7',
    (_p, _q) == (T, 14 * S), '(8+3sqrt7)^14 = %d + %d sqrt7' % (_p, _q))
_least = None
_pp, _qq = 1, 0
for _k in range(1, 40):
    _pp, _qq = FUND[0] * _pp + 7 * FUND[1] * _qq, FUND[1] * _pp + FUND[0] * _qq
    if _pp * _pp - 7 * _qq * _qq == 1 and _qq % 14 == 0:
        _least = _k
        break
chk('14_is_the_LEAST_power_of_that_unit_with_norm_1_and_14_dividing_w',
    _least == FUND_POWER, 'least k = %s = lcm(2,7)' % _least)
chk('1372_is_1_mod_3_and_T_is_1_mod_3_and_S_is_0_mod_3',
    D % 3 == 1 and T % 3 == 1 and S % 3 == 0,
    'so (u,a) -> (uT+1372aS, uS+aT) is the identity map mod 3')


def step(u, a):
    return u * T + D * a * S, u * S + a * T


_u, _a = 37, 1
chk('the_seed_37_1_has_norm_minus_3', _u * _u - D * _a * _a == -3)
_bad = []
for _k in range(20):
    if _u <= 0 or _u % 3 != 1 or ((_u - 1) // 2) % 3 != 0 \
            or _u * _u - D * _a * _a != -3:
        _bad.append(_k)
    _u, _a = step(_u, _a)
chk('from_the_seed_37_1_every_u_is_positive_and_1_mod_3_and_3_divides_x',
    not _bad, '20 iterates, 0 exceptions; x=18 is the k=0 member, '
    'so that branch of the family contains no prime at all')
_u, _a = 37, -1
_chain = []
_signs = []
for _k in range(4):
    _signs.append(1 if _u > 0 else -1)
    _chain.append((abs(_u) - 1) // 2)
    _u, _a = step(_u, _a)
chk('from_the_conjugate_seed_37_minus_1_u_is_negative_from_the_second_term_on',
    _signs == [1, -1, -1, -1], 'signs of u: ' + ','.join(map(str, _signs)))
chk('there_x_equals_abs_u_minus_1_over_2_has_2_15_32_49_digits',
    [len(str(v)) for v in _chain] == [2, 15, 32, 49],
    'digit counts ' + ','.join(str(len(str(v))) for v in _chain))
chk('the_second_iterate_from_37_minus_1_IS_the_witness', _chain[2] == X)
chk('every_term_after_the_seed_there_has_x_equal_2_mod_3',
    all(v % 3 == 2 for v in _chain[1:]), 'x mod 3 = ' +
    ','.join(str(v % 3) for v in _chain[1:]))
note('so the first member of that branch not automatically divisible by 3 has '
     '15 digits (%d), which is why a search over x <= 10^9 could not find one.'
     % _chain[1])

# --- 6 -------------------------------------------------------------------
head('Step 6: controls on the decider, both polarities')
chk('control_positive_343_is_square_full', square_full_small(343), '343 = 7^3, '
    'the least square-full value of x^2+x+1, at the composite x=18 = a(2) of '
    'OEIS A296376; NOT an example recorded in the source paper')
chk('control_positive_7906143973_is_square_full', square_full_small(7906143973),
    '= 1897^2 * 13^3, the composite near-miss at x=88916; cube part b=13, so '
    'the b=7 Pell branch is not the whole family')
chk('control_positive_72_and_4_and_8_and_9_and_1_are_square_full',
    all(square_full_small(v) for v in (72, 4, 8, 9, 1)))
chk('control_negative_7_3_21_57_12_31_are_not_square_full',
    not any(square_full_small(v) for v in (7, 3, 21, 57, 12, 31)))
chk('anti_control_the_perfect_power_misreading_is_refuted_by_7906143973',
    square_full_small(7906143973) and not perfect_power(7906143973),
    'square-full but not a perfect power, so "square-full" is NOT "perfect power"')
chk('anti_control_x_18_is_square_full_but_18_is_NOT_prime',
    square_full_small(343) and not prime_by_trial(18))
chk('anti_control_x_88916_is_square_full_but_88916_is_NOT_prime',
    square_full_small(7906143973) and not prime_by_trial(88916))
_counts = [powerful_count(10 ** k) for k in range(1, 9)]
chk('published_sequence_A118896_reproduced_for_k_1_to_8', _counts == A118896,
    ','.join(map(str, _counts)))
# the census enumerates n = a^2 b^3 with b squarefree; that this parameterisation
# reaches EVERY square-full number is what makes the census exhaustive, so it is
# checked against complete factorisation rather than assumed.
_EQ = 10 ** 5
_adm5, _sqf5 = sieve_flags(isqrt(_EQ) + 1)
_enum = set()
_b = 1
while _b ** 3 <= _EQ:
    if _b == 1 or _sqf5[_b]:
        _a = 1
        while _a * _a * _b ** 3 <= _EQ:
            _enum.add(_a * _a * _b ** 3)
            _a += 1
    _b += 1
_brute5 = set(v for v in range(1, _EQ + 1) if square_full_small(v))
chk('the_a_squared_b_cubed_enumeration_reproduces_every_square_full_number_to_10^5',
    _enum == _brute5, '%d numbers, both ways' % len(_enum))
chk('b_equals_1_contributes_no_census_hit_since_n_is_never_a_perfect_square',
    all((v * v + v + 1) - v * v == v + 1 and (v + 1) ** 2 - (v * v + v + 1) == v
        for v in range(1, 1001)),
    'n - x^2 = x+1 > 0 and (x+1)^2 - n = x > 0 identically, so x^2 < n < (x+1)^2')

# --- 7 -------------------------------------------------------------------
head('Step 7: the census -- no prime witness with x <= 10^8')
_brute = sorted(x for x in range(1, BRUTE_XMAX + 1)
                if square_full_small(x * x + x + 1))
chk('brute_force_factorisation_up_to_x_3000_finds_exactly_x_18', _brute == [18],
    'hits %s' % _brute)
_c0 = census(BRUTE_XMAX)
chk('the_pruned_census_agrees_with_brute_force_up_to_x_3000',
    [h[0] for h in _c0[0]] == _brute, 'hits %s' % [h[0] for h in _c0[0]])
_c1, _pairs1, _na1, _nb1 = census(ANTIPRUNE_XMAX, restrict=True)
_c2, _pairs2, _na2, _nb2 = census(ANTIPRUNE_XMAX, restrict=False)
chk('anti_prune_control_pruned_and_unpruned_censuses_agree_up_to_x_10^6',
    [h[0] for h in _c1] == [h[0] for h in _c2],
    'pruned %d pairs vs unpruned %d pairs, identical hit set %s'
    % (_pairs1, _pairs2, [h[0] for h in _c1]))
_c3, _pairs3, _na3, _nb3 = census(CENSUS_XMAX)
note('census to x <= %d: %d admissible a, %d admissible squarefree b, '
     '%d (a,b) pairs examined' % (CENSUS_XMAX, _na3, _nb3, _pairs3))
for _h in _c3:
    note('  x = %d  n = %d = %d^2 * %d^3' % (_h[0], _h[0] ** 2 + _h[0] + 1, _h[1], _h[2]))
chk('census_to_10^8_finds_exactly_the_two_composite_near_misses_18_and_88916',
    [h[0] for h in _c3] == NEAR_MISSES, 'hits %s' % [h[0] for h in _c3])
chk('neither_of_them_is_prime_so_no_prime_witness_has_x_at_most_10^8',
    not any(prime_by_trial(h[0]) for h in _c3))
chk('the_witness_is_larger_than_the_census_bound', X > CENSUS_XMAX,
    'x has 32 digits; the census decides only the first 8')

# --- 8 -------------------------------------------------------------------
head('Step 8: the answer')
_answer = (prove_prime(X)
           and all(e >= 2 for _, e in N_FACTORS)
           and N % (_rad * _rad) == 0
           and N == X * X + X + 1)
chk('x_is_prime_AND_x^2+x+1_is_square_full_so_Question_sqfull_answers_YES', _answer)

print()
note('SCOPE -- what this program does NOT cover. (i) Nothing bibliographic is '
     're-derived: the locator of the question (restrained.tex line 852 of the '
     'arXiv:2607.19746v1 e-print), the provenance of the witness integer as '
     'a(27) of OEIS A296376 and line 27 of its b-file b296376.txt, and the 1902 '
     'Majol reference are transcribed from those sources and are not fetched or '
     'checked here, and neither is the (negative) fact that x=18 is NOT an '
     'example recorded in that e-print.  (ii) NOT RE-RUN: any x > 10^8. '
     'The census in Step 7 decides '
     'every x with 1 <= x <= 10^8 and finds no prime witness there, so no claim '
     'of minimality is made or supported for the exhibited 32-digit x; a smaller '
     'prime witness with 9 to 31 digits is not excluded by anything here.  '
     '(iii) NOT RE-RUN: whether infinitely many primes x have x^2+x+1 '
     'square-full, and whether the two branches exhibited in Step 5 exhaust the '
     'solutions of u^2 - 1372 a^2 = -3 -- Step 5 checks 20 members of one branch '
     'and 4 of the other, which supports the remark about why the object is hard '
     'to find and proves no structure theorem.  (iv) The Lucas certificates '
     'prove primality of x and of the factors of n; they say nothing about '
     'whether the factorisation of n printed in the paper is the only one, '
     'beyond the re-multiplication check in Step 1.')

print()
if _state['fail'] == 0:
    print('VERDICT: ALL %d CHECKS PASS' % _state['pass'])
    sys.exit(0)
print('VERDICT: %d of %d CHECKS FAILED' % (_state['fail'],
                                           _state['pass'] + _state['fail']))
sys.exit(1)
