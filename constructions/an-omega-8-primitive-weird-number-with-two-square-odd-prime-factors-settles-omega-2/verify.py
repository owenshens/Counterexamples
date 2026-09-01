#!/usr/bin/env python3
"""Independent verification of the note

    "A Primitive Weird Number with Two Square Odd Prime Factors and Omega = 8"

Python 3.9+, STANDARD LIBRARY ONLY (no third-party module, no external data file), exact
integer / Fraction arithmetic throughout: no floating-point value ever decides a check.

The inputs taken from the paper are the strings it prints -- the factorisation, the value of m,
the product form of sigma(m), and the five-row table of maximal divisors. One further input, the
index sequence of Step 7, is NOT printed by the paper: the submitted version contains no Remark,
no index sequence, no prefix ratio and no interval of 167 primes, so Step 7 corresponds to nothing
in the paper and no claim of the paper rests on it. Everything else is re-derived here from the
factorisation alone and compared against what the paper prints.

Each check prints one `PASS <name> [detail]` line; the program closes with a single
`VERDICT: ALL <n> CHECKS PASS` line and exits 0 if and only if every check passed.
"""

import sys
from fractions import Fraction
from math import isqrt

# ----------------------------------------------------------------------------------------------
# THE OBJECT, EXACTLY AS THE PAPER PRINTS IT
# ----------------------------------------------------------------------------------------------
FACTORIZATION_AS_PRINTED = "2^2 * 13^2 * 19^2 * 46219 * 1108619"
M_AS_PRINTED = 12504224434300196
SIGMA_AS_PRINTED = 25008448868600400
SIGMA_PRODUCT_FORM_AS_PRINTED = "7 * 183 * 381 * 46220 * 1108620"
TWO_M_AS_PRINTED = 25008448868600392
DELTA_AS_PRINTED = 8
DIVISORS_AT_MOST_DELTA_AS_PRINTED = [1, 2, 4]
SUM_OF_THOSE_AS_PRINTED = 7
SQUARE_ODD_PRIMES_AS_PRINTED = [13, 19]
OMEGA_AS_PRINTED = 8
SMALL_OMEGA_AS_PRINTED = 5

# The table of section 2, item (iv): x, sigma(x), 2x -- transcribed digit for digit.
MAXIMAL_DIVISOR_TABLE_AS_PRINTED = [
    ("m/2", 2, 10717906657971600, 12504224434300196),
    ("m/13", 13, 1913214667543200, 1923726836046184),
    ("m/19", 19, 1312779468168000, 1316234150978968),
    ("m/46219", 46219, 541074185820, 541085892568),
    ("m/1108619", 1108619, 22558179420, 22558199768),
]

# An index sequence for m, in the terminology of Amato-Hasler-Melfi-Parton, written
# [1^2, 2^2, 1^2, 167, -1]; as (index, exponent) pairs in increasing prime order that is the list
# below. The SUBMITTED PAPER PRINTS NONE OF THIS: it has no Remark, no index sequence, no prefix
# ratio 488061/11 and no interval of 167 primes. The "AS_PRINTED" suffixes on the next five names
# are inherited from an earlier draft that did print them; the values are transcribed, not read
# off the submitted paper, and Step 7 below checks material the paper does not claim.
INDEX_SEQUENCE_AS_PRINTED = [(1, 2), (2, 2), (1, 2), (167, 1), (-1, 1)]
PREFIX_C_NUMERATOR_AS_PRINTED = 488061
PREFIX_C_DENOMINATOR_AS_PRINTED = 11
PRIMES_IN_INTERVAL_AS_PRINTED = 167
FIRST_PRIME_IN_INTERVAL_AS_PRINTED = 44371
LAST_PRIME_IN_INTERVAL_AS_PRINTED = 46219

# Bounds that the paper's section 1 quotes from Open Question 1 of Amato-Hasler-Melfi-Parton, as
# that question stands in arXiv:1802.07178v2; the paper states that the published version was not
# compared with v2. So these are the ARXIV v2 values, TRANSCRIBED and not re-derived here; see the
# SCOPE lines at the end. The two identifiers keep the word PUBLISHED for continuity with the
# recorded run, but what is checked against them is the arXiv v2 text, not a published galley.
PUBLISHED_LOWER_BOUND_OMEGA_2 = 8
PUBLISHED_UPPER_BOUND_OMEGA_2 = 12

# ----------------------------------------------------------------------------------------------
# CHECK HARNESS
# ----------------------------------------------------------------------------------------------
_passed = 0
_failed = 0


def check(name, ok, detail=""):
    global _passed, _failed
    if ok:
        _passed += 1
        print("PASS %s%s" % (name, (" [%s]" % detail) if detail else ""))
    else:
        _failed += 1
        print("FAIL %s%s" % (name, (" [%s]" % detail) if detail else ""))


def note(text):
    print("NOTE " + text)


def step(text):
    print()
    print("=== " + text)


# ----------------------------------------------------------------------------------------------
# ARITHMETIC (exact, standard library only)
# ----------------------------------------------------------------------------------------------
def is_prime(n):
    """Deterministic trial division. n < 2**63 here, so this is exact and fast enough."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def sigma_pe(p, e):
    """sigma(p**e) as an exact integer, by summation (no division, no inverse)."""
    s = 0
    t = 1
    for _ in range(e + 1):
        s += t
        t *= p
    return s


def parse_factorization(text):
    """'2^2 * 13^2 * 19 * 5' -> [(2, 2), (13, 2), (19, 1), (5, 1)]. The paper's own notation."""
    out = []
    for tok in text.split("*"):
        tok = tok.strip()
        if "^" in tok:
            b, e = tok.split("^")
            out.append((int(b.strip()), int(e.strip())))
        else:
            out.append((int(tok), 1))
    return out


def divisor_exponent_tuples(fac):
    """Every divisor of prod p**e, as an exponent tuple, in no particular order."""
    tuples = [()]
    for _p, e in fac:
        tuples = [t + (i,) for t in tuples for i in range(e + 1)]
    return tuples


def value_of(fac, expo):
    v = 1
    for (p, _e), i in zip(fac, expo):
        v *= p ** i
    return v


def sigma_of(fac, expo):
    s = 1
    for (p, _e), i in zip(fac, expo):
        s *= sigma_pe(p, i)
    return s


def sieve_primes(limit):
    """All primes <= limit."""
    if limit < 2:
        return []
    flags = bytearray([1]) * (limit + 1)
    flags[0] = flags[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if flags[i]:
            flags[i * i::i] = bytearray(len(flags[i * i::i]))
    return [i for i in range(limit + 1) if flags[i]]


def sigma_full(n):
    """sigma(n) for a small n, by trial division. Used only on the control integers."""
    s = 0
    d = 1
    while d * d <= n:
        if n % d == 0:
            s += d
            q = n // d
            if q != d:
                s += q
        d += 1
    return s


def divisors_full(n):
    ds = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            ds.append(d)
            q = n // d
            if q != d:
                ds.append(q)
        d += 1
    return sorted(ds)


def representable(target, values):
    """True iff some subset of the distinct integers `values` sums exactly to `target`.
    Exact subset-sum by dynamic programming over a set of reachable sums."""
    reach = {0}
    for v in values:
        if v > target:
            continue
        reach |= {r + v for r in reach if r + v <= target}
        if target in reach:
            return True
    return target in reach


def weird_by_criterion(n):
    """(is_abundant, is_weird, Delta, witness_or_None) using Lemma 2 of the paper: an abundant n
    is weird iff Delta(n) is not a sum of distinct proper divisors of n. Only divisors <= Delta
    can occur in such a sum, so the search is over those alone."""
    delta = sigma_full(n) - 2 * n
    if delta <= 0:
        return (False, False, delta, None)
    eligible = [d for d in divisors_full(n) if d <= delta and d != n]
    return (True, not representable(delta, eligible), delta, tuple(eligible))


# ----------------------------------------------------------------------------------------------
print("verification of the note: a primitive weird number with two square odd prime")
print("factors and Omega = 8  --  m = %s" % FACTORIZATION_AS_PRINTED)
print("python %s, exact integer and Fraction arithmetic only" % sys.version.split()[0])

# ----------------------------------------------------------------------------------------------
step("Step 1: the exhibited object, read from the paper's own notation")
# ----------------------------------------------------------------------------------------------
fac = parse_factorization(FACTORIZATION_AS_PRINTED)
note("parsed factorisation: %s" % (fac,))
check("factorization_string_parses_to_five_prime_powers", len(fac) == 5,
      "%d prime powers" % len(fac))
check("bases_are_in_strictly_increasing_order",
      all(fac[i][0] < fac[i + 1][0] for i in range(len(fac) - 1)),
      "bases %s" % ([p for p, _ in fac],))
check("all_five_bases_are_prime", all(is_prime(p) for p, _ in fac),
      "2, 13, 19, 46219, 1108619 all prime by trial division")
check("cofactor_46219_is_prime_no_factor_up_to_isqrt",
      is_prime(46219) and all(46219 % d for d in range(2, isqrt(46219) + 1)),
      "isqrt(46219) = %d" % isqrt(46219))
check("cofactor_1108619_is_prime_no_factor_up_to_isqrt",
      is_prime(1108619) and all(1108619 % d for d in range(2, isqrt(1108619) + 1)),
      "isqrt(1108619) = %d" % isqrt(1108619))

m = 1
for p, e in fac:
    m *= p ** e
check("product_of_prime_powers_equals_the_m_printed_in_the_paper", m == M_AS_PRINTED,
      "m = %d" % m)

omega_big = sum(e for _p, e in fac)
omega_small = len(fac)
check("Omega_of_m_is_8", omega_big == OMEGA_AS_PRINTED, "2+2+2+1+1 = %d" % omega_big)
check("omega_of_m_is_5", omega_small == SMALL_OMEGA_AS_PRINTED,
      "%d distinct primes" % omega_small)

sq_ge2 = sorted(p for p, e in fac if p % 2 == 1 and e >= 2)
sq_eq2 = sorted(p for p, e in fac if p % 2 == 1 and e == 2)
check("square_odd_prime_factors_are_exactly_13_and_19", sq_ge2 == SQUARE_ODD_PRIMES_AS_PRINTED,
      "odd p with p^2 | m: %s" % (sq_ge2,))
check("there_are_exactly_two_square_odd_prime_factors", len(sq_ge2) == 2,
      "n = 2 line of Open Question 1")
check("the_two_readings_of_square_odd_prime_factor_agree_here", sq_ge2 == sq_eq2,
      "exponent >= 2 and exponent exactly 2 give the same set %s" % (sq_eq2,))
check("no_odd_prime_of_m_has_exponent_three_or_more",
      all(e <= 2 for p, e in fac if p % 2 == 1),
      "odd exponents %s" % ([e for p, e in fac if p % 2 == 1],))

# ----------------------------------------------------------------------------------------------
step("Step 2: sigma(m) by three independent routes, and Delta(m)")
# ----------------------------------------------------------------------------------------------
sigma_mult = 1
for p, e in fac:
    sigma_mult *= sigma_pe(p, e)
check("sigma_by_multiplicative_formula_matches_the_paper", sigma_mult == SIGMA_AS_PRINTED,
      "sigma(m) = %d" % sigma_mult)

prod_form = 1
for tok in SIGMA_PRODUCT_FORM_AS_PRINTED.split("*"):
    prod_form *= int(tok.strip())
check("the_papers_printed_product_form_of_sigma_evaluates_to_the_same_integer",
      prod_form == sigma_mult, "%s = %d" % (SIGMA_PRODUCT_FORM_AS_PRINTED, prod_form))

expo_tuples = divisor_exponent_tuples(fac)
all_divisors = sorted(value_of(fac, t) for t in expo_tuples)
check("number_of_divisors_of_m_is_108", len(all_divisors) == 108,
      "3*3*3*2*2 = %d, all distinct" % len(set(all_divisors)))
check("sigma_by_explicit_summation_over_all_108_divisors_matches",
      sum(all_divisors) == sigma_mult, "sum of divisors = %d" % sum(all_divisors))

check("two_m_matches_the_paper", 2 * m == TWO_M_AS_PRINTED, "2m = %d" % (2 * m))
delta = sigma_mult - 2 * m
check("Delta_of_m_is_8", delta == DELTA_AS_PRINTED, "Delta = sigma(m) - 2m = %d" % delta)
check("m_is_abundant", delta > 0, "Delta = %d > 0" % delta)

# ----------------------------------------------------------------------------------------------
step("Step 3: m is weird -- Lemma 2 plus a sum of three integers")
# ----------------------------------------------------------------------------------------------
least_odd_prime = min(p for p, _e in fac if p % 2 == 1)
check("least_odd_prime_factor_of_m_exceeds_Delta", least_odd_prime > delta,
      "least odd prime %d > Delta %d" % (least_odd_prime, delta))
small = [d for d in all_divisors if d <= delta]
check("divisors_of_m_at_most_Delta_are_exactly_1_2_4",
      small == DIVISORS_AT_MOST_DELTA_AS_PRINTED, "%s" % (small,))
check("their_sum_is_7_which_is_less_than_Delta",
      sum(small) == SUM_OF_THOSE_AS_PRINTED and sum(small) < delta,
      "1+2+4 = %d < %d" % (sum(small), delta))
proper = [d for d in all_divisors if d != m]
eligible = [d for d in proper if d <= delta]
check("no_subset_of_the_proper_divisors_of_m_sums_to_Delta",
      not representable(delta, eligible),
      "exhaustive over the %d eligible divisors (all subsets)" % len(eligible))
check("m_is_weird_by_Lemma_2", delta > 0 and not representable(delta, eligible),
      "abundant and Delta not representable")

# both polarities of the decider, on published integers
for n_ctrl, want_weird, label in [(70, True, "A002975"), (836, True, "A002975"),
                                  (4030, True, "A002975"), (5830, True, "A002975"),
                                  (12, False, "semiperfect"), (20, False, "semiperfect"),
                                  (945, False, "semiperfect")]:
    ab, wd, dl, elig = weird_by_criterion(n_ctrl)
    check("control_%d_decided_%s" % (n_ctrl, "WEIRD" if want_weird else "NOT_WEIRD"),
          ab and wd == want_weird,
          "Delta=%d, eligible divisors %s, %s (%s)"
          % (dl, list(elig), "weird" if wd else "semiperfect", label))

# ----------------------------------------------------------------------------------------------
step("Step 4: m is primitive abundant, hence primitive weird")
# ----------------------------------------------------------------------------------------------
table_ok = True
table_detail = []
for label, q, sigma_printed, twox_printed in MAXIMAL_DIVISOR_TABLE_AS_PRINTED:
    x = m // q
    ex = [(i - 1 if p == q else i) for (p, _e), i in zip(fac, [e for _p, e in fac])]
    sx = sigma_of(fac, tuple(ex))
    if not (sx == sigma_printed and 2 * x == twox_printed):
        table_ok = False
        table_detail.append("%s recomputed sigma=%d 2x=%d" % (label, sx, 2 * x))
check("the_five_row_table_of_maximal_divisors_matches_the_paper_digit_for_digit", table_ok,
      "; ".join(table_detail) if table_detail else "5 of 5 rows, sigma(x) and 2x both")

max_def = all(sigma_of(fac, tuple((i - 1 if p == q else i)
                                  for (p, _e), i in zip(fac, [e for _p, e in fac])))
              < 2 * (m // q)
              for _label, q, _s, _t in MAXIMAL_DIVISOR_TABLE_AS_PRINTED)
check("every_one_of_the_five_maximal_divisors_is_strictly_deficient", max_def,
      "sigma(m/q) < 2(m/q) for q in {2, 13, 19, 46219, 1108619}")

deficient_all = True
for t in expo_tuples:
    v = value_of(fac, t)
    if v == m:
        continue
    if sigma_of(fac, t) >= 2 * v:
        deficient_all = False
        break
check("all_107_proper_divisors_of_m_are_strictly_deficient", deficient_all,
      "checked one by one, exact integers")

# the monotonicity the proof of item (iv) leans on, verified on this divisor lattice
mono = True
abund = {t: Fraction(sigma_of(fac, t), value_of(fac, t)) for t in expo_tuples}
for a in expo_tuples:
    for b in expo_tuples:
        if all(x <= y for x, y in zip(a, b)) and abund[a] > abund[b]:
            mono = False
            break
    if not mono:
        break
check("abundancy_is_nondecreasing_along_divisibility_on_the_divisor_lattice_of_m", mono,
      "all %d ordered pairs of the 108 divisors, exact Fractions" % (len(expo_tuples) ** 2))

check("m_is_primitive_abundant", delta > 0 and deficient_all,
      "abundant, and every proper divisor deficient")
check("m_is_a_primitive_weird_number",
      delta > 0 and not representable(delta, eligible) and deficient_all,
      "weird, and no proper divisor is even abundant, so none is weird")

# ----------------------------------------------------------------------------------------------
step("Step 5: the closure Omega_2 = 8, conditional on the quoted lower bound Omega_2 >= 8")
# ----------------------------------------------------------------------------------------------
check("the_witness_meets_the_lower_bound_for_Omega_2_quoted_from_arXiv_1802_07178v2",
      omega_big == PUBLISHED_LOWER_BOUND_OMEGA_2,
      "Omega(m) = %d, lower bound %d quoted from arXiv:1802.07178v2 and not reproved here; the "
      "paper does not compare that source with its published version"
      % (omega_big, PUBLISHED_LOWER_BOUND_OMEGA_2))
check("the_witness_is_four_below_the_transcribed_upper_bound_12_not_a_claim_of_the_paper",
      PUBLISHED_UPPER_BOUND_OMEGA_2 - omega_big == 4,
      "12 - 8 = %d; the arithmetic is correct, but the submitted paper makes no such comparison, "
      "so this check now corresponds to no claim of the paper"
      % (PUBLISHED_UPPER_BOUND_OMEGA_2 - omega_big))
check("Omega_2_equals_8",
      omega_big <= PUBLISHED_UPPER_BOUND_OMEGA_2
      and omega_big == PUBLISHED_LOWER_BOUND_OMEGA_2
      and len(sq_ge2) == 2 and delta > 0 and not representable(delta, eligible)
      and deficient_all,
      "a PWN with exactly 2 square odd primes and Omega = 8 gives Omega_2 <= 8; the source's "
      "own theorem gives Omega_2 >= 8")
check("the_general_bound_Omega_n_ge_8_for_n_ge_2_is_attained_at_n_equals_2",
      omega_big == 8 and len(sq_ge2) == 2,
      "attained by m -- granting, as the paper's corollary does, the quoted bound Omega_n >= 8")

# ----------------------------------------------------------------------------------------------
step("Step 6: consistency with the third-party record of the same integer")
# ----------------------------------------------------------------------------------------------
check("the_defining_relation_of_OEIS_A063788_holds_for_m",
      sigma_mult == 2 * m + omega_big,
      "sigma(m) = 2m + Omega(m): %d = %d + %d" % (sigma_mult, 2 * m, omega_big))

# ----------------------------------------------------------------------------------------------
step("Step 7: an index sequence for m -- the submitted paper prints no such sequence, so these "
     "six checks correspond to nothing in it and establish more than it claims")
# ----------------------------------------------------------------------------------------------
primes = sieve_primes(1108620)
prime_set = set(primes)


def count_primes_in(lo_exclusive, hi_inclusive):
    """#{ primes q : lo_exclusive < q <= hi_inclusive }, both bounds exact Fractions/ints."""
    return sum(1 for q in primes if Fraction(q) > lo_exclusive and Fraction(q) <= hi_inclusive)


idx = []
n_pref = 1
s_pref = 1
c_at_third_prefix = None
for p, e in fac:
    d_pref = 2 * n_pref - s_pref
    assert d_pref > 0, "prefix not deficient"
    c = Fraction(s_pref, d_pref)
    if n_pref == 244036:
        c_at_third_prefix = c
    if Fraction(p) > c:
        k = sum(1 for q in primes if Fraction(q) > c and q <= p)
        idx.append((k, e))
    else:
        k = sum(1 for q in primes if q >= p and Fraction(q) <= c)
        idx.append((-k, e))
    n_pref *= p ** e
    s_pref *= sigma_pe(p, e)

check("the_index_sequence_re_derived_here_matches_the_transcribed_value",
      idx == INDEX_SEQUENCE_AS_PRINTED,
      "%s -- transcribed value, not printed by the paper" % (idx,))
check("the_prefix_ratio_of_the_transcribed_index_sequence_is_correct",
      c_at_third_prefix == Fraction(PREFIX_C_NUMERATOR_AS_PRINTED,
                                    PREFIX_C_DENOMINATOR_AS_PRINTED),
      "sigma(n)/(2n-sigma(n)) = %s for n = 2^2*13^2*19^2 = 244036; the paper prints no such ratio"
      % (c_at_third_prefix,))
in_interval = [q for q in primes
               if Fraction(q) > c_at_third_prefix and q <= 46219]
check("the_interval_contains_exactly_167_primes",
      len(in_interval) == PRIMES_IN_INTERVAL_AS_PRINTED,
      "%d primes in (488061/11, 46219]; the paper prints no such interval" % len(in_interval))
check("the_first_and_last_primes_of_that_interval_match_the_transcribed_sequence",
      in_interval[0] == FIRST_PRIME_IN_INTERVAL_AS_PRINTED
      and in_interval[-1] == LAST_PRIME_IN_INTERVAL_AS_PRINTED,
      "first %d, last %d" % (in_interval[0], in_interval[-1]))

# Delta re-derived through the terminal prime, a route that never forms sigma(m):
n4 = 1
s4 = 1
for p, e in fac[:-1]:
    n4 *= p ** e
    s4 *= sigma_pe(p, e)
r = fac[-1][0]
d4 = 2 * n4 - s4
check("Delta_re_derived_by_the_terminal_prime_recursion", s4 - r * d4 == delta,
      "sigma(n) - r*(2n - sigma(n)) = %d - %d*%d = %d" % (s4, r, d4, s4 - r * d4))
check("the_terminal_prime_1108619_is_in_the_prime_sieve_used_above", r in prime_set,
      "sieve limit 1108620")

# ----------------------------------------------------------------------------------------------
print()
note("SCOPE. Re-derived here, from the factorisation printed in the paper and nothing else: "
     "the value of m, Omega(m) = 8, omega(m) = 5, the two square odd prime factors, sigma(m) by "
     "three routes, Delta(m) = 8, weirdness via Lemma 2 with an exhaustive subset-sum over the "
     "eligible divisors, primitive abundance by checking ALL 107 proper divisors (not only the "
     "five maximal ones), the five-row table digit for digit, the A063788 relation, and the "
     "index sequence of Step 7 -- which the paper does not print.")
note("NOT CLAIMED BY THE PAPER. Step 7's six checks re-derive an index sequence, a prefix ratio "
     "488061/11 and an interval of 167 primes, none of which appears in the submitted paper: it "
     "contains no Remark and prints no such quantities. Step 5's check that the witness is four "
     "below the upper bound 12 is likewise a comparison the paper no longer makes. Those checks "
     "are correct and are retained rather than deleted, but nothing in the paper rests on them; a "
     "transcript establishing more than its paper claims is not a defect in the paper.")
note("NOT RE-RUN, and this program makes no claim about any of it: the lower bound Omega_2 >= 8 "
     "is TRANSCRIBED from Theorem thm:patterns of Amato-Hasler-Melfi-Parton as it stands in "
     "arXiv:1802.07178v2 and is not proved here, so the equality Omega_2 = 8 above is conditional "
     "on that quoted theorem; the upper bound 12 is likewise transcribed from their table in the "
     "same arXiv v2 source. Neither bound is checked against the published journal version, which "
     "the paper states was not compared with v2, so 'published' would overstate what was checked. "
     "The exhaustive census that FOUND "
     "m is not re-run here and none of its by-products (that the cell Omega = 8 with at least "
     "two squared odd primes contains exactly one PWN, hence Omega_3 >= 9) is checked -- the "
     "paper does not claim them either. No value of Omega_n for n >= 3 is touched. The "
     "definition of the index sequence used in Step 7 is transcribed from the source's "
     "terminology; it is descriptive, and no claim of the paper rests on it.")

print()
if _failed:
    print("VERDICT: %d CHECK(S) FAILED out of %d" % (_failed, _passed + _failed))
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % _passed)
sys.exit(0)
