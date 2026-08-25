#!/usr/bin/env python3
"""Verification of counterexamples to a sorting conjecture for reduced
(stable) Kronecker coefficients.

TAKEN FROM THE PAPER (inputs, not re-derived here)
--------------------------------------------------
  * Definition of the reduced Kronecker coefficient as the stable limit
    gbar(a,b;c) = lim_N g(a[N], b[N]; c[N]) with a[N] = (N-|a|, a).
  * The sorting maps: (d_1,d_2,...) = decreasing rearrangement of all parts
    of lambda and mu together; sort1 = (d_1,d_3,...), sort2 = (d_2,d_4,...).
  * The conjecture under attack: gbar(sort1(lambda,mu), sort2(lambda,mu); nu)
    >= gbar(lambda, mu; nu) for every partition nu.
  * The exhibited counterexample family: lambda = (k), mu = (k-1,1), nu = (1),
    for k >= 2, whose sorted pair is ((k,1), (k-1)); the two asserted
    coefficient values are 0 on the left and 1 on the right.
  * The closed formula asserted for nu = (1):
    gbar(alpha,beta;(1)) = |R(alpha) cap R(beta)| + 1[alpha in R(beta)]
    + 1[beta in R(alpha)], where R(.) deletes one removable corner.
  * The branching identity g(sigma,(N-1,1);rho) = |R(sigma) cap R(rho)|
    - delta(sigma,rho) used inside the proof of that formula.
  * The claimed minimality statement: up to swapping the two factors, the
    k = 2 instance is the only counterexample among pairs of nonempty
    partitions with at most four boxes in total, the only two pairs of that
    size moved by sorting being ((1),(1,1,1)) and ((2),(1,1)).
  * The claimed multipartition extension being refuted by the same family,
    AND, in particular, the identification of the two partitions written
    (lambda cup mu)^{[1,2]} and (lambda cup mu)^{[2,2]} in the cited
    multipartition conjecture with the pair (sort1, sort2).  That
    identification is a reading of the notation of an external paper.  It is
    TAKEN ON TRUST here and no check below tests it: coding the r-block
    index-residue split and comparing its r=2 output with sort1, sort2 would
    only compare one transcription of a single definition with another, so no
    such check is reported.  The program's closing paragraph says this in its
    own output.

DERIVED HERE (computed, nothing asserted)
-----------------------------------------
  * Ordinary Kronecker coefficients g(a,b;c) computed from scratch by
    characters of the symmetric group (integer arithmetic only), via
    g = (1/n!) sum_rho |C_rho| chi_a(rho) chi_b(rho) chi_c(rho).
  * Reduced Kronecker coefficients obtained only after the sequence
    g(a[N],b[N];c[N]) is observed constant over STAB_RUNS_LOAD consecutive
    admissible N; a non-constant window is an error, not a silent average.
    A longer window of STAB_RUNS_WINDOW consecutive admissible N is exhibited
    only for the members k in K_STAB_WINDOW_RANGE, not for the whole
    character range; the closing paragraph states which range gets which.
  * The sorting maps recomputed from their definition: the decreasing
    rearrangement of the parts of (k) and (k-1,1) is (k,k-1,1), whence
    sort1 = (k,1) and sort2 = (k-1).
  * The strict violation gbar(sort1,sort2;(1)) < gbar(lambda,mu;(1)), both
    sides computed from characters, for every k in the character range, and
    from the independently validated closed formula for a much longer range.
  * The closed formula and the branching identity checked against the
    character computation on every small pair.
  * An exhaustive census over all pairs of nonempty partitions with total
    size <= 4 (the minimality claim) and a wider sweep to total size 6,
    together with evidence that the truncation of the nu range is safe.
"""

import sys
from fractions import Fraction

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)


def finish():
    n = len(CHECKS)
    bad = [c for c in CHECKS if not c[1]]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


# ----------------------------------------------------------------------
# Partitions
# ----------------------------------------------------------------------

_PART_CACHE = {}


def partitions(n, cap=None):
    """All partitions of n as weakly decreasing tuples, parts <= cap."""
    if cap is None:
        cap = n
    key = (n, cap)
    if key in _PART_CACHE:
        return _PART_CACHE[key]
    if n == 0:
        out = [()]
    else:
        out = []
        for first in range(min(n, cap), 0, -1):
            for rest in partitions(n - first, first):
                out.append((first,) + rest)
    _PART_CACHE[key] = out
    return out


def is_partition(p):
    return all(isinstance(x, int) and x > 0 for x in p) and all(
        p[i] >= p[i + 1] for i in range(len(p) - 1))


def norm(p):
    """Drop trailing zeros; used after arithmetic on part lists."""
    q = [x for x in p if x != 0]
    return tuple(q)


def removable_corners(p):
    """R(p): partitions obtained by deleting one removable corner box."""
    out = set()
    for i in range(len(p)):
        if i + 1 == len(p) or p[i] > p[i + 1]:
            q = list(p)
            q[i] -= 1
            out.add(norm(q))
    return out


# ----------------------------------------------------------------------
# Characters of the symmetric group (Murnaghan-Nakayama, exact integers)
# ----------------------------------------------------------------------

def beta_set(lam):
    """First-column hook lengths: distinct integers lam_i + (m-1-i)."""
    m = len(lam)
    return [lam[i] + (m - 1 - i) for i in range(m)]


def from_beta(bs):
    """Inverse of beta_set on a set of distinct non-negative integers."""
    b = sorted(bs, reverse=True)
    m = len(b)
    lam = [b[i] - (m - 1 - i) for i in range(m)]
    return norm(lam)


def strip_removals(lam, r):
    """All (mu, sign) with mu obtained from lam by removing a rim hook of
    size r; sign = (-1)^(rows spanned - 1)."""
    bs = beta_set(lam)
    bset = set(bs)
    out = []
    for b in bs:
        t = b - r
        if t < 0 or t in bset:
            continue
        crossed = sum(1 for x in bset if t < x < b)
        newset = set(bset)
        newset.discard(b)
        newset.add(t)
        out.append((from_beta(newset), -1 if crossed % 2 else 1))
    return out


_CHI = {}


def chi(lam, rho):
    """Irreducible character chi_lam evaluated on cycle type rho."""
    lam = norm(lam)
    rho = tuple(sorted((x for x in rho if x), reverse=True))
    if not is_partition(lam):
        raise ValueError("not a partition: %s" % (lam,))
    if sum(lam) != sum(rho):
        raise ValueError("size mismatch")
    if not rho:
        return 1
    key = (lam, rho)
    if key in _CHI:
        return _CHI[key]
    r, rest = rho[0], rho[1:]
    total = 0
    for mu, sgn in strip_removals(lam, r):
        total += sgn * chi(mu, rest)
    _CHI[key] = total
    return total


# ----------------------------------------------------------------------
# Ordinary and reduced Kronecker coefficients
# ----------------------------------------------------------------------

def centraliser_order(rho):
    """z_rho = prod_i i^{m_i} m_i!."""
    z = 1
    seen = {}
    for x in rho:
        seen[x] = seen.get(x, 0) + 1
    for part, mult in seen.items():
        z *= part ** mult
        for j in range(2, mult + 1):
            z *= j
    return z


_KRON = {}


def kron(a, b, c):
    """g_{a,b}^{c} = sum_rho chi_a(rho) chi_b(rho) chi_c(rho) / z_rho."""
    n = sum(a)
    if sum(b) != n or sum(c) != n:
        raise ValueError("Kronecker coefficient needs equal sizes")
    key = tuple(sorted([norm(a), norm(b), norm(c)]))
    if key in _KRON:
        return _KRON[key]
    acc = Fraction(0)
    for rho in partitions(n):
        ca = chi(a, rho)
        if ca == 0:
            continue
        cb = chi(b, rho)
        if cb == 0:
            continue
        cc = chi(c, rho)
        if cc == 0:
            continue
        acc += Fraction(ca * cb * cc, centraliser_order(rho))
    if acc.denominator != 1:
        raise ArithmeticError("non-integral Kronecker coefficient")
    val = acc.numerator
    _KRON[key] = val
    return val


def pad(alpha, N):
    """alpha[N] = (N - |alpha|, alpha), or None if not a partition."""
    first = N - sum(alpha)
    if alpha and first < alpha[0]:
        return None
    if first < 0:
        return None
    return norm((first,) + tuple(alpha))


_RED = {}

# Length of the constant window demanded of the ordinary sequence before its
# limit is accepted.  STAB_RUNS_LOAD is what every value in this program
# (including the two load-bearing ones) is built on; STAB_RUNS_WINDOW is the
# longer window exhibited by check_murnaghan_stability_observed, and only for
# the members listed in K_STAB_WINDOW_RANGE below.  Both numbers are printed
# from these constants, never spelled out in prose.
STAB_RUNS_LOAD = 3
STAB_RUNS_WINDOW = 5


def reduced_kron(alpha, beta, gamma, runs=STAB_RUNS_LOAD, extra=6):
    """Stable value of g_{alpha[N],beta[N]}^{gamma[N]}.

    Returns (value, sequence_of_(N,value)).  Stability is OBSERVED: the
    value must repeat `runs` times in a row.  Raises if it does not.
    """
    alpha, beta, gamma = norm(alpha), norm(beta), norm(gamma)
    key = (alpha, beta, gamma, runs)
    if key in _RED:
        return _RED[key]
    base = max(sum(alpha) + (alpha[0] if alpha else 0),
               sum(beta) + (beta[0] if beta else 0),
               sum(gamma) + (gamma[0] if gamma else 0), 2)
    seq = []
    N = base
    limit = base + sum(alpha) + sum(beta) + sum(gamma) + extra
    while N <= limit:
        a, b, c = pad(alpha, N), pad(beta, N), pad(gamma, N)
        if a is None or b is None or c is None:
            N += 1
            continue
        v = kron(a, b, c)
        seq.append((N, v))
        tail = seq[-runs:]
        if len(tail) == runs and all(x == v for _, x in tail):
            _RED[key] = (v, seq)
            return _RED[key]
        N += 1
    raise ArithmeticError("no stabilisation for %s %s %s: %s"
                          % (alpha, beta, gamma, seq))


# ----------------------------------------------------------------------
# VALUES TAKEN FROM THE PAPER.  Corrupting any of these must make a
# check below report FAIL; nothing here is compared only to itself.
# ----------------------------------------------------------------------

PAPER_NU = (1,)                      # the partition nu at which (eq) fails
PAPER_LHS_VALUE = 0                  # gbar_{(k,1),(k-1)}^{(1)}
PAPER_RHS_VALUE = 1                  # gbar_{(k),(k-1,1)}^{(1)}
PAPER_K_MIN = 2                      # family runs over k >= 2
PAPER_MINIMAL_PAIR = ((2,), (1, 1))  # the claimed unique small failure
PAPER_UNIQUENESS_TOTAL = 4           # ... among |lambda|+|mu| <= 4
PAPER_MOVED_PAIRS_TOTAL4 = [((1,), (1, 1, 1)), ((2,), (1, 1))]
PAPER_COLUMNS_PAIR = ((1,), (1, 1, 1))   # covered by Gui's Theorem 5.6


def family_lambda(k):
    return (k,)


def family_mu(k):
    return norm((k - 1, 1))


def family_sort1(k):
    return (k, 1)


def family_sort2(k):
    return norm((k - 1,))


# ----------------------------------------------------------------------
# The sorting maps of the conjecture
# ----------------------------------------------------------------------

def delta_sequence(lam, mu):
    """Decreasing rearrangement of all parts of lam and mu together."""
    return tuple(sorted(list(lam) + list(mu), reverse=True))


def sort1(lam, mu):
    return norm(delta_sequence(lam, mu)[0::2])


def sort2(lam, mu):
    return norm(delta_sequence(lam, mu)[1::2])



# NOTE.  An earlier revision of this program also carried a check named
# "multipartition_two_block_case_is_the_same_map", which built the r-block
# index-residue split and compared its r=2 output with (sort1, sort2).  That
# check was VACUOUS: the r=2 slices of the sorted concatenation are the very
# expressions sort1 and sort2 evaluate, so the comparison could not fail and
# it inflated the check count by one.  It has been deleted.  The step it
# pretended to establish -- that the notation of the cited multipartition
# conjecture means the index-residue split -- is a reading of an external
# paper that no program can settle; it is disclosed as taken on trust in the
# closing "NOT RE-RUN HERE" paragraph printed by main().


# ----------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------

DIM_SELF_TEST_N_MAX = 8    # sum of squared dimensions checked for n <= this
SIGN_SELF_TEST_N_MAX = 7   # trivial/sign characters checked for n <= this


def check_character_machinery():
    """Independent validation of the character routine used everywhere."""
    ok = True
    facts = []
    fact = 1
    for n in range(1, DIM_SELF_TEST_N_MAX + 1):
        fact *= n
        s = sum(chi(p, (1,) * n) ** 2 for p in partitions(n))
        facts.append("%d:%s" % (n, "ok" if s == fact else "BAD"))
        ok = ok and s == fact
    for n in range(1, SIGN_SELF_TEST_N_MAX + 1):
        for rho in partitions(n):
            if chi((n,), rho) != 1:
                ok = False
            sign = 1
            for part in rho:
                if part % 2 == 0:
                    sign = -sign
            if chi(tuple([1] * n), rho) != sign:
                ok = False
    if chi((2, 1), (3,)) != -1 or chi((2, 1), (2, 1)) != 0:
        ok = False
    return ck("character_machinery_self_test", ok,
              "sum of squared dimensions = n! for n=1..%d (%s); trivial and "
              "sign characters reproduced for n=1..%d"
              % (DIM_SELF_TEST_N_MAX, ",".join(facts), SIGN_SELF_TEST_N_MAX))


def check_padding_definition():
    """The stable padding used to define the reduced coefficients must be
    exactly alpha[N] = (N - |alpha|, alpha), admissible precisely when
    N - |alpha| >= alpha_1."""
    ok = True
    samples = []
    for n in range(0, 6):
        for alpha in partitions(n):
            for N in range(0, 16):
                got = pad(alpha, N)
                first = N - n
                legal = first >= 0 and (not alpha or first >= alpha[0])
                if not legal:
                    if got is not None:
                        ok = False
                    continue
                want = norm((first,) + tuple(alpha))
                if got != want or sum(got) != N:
                    ok = False
                if first > 0 and tuple(got[1:]) != tuple(alpha):
                    ok = False
                samples.append((alpha, N, got))
    if pad((2, 1), 7) != (4, 2, 1) or pad((2, 1), 4) is not None:
        ok = False
    return ck("stable_padding_definition", ok,
              "%d instances; e.g. (2,1)[7]=%s and (2,1)[4] inadmissible"
              % (len(samples), pad((2, 1), 7)))


# equation (branching) checked on all sigma,rho of size at most this
BRANCHING_N_MAX = 8


def check_branching_identity():
    """Equation (branching) of the paper: for sigma, rho of size N >= 2,
    g_{sigma,(N-1,1)}^{rho} = |R(sigma) cap R(rho)| - delta_{sigma,rho}."""
    ok = True
    tested = 0
    for N in range(2, BRANCHING_N_MAX + 1):
        second = norm((N - 1, 1))
        for sig in partitions(N):
            for rho in partitions(N):
                lhs = kron(sig, second, rho)
                rhs = len(removable_corners(sig) & removable_corners(rho))
                if sig == rho:
                    rhs -= 1
                tested += 1
                if lhs != rhs:
                    ok = False
    return ck("branching_identity_eq_branching", ok,
              "%d ordered pairs sigma,rho with |sigma|=2..%d"
              % (tested, BRANCHING_N_MAX))


K_CHAR_RANGE = list(range(PAPER_K_MIN, 13))   # k verified from characters
K_LEMMA_RANGE = list(range(PAPER_K_MIN, 401))  # k verified from the lemma
# The longer STAB_RUNS_WINDOW-term stability window is exhibited only for these
# members; k outside it is accepted on a STAB_RUNS_LOAD-term window.
K_STAB_WINDOW_RANGE = K_CHAR_RANGE[:6]


def check_exhibited_object():
    """Decode, count and print back the exhibited counterexample family."""
    ok = True
    lines = []
    for k in K_CHAR_RANGE[:6]:
        lam, mu = family_lambda(k), family_mu(k)
        s1, s2 = family_sort1(k), family_sort2(k)
        lines.append("k=%d lambda=%s mu=%s sort1=%s sort2=%s" %
                     (k, lam, mu, s1, s2))
    for k in K_LEMMA_RANGE:
        lam, mu = family_lambda(k), family_mu(k)
        s1, s2 = family_sort1(k), family_sort2(k)
        for p in (lam, mu, s1, s2, PAPER_NU):
            if not is_partition(p) or len(p) == 0:
                ok = False
        if sum(lam) != k or sum(mu) != k:
            ok = False
        if sum(s1) != k + 1 or sum(s2) != k - 1:
            ok = False
        if len(lam) != 1 or len(mu) != (2 if k >= 2 else 1):
            ok = False
    print("  exhibited: " + " | ".join(lines))
    return ck("exhibited_object_wellformed", ok,
              "%d members k=%d..%d; each of lambda,mu,sort1,sort2 a nonempty "
              "partition, |lambda|=|mu|=k, |sort1|+|sort2|=2k"
              % (len(K_LEMMA_RANGE), K_LEMMA_RANGE[0], K_LEMMA_RANGE[-1]))


def check_sorting_maps():
    """Recompute sort1, sort2 from the definition; they must equal the
    partitions (k,1) and (k-1) named in the paper."""
    ok = True
    for k in K_LEMMA_RANGE:
        lam, mu = family_lambda(k), family_mu(k)
        d = delta_sequence(lam, mu)
        want_d = tuple(sorted([k, k - 1, 1], reverse=True))
        if d != want_d:
            ok = False
        if sort1(lam, mu) != family_sort1(k):
            ok = False
        if sort2(lam, mu) != family_sort2(k):
            ok = False
    return ck("sorting_maps_recomputed", ok,
              "delta=(k,k-1,1), sort1=(k,1), sort2=(k-1) for k=%d..%d"
              % (K_LEMMA_RANGE[0], K_LEMMA_RANGE[-1]))


def check_hypotheses():
    """Every hypothesis of the conjecture holds for the exhibited data:
    lambda, mu are partitions, nu is a partition, the sorting maps consume
    exactly the multiset of parts, and the pair is genuinely moved (so the
    inequality is not the trivial one)."""
    ok = True
    moved = 0
    for k in K_LEMMA_RANGE:
        lam, mu = family_lambda(k), family_mu(k)
        s1, s2 = sort1(lam, mu), sort2(lam, mu)
        multiset_in = sorted(list(lam) + list(mu))
        multiset_out = sorted(list(s1) + list(s2))
        if multiset_in != multiset_out:
            ok = False
        if sum(s1) + sum(s2) != sum(lam) + sum(mu):
            ok = False
        if not is_partition(s1) or not is_partition(s2):
            ok = False
        if sorted([s1, s2]) != sorted([lam, mu]):
            moved += 1
    if moved != len(K_LEMMA_RANGE):
        ok = False
    return ck("hypotheses_of_the_conjecture_hold", ok,
              "part multiset preserved and {sort1,sort2} != {lambda,mu} for "
              "all %d members (no trivial fixed point)" % moved)


def check_reduced_values_from_characters():
    """The two load-bearing numbers, computed from S_N characters through
    the stable limit -- not read off the paper."""
    ok = True
    rows = []
    for k in K_CHAR_RANGE:
        lo, lseq = reduced_kron(family_sort1(k), family_sort2(k), PAPER_NU,
                                runs=STAB_RUNS_LOAD)
        hi, hseq = reduced_kron(family_lambda(k), family_mu(k), PAPER_NU,
                                runs=STAB_RUNS_LOAD)
        rows.append("k=%d:%d/%d" % (k, lo, hi))
        if lo != PAPER_LHS_VALUE or hi != PAPER_RHS_VALUE:
            ok = False
        # The window each value rests on must really be there: at least
        # STAB_RUNS_LOAD terms, the last STAB_RUNS_LOAD of them equal to the
        # accepted value, and at CONSECUTIVE admissible N (no gap glossed over).
        for seq, val in ((lseq, lo), (hseq, hi)):
            if len(seq) < STAB_RUNS_LOAD:
                ok = False
                continue
            tail = seq[-STAB_RUNS_LOAD:]
            if any(v != val for _, v in tail):
                ok = False
            ns = [n for n, _ in tail]
            if ns != list(range(ns[0], ns[0] + STAB_RUNS_LOAD)):
                ok = False
    print("  gbar(sort1,sort2;nu)/gbar(lambda,mu;nu): " + " ".join(rows))
    return ck("reduced_coefficients_from_characters", ok,
              "gbar(sort1,sort2;(1))=%d and gbar(lambda,mu;(1))=%d for k=%d..%d,"
              " each accepted only after %d consecutive admissible N agree"
              % (PAPER_LHS_VALUE, PAPER_RHS_VALUE,
                 K_CHAR_RANGE[0], K_CHAR_RANGE[-1], STAB_RUNS_LOAD))


def check_murnaghan_stability_observed():
    """The reduced coefficients are limits.  Confirm the ordinary sequence
    g_{alpha[N],beta[N]}^{nu[N]} is genuinely constant over STAB_RUNS_WINDOW
    consecutive admissible N, rather than sampled once.  This longer window is
    exhibited only for the members k in K_STAB_WINDOW_RANGE, which is SHORTER
    than K_CHAR_RANGE; every k in K_CHAR_RANGE outside it rests on the
    STAB_RUNS_LOAD-term window of check_reduced_values_from_characters.  The
    shortfall is stated in the closing NOT RE-RUN paragraph."""
    ok = True
    windows = []
    for k in K_STAB_WINDOW_RANGE:
        for (a, b) in [(family_sort1(k), family_sort2(k)),
                       (family_lambda(k), family_mu(k))]:
            val, seq = reduced_kron(a, b, PAPER_NU, runs=STAB_RUNS_WINDOW)
            tail = seq[-STAB_RUNS_WINDOW:]
            if len(tail) < STAB_RUNS_WINDOW or any(v != val for _, v in tail):
                ok = False
            ns = [n for n, _ in tail]
            if ns != list(range(ns[0], ns[0] + len(ns))):
                ok = False   # the window must be consecutive N, not a sample
            windows.append("k=%d N=%d..%d" % (k, tail[0][0], tail[-1][0]))
    if len(windows) != 2 * len(K_STAB_WINDOW_RANGE):
        ok = False
    print("  stability windows (%d, k=%d..%d): %s"
          % (len(windows), K_STAB_WINDOW_RANGE[0], K_STAB_WINDOW_RANGE[-1],
             "; ".join(windows)))
    return ck("murnaghan_stability_observed", ok,
              "constant over %d consecutive N in each of %d windows, k=%d..%d "
              "only; k=%d..%d rests on the %d-term window instead"
              % (STAB_RUNS_WINDOW, len(windows), K_STAB_WINDOW_RANGE[0],
                 K_STAB_WINDOW_RANGE[-1], K_STAB_WINDOW_RANGE[-1] + 1,
                 K_CHAR_RANGE[-1], STAB_RUNS_LOAD))


def check_violation_of_conclusion():
    """LOAD-BEARING.  Both sides are computed; the conjectured inequality
    gbar(sort1,sort2;nu) >= gbar(lambda,mu;nu) must be strictly violated.
    The detail line reports the pairs of values ACTUALLY COMPUTED here, so
    that it cannot echo a value quoted by the paper but never obtained."""
    ok = True
    strict = 0
    observed = set()
    for k in K_CHAR_RANGE:
        lo = reduced_kron(family_sort1(k), family_sort2(k), PAPER_NU)[0]
        hi = reduced_kron(family_lambda(k), family_mu(k), PAPER_NU)[0]
        observed.add((lo, hi))
        if lo < hi:
            strict += 1
        else:
            ok = False
    if strict != len(K_CHAR_RANGE):
        ok = False
    shown = ", ".join("(%d,%d)" % (a, b) for (a, b) in sorted(observed))
    return ck("conjectured_inequality_violated", ok,
              "computed (gbar_sorted, gbar_unsorted) in {%s} at nu=%s; strict "
              "violation on %d of %d members k=%d..%d"
              % (shown, PAPER_NU, strict, len(K_CHAR_RANGE),
                 K_CHAR_RANGE[0], K_CHAR_RANGE[-1]))


def lemma_formula(alpha, beta):
    ra, rb = removable_corners(alpha), removable_corners(beta)
    return (len(ra & rb) + (1 if norm(alpha) in rb else 0)
            + (1 if norm(beta) in ra else 0))


LEMMA_SIZE_CAP = 5


def check_lemma_one():
    """The paper's Lemma: gbar_{alpha,beta}^{(1)} equals
    |R(alpha) cap R(beta)| + 1[alpha in R(beta)] + 1[beta in R(alpha)].
    Compared against the character computation on every pair of partitions
    of size at most LEMMA_SIZE_CAP."""
    ok = True
    tested = 0
    worst = None
    for na in range(0, LEMMA_SIZE_CAP + 1):
        for a in partitions(na):
            for nb in range(0, LEMMA_SIZE_CAP + 1):
                for b in partitions(nb):
                    got = reduced_kron(a, b, PAPER_NU)[0]
                    want = lemma_formula(a, b)
                    tested += 1
                    if got != want:
                        ok = False
                        worst = (a, b, want, got)
    return ck("lemma_closed_formula_matches_characters", ok,
              "%d ordered pairs with |alpha|,|beta| <= %d%s"
              % (tested, LEMMA_SIZE_CAP,
                 "" if ok else "; mismatch at %s" % (worst,)))


def check_lemma_proof_steps():
    """The two combinatorial statements the paper's proof of the Theorem
    rests on, recomputed: (k) and (k-1,1) have the single common corner
    predecessor (k-1); and (k,1), (k-1) differ in size by 2 so every term
    of the lemma vanishes."""
    ok = True
    for k in K_LEMMA_RANGE:
        lam, mu = family_lambda(k), family_mu(k)
        common = removable_corners(lam) & removable_corners(mu)
        if common != {norm((k - 1,))}:
            ok = False
        if lam in removable_corners(mu) or mu in removable_corners(lam):
            ok = False
        s1, s2 = family_sort1(k), family_sort2(k)
        if sum(s1) - sum(s2) != 2:
            ok = False
        if removable_corners(s1) & removable_corners(s2):
            ok = False
        if s1 in removable_corners(s2) or s2 in removable_corners(s1):
            ok = False
    return ck("lemma_proof_steps_recomputed", ok,
              "unique common corner (k-1) on the right, empty on the left, "
              "for k=%d..%d" % (K_LEMMA_RANGE[0], K_LEMMA_RANGE[-1]))


def check_infinite_family_via_lemma():
    """Extend the theorem past the character range using the closed formula
    validated above: the family violates the inequality for every k in
    K_LEMMA_RANGE."""
    ok = True
    bad = []
    for k in K_LEMMA_RANGE:
        lo = lemma_formula(family_sort1(k), family_sort2(k))
        hi = lemma_formula(family_lambda(k), family_mu(k))
        if lo != PAPER_LHS_VALUE or hi != PAPER_RHS_VALUE or not lo < hi:
            ok = False
            bad.append(k)
    return ck("infinite_family_via_validated_lemma", ok,
              "%d values of k in [%d,%d] all give %d < %d%s"
              % (len(K_LEMMA_RANGE), K_LEMMA_RANGE[0], K_LEMMA_RANGE[-1],
                 PAPER_LHS_VALUE, PAPER_RHS_VALUE,
                 "" if ok else "; failures at %s" % bad[:5]))


def unordered_pairs(total):
    """All unordered pairs of nonempty partitions with |lam|+|mu| <= total."""
    out = []
    seen = set()
    for na in range(1, total):
        for a in partitions(na):
            for nb in range(1, total - na + 1):
                for b in partitions(nb):
                    key = tuple(sorted([a, b]))
                    if key in seen:
                        continue
                    seen.add(key)
                    out.append((a, b))
    return out


def failures_for_pair(lam, mu, nu_cap):
    """All nu with |nu| <= nu_cap at which the conjectured inequality fails."""
    s1, s2 = sort1(lam, mu), sort2(lam, mu)
    out = []
    for n in range(0, nu_cap + 1):
        for nu in partitions(n):
            lo = reduced_kron(s1, s2, nu)[0]
            hi = reduced_kron(lam, mu, nu)[0]
            if lo < hi:
                out.append((nu, lo, hi))
    return out


SWEEP_TOTAL = 6      # widest exhaustive census run here
NU_OVERSHOOT = 2     # how far past |lambda|+|mu| the nu window is probed


def check_nu_truncation():
    """The census below ranges nu over |nu| <= |lambda|+|mu|.  Confirm no
    coefficient in the next NU_OVERSHOOT sizes above that range is nonzero,
    for every pair and both sides.  This is EVIDENCE that the truncation loses
    no potential failure, not a proof that it cannot: nu beyond
    total+NU_OVERSHOOT is untested, and the closing paragraph says so."""
    ok = True
    tested = 0
    for (lam, mu) in unordered_pairs(SWEEP_TOTAL):
        tot = sum(lam) + sum(mu)
        s1, s2 = sort1(lam, mu), sort2(lam, mu)
        for n in range(tot + 1, tot + NU_OVERSHOOT + 1):
            for nu in partitions(n):
                for (a, b) in ((s1, s2), (lam, mu)):
                    tested += 1
                    if reduced_kron(a, b, nu)[0] != 0:
                        ok = False
    return ck("nu_range_truncation_is_safe", ok,
              "%d coefficients with |nu| in [total+1,total+%d] all vanish; "
              "|nu| > total+%d not probed"
              % (tested, NU_OVERSHOOT, NU_OVERSHOOT))


def check_uniqueness_census():
    """Exhaustive census: over all unordered pairs of nonempty partitions
    with |lambda|+|mu| <= 4 and all nu with |nu| <= 4, exactly one pair
    admits a failure, and it is the pair named in the paper."""
    pairs = unordered_pairs(PAPER_UNIQUENESS_TOTAL)
    failing = []
    for (lam, mu) in pairs:
        fs = failures_for_pair(lam, mu, PAPER_UNIQUENESS_TOTAL)
        if fs:
            failing.append((tuple(sorted([lam, mu])), fs))
    want = tuple(sorted([norm(PAPER_MINIMAL_PAIR[0]),
                         norm(PAPER_MINIMAL_PAIR[1])]))
    ok = len(failing) == 1 and failing and failing[0][0] == want
    detail = "; ".join("%s at nu=%s (%d<%d)" % (p, f[0][0], f[0][1], f[0][2])
                       for p, f in failing)
    return ck("uniqueness_census_total_at_most_%d" % PAPER_UNIQUENESS_TOTAL, ok,
              "%d unordered pairs searched, %d failing: %s"
              % (len(pairs), len(failing), detail))


def check_corollary_proof_steps():
    """The paper's proof of minimality: for total size at most 3 sorting
    fixes the unordered pair (so the inequality is an equality), and at
    total size exactly 4 the only pairs sorting moves are the two named."""
    ok = True
    below = PAPER_UNIQUENESS_TOTAL - 1
    for (lam, mu) in unordered_pairs(PAPER_UNIQUENESS_TOTAL):
        s1, s2 = sort1(lam, mu), sort2(lam, mu)
        fixed = sorted([s1, s2]) == sorted([lam, mu])
        if sum(lam) + sum(mu) <= below and not fixed:
            ok = False
    moved = sorted(tuple(sorted([lam, mu]))
                   for (lam, mu) in unordered_pairs(PAPER_UNIQUENESS_TOTAL)
                   if sorted([sort1(lam, mu), sort2(lam, mu)])
                   != sorted([lam, mu]))
    want = sorted(tuple(sorted([norm(a), norm(b)]))
                  for (a, b) in PAPER_MOVED_PAIRS_TOTAL4)
    if moved != want:
        ok = False
    return ck("minimality_proof_steps_recomputed", ok,
              "size<=%d pairs all fixed by sorting; moved pairs at size %d = %s"
              % (below, PAPER_UNIQUENESS_TOTAL, moved))


def check_columns_pair_holds():
    """The other pair sorting moves at total size 4 must satisfy the
    inequality (the paper attributes it to Gui's column theorem); computed
    here for every nu with |nu| <= 4."""
    lam, mu = norm(PAPER_COLUMNS_PAIR[0]), norm(PAPER_COLUMNS_PAIR[1])
    s1, s2 = sort1(lam, mu), sort2(lam, mu)
    rows = []
    ok = True
    for n in range(0, PAPER_UNIQUENESS_TOTAL + 1):
        for nu in partitions(n):
            lo = reduced_kron(s1, s2, nu)[0]
            hi = reduced_kron(lam, mu, nu)[0]
            rows.append("%s:%d>=%d" % (nu, lo, hi))
            if lo < hi:
                ok = False
    return ck("columns_pair_satisfies_inequality", ok,
              "pair %s/%s sorts to %s/%s; %d values of nu with |nu|<=%d "
              "checked: %s"
              % (lam, mu, s1, s2, len(rows), PAPER_UNIQUENESS_TOTAL,
                 " ".join(rows[:6]) + " ..."))


def check_minimal_pair_is_family_member():
    """The unique small failure must be the k=2 member of the family, up to
    interchanging the two factors."""
    got = tuple(sorted([norm(PAPER_MINIMAL_PAIR[0]),
                        norm(PAPER_MINIMAL_PAIR[1])]))
    fam = tuple(sorted([family_lambda(PAPER_K_MIN), family_mu(PAPER_K_MIN)]))
    ok = got == fam
    # "up to interchanging the factors": both orderings must fail, and the
    # sorting maps must be insensitive to the order of the two arguments.
    a, b = PAPER_MINIMAL_PAIR
    for (x, y) in ((a, b), (b, a)):
        if sort1(x, y) != sort1(a, b) or sort2(x, y) != sort2(a, b):
            ok = False
        if not failures_for_pair(x, y, PAPER_UNIQUENESS_TOTAL):
            ok = False
    # symmetry of the reduced coefficient in its two lower indices
    for nu in partitions(2) + partitions(3):
        if reduced_kron(a, b, nu)[0] != reduced_kron(b, a, nu)[0]:
            ok = False
    return ck("minimal_pair_is_the_k_equals_%d_member" % PAPER_K_MIN, ok,
              "census pair %s equals {lambda,mu} at k=%d = %s, and both "
              "orderings fail" % (got, PAPER_K_MIN, fam))


def check_extended_sweep():
    """Wider exhaustive census, |lambda|+|mu| <= SWEEP_TOTAL: list every
    failing unordered pair, confirm none has total size <= 3, confirm the
    total-4 answer is unchanged, and confirm the k=3 member appears."""
    failing = []
    for (lam, mu) in unordered_pairs(SWEEP_TOTAL):
        fs = failures_for_pair(lam, mu, sum(lam) + sum(mu))
        if fs:
            failing.append((tuple(sorted([lam, mu])), fs[0]))
    keys = [p for p, _ in failing]
    ok = True
    below = PAPER_UNIQUENESS_TOTAL - 1
    for p in keys:
        if sum(p[0]) + sum(p[1]) <= below:
            ok = False
    small = [p for p in keys if sum(p[0]) + sum(p[1]) <= PAPER_UNIQUENESS_TOTAL]
    if small != [tuple(sorted([norm(PAPER_MINIMAL_PAIR[0]),
                               norm(PAPER_MINIMAL_PAIR[1])]))]:
        ok = False
    kwit = PAPER_K_MIN + 1
    fam3 = tuple(sorted([family_lambda(kwit), family_mu(kwit)]))
    if fam3 not in keys:
        ok = False
    print("  failing pairs up to total %d: %s" % (SWEEP_TOTAL, keys))
    return ck("extended_census_total_at_most_%d" % SWEEP_TOTAL, ok,
              "%d unordered pairs searched, %d failing, none of total size "
              "<= %d, k=%d member %s present"
              % (len(unordered_pairs(SWEEP_TOTAL)), len(keys), below, kwit,
                 fam3))


def main():
    print("Counterexamples to a sorting conjecture for reduced Kronecker "
          "coefficients: verification")
    print("All coefficients are computed from symmetric-group characters "
          "with exact integer/rational arithmetic.")
    plan = [
        ("character_machinery_self_test", check_character_machinery),
        ("stable_padding_definition", check_padding_definition),
        ("branching_identity_eq_branching", check_branching_identity),
        ("exhibited_object_wellformed", check_exhibited_object),
        ("sorting_maps_recomputed", check_sorting_maps),
        ("hypotheses_of_the_conjecture_hold", check_hypotheses),
        ("reduced_coefficients_from_characters",
         check_reduced_values_from_characters),
        ("murnaghan_stability_observed", check_murnaghan_stability_observed),
        ("conjectured_inequality_violated", check_violation_of_conclusion),
        ("lemma_closed_formula_matches_characters", check_lemma_one),
        ("lemma_proof_steps_recomputed", check_lemma_proof_steps),
        ("infinite_family_via_validated_lemma",
         check_infinite_family_via_lemma),
        ("nu_range_truncation_is_safe", check_nu_truncation),
        ("uniqueness_census_total_at_most_%d" % PAPER_UNIQUENESS_TOTAL,
         check_uniqueness_census),
        ("minimality_proof_steps_recomputed", check_corollary_proof_steps),
        ("columns_pair_satisfies_inequality", check_columns_pair_holds),
        ("minimal_pair_is_the_k_equals_%d_member" % PAPER_K_MIN,
         check_minimal_pair_is_family_member),
        ("extended_census_total_at_most_%d" % SWEEP_TOTAL,
         check_extended_sweep),
    ]
    for name, fn in plan:
        before = len(CHECKS)
        try:
            fn()
        except Exception as exc:            # a corrupted input must not
            if len(CHECKS) == before:       # silence the verdict line
                ck(name, False, "raised %s: %s"
                   % (type(exc).__name__, str(exc)[:120]))
    print("NOT RE-RUN HERE: (i) the reduced Kronecker coefficients are "
          "limits, and stability is OBSERVED over finitely many N rather than "
          "proved -- the two load-bearing values are accepted for k=%d..%d "
          "once %d consecutive admissible N agree, and the longer %d-term "
          "constant window is exhibited only for k=%d..%d, so no window "
          "longer than %d terms is shown for k=%d..%d; (ii) the "
          "character-based confirmation of the family covers "
          "k=%d..%d, and k=%d..%d rests on the closed formula for nu=(1), "
          "which is itself validated against characters on all pairs of "
          "partitions of size <= %d; (iii) the exhaustive census over pairs "
          "runs to total size %d, beyond the total size %d needed for the "
          "minimality claim, but its nu range is TRUNCATED at "
          "|nu| <= |lambda|+|mu|: that truncation is supported only by the "
          "evidence that every coefficient with |nu| in [total+1,total+%d] "
          "vanishes, and |nu| > total+%d is untested; (iv) no external "
          "catalogue or table of "
          "Kronecker coefficients is consulted; (v) the multipartition "
          "corollary is NOT verified here beyond its two-factor content: it "
          "rests on identifying the two partitions written (lambda cup "
          "mu)^[1,2] and (lambda cup mu)^[2,2] in the cited multipartition "
          "conjecture with (sort1,sort2), which is a reading of an external "
          "paper's notation and is TAKEN ON TRUST -- no check above tests "
          "it, because coding the r-block index-residue split and comparing "
          "its r=2 output with sort1,sort2 would only compare one "
          "transcription of a single definition with another and could not "
          "fail; (vi) nothing here consults the cited paper itself (no "
          "network access and no external file is read), so the numbering "
          "and wording of the external statements quoted as Conjecture 5.5, "
          "Theorem 5.6 and Conjecture 5.7 of its published version are TAKEN "
          "ON TRUST, as is the attribution of the pair %s/%s to that "
          "Theorem 5.6 -- the inequality for that pair is nevertheless "
          "computed here, only its attribution is not."
          % (K_CHAR_RANGE[0], K_CHAR_RANGE[-1], STAB_RUNS_LOAD,
             STAB_RUNS_WINDOW, K_STAB_WINDOW_RANGE[0], K_STAB_WINDOW_RANGE[-1],
             STAB_RUNS_LOAD, K_STAB_WINDOW_RANGE[-1] + 1, K_CHAR_RANGE[-1],
             K_CHAR_RANGE[0], K_CHAR_RANGE[-1], K_CHAR_RANGE[-1] + 1,
             K_LEMMA_RANGE[-1], LEMMA_SIZE_CAP, SWEEP_TOTAL,
             PAPER_UNIQUENESS_TOTAL, NU_OVERSHOOT, NU_OVERSHOOT,
             norm(PAPER_COLUMNS_PAIR[0]), norm(PAPER_COLUMNS_PAIR[1])))
    finish()


if __name__ == "__main__":
    sys.setrecursionlimit(20000)
    main()
