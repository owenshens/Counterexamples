#!/usr/bin/env python3
"""
Verification of the claim ccn(chi_(3,3); S_6) = ccn(chi_(2,2,2); S_6) = 5,
and of the consequent failure of two published bounds at n = 6.

Here ccn(theta; G) is the least positive integer k with c(theta^k) = Irr(G),
where theta^k is the pointwise k-th power (the character of the k-fold tensor
power) and c(.) is the set of irreducible constituents.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs; nothing here is computed by the paper's
own machinery, every one of them is re-derived below and compared):

  * n = 6, and the two partitions (3,3) and (2,2,2) = (3,3)'.
  * The exhibited class-function table, transcribed verbatim, in the paper's
    column order 1^6, 21^4, 2^2 1^2, 2^3, 31^3, 321, 3^2, 41^2, 42, 51, 6:
        |C_mu|            1 15 45 15 40 120 40 90 90 144 120
        chi_(3,3)(mu)     5  1  1 -3 -1   1  2 -1 -1   0   0
        eps(mu)           1 -1  1 -1  1  -1  1 -1  1   1  -1
  * The claimed Young-rule decompositions
        pi_3 = chi_(6)+chi_(5,1)+chi_(4,2)+chi_(3,3),
        pi_2 = chi_(6)+chi_(5,1)+chi_(4,2),
    where pi_r is the permutation character on r-subsets, and the claimed
    consequence chi_(3,3) = pi_3 - pi_2.
  * The displayed arithmetic
        720*<chi^4, eps> = 625-15+45-1215+40-120+640-90+90 = 0.
  * The claimed value 5 for both covering numbers.
  * The two statements acted on:
      Conjecture 36(3), upper half: for even n >= 5,
          ccn(chi_(n/2,n/2); S_n) <= ceil(log_2 n) + 1.
      Theorem 1: for n >= 5 and lambda not in
          {(n), (1^n), (n-1,1), (2,1^{n-2})},
          ccn(chi_lambda; S_n) <= ceil(2(n-1)/3).
  * The cited external input (Miller): for n > 4 every nonlinear
    psi in Irr(S_n) has c(psi^{n-1}) = Irr(S_n).  Used by the paper only for
    the upper bound; it is re-verified here at n = 6 by direct computation,
    so no check below depends on trusting it.

DERIVED HERE (computed from scratch, exact integer arithmetic only):

  * The 11 partitions of 6, their conjugacy class sizes, and the full
    11 x 11 character table of S_6, built by the Murnaghan-Nakayama rule
    and independently rebuilt by Jacobi-Trudi applied to the fixed-point
    counts of the actions on ordered set partitions; both are checked
    against row and column orthogonality.
  * The permutation characters pi_2, pi_3 from the generating function
    pi_r(mu) = [x^r] prod_i (1 + x^{mu_i}), their decompositions, and
    pi_3 - pi_2.
  * All multiplicities <chi^k, psi> for k = 1..5 and the exact sets of
    irreducible constituents of chi^k, hence the covering numbers
    themselves; the irreducibles MISSING from chi^k for k <= 4 are named.
  * ceil(log_2 6) + 1 and ceil(2*5/3), by integer arithmetic.
  * A census: ccn(chi_lambda; S_n) for every partition lambda of every
    n in 5..11, and the complete list of pairs (n, lambda) violating the
    Theorem 1 bound over that range.
--------------------------------------------------------------------------
"""

import sys
from itertools import permutations
from fractions import Fraction

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return bool(ok)


# ----------------------------------------------------------------------
# partitions, classes
# ----------------------------------------------------------------------

def partitions(n, cap=None):
    """All partitions of n as weakly decreasing tuples, reverse-lex order."""
    if cap is None:
        cap = n
    if n == 0:
        return [()]
    out = []
    for first in range(min(n, cap), 0, -1):
        for rest in partitions(n - first, first):
            out.append((first,) + rest)
    return out


def factorial(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def class_size(mu, n):
    """|C_mu| = n! / (prod mu_i * prod_j m_j!) with m_j the part multiplicities."""
    den = 1
    for p in mu:
        den *= p
    mult = {}
    for p in mu:
        mult[p] = mult.get(p, 0) + 1
    for m in mult.values():
        den *= factorial(m)
    num = factorial(n)
    assert num % den == 0
    return num // den


def conjugate(lam):
    if not lam:
        return ()
    return tuple(sum(1 for p in lam if p > j) for j in range(lam[0]))


def parse_label(s):
    """Cycle-type label -> partition, e.g. '2^21^2' -> (2,2,1,1), '21^4' ->
    (2,1,1,1,1).  Single-digit bases and exponents, which is all that occurs
    for n = 6; a misreading would be caught by the partition-set check."""
    out = []
    i = 0
    while i < len(s):
        base = int(s[i])
        i += 1
        exp = 1
        if i < len(s) and s[i] == "^":
            exp = int(s[i + 1])
            i += 2
        out.extend([base] * exp)
    return tuple(sorted(out, reverse=True))


# ----------------------------------------------------------------------
# character table by Murnaghan-Nakayama (rim hooks via first-column hooks)
# ----------------------------------------------------------------------

_MN = {}


def mn_char(lam, mu):
    """chi_lam(mu) by the Murnaghan-Nakayama rule; exact integer."""
    key = (lam, mu)
    if key in _MN:
        return _MN[key]
    if not mu:
        val = 1 if not lam else 0
        _MN[key] = val
        return val
    r = mu[0]
    rest = mu[1:]
    total = 0
    ell = len(lam)
    beta = [lam[i] + (ell - 1 - i) for i in range(ell)]
    bset = set(beta)
    for i in range(ell):
        new = beta[i] - r
        if new < 0 or new in bset:
            continue
        height = sum(1 for b in beta if new < b < beta[i])
        nb = sorted([new if k == i else beta[k] for k in range(ell)],
                    reverse=True)
        nu = tuple(nb[k] - (ell - 1 - k) for k in range(ell))
        while nu and nu[-1] == 0:
            nu = nu[:-1]
        total += (-1) ** height * mn_char(nu, rest)
    _MN[key] = total
    return total


def char_table(n):
    """(parts, sizes, table) with table[a][b] = chi_{parts[a]}(parts[b])."""
    parts = partitions(n)
    sizes = [class_size(mu, n) for mu in parts]
    table = [[mn_char(lam, mu) for mu in parts] for lam in parts]
    return parts, sizes, table


def inner(f, g, sizes, order):
    """<f,g> = (1/|G|) sum_mu |C_mu| f(mu) g(mu); returns an exact Fraction."""
    s = 0
    for i, sz in enumerate(sizes):
        s += sz * f[i] * g[i]
    return Fraction(s, order)


# ----------------------------------------------------------------------
# independent character table: Jacobi-Trudi over induced-trivial characters
# ----------------------------------------------------------------------

_H = {}


def h_value(comp, mu):
    """Value at class mu of Ind_{S_comp}^{S_n} 1, i.e. the number of ordered
    set partitions of [n] with block sizes comp fixed by a permutation of
    cycle type mu: the number of ways to assign the cycles of mu to the
    blocks so that block i receives cycles of total length comp_i."""
    if any(c < 0 for c in comp):
        return 0
    caps = tuple(c for c in comp if c > 0)
    if sum(caps) != sum(mu):
        return 0
    key = (tuple(sorted(caps)), mu)
    if key in _H:
        return _H[key]

    def rec(idx, rem):
        if idx == len(mu):
            return 1 if all(r == 0 for r in rem) else 0
        cyc = mu[idx]
        tot = 0
        seen = set()
        for j, r in enumerate(rem):
            if r >= cyc and (r, j) not in seen:
                seen.add((r, j))
                nxt = list(rem)
                nxt[j] = r - cyc
                tot += rec(idx + 1, tuple(nxt))
        return tot

    val = rec(0, caps)
    _H[key] = val
    return val


def jt_char_row(lam, classes):
    """chi_lam by Jacobi-Trudi: s_lam = det(h_{lam_i - i + j})."""
    ell = len(lam)
    row = []
    for mu in classes:
        tot = 0
        for sigma in permutations(range(ell)):
            sign = 1
            for a in range(ell):
                for b in range(a + 1, ell):
                    if sigma[a] > sigma[b]:
                        sign = -sign
            comp = tuple(lam[i] - (i + 1) + (sigma[i] + 1) for i in range(ell))
            if any(c < 0 for c in comp):
                continue
            tot += sign * h_value(comp, mu)
        row.append(tot)
    return row


# ----------------------------------------------------------------------
# constituents of tensor powers and covering numbers
# ----------------------------------------------------------------------

def multiplicities(chi, k, sizes, table, order):
    """Exact list of <chi^k, chi_lam> over all lam, as Fractions."""
    pw = [c ** k for c in chi]
    return [inner(pw, row, sizes, order) for row in table]


def constituents(chi, k, sizes, table, order):
    """(frozenset of indices of irreducible constituents of chi^k,
        True iff every multiplicity came out a nonnegative integer)."""
    mult = multiplicities(chi, k, sizes, table, order)
    integral = all(m.denominator == 1 and m >= 0 for m in mult)
    return frozenset(i for i, m in enumerate(mult) if m != 0), integral


def covering_number(chi, sizes, table, order, kmax):
    """Least k <= kmax with c(chi^k) = Irr, else None; plus per-k missing sets."""
    full = frozenset(range(len(table)))
    missing = []
    found = None
    for k in range(1, kmax + 1):
        cons, integral = constituents(chi, k, sizes, table, order)
        if not integral:
            return ("NONINTEGRAL", missing)
        missing.append(full - cons)
        if cons == full and found is None:
            found = k
    return (found, missing)


def ceil_log2(n):
    """ceil(log_2 n) for n >= 1 by integer arithmetic only."""
    k = 0
    v = 1
    while v < n:
        v *= 2
        k += 1
    return k


def ceil_div(a, b):
    return -((-a) // b)


def covers_at(chi, k, sizes, table, order):
    """Fast integer path: (covers_all, all_multiplicities_nonneg_integers)."""
    pw = [c ** k for c in chi]
    w = [sz * p for sz, p in zip(sizes, pw)]
    ok = True
    nzero = 0
    for row in table:
        s = 0
        for a, wa in enumerate(w):
            if wa:
                s += wa * row[a]
        if s % order != 0 or s // order < 0:
            ok = False
        if s != 0:
            nzero += 1
    return (nzero == len(table), ok)


def ccn_fast(chi, sizes, table, order, kmax):
    for k in range(1, kmax + 1):
        full, ok = covers_at(chi, k, sizes, table, order)
        if not ok:
            return None
        if full:
            return k
    return None


# ----------------------------------------------------------------------
# values transcribed from the paper
# ----------------------------------------------------------------------

N = 6
ORDER = 720
PAPER_LABELS = ["1^6", "21^4", "2^21^2", "2^3", "31^3", "321", "3^2",
                "41^2", "42", "51", "6"]
PAPER_SIZES = [1, 15, 45, 15, 40, 120, 40, 90, 90, 144, 120]
PAPER_CHI33 = [5, 1, 1, -3, -1, 1, 2, -1, -1, 0, 0]
PAPER_EPS = [1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1]
PAPER_PI3 = [(6,), (5, 1), (4, 2), (3, 3)]
PAPER_PI2 = [(6,), (5, 1), (4, 2)]
PAPER_TERMS_X4EPS = [625, -15, 45, -1215, 40, -120, 640, -90, 90]
PAPER_TERM_CLASSES = ["1^6", "21^4", "2^21^2", "2^3", "31^3", "321", "3^2",
                      "41^2", "42"]
PAPER_CCN = 5
PAPER_CONJ363_RHS = 4
PAPER_THM1_RHS = 4
PAPER_THM1_EXCLUDED = [(6,), (1, 1, 1, 1, 1, 1), (5, 1), (2, 1, 1, 1, 1)]


# ----------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------

def check_exhibited_object(parts, sizes, table):
    """(1) decode the exhibited table, count it, print it back."""
    dec = [parse_label(s) for s in PAPER_LABELS]
    ok = (len(dec) == 11 and len(set(dec)) == 11
          and all(sum(d) == N for d in dec)
          and set(dec) == set(parts))
    ck("exhibited_class_labels_decode_to_the_11_partitions_of_6", ok,
       "decoded: " + " ".join("+".join(str(x) for x in d) for d in dec))
    col = [parts.index(d) for d in dec]

    got = [sizes[c] for c in col]
    ck("exhibited_class_sizes_match_computed_and_sum_to_720",
       got == PAPER_SIZES and sum(got) == ORDER,
       "computed " + " ".join(map(str, got)))

    eps_idx = parts.index(tuple([1] * N))
    eps = [table[eps_idx][c] for c in col]
    par = [(-1) ** (N - len(dec[i])) for i in range(11)]
    ck("exhibited_sign_row_is_the_sign_character",
       eps == PAPER_EPS and par == PAPER_EPS,
       "computed " + " ".join("%+d" % v for v in eps))

    i33 = parts.index((3, 3))
    chi = [table[i33][c] for c in col]
    ck("exhibited_chi_33_row_matches_computed_character_table",
       chi == PAPER_CHI33,
       "computed " + " ".join("%+d" % v for v in chi))
    return col


def check_table_sound(parts, sizes, table):
    """The computed table really is the character table of S_6."""
    ok_rows = True
    for a in range(len(parts)):
        for b in range(len(parts)):
            v = inner(table[a], table[b], sizes, ORDER)
            if v != (1 if a == b else 0):
                ok_rows = False
    ck("computed_table_rows_are_orthonormal", ok_rows,
       "%dx%d inner products" % (len(parts), len(parts)))

    ok_cols = True
    for a in range(len(parts)):
        for b in range(len(parts)):
            s = sum(table[r][a] * table[r][b] for r in range(len(parts)))
            want = ORDER // sizes[a] if a == b else 0
            if s != want:
                ok_cols = False
    ck("computed_table_columns_are_orthogonal", ok_cols,
       "column norms equal |G|/|C_mu|")

    rebuilt = [jt_char_row(lam, parts) for lam in parts]
    ck("second_independent_construction_reproduces_the_table",
       rebuilt == table,
       "Jacobi-Trudi over fixed-point counts vs Murnaghan-Nakayama")

    ident = parts.index(tuple([1] * N))
    degs = [table[a][ident] for a in range(len(parts))]
    i33, i222 = parts.index((3, 3)), parts.index((2, 2, 2))
    ck("degrees_square_sum_to_720_and_both_degrees_are_5",
       sum(d * d for d in degs) == ORDER and degs[i33] == 5
       and degs[i222] == 5,
       "degrees " + " ".join(map(str, degs)))


def check_youngs_rule(parts, sizes, table, col):
    """The paper's route to chi_(3,3): pi_3 - pi_2 with Young's rule."""
    def gen_coeff(mu, r):
        poly = [1]
        for m in mu:
            new = [0] * (len(poly) + m)
            for d, c in enumerate(poly):
                new[d] += c
                new[d + m] += c
            poly = new
        return poly[r] if r < len(poly) else 0

    pi = {}
    ok_gf = True
    for r in (2, 3):
        by_gf = [gen_coeff(mu, r) for mu in parts]
        by_ind = [h_value((r, N - r), mu) for mu in parts]
        if by_gf != by_ind:
            ok_gf = False
        pi[r] = by_gf
    ck("pi_r_generating_function_equals_the_r_subset_permutation_character",
       ok_gf,
       "pi_3 = " + " ".join(str(pi[3][c]) for c in col))

    ok_dec = True
    detail = []
    for r, claimed in ((3, PAPER_PI3), (2, PAPER_PI2)):
        mult = [inner(pi[r], row, sizes, ORDER) for row in table]
        got = sorted([parts[i] for i, m in enumerate(mult) if m != 0])
        if got != sorted(claimed) or any(m not in (0, 1) for m in mult):
            ok_dec = False
        detail.append("pi_%d has %d constituents" % (r, len(got)))
    ck("youngs_rule_decompositions_of_pi_3_and_pi_2_are_as_stated", ok_dec,
       "; ".join(detail))

    diff = [pi[3][i] - pi[2][i] for i in range(len(parts))]
    i33 = parts.index((3, 3))
    ck("pi_3_minus_pi_2_equals_chi_33_on_every_class",
       diff == table[i33] and [diff[c] for c in col] == PAPER_CHI33,
       "difference " + " ".join("%+d" % diff[c] for c in col))


def check_paper_arithmetic(parts, sizes, col):
    """The paper's displayed vanishing inner products, from its OWN row."""
    inv = {}
    for k, c in enumerate(col):
        inv[c] = k
    chi = [PAPER_CHI33[inv[i]] for i in range(len(parts))]
    eps = [PAPER_EPS[inv[i]] for i in range(len(parts))]
    one = [1] * len(parts)

    ck("transcribed_row_has_norm_one_so_it_is_an_irreducible_character",
       inner(chi, chi, sizes, ORDER) == 1
       and chi[parts.index(tuple([1] * N))] == 5,
       "<chi,chi> = %s, degree %d"
       % (inner(chi, chi, sizes, ORDER), chi[parts.index(tuple([1] * N))]))

    vals = [inner(chi, one, sizes, ORDER),
            inner([c * c for c in chi], eps, sizes, ORDER),
            inner([c ** 3 for c in chi], one, sizes, ORDER),
            inner([c ** 4 for c in chi], eps, sizes, ORDER)]
    ck("four_displayed_inner_products_all_vanish",
       all(v == 0 for v in vals),
       "<chi,1>=%s <chi^2,eps>=%s <chi^3,1>=%s <chi^4,eps>=%s"
       % tuple(str(v) for v in vals))

    terms = []
    for lab in PAPER_TERM_CLASSES:
        k = PAPER_LABELS.index(lab)
        terms.append(PAPER_SIZES[k] * PAPER_CHI33[k] ** 4 * PAPER_EPS[k])
    rest = [PAPER_SIZES[k] * PAPER_CHI33[k] ** 4 * PAPER_EPS[k]
            for k in range(11) if PAPER_LABELS[k] not in PAPER_TERM_CLASSES]
    ck("displayed_expansion_of_720_times_chi4_eps_is_term_by_term_correct",
       terms == PAPER_TERMS_X4EPS and all(t == 0 for t in rest)
       and sum(terms) == 0,
       "terms " + " ".join("%+d" % t for t in terms) + " sum 0; "
       + "omitted classes contribute " + " ".join("%+d" % t for t in rest))


def check_covering_numbers(parts, sizes, table):
    """The theorem itself: ccn = 5 for BOTH (3,3) and (2,2,2), derived from the
    computed table, with the missing irreducibles named for every k <= 4."""
    full = frozenset(range(len(parts)))
    results = {}
    for lam in ((3, 3), (2, 2, 2)):
        idx = parts.index(lam)
        chi = table[idx]

        found, missing = covering_number(chi, sizes, table, ORDER, PAPER_CCN + 1)
        ck("all_multiplicities_of_chi_%s_powers_are_nonneg_integers"
           % "".join(map(str, lam)),
           found != "NONINTEGRAL",
           "exact Fraction inner products, k = 1..%d" % (PAPER_CCN + 1))
        if found == "NONINTEGRAL":
            continue

        for k in range(1, PAPER_CCN):
            miss = missing[k - 1]
            names = " ".join("(" + ",".join(map(str, parts[i])) + ")"
                             for i in sorted(miss))
            ck("chi_%s_power_%d_omits_at_least_one_irreducible"
               % ("".join(map(str, lam)), k),
               len(miss) > 0,
               "%d missing: %s" % (len(miss), names))

        miss5 = missing[PAPER_CCN - 1]
        ck("chi_%s_power_%d_contains_every_irreducible"
           % ("".join(map(str, lam)), PAPER_CCN),
           len(miss5) == 0 and (full - miss5) == full,
           "all %d irreducibles of S_6 occur" % len(parts))

        fast = ccn_fast(chi, sizes, table, ORDER, PAPER_CCN + 1)
        ck("ccn_chi_%s_equals_the_paper_value_by_two_independent_code_paths"
           % "".join(map(str, lam)),
           found == PAPER_CCN and fast == PAPER_CCN,
           "Fraction path %s, integer path %s, paper says %d"
           % (found, fast, PAPER_CCN))
        results[lam] = found

    ck("the_two_covering_numbers_are_equal_as_the_theorem_asserts",
       len(results) == 2 and results[(3, 3)] == results[(2, 2, 2)],
       "ccn(chi_33) = %s, ccn(chi_222) = %s"
       % (results.get((3, 3)), results.get((2, 2, 2))))


def check_sign_twist(parts, sizes, table):
    """The paper's last paragraph: (3,3)' = (2,2,2), chi_222 = chi_33 * eps,
    (chi eps)^k = chi^k eps^k, and tensoring by eps permutes Irr(S_6)."""
    eps_idx = parts.index(tuple([1] * N))
    eps = table[eps_idx]

    ck("conjugate_of_33_is_222",
       conjugate((3, 3)) == (2, 2, 2) and conjugate((2, 2, 2)) == (3, 3),
       "(3,3)' = " + str(conjugate((3, 3))))

    i33, i222 = parts.index((3, 3)), parts.index((2, 2, 2))
    twist = [a * b for a, b in zip(table[i33], eps)]
    ck("chi_222_equals_chi_33_times_the_sign_character",
       twist == table[i222],
       "on all %d classes" % len(parts))

    perm = {}
    ok_perm = True
    for a, lam in enumerate(parts):
        t = [v * e for v, e in zip(table[a], eps)]
        hits = [b for b in range(len(parts)) if table[b] == t]
        if len(hits) != 1:
            ok_perm = False
            continue
        perm[a] = hits[0]
        if hits[0] != parts.index(conjugate(lam)):
            ok_perm = False
    ck("tensoring_by_eps_permutes_Irr_S6_and_is_conjugation_of_partitions",
       ok_perm and sorted(perm.values()) == list(range(len(parts))),
       "involution on %d irreducibles" % len(perm))

    # (chi eps)^k = chi^k eps^k is a ring identity and holds for ANY inputs, so
    # on its own it is vacuous.  The paper's actual claim is the consequence:
    # tensoring by eps^k carries c(chi^k) bijectively onto c((chi eps)^k), i.e.
    # onto the conjugate partitions for odd k and onto itself for even k.  That
    # is computed here class by class, so a corrupted table breaks it.
    ok_pw = True
    detail_pw = []
    for k in range(1, PAPER_CCN + 1):
        lhs = [(v * e) ** k for v, e in zip(table[i33], eps)]
        rhs = [(v ** k) * (e ** k) for v, e in zip(table[i33], eps)]
        if lhs != rhs:
            ok_pw = False
        c33, int33 = constituents(table[i33], k, sizes, table, ORDER)
        c222, int222 = constituents(table[i222], k, sizes, table, ORDER)
        if not (int33 and int222):
            ok_pw = False
        if k % 2 == 1:
            want = frozenset(parts.index(conjugate(parts[a])) for a in c33)
        else:
            want = c33
        if c222 != want:
            ok_pw = False
        detail_pw.append("k=%d: %d vs %d constituents"
                         % (k, len(c33), len(c222)))
    ck("constituents_of_chi_eps_power_k_are_exactly_the_eps_twists_of_those_"
       "of_chi_power_k", ok_pw,
       "k = 1..%d; " % PAPER_CCN + "; ".join(detail_pw))


def check_miller_input_at_n6(parts, sizes, table):
    """The paper's only external input, re-derived at n = 6 so that no check
    here depends on trusting Miller: every nonlinear psi in Irr(S_6) has
    c(psi^5) = Irr(S_6)."""
    ident = parts.index(tuple([1] * N))
    full = frozenset(range(len(parts)))
    nonlinear = 0
    ok = True
    for a in range(len(parts)):
        if table[a][ident] == 1:
            continue
        nonlinear += 1
        cons, integral = constituents(table[a], N - 1, sizes, table, ORDER)
        if not integral or cons != full:
            ok = False
    ck("miller_bound_reverified_directly_every_nonlinear_irr_covers_at_power_5",
       ok and nonlinear == len(parts) - 2,
       "%d nonlinear irreducibles of S_6, all cover at k = %d"
       % (nonlinear, N - 1))


def check_published_bounds():
    """The corollary's arithmetic and its exclusion clause, at n = 6."""
    conj = ceil_log2(N) + 1
    ck("conjecture_36_3_upper_rhs_at_n6_is_4_by_integer_arithmetic",
       conj == PAPER_CONJ363_RHS,
       "ceil(log_2 %d) + 1 = %d" % (N, conj))

    thm = ceil_div(2 * (N - 1), 3)
    ck("theorem_1_rhs_at_n6_is_4_by_integer_arithmetic",
       thm == PAPER_THM1_RHS,
       "ceil(2*%d/3) = %d" % (N - 1, thm))

    excluded = [(N,), tuple([1] * N), (N - 1, 1), (2,) + tuple([1] * (N - 2))]
    ck("theorem_1_exclusion_set_instantiated_at_n6_is_the_papers_list",
       excluded == PAPER_THM1_EXCLUDED,
       " ".join("(" + ",".join(map(str, e)) + ")" for e in excluded))
    ck("neither_33_nor_222_is_excluded_by_theorem_1",
       (3, 3) not in excluded and (2, 2, 2) not in excluded,
       "both partitions of 6 are inside the theorem's scope")

    ck("the_derived_value_5_strictly_exceeds_both_published_bounds",
       PAPER_CCN > conj and PAPER_CCN > thm,
       "5 > %d (Conjecture 36(3) upper) and 5 > %d (Theorem 1)"
       % (conj, thm))


def check_census_5_to_11():
    """The remark's claim that a direct verification over n <= 11 cannot be
    correct: compute ccn for every partition of every n in 5..11 and list the
    complete set of Theorem 1 violations over that range."""
    viol = []
    linear_never_cover = True
    all_integral = True
    tables_orthonormal = True
    per_n = []
    for n in range(5, 12):
        parts, sizes, table = char_table(n)
        order = factorial(n)
        ident = parts.index(tuple([1] * n))
        bound = ceil_div(2 * (n - 1), 3)
        excluded = set([(n,), tuple([1] * n), (n - 1, 1),
                        (2,) + tuple([1] * (n - 2))])
        # The census tables carry the "only two violations" claim, so they are
        # verified to be genuine character tables rather than assumed to be.
        for ra in range(len(parts)):
            for rb in range(ra, len(parts)):
                acc = 0
                for cc in range(len(parts)):
                    acc += sizes[cc] * table[ra][cc] * table[rb][cc]
                if acc != (order if ra == rb else 0):
                    tables_orthonormal = False
        nv = 0
        for a, lam in enumerate(parts):
            for k in range(1, n):
                _, ok_int = covers_at(table[a], k, sizes, table, order)
                if not ok_int:
                    all_integral = False
            c = ccn_fast(table[a], sizes, table, order, n - 1)
            if table[a][ident] == 1:
                if c is not None:
                    linear_never_cover = False
                continue
            if c is None or c > bound:
                viol.append((n, lam, c, bound, lam in excluded))
                nv += 1
        per_n.append("n=%d bound=%d violations=%d" % (n, bound, nv))

    ck("census_character_tables_for_n_5_to_11_are_orthonormal",
       tables_orthonormal,
       "every census table verified against row orthonormality before use")
    ck("census_all_tensor_power_multiplicities_are_nonneg_integers",
       all_integral, "n = 5..11, k = 1..n-1")
    ck("census_linear_characters_never_cover_so_ccn_is_undefined_for_them",
       linear_never_cover,
       "trivial and sign character omitted from the census as in Theorem 1")

    unexcluded = sorted((n, lam) for n, lam, c, b, ex in viol if not ex)
    ck("census_over_n_5_to_11_finds_violations_of_the_theorem_1_bound",
       len(viol) > 0,
       "; ".join(per_n))
    # Stated separately from the exactness claim below: the paper asserts only
    # that (3,3) and (2,2,2) violate the bound at n = 6, never that they are the
    # only violations over 5..11.  If the stronger claim below ever fails, this
    # one localises the failure and shows the paper's corollary still stands.
    ck("the_census_reproduces_the_corollarys_two_counterexamples_at_n_6",
       (N, (3, 3)) in unexcluded and (N, (2, 2, 2)) in unexcluded,
       "both partitions of 6 violate ceil(2(n-1)/3) = %d"
       % ceil_div(2 * (N - 1), 3))
    ck("the_only_unexcluded_theorem_1_violations_in_n_5_to_11_are_the_papers_two",
       unexcluded == [(6, (2, 2, 2)), (6, (3, 3))],
       "unexcluded: " + " ".join("n=%d (%s)" % (n, ",".join(map(str, l)))
                                 for n, l in unexcluded))
    print("     census detail, every (n, lambda) with ccn > ceil(2(n-1)/3):")
    for n, lam, c, b, ex in viol:
        print("       n=%2d lambda=(%-22s) ccn=%s bound=%d  %s"
              % (n, ",".join(map(str, lam)), c, b,
                 "EXCLUDED by Theorem 1" if ex
                 else "NOT excluded -> counterexample"))


# Assertions of the paper that no check above tests.  Each entry is
# (one-clause tag, block of lines).  The SCOPE block below prints the blocks,
# the closing NOT RE-RUN line is assembled from the tags, and the count in
# both is len() of this list, so the two cannot drift apart.
SCOPE_NOT_TESTED = [
    ("the abstract's priority claim to Sun-Zhang-Zhu, a statement about the "
     "literature",
     ["The abstract's priority claim, that the value 5 is already",
      "implicit in the power formulas of Sun-Zhang-Zhu, is a",
      "statement about the literature and is not tested here."]),
    ("the remark's claim that the LOWER inequality of Conjecture 36(3) still "
     "holds at n = 6, which the paper never reproduces",
     ["The remark's claim that the LOWER inequality of Conjecture",
      "36(3) remains valid at n = 6 is not tested: the paper does",
      "not reproduce the lower bound, so there is no expression to",
      "evaluate.  Nothing above depends on it."]),
    ("anything beyond n = 11, the census range, so the remark's even n >= 8 "
     "disclaimer is untested",
     ["The remark's disclaimer that the counterexample does not",
      "address even n >= 8 is a non-claim; the census above stops",
      "at n = 11 and asserts nothing beyond it."]),
    ("whether Kundu and Velmurugan's SageMath verification was in fact "
     "carried out, which is unobservable from outside their paper",
     ["That Kundu and Velmurugan's SageMath run in fact produced",
      "the bound for n <= 11 cannot be checked from outside their",
      "paper.  The census shows only that the bound itself is",
      "false at n = 6, hence that any verification affirming it",
      "over n <= 11 must be in error."]),
    ("Miller's general n > 4 theorem, re-derived here only at n = 6, the "
     "sole case any check uses",
     ["Miller's theorem is used here only at n = 6, where it is",
      "re-derived by direct computation; the general n > 4",
      "statement is not tested and no check relies on it."]),
    ("the faithfulness of the quotations of Conjecture 36(3) and of Kundu "
     "and Velmurugan's Theorem 1, which are transcribed here, not read from "
     "the source",
     ["The WORDING of the two statements acted on -- Conjecture",
      "36(3) and Theorem 1 of Kundu and Velmurugan, including the",
      "latter's exclusion set -- is transcribed from their paper.",
      "This program reads no literature, so it instantiates those",
      "quotations at n = 6 and compares them with the derived",
      "value 5; it cannot confirm that the quotations are faithful.",
      "Everything else in the corollary is derived here."]),
]


def print_scope():
    """State plainly what the paper asserts that this program does NOT test,
    so the coverage is not silently narrower than the paper."""
    print("")
    print("SCOPE: the %d assertions of the paper NOT covered by the checks "
          "above" % len(SCOPE_NOT_TESTED))
    for num, (_, para) in enumerate(SCOPE_NOT_TESTED, 1):
        print("  %d. %s" % (num, para[0]))
        for line in para[1:]:
            print("     " + line)
    print("")
    gaps = "; ".join("(%d) %s" % (num, tag) for num, (tag, _)
                     in enumerate(SCOPE_NOT_TESTED, 1))
    print("NOT RE-RUN: the paper's full text.  What is re-derived from "
          "scratch, with no paper value taken on trust, is Theorem 1 "
          "-- both covering numbers, from an S_6 character table "
          "built twice by independent routes (Murnaghan-Nakayama and "
          "Jacobi-Trudi) and checked against row and column orthogonality -- "
          "every arithmetic step its proof displays, and, for Corollary 2, "
          "the two right-hand sides at n = 6, the exclusion set, and the "
          "strict inequality 5 > 4.  The %d assertions listed under SCOPE "
          "above lie outside that: %s."
          % (len(SCOPE_NOT_TESTED), gaps))


def main():
    parts, sizes, table = char_table(N)
    ck("computed_group_order_is_720", factorial(N) == ORDER,
       "6! = %d" % factorial(N))
    col = check_exhibited_object(parts, sizes, table)
    check_table_sound(parts, sizes, table)
    check_youngs_rule(parts, sizes, table, col)
    check_paper_arithmetic(parts, sizes, col)
    check_covering_numbers(parts, sizes, table)
    check_sign_twist(parts, sizes, table)
    check_miller_input_at_n6(parts, sizes, table)
    check_published_bounds()
    check_census_5_to_11()
    print_scope()
    print("")


# Every check function above must run.  A run that executes fewer than this
# many checks has skipped a section, and "no failures" would then be
# meaningless rather than reassuring, so it is reported as a failure.
MIN_CHECKS = 40


if __name__ == "__main__":
    main()
    nfail = sum(1 for _, ok in CHECKS if not ok)
    ntot = len(CHECKS)
    if ntot < MIN_CHECKS:
        print("VERDICT: ONLY %d CHECKS RAN, EXPECTED AT LEAST %d "
              "- the run is incomplete and proves nothing"
              % (ntot, MIN_CHECKS))
        sys.exit(1)
    if nfail == 0:
        print("VERDICT: ALL %d CHECKS PASS" % ntot)
    else:
        print("VERDICT: %d OF %d CHECKS FAILED" % (nfail, ntot))
    sys.exit(0 if nfail == 0 else 1)
