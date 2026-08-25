#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verifier for the accompanying paper: a five-variable self-adhesivity
counterexample, i.e. a supermodular integer set function f on
N = {0,1,2,3,4} whose induced structural semi-graphoid M is
self-adhesive-at-every-overlap yet is NOT induced by a self-adhesive
supermodular function (Lemmas 2/4/5 and Remark 3 of the paper).

Python 3.9 compatible.  STANDARD LIBRARY ONLY.  Exact arithmetic only
(int / fractions.Fraction); no floating point anywhere.

=====================  WHAT IS TAKEN FROM THE PAPER  =====================
  W1  the 20 elementary statements of the model M          (witness)
  W2  the 32 values of the standardized set function f      (witness)
  W3  the linear certificate of Lemma 5: 10 LHS terms with
      coefficient +1, 11 RHS terms with the printed signed
      integer coefficients (+1 everywhere except 24|0 -> -1
      and 35|012 -> +2)                                    (witness)
  W4  the choice of overlap L = {0,1,2} and the copy
      b(0)=0 b(1)=1 b(2)=2 b(3)=5 b(4)=6                   (a choice)
  DEF the definitions: elementary statement, elementary difference,
      supermodularity, standardization, the block set [A,B|C], the
      canonical copy b_L, the adhesion set A_L, the elementary
      semi-graphoid rule and its closure sg(.), the marginal S|N.
  EXP expected-value constants used ONLY as comparison targets, never
      as inputs: the difference histogram (20,20,25,14,1); the 20 rows
      (|A_L|,|sg(A_L)|) of the Lemma-4 table; rank 19; ambient dim 26;
      face dim 7; Delta_f(04|3) = 1; |M| = 20.
  AUD six (|A_L|,|sg(A_L)|) pairs for overlaps the paper does NOT
      tabulate, as reported by a second, independent closure
      implementation that is NOT distributed with this bundle, so a
      reader cannot inspect it.  These are NOT paper data and are NOT
      part of any pass/fail decision: the comparison is printed as an
      unauditable implementation-vs-implementation remark only, and the
      only closure computation a reader can check is the one below.

=====================  WHAT THE PROGRAM DERIVES  ========================
  * all 80 elementary differences Delta_f(s) from the f table alone,
    their histogram, and min >= 0 (supermodularity)
  * the zero set M_f, and its two-sided set equality with the
    transcribed M  (so W1 is not merely trusted)
  * sg(M) = M on N  (M is semi-graphoid closed)
  * the dimension of the standardized space, by enumeration
  * the 20 x 26 integer matrix of active equalities, its exact rank
    over Q by Fraction Gaussian elimination, and the face dimension
  * strict positivity of all 60 inactive differences
  * for ALL 32 subsets L of N (not only the 20 the paper tabulates):
    b_L, b_L(M), the bridge block set and its cardinality, |A_L|,
    |sg(A_L)| by fixed-point closure, and sg(A_L)|N == M
  * |F| = 54 for the Lemma-5 forced set and F == A_{012} built by the
    independent Lemma-4 code path
  * the coefficientwise residual of the certificate over all 128
    coordinates of 2^{0..6}, and that it is the zero vector
  * an ANTI-VACUITY probe on the residual machinery: flipping the
    coefficient of 24|0 must move the residual to the independently
    derived value -(1-c)*Delta_{24|0}.  (The old claim that "every
    single-coefficient perturbation breaks the identity" is a theorem
    about every certificate -- true of a wrong one too -- and is now
    printed for information only, not gated.)
  * that the contradiction target 04|3 REALLY IS a term of the
    certificate LHS with a positive coefficient, that its support lies
    in N, and that deleting it from the LHS leaves exactly
    -c*Delta_{04|3} -- the link Lemma 5 needs and that nothing in the
    original file asserted
  * group membership (M / b(M) / bridge) of the 11 RHS terms, counted
    from DERIVED membership (5/3/3), and that the three groups
    partition the RHS
  * that every statement of F is a legal elementary statement on {0..6}
  * that assertions are live (a -O run inertifies ~20 internal guards)
  * 04|3 not in M and not in F; no LHS term in F; Delta_f(04|3)

===============  DECLARED, UNVERIFIABLE ASSUMPTIONS  ====================
printed by the program; these are citations to [BBS], not checked here:
  (a) one canonical copy suffices [BBS Def 4.3 Rmk 1, Lem 4.5]
  (b) sg^diamond(M|L) = sg(A_L)|N
  (c) only |L| in {2,3} need computation [BBS Cor 4.12, Lem 4.15-4.16]
      -- REMOVED by this verifier: all 32 overlaps are computed
  (d) marginal equality + the adhesion condition force delta_s = 0
      for every s in F, for every supermodular g
  (e) "f is not a coatom" / not among the coatoms of [BBS Sec 7.6]
  (f) minimality of five variables [BBS Cor 4.17, Sec 7.4]
"""

import sys
from itertools import combinations
from collections import Counter, defaultdict
from fractions import Fraction

# ----------------------------------------------------------------- ground set
N = frozenset({0, 1, 2, 3, 4})
NLIST = (0, 1, 2, 3, 4)

# ------------------------------------------ W2: the 32 values of f (verbatim)
# printed by the paper as  value -> list of subsets ; transcribed as such so
# that a duplicated or missing subset in the transcription is detectable.
F_TABLE_BY_VALUE = [
    (0, ["", "0", "1", "2", "3", "4", "01", "12", "13", "14", "24"]),
    (1, ["02", "03", "23", "34"]),
    (2, ["04", "014", "234"]),
    (3, ["012", "013", "123", "024", "124", "034", "134"]),
    (4, ["023"]),
    (6, ["0123", "0124"]),
    (7, ["0134", "0234", "1234"]),
    (12, ["01234"]),
]

# ------------------------------------------- W1: the 20 statements of M, as
# printed by the paper ("ij|K"; empty K written as the empty string)
M_WITNESS_STRINGS = [
    "01|", "01|4", "01|23", "01|24", "01|234",
    "02|13", "03|4", "03|12", "04|12",
    "12|", "12|03", "13|", "13|02", "14|", "14|0", "23|01",
    "24|", "24|0", "24|3", "34|0",
]

# ------------------- W3: the linear certificate (eq:certificate), verbatim ---
CERT_LHS = ["04|3", "01|235", "01|245", "05|234", "15|0234",
            "34|25", "35|02", "35|12", "45|02", "45|12"]
CERT_RHS = [  # (statement, signed integer coefficient, paper's group label)
    ("01|23", 1, "M"), ("04|12", 1, "M"), ("24|0", -1, "M"),
    ("24|3", 1, "M"), ("34|0", 1, "M"),
    ("01|25", 1, "b(M)"), ("05|12", 1, "b(M)"), ("15|02", 1, "b(M)"),
    ("35|012", 2, "bridge"), ("45|012", 1, "bridge"), ("45|0123", 1, "bridge"),
]

# ------------- W4/EXP: the Lemma-4 table  L -> (|A_L|, |sg(A_L)|) -----------
PAPER_LEMMA4_TABLE = {
    "01": (183, 430), "02": (184, 494), "03": (184, 437), "04": (184, 577),
    "12": (183, 419), "13": (183, 374), "14": (183, 346), "23": (184, 405),
    "24": (183, 248), "34": (184, 325),
    "012": (54, 172), "013": (54, 154), "014": (52, 52), "023": (56, 169),
    "024": (54, 81), "034": (54, 67), "123": (54, 165), "124": (53, 108),
    "134": (54, 67), "234": (54, 60),
}

# ------------------------------ EXP: other reported numbers -----------------
EXP_HISTOGRAM = (20, 20, 25, 14, 1)   # counts of Delta_f = r for r = 0..4
EXP_RANK = 19
EXP_AMBIENT_DIM = 26
EXP_FACE_DIM = 7
EXP_DELTA_04_3 = 1
EXP_M_SIZE = 20
EXP_F_SIZE = 54

# Sizes for six of the twelve overlaps the paper does NOT tabulate, as reported
# by a SECOND, INDEPENDENT closure implementation that is not distributed with
# this bundle, and not by the paper.  Because that implementation is not shipped
# here, the agreement below is not something a reader can audit; it is printed
# only as an informational cross-check between two implementations, on rows the
# paper never tabulates, and a mismatch is printed loudly but does not by itself
# impeach the paper.  No pass/fail verdict depends on these numbers.
AUDIT_UNTABULATED = {
    "": (6440, 7680),
    "0": (1064, 2047), "1": (1064, 2147), "2": (1064, 1881),
    "3": (1064, 1760), "4": (1064, 1855),
}


# ===========================================================================
#  basic representation helpers
#  a subset is a frozenset of ints; an elementary statement is a pair
#  (frozenset({i,j}), frozenset(K)) with |{i,j}| = 2 and K disjoint from {i,j}
# ===========================================================================

def parse_set(text):
    """'0134' -> frozenset({0,1,3,4}); '' -> frozenset(). Single digits only."""
    out = set()
    for ch in text.strip():
        if ch in " ,{}":
            continue
        assert ch.isdigit(), "bad digit %r in %r" % (ch, text)
        out.add(int(ch))
    return frozenset(out)


def parse_statement(text):
    """'01|234' -> (frozenset({0,1}), frozenset({2,3,4}))."""
    assert "|" in text, "statement %r lacks '|'" % (text,)
    left, right = text.split("|", 1)
    pair = parse_set(left)
    cond = parse_set(right)
    assert len(pair) == 2, "pair part of %r is not of size 2" % (text,)
    assert not (pair & cond), "K overlaps the pair in %r" % (text,)
    return (pair, cond)


def show_set(S):
    return "".join(str(x) for x in sorted(S)) if S else "-"


def show_statement(s):
    P, K = s
    a, b = sorted(P)
    return "%d%d|%s" % (a, b, "" if not K else "".join(str(x) for x in sorted(K)))


def subsets(elements):
    """Yield every subset of the iterable `elements` as a frozenset."""
    els = tuple(sorted(elements))
    for r in range(len(els) + 1):
        for c in combinations(els, r):
            yield frozenset(c)


def elementary_statements(ground):
    """All elementary statements on the ground set `ground`."""
    els = tuple(sorted(ground))
    out = []
    for i, j in combinations(els, 2):
        rest = [x for x in els if x != i and x != j]
        for K in subsets(rest):
            out.append((frozenset((i, j)), K))
    return out


def support(s):
    """P union K -- the set of variables a statement mentions."""
    return s[0] | s[1]


# ===========================================================================
#  check harness
# ===========================================================================

RESULTS = []          # list of (name, ok, detail)


def record(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print("%s %s%s" % ("PASS" if ok else "FAIL", name,
                       ("  -- " + detail) if detail else ""))
    return bool(ok)


def note(text):
    print("      " + text)


# ===========================================================================
#  CHECK f_table_well_formed  --  parse W2 into a dict on all 32 subsets
# ===========================================================================

def build_f():
    """Return (f, ok, detail).  f maps frozenset -> int for all 32 subsets."""
    f = {}
    problems = []
    for value, keys in F_TABLE_BY_VALUE:
        for key in keys:
            S = parse_set(key)
            if S in f:
                problems.append("duplicate subset %s (values %d and %d)"
                                % (show_set(S), f[S], value))
                continue
            assert not (S - N), "subset %s leaves N" % (show_set(S),)
            f[S] = int(value)
    all_subsets = set(subsets(N))
    missing = sorted(all_subsets - set(f), key=lambda S: (len(S), sorted(S)))
    extra = sorted(set(f) - all_subsets, key=lambda S: (len(S), sorted(S)))
    if missing:
        problems.append("missing subsets: " + ",".join(show_set(S) for S in missing))
    if extra:
        problems.append("extra subsets: " + ",".join(show_set(S) for S in extra))
    if len(f) != 32:
        problems.append("table has %d distinct subsets, expected 32" % len(f))
    # standardization
    if f.get(frozenset(), None) != 0:
        problems.append("f(empty) = %r, expected 0" % (f.get(frozenset()),))
    for i in NLIST:
        if f.get(frozenset((i,)), None) != 0:
            problems.append("f({%d}) = %r, expected 0" % (i, f.get(frozenset((i,)))))
    for S, v in f.items():
        if not isinstance(v, int) or isinstance(v, bool):
            problems.append("f(%s) is not a plain int: %r" % (show_set(S), v))
    ok = not problems
    detail = ("32 subsets, all distinct, standardized"
              if ok else "; ".join(problems))
    return f, ok, detail


# ===========================================================================
#  elementary difference Delta_f(ij|K) = f(Kij) + f(K) - f(Ki) - f(Kj)
# ===========================================================================

def delta(f, s):
    P, K = s
    i, j = sorted(P)
    v = f[K | P] + f[K] - f[K | frozenset((i,))] - f[K | frozenset((j,))]
    assert isinstance(v, int) and not isinstance(v, bool), \
        "non-int difference at %s: %r" % (show_statement(s), v)
    return v


def all_differences(f):
    """dict statement -> Delta_f(statement) over all 80 statements on N."""
    stmts = elementary_statements(N)
    assert len(stmts) == len(set(stmts)), "duplicate statement in enumeration"
    return dict((s, delta(f, s)) for s in stmts), stmts


def print_differences(diffs):
    """Report the full list of 80 values, grouped by conditioning-set size."""
    items = sorted(diffs.items(),
                   key=lambda kv: (sorted(kv[0][0]), len(kv[0][1]), sorted(kv[0][1])))
    line = []
    for s, v in items:
        line.append("%s=%d" % (show_statement(s), v))
        if len(line) == 8:
            note(" ".join(line))
            line = []
    if line:
        note(" ".join(line))


# ===========================================================================
#  standardized coordinate space and the active-equality matrix
# ===========================================================================

def standardized_coordinates():
    """The subsets S of N with |S| >= 2, in a fixed order.  Count is derived."""
    coords = [S for S in subsets(N) if len(S) >= 2]
    coords.sort(key=lambda S: (len(S), sorted(S)))
    return coords


def statement_row(s, index):
    """Row vector (list of ints) of the functional ij|K in standardized space.
    +1 at Kij, +1 at K, -1 at Ki, -1 at Kj; any term whose set has size <= 1
    is DROPPED (standardization makes it 0).  Terms accumulate."""
    P, K = s
    i, j = sorted(P)
    row = [0] * len(index)
    for coeff, S in ((1, K | P), (1, K),
                     (-1, K | frozenset((i,))), (-1, K | frozenset((j,)))):
        if len(S) <= 1:
            continue
        row[index[S]] += coeff
    for x in row:
        assert isinstance(x, int) and not isinstance(x, bool)
    return row


def exact_rank(rows):
    """Rank over Q by Gaussian elimination in Fraction.  No float ever."""
    mat = [[Fraction(x) for x in r] for r in rows]
    for r in mat:
        for x in r:
            assert isinstance(x, Fraction), "non-Fraction matrix entry"
    ncols = len(mat[0]) if mat else 0
    rank = 0
    row_at = 0
    for col in range(ncols):
        piv = None
        for r in range(row_at, len(mat)):
            if mat[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        mat[row_at], mat[piv] = mat[piv], mat[row_at]
        pv = mat[row_at][col]
        mat[row_at] = [x / pv for x in mat[row_at]]
        for r in range(len(mat)):
            if r != row_at and mat[r][col] != 0:
                fac = mat[r][col]
                mat[r] = [a - fac * b for a, b in zip(mat[r], mat[row_at])]
        row_at += 1
        rank += 1
        if row_at == len(mat):
            break
    for r in mat:                       # nothing leaked to float during pivoting
        for x in r:
            assert isinstance(x, Fraction), "elimination produced %r" % (x,)
    assert isinstance(rank, int) and not isinstance(rank, bool)
    return rank


# ===========================================================================
#  ELEMENTARY SEMI-GRAPHOID CLOSURE
#  rule: whenever  ij|K  and  il|jK  are both present, add  il|K  and  ij|lK
#  (i is the shared pivot).  Iterating over ORDERED triples makes the reverse
#  half of the iff redundant: swapping the roles of j and l maps one instance
#  of the rule onto the other, so this single sweep is the full closure.
#  Driven by the present statements rather than by all (i,j,l,K) tuples --
#  mathematically identical, far cheaper on a 10-element ground set.
# ===========================================================================

def sg_closure(seed, ground, sweep_counter=None):
    """Smallest superset of `seed` closed under the elementary rule.
    Repeats full sweeps over the CURRENT set until a sweep adds nothing;
    because a newly added statement may pair with any older one, only a
    no-growth full sweep certifies the fixed point."""
    ground = frozenset(ground)
    cur = set(seed)
    for P, K in cur:
        assert len(P) == 2 and not (P & K), "malformed seed statement"
        assert not ((P | K) - ground), "seed statement leaves the ground set"
    sweeps = 0
    while True:
        sweeps += 1
        grew = False
        for P, K in list(cur):
            i0, j0 = sorted(P)
            for (i, j) in ((i0, j0), (j0, i0)):
                used = K | P
                Kj = K | frozenset((j,))
                for l in ground:
                    if l in used:
                        continue
                    if (frozenset((i, l)), Kj) not in cur:
                        continue
                    for new in ((frozenset((i, l)), K),
                                (P, K | frozenset((l,)))):
                        if new not in cur:
                            cur.add(new)
                            grew = True
        if not grew:
            break
    if sweep_counter is not None:
        sweep_counter.append(sweeps)
    return cur


# ===========================================================================
#  BLOCK SET  [A,B|C] = { ij|K : i in A, j in B, C subset K subset ABC\{i,j} }
# ===========================================================================

def block_set(A, B, C):
    A, B, C = frozenset(A), frozenset(B), frozenset(C)
    assert not (A & B) and not (A & C) and not (B & C), "A,B,C not disjoint"
    ABC = A | B | C
    out = set()
    for i in sorted(A):
        for j in sorted(B):
            free = ABC - C - frozenset((i, j))
            for extra in subsets(free):
                out.add((frozenset((i, j)), C | extra))
    return out


def block_set_closed_form(A, B):
    """|A|*|B| * 2^(|A|+|B|-2), valid when C is disjoint from A and B.
    INTEGER ONLY: when |A| or |B| is 0 the block set is empty, so return 0
    rather than evaluating 2 ** (negative), which in Python yields a FLOAT
    (0.25 or 0.5) and would smuggle floating point into an exact verifier.
    The source-level '/' probe cannot see a float produced by '**'."""
    a, b = len(A), len(B)
    if a == 0 or b == 0:
        return 0
    exponent = a + b - 2
    assert exponent >= 0, "unreachable: a,b >= 1 forces a+b-2 >= 0"
    out = a * b * (2 ** exponent)
    assert isinstance(out, int) and not isinstance(out, bool), \
        "block-set closed form left the integers: %r" % (out,)
    return out


# ===========================================================================
#  CANONICAL COPY b_L  (fixes L pointwise, sends N\L to fresh labels 5,6,...)
# ===========================================================================

def canonical_copy(L):
    """Return (b, N_L, G_L) with b a dict N -> N_L."""
    L = frozenset(L)
    assert not (L - N), "L is not a subset of N"
    b = dict((x, x) for x in sorted(L))
    for t, x in enumerate(sorted(N - L), start=1):
        b[x] = 4 + t                      # fresh labels 5, 6, ...
    assert set(b.keys()) == set(N)
    N_L = frozenset(b.values())
    assert len(N_L) == 5, "b_L is not injective"
    assert (N & N_L) == L, "N n N_L = %s, expected %s" % (show_set(N & N_L), show_set(L))
    G_L = N | N_L
    assert len(G_L) == 10 - len(L)
    return b, N_L, G_L


def apply_copy(b, S):
    """Image of a set of statements under the relabelling b."""
    out = set()
    for P, K in S:
        out.add((frozenset(b[x] for x in P), frozenset(b[x] for x in K)))
    return out


def adhesion_set(M, L):
    """A_L = M u b_L(M) u [ N\\L , b_L(N\\L) | L ], plus derived diagnostics."""
    L = frozenset(L)
    b, N_L, G_L = canonical_copy(L)
    bM = apply_copy(b, M)
    A_part = N - L                      # = N \ N_L
    B_part = frozenset(b[x] for x in A_part)   # = N_L \ N
    bridge = block_set(A_part, B_part, L) if (A_part and B_part) else set()
    A_L = set(M) | bM | bridge
    inside_L = set(s for s in M if not (support(s) - L))
    return {"L": L, "b": b, "N_L": N_L, "G_L": G_L, "bM": bM,
            "bridge": bridge, "A_L": A_L, "inside_L": inside_L,
            "A": A_part, "B": B_part}


def marginal(S, ground):
    ground = frozenset(ground)
    return set(s for s in S if not (support(s) - ground))


# ===========================================================================
#  CERTIFICATE:  delta_{ij|K} as a formal integer vector on all 2^{0..6}
#  coordinates (NO standardization, nothing assumed about g)
# ===========================================================================

def add_delta(acc, s, coeff):
    """acc[S] += coeff * (formal expansion of delta_s).  Integers only."""
    P, K = s
    i, j = sorted(P)
    assert i != j, "degenerate pair"
    assert not (P & K), "K meets the pair"
    for c, S in ((1, K | P), (1, K),
                 (-1, K | frozenset((i,))), (-1, K | frozenset((j,)))):
        acc[S] += c * coeff
        assert isinstance(acc[S], int) and not isinstance(acc[S], bool)


def certificate_residual(lhs_terms, rhs_terms):
    """(sum of LHS) - (sum of coeff * RHS) as a dict frozenset -> int.
    Only the nonzero entries are returned."""
    acc = defaultdict(int)
    for text in lhs_terms:
        add_delta(acc, parse_statement(text), 1)
    for text, coeff, _group in rhs_terms:
        add_delta(acc, parse_statement(text), -int(coeff))
    return dict((S, v) for S, v in acc.items() if v != 0)


def show_residual(res):
    if not res:
        return "zero vector on all 128 coordinates"
    items = sorted(res.items(), key=lambda kv: (len(kv[0]), sorted(kv[0])))
    return "{" + ", ".join("%s:%+d" % (show_set(S), v) for S, v in items) + "}"


# ===========================================================================
#  SECTION 1  --  Lemma 2: f is supermodular and its zero set is M
# ===========================================================================

def section_paper_lemma2():
    print("--- SECTION 1: Lemma 2 (f supermodular, M_f = M) " + "-" * 22)

    f, ok, detail = build_f()
    record("f_table_well_formed", ok, detail)
    if set(f) != set(subsets(N)):
        # every later check reads f on all 32 subsets; continuing would raise a
        # KeyError instead of reporting a verdict, so stop cleanly here.
        print("FATAL: f is not defined on all 32 subsets; remaining checks "
              "cannot run.")
        print("VERDICT: %d OF %d CHECKS FAILED"
              % (len([1 for _n, o, _d in RESULTS if not o]), len(RESULTS)))
        sys.exit(1)

    stmts = elementary_statements(N)
    n = len(N)
    expect_count = (n * (n - 1) // 2) * (2 ** (n - 2))
    ok = (len(stmts) == 80 == expect_count and len(set(stmts)) == 80)
    record("eighty_elementary_statements", ok,
           "enumerated %d statements, distinct %d, C(5,2)*2^3 = %d"
           % (len(stmts), len(set(stmts)), expect_count))

    diffs, _ = all_differences(f)
    note("all 80 elementary differences Delta_f(s):")
    print_differences(diffs)
    values = list(diffs.values())
    bad = sorted([s for s in diffs if diffs[s] < 0],
                 key=lambda s: show_statement(s))
    ok = (not bad) and all(isinstance(v, int) for v in values)
    record("supermodularity_of_f", ok,
           "min Delta_f = %d over %d statements%s"
           % (min(values), len(values),
              "" if ok else "; NEGATIVE at " + ",".join(show_statement(s) for s in bad)))

    cnt = Counter(values)
    derived_hist = tuple(cnt.get(r, 0) for r in range(5))
    outside = sorted(k for k in cnt if k < 0 or k > 4)
    ok = (derived_hist == EXP_HISTOGRAM and not outside
          and sum(cnt.values()) == 80)
    record("difference_histogram", ok,
           "derived (r=0..4) = %s, paper %s, values outside 0..4: %s, total %d"
           % (str(derived_hist), str(EXP_HISTOGRAM),
              ("none" if not outside else str(outside)), sum(cnt.values())))

    M_f = set(s for s in diffs if diffs[s] == 0)
    witness = set(parse_statement(t) for t in M_WITNESS_STRINGS)
    dup_free = (len(witness) == len(M_WITNESS_STRINGS) == EXP_M_SIZE)
    only_derived = sorted(M_f - witness, key=show_statement)
    only_witness = sorted(witness - M_f, key=show_statement)
    ok = dup_free and (M_f == witness)
    note("M_f (derived zero set, %d statements): %s"
         % (len(M_f), " ".join(sorted((show_statement(s) for s in M_f)))))
    note("in M_f but not in the transcribed M: %s"
         % (" ".join(show_statement(s) for s in only_derived) or "none"))
    note("in the transcribed M but not in M_f: %s"
         % (" ".join(show_statement(s) for s in only_witness) or "none"))
    record("M_f_equals_M", ok,
           "|transcribed M| = %d (of %d printed strings), |M_f| = %d, equal = %s"
           % (len(witness), len(M_WITNESS_STRINGS), len(M_f), M_f == witness))

    M = witness if ok else M_f
    sweeps = []
    closure_M = sg_closure(M, N, sweeps)
    added = sorted(closure_M - set(M), key=show_statement)
    record("M_is_semi_graphoid_closed", closure_M == set(M) and len(closure_M) == 20,
           "|sg(M)| = %d (fixed point after %d sweep(s)); added: %s"
           % (len(closure_M), sweeps[0] if sweeps else -1,
              " ".join(show_statement(s) for s in added) or "nothing"))

    return f, diffs, set(M)


# ===========================================================================
#  SECTION 2  --  Remark 3: rank 19, ambient dim 26, face dim 7, interior
# ===========================================================================

def section_paper_remark3(diffs, M):
    print("--- SECTION 2: Remark 3 (rank / face dimension) " + "-" * 24)

    coords = standardized_coordinates()
    derived_dim = len(coords)
    formula_dim = 2 ** len(N) - len(N) - 1
    by_size = Counter(len(S) for S in coords)
    ok = (derived_dim == formula_dim == EXP_AMBIENT_DIM)
    record("ambient_dimension_26", ok,
           "derived dim = %d (2^5-5-1 = %d), paper %d; sizes %s"
           % (derived_dim, formula_dim, EXP_AMBIENT_DIM,
              ", ".join("|S|=%d:%d" % (k, by_size[k]) for k in sorted(by_size))))

    index = dict((S, k) for k, S in enumerate(coords))
    ordered_M = sorted(M, key=show_statement)
    rows = [statement_row(s, index) for s in ordered_M]
    note("active-equality matrix: %d rows x %d columns" % (len(rows), derived_dim))
    for s, r in zip(ordered_M, rows):
        nz = [(show_set(coords[k]), r[k]) for k in range(derived_dim) if r[k]]
        note("  %-8s %s" % (show_statement(s),
                            " ".join("%s%+d" % (a, c) for a, c in nz)))
    rank = exact_rank(rows)
    face_dim = derived_dim - rank
    ok = (rank == EXP_RANK and face_dim == EXP_FACE_DIM and rank < len(rows))
    record("rank_of_active_equalities_is_19", ok,
           "derived rank = %d (paper %d), face dim = %d - %d = %d (paper %d), "
           "rank < %d rows so a dependency exists: %s"
           % (rank, EXP_RANK, derived_dim, rank, face_dim, EXP_FACE_DIM,
              len(rows), rank < len(rows)))

    inactive = dict((s, v) for s, v in diffs.items() if s not in M)
    nonpos = sorted([s for s in inactive if inactive[s] <= 0], key=show_statement)
    ok = (len(inactive) == 60 and not nonpos and (80 - len(M) == 60))
    record("f_is_in_the_relative_interior_of_that_face", ok,
           "%d inactive statements (80 - %d = %d), min slack = %d, non-positive: %s"
           % (len(inactive), len(M), 80 - len(M),
              min(inactive.values()) if inactive else 0,
              " ".join(show_statement(s) for s in nonpos) or "none"))
    return rank, face_dim


# ===========================================================================
#  SECTION 3  --  Lemma 4: A_L, sg(A_L) and the N-marginal, for ALL 32 L
# ===========================================================================

def analyse_overlap(M, L):
    """Derive every Lemma-4 quantity for one overlap L.  Returns a dict."""
    d = adhesion_set(M, L)
    bridge_n = len(d["bridge"])
    closed = block_set_closed_form(d["A"], d["B"]) if (d["A"] and d["B"]) else 0
    inside = len(d["inside_L"])
    predicted = 2 * len(M) - inside + bridge_n
    sweeps = []
    closure = sg_closure(d["A_L"], d["G_L"], sweeps)
    marg = marginal(closure, N)
    d.update({
        "bridge_n": bridge_n,
        "bridge_closed_form": closed,
        "bridge_ok": bridge_n == closed,
        "inside": inside,
        "size_A": len(d["A_L"]),
        "size_predicted": predicted,
        "size_formula_ok": len(d["A_L"]) == predicted,
        "size_sg": len(closure),
        "sweeps": sweeps[0] if sweeps else -1,
        "closure": closure,
        "marginal": marg,
        "marginal_is_M": marg == set(M),
        "marginal_extra": sorted(marg - set(M), key=show_statement),
        "marginal_missing": sorted(set(M) - marg, key=show_statement),
        "ground_size": len(d["G_L"]),
    })
    return d


def overlap_line(key, d, tabulated):
    tag = ""
    if tabulated is not None:
        tag = " paper (%d,%d)" % tabulated
    else:
        tag = " [NEW derived data -- not tabulated in the paper]"
    return ("L=%-5s |G_L|=%2d  M-inside-L=%d  bridge=%-4d  |A_L|=%-5d  "
            "|sg(A_L)|=%-5d  sweeps=%d  marginal==M:%s%s"
            % (key, d["ground_size"], d["inside"], d["bridge_n"], d["size_A"],
               d["size_sg"], d["sweeps"], d["marginal_is_M"], tag))


def section_paper_lemma4(M):
    print("--- SECTION 3: Lemma 4 (all 32 overlaps L) " + "-" * 29)
    all_L = sorted(subsets(N), key=lambda S: (len(S), sorted(S)))
    data = {}
    for L in all_L:
        key = show_set(L) if L else ""
        data[key if key != "-" else ""] = analyse_overlap(M, L)

    def dkey(L):
        return show_set(L) if L else ""

    # ---- tabulated rows: |A_L| ------------------------------------------
    size_bad, sg_bad, bridge_bad, formula_bad = [], [], [], []
    for key in sorted(PAPER_LEMMA4_TABLE, key=lambda k: (len(k), k)):
        d = data[key]
        exp_A, exp_sg = PAPER_LEMMA4_TABLE[key]
        note(overlap_line(key, d, (exp_A, exp_sg)))
        if d["size_A"] != exp_A:
            size_bad.append("%s: got %d want %d" % (key, d["size_A"], exp_A))
        if d["size_sg"] != exp_sg:
            sg_bad.append("%s: got %d want %d" % (key, d["size_sg"], exp_sg))
        if not d["bridge_ok"]:
            bridge_bad.append("%s: %d vs closed form %d"
                              % (key, d["bridge_n"], d["bridge_closed_form"]))
        if not d["size_formula_ok"]:
            formula_bad.append("%s: |A_L|=%d vs 40-%d+%d=%d"
                               % (key, d["size_A"], d["inside"],
                                  d["bridge_n"], d["size_predicted"]))
    ok = not (size_bad or bridge_bad or formula_bad)
    record("A_L_construction_sizes", ok,
           ("all 20 tabulated |A_L| reproduce; bridge cardinality matches "
            "|A|*|B|*2^(|A|+|B|-2) (144 for |L|=2, 16 for |L|=3); "
            "|A_L| = 40 - #(M inside L) + bridge holds in all 20 rows")
           if ok else "; ".join(size_bad + bridge_bad + formula_bad))

    fixed_pt = data["014"]
    ok = (not sg_bad)
    record("closure_sizes_sg_A_L", ok,
           ("all 20 tabulated |sg(A_L)| reproduce; L=014 is a fixed point "
            "(|A_L|=%d = |sg(A_L)|=%d), so the closure adds nothing spurious"
            % (fixed_pt["size_A"], fixed_pt["size_sg"]))
           if ok else "; ".join(sg_bad))

    # ---- tabulated rows: the actual content, marginal == M ---------------
    marg_bad = []
    for key in sorted(PAPER_LEMMA4_TABLE, key=lambda k: (len(k), k)):
        d = data[key]
        if not d["marginal_is_M"]:
            marg_bad.append("%s: extra %s missing %s"
                            % (key,
                               " ".join(show_statement(s) for s in d["marginal_extra"]) or "none",
                               " ".join(show_statement(s) for s in d["marginal_missing"]) or "none"))
    # NOTE ON WHERE THE CONTENT IS: the inclusion marg superset M is STRUCTURAL
    # (M subset A_L subset sg(A_L) and every statement of M is supported on N),
    # so `marginal_missing` can never be non-empty and that half of the equality
    # is a tautology.  All the content lies in `marginal_extra` being empty,
    # i.e. the closure adds no new statement supported on N.
    always_true = all(not (set(M) - data[k]["marginal"]) for k in PAPER_LEMMA4_TABLE)
    note("structural half (marg superset M) held in all 20 rows as it must: %s; "
         "the informative half is that no row has an EXTRA N-statement"
         % always_true)
    record("marginal_of_each_closure_is_M", not marg_bad,
           "sg(A_L)|N == M in all 20 tabulated rows (informative half: no "
           "closure produced a statement on N outside M)"
           if not marg_bad else "; ".join(marg_bad))

    # ---- all 32, including the 12 the paper does not tabulate ------------
    note("the 12 overlaps the paper does NOT tabulate (derived here, new data):")
    untab_bad = []
    for L in all_L:
        key = dkey(L)
        if key in PAPER_LEMMA4_TABLE:
            continue
        d = data[key]
        note(overlap_line(key if key else "(empty)", d, None))
        if not d["marginal_is_M"]:
            untab_bad.append("%s: extra %s" % (key,
                             " ".join(show_statement(s) for s in d["marginal_extra"])))
        # the arithmetic cross-checks were previously inspected ONLY on the 20
        # rows where the paper already supplies |A_L|; run them on the other 12
        # too, so that every one of the 32 constructions is self-consistent.
        if not d["bridge_ok"]:
            untab_bad.append("%s: bridge %d vs closed form %d"
                             % (key or "(empty)", d["bridge_n"],
                                d["bridge_closed_form"]))
        if not d["size_formula_ok"]:
            untab_bad.append("%s: |A_L|=%d vs 40-%d+%d=%d"
                             % (key or "(empty)", d["size_A"], d["inside"],
                                d["bridge_n"], d["size_predicted"]))
    # informational: agree with a second, independent implementation?
    agree, disagree = 0, []
    for key, (aA, aSG) in sorted(AUDIT_UNTABULATED.items(), key=lambda kv: kv[0]):
        d = data[key]
        if (d["size_A"], d["size_sg"]) == (aA, aSG):
            agree += 1
        else:
            disagree.append("L=%s this run (%d,%d) vs independent audit (%d,%d)"
                            % (key or "(empty)", d["size_A"], d["size_sg"], aA, aSG))
    note("CROSS-CHECK against a second implementation that is NOT shipped with "
         "this bundle, on %d of the 12 untabulated rows (these are NOT paper "
         "values, the other implementation cannot be inspected by a reader, and "
         "no verdict depends on this line): %d agree, %d differ%s"
         % (len(AUDIT_UNTABULATED), agree, len(disagree),
            "" if not disagree else " -- MISMATCH: " + "; ".join(disagree)))
    every = [dkey(L) for L in all_L]
    bad32 = [k for k in every if not data[k]["marginal_is_M"]]
    n_untab = len(every) - len(PAPER_LEMMA4_TABLE)
    record("all_32_overlaps_not_just_20",
           len(every) == 32 and n_untab == 12 and not bad32 and not untab_bad,
           "sg(A_L)|N == M for all %d subsets L of N, including the %d rows the "
           "paper does not tabulate (this REMOVES the paper's citation-only "
           "reduction to |L| in {2,3})" % (len(every), n_untab)
           if not (bad32 or untab_bad)
           else "failed for L in %s; untabulated failures: %s"
                % (bad32, "; ".join(untab_bad) or "none"))
    return data


# ===========================================================================
#  SECTION 4  --  Lemma 5: the forced set F and the linear certificate
# ===========================================================================

def section_paper_lemma5(f, M, lemma4_data):
    print("--- SECTION 4: Lemma 5 (forced set F, certificate) " + "-" * 21)

    # W4: the paper's copy b for L = {0,1,2}
    b_paper = {0: 0, 1: 1, 2: 2, 3: 5, 4: 6}
    b_canon, _, G_L = canonical_copy(frozenset({0, 1, 2}))
    note("paper's b = %s ; canonical b_{012} = %s ; identical: %s"
         % (b_paper, b_canon, b_paper == b_canon))
    bM = apply_copy(b_paper, M)
    bridge = block_set({3, 4}, {5, 6}, {0, 1, 2})
    F = set(M) | bM | bridge
    A_012 = lemma4_data["012"]["A_L"]
    parts_ok = (len(M) == 20 and len(bM) == 20 and len(bridge) == 16
                and len(F) == EXP_F_SIZE)
    same = (F == A_012)
    ground7 = frozenset(range(7))
    on_ground = all(not (support(s) - ground7) for s in F)
    ground_ok = (G_L == ground7)
    # previously only PRINTED; the equality of the paper's b with the canonical
    # b_{012} is gated here rather than left to the indirect F == A_012 test.
    copy_ok = (b_paper == b_canon)
    if not copy_ok:
        note("MISMATCH: the paper's copy b is not the canonical b_{012}")
    record("forced_set_F_construction",
           parts_ok and same and on_ground and ground_ok and copy_ok,
           "|M|=%d |b(M)|=%d |bridge|=%d |F|=%d (paper %d); F == A_{012} from the "
           "Lemma-4 code path: %s; G_{012} == {0..6}: %s; all of F lives on "
           "{0..6}: %s; paper's b == canonical b_{012}: %s"
           % (len(M), len(bM), len(bridge), len(F), EXP_F_SIZE, same,
              ground_ok, on_ground, copy_ok))

    # ---- the identity, coefficientwise over all 128 coordinates ----------
    res = certificate_residual(CERT_LHS, CERT_RHS)
    ncoords = 2 ** 7
    note("certificate LHS: %s" % " + ".join(CERT_LHS))
    note("certificate RHS: %s"
         % " ".join("%+d*%s[%s]" % (c, t, g) for t, c, g in CERT_RHS))
    note("residual (LHS - RHS) over the %d coordinates of 2^{0..6}: %s"
         % (ncoords, show_residual(res)))
    record("certificate_identity_coefficientwise", not res,
           "residual is the zero vector on all %d coordinates, so the identity "
           "holds for EVERY set function g (no supermodularity assumed)" % ncoords
           if not res else "nonzero residual " + show_residual(res))

    # ---- ANTI-VACUITY PROBE on certificate_residual ----------------------
    # The old framing of this block ("the printed minus sign is load bearing")
    # was decoration: perturbing one coefficient by eps changes the residual by
    # -eps * Delta_t, and a single Delta functional is NEVER the zero vector
    # (its four sets K+P, K, K+i, K+j are pairwise distinct because i != j and
    # neither i nor j lies in K).  So "every perturbation breaks it" is a
    # theorem about every possible certificate, true even of a wrong one, and
    # proves nothing about this one.
    # What IS worth checking is that certificate_residual is not vacuously
    # returning the empty dict: flipping the coefficient of 24|0 from -1 to +1
    # must move the residual by exactly -(1 - (-1)) * Delta_{24|0}, a value
    # DERIVED here rather than taken from anybody's run.
    flipped = [(t, (1 if t == "24|0" else c), g) for t, c, g in CERT_RHS]
    res_flip = certificate_residual(CERT_LHS, flipped)
    old_coeff = [c for t, c, _g in CERT_RHS if t == "24|0"]
    pred = defaultdict(int)
    for c in old_coeff:                      # normally exactly one, c = -1
        add_delta(pred, parse_statement("24|0"), -(1 - c))
    pred_nz = dict((S, v) for S, v in pred.items() if v != 0)
    flip_matches = (res_flip == pred_nz) and bool(pred_nz)
    note("printed coefficient of 24|0: %s; flipping it to +1 moves the residual "
         "to %s, derived prediction -(1-c)*Delta_{24|0} = %s; agree: %s"
         % (old_coeff, show_residual(res_flip), show_residual(pred_nz),
            flip_matches))
    slack = []
    for idx in range(len(CERT_LHS)):
        for eps in (1, -1):
            # perturbing the LHS coefficient of term t by eps is the same as
            # appending (t, -eps) to the RHS list (RHS enters with a minus).
            pert_rhs = list(CERT_RHS) + [(CERT_LHS[idx], -eps, "perturbation")]
            if not certificate_residual(CERT_LHS, pert_rhs):
                slack.append("LHS term %s tolerates %+d" % (CERT_LHS[idx], eps))
    for k in range(len(CERT_RHS)):
        for eps in (1, -1):
            pert = list(CERT_RHS)
            t, c, g = pert[k]
            pert[k] = (t, c + eps, g)
            if not certificate_residual(CERT_LHS, pert):
                slack.append("RHS term %s tolerates %+d" % (t, eps))
    total_pert = 2 * (len(CERT_LHS) + len(CERT_RHS))
    note("all %d single-coefficient +/-1 perturbations leave a nonzero residual "
         "(%d exceptions) -- STRUCTURAL, true of every certificate, reported for "
         "information only and NOT part of this check's verdict"
         % (total_pert, len(slack)))
    ok = flip_matches
    record("residual_probe_is_not_vacuous", ok,
           "certificate_residual demonstrably responds to its input: perturbing "
           "the coefficient of 24|0 moves the residual to exactly the "
           "independently derived value %s, so the zero residual reported above "
           "is a computed fact and not an empty accumulator"
           % show_residual(res_flip)
           if ok else "residual after flipping 24|0 is %s but the derived "
                      "prediction is %s"
                      % (show_residual(res_flip), show_residual(pred_nz)))

    section_cert_hygiene(F, M, bM, bridge)
    section_target_is_on_the_lhs(M, F)
    section_contradiction(f, M, F)
    return F, bM, bridge


# ===========================================================================
#  CHECK certificate_term_hygiene
# ===========================================================================

def section_cert_hygiene(F, M, bM, bridge):
    ground7 = frozenset(range(7))
    problems = []
    for text in CERT_LHS + [t for t, _c, _g in CERT_RHS]:
        P, K = parse_statement(text)          # parse_statement asserts i!=j, K n P = 0
        if len(P) != 2 or (P & K):
            problems.append("%s is not a valid elementary statement" % text)
        if support((P, K)) - ground7:
            problems.append("%s leaves the ground set {0..6}" % text)
    # every statement of the forced set must itself be a legal elementary
    # statement on the adhesion ground set, or "delta_s = 0 for s in F" is not
    # even a well-formed hypothesis
    for s in F:
        P, K = s
        if len(P) != 2 or (P & K) or (support(s) - ground7):
            problems.append("F contains the malformed statement %s"
                            % show_statement(s))
    groups = {"M": set(M), "b(M)": set(bM), "bridge": set(bridge)}
    # counts are DERIVED from membership in the three constructed sets, not
    # tallied from the paper's own group labels (that would be a tautology on
    # the transcription).
    counts = Counter()
    for text, coeff, claimed in CERT_RHS:
        s = parse_statement(text)
        member = sorted(g for g, S in groups.items() if s in S)
        note("RHS %-8s coeff %+d  claimed %-6s  actually in %s"
             % (text, coeff, claimed, ",".join(member) or "NONE"))
        if claimed not in member:
            problems.append("%s claimed %s but lies in %s"
                            % (text, claimed, member or "no group"))
        if s not in F:
            problems.append("%s is not in F at all" % text)
        if len(member) != 1:
            problems.append("%s lies in %d of the three groups (%s), so the "
                            "paper's M / b(M) / bridge grouping is not a "
                            "partition of the RHS"
                            % (text, len(member), ",".join(member) or "none"))
        for g in member:
            counts[g] += 1
    # every LHS coefficient is +1: DERIVED by checking the LHS term list has no
    # repetition (a repeated term would silently carry coefficient 2 or more)
    lhs_parsed = [parse_statement(t) for t in CERT_LHS]
    lhs_multi = sorted(t for t, c in Counter(CERT_LHS).items() if c > 1)
    lhs_coeffs_ok = (len(set(lhs_parsed)) == len(lhs_parsed) == 10
                     and not lhs_multi)
    if not lhs_coeffs_ok:
        problems.append("LHS is not 10 distinct terms of coefficient +1 "
                        "(repeats: %s)" % (lhs_multi or "none as text but "
                                           "duplicated after parsing"))
    group_counts = (counts["M"], counts["b(M)"], counts["bridge"])
    if group_counts != (5, 3, 3):
        problems.append("RHS group sizes %s, paper says (5,3,3)" % (group_counts,))
    record("certificate_term_hygiene", not problems and lhs_coeffs_ok,
           "all %d terms are valid elementary statements on {0..6}; every LHS "
           "coefficient is +1; RHS group sizes M/b(M)/bridge = %s as the paper "
           "attributes them" % (len(CERT_LHS) + len(CERT_RHS), group_counts)
           if not problems else "; ".join(problems))


# ===========================================================================
#  CHECK contradiction_is_genuine
# ===========================================================================

TARGET_TEXT = "04|3"        # the statement Lemma 5 forces to be zero


def section_target_is_on_the_lhs(M, F):
    """CHECK certificate_forces_the_target.

    The rest of the file proves (i) the 21-term identity holds coefficientwise,
    (ii) every RHS term lies in F so its delta vanishes, (iii) 04|3 is not in M.
    NONE of that connects (i)-(ii) to (iii).  The missing links, both of which
    Lemma 5 needs and neither of which was asserted anywhere, are:

      L1  the target is a term of the LHS, carrying a POSITIVE coefficient
          (otherwise 'all LHS terms vanish' says nothing about the target);
      L2  the target is a statement on N, i.e. support(04|3) subset N
          (otherwise it is not in the N-marginal and there is no contradiction).

    L1 is derived, not asserted: the coefficient of the target in the LHS is
    recomputed from CERT_LHS by counting occurrences, and the identity is
    re-evaluated with the target deleted from the LHS -- the residual must then
    be exactly -c * Delta_target, which pins the coefficient c independently of
    the occurrence count."""
    target = parse_statement(TARGET_TEXT)
    problems = []

    # ---- L1a: occurrence count of the target among the LHS terms -----------
    lhs_parsed = [parse_statement(t) for t in CERT_LHS]
    coeff = sum(1 for s in lhs_parsed if s == target)
    if coeff <= 0:
        problems.append("%s does not occur on the certificate LHS at all, so "
                        "the identity forces nothing about it" % TARGET_TEXT)

    # ---- L1b: independent pin on that coefficient --------------------------
    # deleting every copy of the target from the LHS must leave exactly
    # -coeff * Delta_target as the residual.
    reduced_lhs = [t for t in CERT_LHS if parse_statement(t) != target]
    res_reduced = certificate_residual(reduced_lhs, CERT_RHS)
    expect = defaultdict(int)
    add_delta(expect, target, -coeff)
    expect_nz = dict((S, v) for S, v in expect.items() if v != 0)
    if res_reduced != expect_nz:
        problems.append("dropping %s from the LHS leaves %s, expected %s"
                        % (TARGET_TEXT, show_residual(res_reduced),
                           show_residual(expect_nz)))
    if coeff > 0 and not res_reduced:
        problems.append("dropping %s from the LHS left the identity intact, so "
                        "its coefficient is really 0" % TARGET_TEXT)

    # ---- L2: the target is an elementary statement on N --------------------
    on_N = not (support(target) - N)
    if not on_N:
        problems.append("support(%s) = %s leaves N = %s, so it cannot appear "
                        "in the N-marginal"
                        % (TARGET_TEXT, show_set(support(target)), show_set(N)))

    # ---- and the target must NOT be one of the forced-zero statements ------
    if target in F:
        problems.append("%s lies in F, so delta = 0 is a hypothesis, not a "
                        "consequence" % TARGET_TEXT)
    if target in set(M):
        problems.append("%s lies in M, so there is no contradiction" % TARGET_TEXT)

    # ---- every LHS coefficient must be strictly positive -------------------
    lhs_counts = Counter(lhs_parsed)
    nonpos = [show_statement(s) for s, c in lhs_counts.items() if c <= 0]
    if nonpos:
        problems.append("LHS terms with non-positive coefficient: %s"
                        % ",".join(nonpos))

    note("target %s: LHS coefficient derived as %+d; support %s subset N: %s; "
         "in F: %s; in M: %s"
         % (TARGET_TEXT, coeff, show_set(support(target)), on_N,
            target in F, target in set(M)))
    note("residual after deleting %s from the LHS: %s (must equal %s)"
         % (TARGET_TEXT, show_residual(res_reduced), show_residual(expect_nz)))
    record("certificate_forces_the_target", not problems,
           "%s occurs on the certificate LHS with derived coefficient %+d, is a "
           "statement on N, is outside F and outside M, and every LHS "
           "coefficient is positive -- so 'sum of LHS = 0' really does force "
           "delta_%s = 0 and really does contradict the required marginal"
           % (TARGET_TEXT, coeff, TARGET_TEXT)
           if not problems else "; ".join(problems))


def section_contradiction(f, M, F):
    target = parse_statement(TARGET_TEXT)
    in_M = target in set(M)
    in_F = target in F
    # independent hand computation from the derived f table
    d = f[frozenset({0, 3, 4})] + f[frozenset({3})] \
        - f[frozenset({0, 3})] - f[frozenset({3, 4})]
    d2 = delta(f, target)
    lhs_in_F = sorted([t for t in CERT_LHS if parse_statement(t) in F])
    ok = ((not in_M) and (not in_F) and d == EXP_DELTA_04_3 == d2 and d > 0
          and not lhs_in_F)
    note("Delta_f(04|3) = f(034)+f(3)-f(03)-f(34) = %d+%d-%d-%d = %d (paper %d)"
         % (f[frozenset({0, 3, 4})], f[frozenset({3})], f[frozenset({0, 3})],
            f[frozenset({3, 4})], d, EXP_DELTA_04_3))
    note("LHS terms of the certificate that lie in F: %s"
         % (" ".join(lhs_in_F) or "none (all 10 are outside F, so none is "
            "trivially zero by hypothesis)"))
    record("contradiction_is_genuine", ok,
           "04|3 in M: %s, in F: %s; Delta_f(04|3) = %d > 0 (two independent "
           "routes agree: %d and %d); none of the 10 LHS terms lies in F"
           % (in_M, in_F, d, d, d2))


# ===========================================================================
#  CHECK exact_arithmetic_only
# ===========================================================================

ALLOWED_MODULES = frozenset(["sys", "itertools", "collections", "fractions"])
THIRD_PARTY_MARKERS = ("nu" "mpy", "sc" "ipy", "sy" "mpy", "netw" "orkx",
                       "pan" "das", "mp" "math")


def code_lines(lines):
    """Yield (lineno, code) for source lines OUTSIDE triple-quoted blocks,
    with comments and single-line quoted spans removed.  Used only by the
    exact-arithmetic self-probe, so that slashes appearing inside the module
    docstring or inside message strings are not mistaken for division."""
    in_doc = False
    for k, raw in enumerate(lines, start=1):
        marks = raw.count('"""') + raw.count("'''")
        if in_doc:
            if marks % 2 == 1:
                in_doc = False
            continue
        if marks % 2 == 1:
            in_doc = True
            continue
        code = raw.split("#")[0]
        out, quote = [], None
        for ch in code:
            if quote is None:
                if ch in ("'", '"'):
                    quote = ch
                else:
                    out.append(ch)
            else:
                if ch == quote:
                    quote = None
        yield k, "".join(out)


def assertions_are_enabled():
    """True iff the interpreter is NOT running with -O / PYTHONOPTIMIZE.

    Roughly twenty of this file's internal guards (non-int values, malformed
    seed statements, injectivity of b_L, 'no Fraction leaked to float') are
    plain `assert` statements.  Under -O every one of them disappears and the
    verifier would happily print ALL CHECKS PASS with those guards absent, so
    whether assertions are live is itself a check."""
    try:
        assert False
        return False
    except AssertionError:
        return True


def section_exact_arithmetic(f, diffs, rank):
    print("--- SECTION 5: exact arithmetic and stdlib-only self-check " + "-" * 13)
    problems = []
    live = assertions_are_enabled()
    note("assertions live (not running under -O): %s" % live)
    if not live:
        problems.append("assertions are DISABLED (-O or PYTHONOPTIMIZE): every "
                        "assert-based guard in this file is inert, so no PASS "
                        "verdict from this run is trustworthy")
    for name in THIRD_PARTY_MARKERS:
        if name in sys.modules:
            problems.append("third-party module %s is loaded" % name)
    for S, v in f.items():
        if isinstance(v, float):
            problems.append("f(%s) is a float" % show_set(S))
    for s, v in diffs.items():
        if not isinstance(v, int) or isinstance(v, bool):
            problems.append("Delta_f(%s) is %r, not an int"
                            % (show_statement(s), v))
    if not isinstance(rank, int) or isinstance(rank, bool):
        problems.append("rank is %r, not an int" % (rank,))
    # a live sanity probe that Fraction elimination really is exact: a matrix
    # whose rank floating point would get wrong is not needed -- instead check
    # that 1/3 round-trips exactly, which no float can do.
    third = Fraction(1, 3)
    if third * 3 != 1:
        problems.append("Fraction arithmetic is not exact")
    if isinstance(third + 0, float):
        problems.append("Fraction arithmetic leaked to float")
    # source-level probe: every import is from the allowed stdlib list, and
    # every true division ('/' not '//') occurs on a line that divides
    # Fractions inside exact_rank.
    try:
        with open(__file__, "r") as fh:
            lines = fh.read().split("\n")
    except (IOError, OSError, NameError):
        lines = None
        note("source self-probe COULD NOT RUN (script path unavailable)")
        # previously this silently degraded the check to four isinstance()
        # tests and still reported PASS; a probe that did not run must not be
        # reported as a satisfied probe.
        problems.append("the source-level division/import probe could not read "
                        "its own source, so 'no floating point anywhere' is "
                        "NOT verified by this run")
    if lines is not None:
        imported, div_lines = [], []
        for k, code in code_lines(lines):
            t = code.strip()
            if t.startswith("import ") or t.startswith("from "):
                root = t.split()[1].split(".")[0]
                imported.append(root)
                if root not in ALLOWED_MODULES:
                    problems.append("disallowed import %r on line %d" % (root, k))
            if "/" in code.replace("//", ""):
                div_lines.append((k, t))
        note("imports found in the source: %s" % ", ".join(sorted(set(imported))))
        note("true-division sites: %s"
             % (", ".join("line %d (%s)" % (k, t) for k, t in div_lines) or "none"))
        for k, t in div_lines:
            # the only true divisions in the program are the two pivot
            # normalisations in exact_rank, where both operands are Fractions
            if "pv" not in t:
                problems.append("line %d divides outside the Fraction pivot "
                                "normalisation: %s" % (k, t))
    record("exact_arithmetic_only", not problems,
           "all stored numbers are int or Fraction; Fraction arithmetic exact; "
           "imports limited to %s; no third-party module loaded"
           % ", ".join(sorted(ALLOWED_MODULES))
           if not problems else "; ".join(problems))
    # (There is deliberately no "sentinel" here.  An earlier comment claimed the
    # probe above only scans text preceding this marker; code_lines() scans the
    # WHOLE file and always did, so the claim was false.  Scanning everything is
    # the stronger behaviour, hence the claim is removed rather than implemented.)


def print_assumptions():
    print("--- DECLARED ASSUMPTIONS (cited to [BBS], NOT verified here) " + "-" * 11)
    for line in [
        "(a) one canonical copy suffices [BBS Def 4.3 Rmk 1, Lem 4.5]",
        "(b) sg^diamond(M|L) = sg(A_L)|N  (definitional in [BBS])",
        "(c) only |L| in {2,3} need computation [BBS Cor 4.12, Lem 4.15-4.16]"
        "  -- NOT relied upon: this verifier computes all 32 overlaps",
        "(d) marginal equality plus the adhesion condition force delta_s = 0"
        " for every s in F and every supermodular g inducing the adhesion",
        "(e) 'f is not a coatom' / not among the five-variable coatoms of"
        " [BBS Sec 7.6] -- the rank-19/face-dim-7 arithmetic IS checked here,"
        " but the step from face dimension 7 to 'not a coatom' is not",
        "(f) minimality of five variables [BBS Cor 4.17, Sec 7.4] -- nothing"
        " about |N| = 4 is computed in this paper or in this verifier",
    ]:
        print("      ASSUMED " + line)


# ===========================================================================
#  main
# ===========================================================================

def main():
    print("=" * 78)
    print("VERIFIER -- five-variable self-adhesivity counterexample")
    print("TAKEN from the paper: the 32 values of f; the 20 statements of M; the")
    print("  certificate's 21 terms and integer coefficients; the overlap L={0,1,2}")
    print("  and the copy b; the definitions.  Expected numbers ((20,20,25,14,1),")
    print("  the Lemma-4 table, rank 19, dim 26, face 7, Delta_f(04|3)=1) are used")
    print("  ONLY as comparison targets.")
    print("DERIVED here: every difference, the histogram, M_f, sg(M), the 20x26")
    print("  matrix and its exact rank, all 32 adhesion sets and their closures and")
    print("  N-marginals, |F|, and the certificate's 128-coordinate residual.")
    print("=" * 78)

    f, diffs, M = section_paper_lemma2()
    print("")
    rank, _face = section_paper_remark3(diffs, M)
    print("")
    lemma4 = section_paper_lemma4(M)
    print("")
    section_paper_lemma5(f, M, lemma4)
    print("")
    section_exact_arithmetic(f, diffs, rank)
    print("")
    print_assumptions()
    print("")

    failed = [name for name, ok, _d in RESULTS if not ok]
    total = len(RESULTS)
    if failed:
        print("FAILED CHECKS: " + ", ".join(failed))
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


if __name__ == "__main__":
    sys.exit(main())
