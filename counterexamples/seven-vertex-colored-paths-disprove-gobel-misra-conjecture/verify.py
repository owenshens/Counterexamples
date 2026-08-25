#!/usr/bin/env python3
"""Verification of a seven-vertex counterexample to the odd-order clause of a
conjecture on equalities between determinant polynomials of colored paths.

Standard library only; exact integer / rational arithmetic throughout.  No
floating point enters any decision.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs; transcribed, never used as evidence)
--------------------------------------------------------------------------
  * the exhibited pair: seven vertices all of one color a, edge words
        p = b b c b c c        q = b c c b b c
  * the polynomial the paper says both determinants equal:
        x^7 - 3x^5(B+C) + 2x^3(B^2+3BC+C^2) - 2xBC(B+C),   B=y_b^2, C=y_c^2
  * the matching index lists printed in the proof:
        size two:   13 14 15 16 24 25 26 35 36 46
        size three: 135 136 146 246
  * the auxiliary values A_1 = 3B+3C, A_2 = 2B^2+6BC+2C^2, A_3 = 2B^2C+2BC^2
  * the reverse word the proof exhibits: reverse(p) = c c b c b b
  * the parity-subsequence table:  p odd (b,c,c) even (b,b,c)
                                   q odd (b,c,b) even (c,b,c)
  * the minimality census table, columns
        m : paths, equal pairs, reflections, constructions, residual
        3 :     30,     14,     14,   0, 0
        4 :    470,    248,    232,  16, 0
        5 : 11,562,  5,780,  5,772,   8, 0
        6 : 394,092, 197,496, 197,016, 480, 0
  * the described edge condition of the odd-order construction: each parity
    class of edge positions is monochromatic and the two colors are exchanged
    between the two paths
  * the positive-definiteness region of the remark: x > 2*max(|y_b|,|y_c|)

--------------------------------------------------------------------------
DERIVED HERE (computed from first principles; nothing below is asserted)
--------------------------------------------------------------------------
  * the matchings of the seven-vertex path, enumerated, counted, printed back,
    and compared with the paper's index lists;
  * the two determinants, computed THREE independent ways -- Leibniz sum over
    all 5040 permutations of the 7x7 symmetric tridiagonal matrix, the
    three-term tridiagonal minor recurrence, and the matching expansion --
    then compared with each other and with the paper's closed form;
  * the hypotheses of the clause being refuted: odd order, common vertex
    coloring, equal vertex- and edge-color multiplicities, equal determinants;
  * the failure of the conclusion: the pair is neither the identity nor a
    reflection, and the described edge condition of the odd-order
    construction is evaluated over all eight orientation/interchange choices
    and is false in every one;
  * exact rational leading-principal-minor tests of the positive-definiteness
    remark, with a negative control outside the stated region;
  * the whole minimality census for 3 <= m <= 6, from the paper's own
    definition of the canonical word sets, regenerating all five columns;
  * the parity-monochromatic-exchange structure of every pair the reflection
    test leaves over below seven vertices: on the edge words in the odd rows
    m = 3, 5, which is the condition the paper quotes, and on the vertex words
    in the even rows m = 4, 6, which is THIS PROGRAM'S OWN ANALOGUE and has no
    textual authority in the paper (see the NOT-RE-RUN lines); each row's test
    is run on the other parity as a discriminating control, and the ambient
    selectivity of the even-order analogue is measured so that the reader can
    see what its 16-of-16 and 480-of-480 hit rates are worth;
  * two seven-vertex sub-censuses, which rediscover the exhibited pair;
  * every count printed anywhere in the output, including the counts quoted in
    the NOT-RE-RUN scope note (11,481 canonical vertex colorings, 1,602
    canonical edge words, their product, 42 two-color edge words) and the
    matching-size profile, is either computed here or cross-derived from a
    closed form, never only pinned to a literal: the path counts of every
    census are re-derived as |vertex words| x |edge words|, the matching
    profile as C(n-k+1,k) with Fibonacci total, and the two alternating
    canonical edge words as the two label orders of an alternation.

NOT RE-RUN (stated plainly).  Every item in this list is printed in full at the
end of the run, one "NOT RE-RUN:" line per item, in this order; this summary
enumerates them all so that nothing here is dropped by a reader who counts:
  (a) the full seven-vertex census over all vertex colorings (18,392,562 paths)
      is outside the time budget of a single-process pure-Python run;
  (b) the two cited theorems of the source conjecture are not reproduced from
      their own statements -- only the edge condition the paper quotes, which
      is the odd-order one -- so the 'constructions' and 'residual' columns of
      the census table are reproduced only as their sum, and only against a
      NECESSARY condition;
  (c) the even rows of the census table are therefore tested only against an
      analogue this program invented, which has no textual authority;
  (d) at seven vertices with a uniform (palindromic) vertex word, the quoted
      odd-order edge condition can never be met by a pair the reflection test
      leaves over, so the "0 of 12 explained" figure of the uniform sub-census
      is a structural artefact and carries no evidential weight;
  (e) the cited source itself is never consulted: its authors, its identifier,
      and the numbering it is quoted under are taken from the paper on trust,
      and this program asserts nothing about any of them.

Output contract: one "PASS <name>" / "FAIL <name>" line per check, then a
final VERDICT line; exit status 0 iff every check passed.
"""

import sys
from fractions import Fraction
from itertools import permutations, product

CHECKS = []

# The counts quoted in the first NOT RE-RUN line below.  They are NOT taken on
# trust: each is the size of a set this program actually builds, recorded into
# SCOPE_FIGURES as it is built, and checked against this table by
# check_scope_note_figures.  The note text itself is formatted from this table,
# so the prose and the checked figures cannot drift apart.
NOTE_FIGURES = {
    "canonical_vertex_colorings_7": 11481,
    "canonical_edge_words_6": 1602,
    "two_color_edge_words_6": 42,
    "un_run_product": 11481 * 1602,
}
SCOPE_FIGURES = {}

NOTES = [
    "the full seven-vertex census over all canonical vertex colorings "
    "({:,} x {:,} = {:,} paths) is not run here; the seven-vertex "
    "work is restricted to a uniform vertex coloring over all {:,} canonical "
    "edge words, and to all {:,} canonical vertex colorings against the {:d} "
    "canonical edge words that use at most two colors.".format(
        NOTE_FIGURES["canonical_vertex_colorings_7"],
        NOTE_FIGURES["canonical_edge_words_6"],
        NOTE_FIGURES["un_run_product"],
        NOTE_FIGURES["canonical_edge_words_6"],
        NOTE_FIGURES["canonical_vertex_colorings_7"],
        NOTE_FIGURES["two_color_edge_words_6"]),
    "the two cited theorems of the source conjecture are not reproduced from "
    "their own statements; only the edge condition of the odd-order "
    "construction, as quoted in the paper, is implemented.  That "
    "condition is necessary for the construction, so refuting it refutes the "
    "alternative, but the 'constructions' column of the census table is "
    "checked only as the count of pairs left over after the reflection test, "
    "and the left-over pairs are then tested against a NECESSARY condition "
    "only -- the parity-monochromatic-exchange structure, on the edge words at "
    "odd order and on the vertex words at even order.  The table's zero "
    "residual is therefore corroborated in every row, never proved.",
    "the EVEN rows of the minimality table are not tested against the theorem "
    "the paper cites for them.  For even m the paper cites Theorem 3.6 of the "
    "source and does not reproduce its statement, so this program has no "
    "even-order condition with any textual authority to implement.  The "
    "condition applied to the even rows here -- parity-monochromatic exchange "
    "on the VERTEX words -- is THE PROGRAM'S OWN INVENTED ANALOGUE of the "
    "odd-order edge condition, inferred from the fact that at even order it is "
    "the vertex positions that split into two parity classes of equal size.  "
    "Nothing in the paper asserts that analogue, and this program does not "
    "verify that Theorem 3.6 implies it.  If it is not in fact implied by "
    "Theorem 3.6, then the 16-of-16 and 480-of-480 figures corroborate NOTHING "
    "about the even rows, and for m = 4 and m = 6 the 'constructions' and "
    "'residual' columns are reproduced here only as their sum.  What those "
    "figures do establish, with no theorem involved, is a measured structural "
    "fact about the census output alone: every even-order pair the reflection "
    "test leaves over has both vertex words alternating in two colors, which "
    "is a vanishingly small fraction of the available vertex-word pairs -- the "
    "fraction is derived and printed by the selectivity check.  In "
    "consequence the EVEN-order half of the paper's minimality Proposition -- "
    "that no counterexample has four or six vertices -- is corroborated by "
    "this program and is NOT independently checkable from the material shipped "
    "beside it: the paper's proof tests each remaining even-order pair against "
    "Theorem 3.6 of the source, and nothing here reproduces that theorem's "
    "statement.  The ODD rows m = 3, 5 are not subject to this caveat, because "
    "there the paper does quote the condition and this program implements what "
    "it quotes.",
    "at odd order the quoted edge condition forces the second edge word to be "
    "the reverse of the first, so under a palindromic vertex word every pair "
    "meeting it is already a reflection.  The seven-vertex uniform sub-census "
    "leaves 12 pairs over and none of them can meet the condition for that "
    "structural reason; that zero is computed and reported as structural, and "
    "carries no independent weight.",
    "the cited source is never consulted by this program, at any point and by "
    "any means.  Everything the paper attributes to it -- the matching "
    "expansion, the two construction theorems, and the conjecture being "
    "refuted, together with the numbers they are cited under -- reaches this "
    "program only as the paper's own quotation of them, and the author names "
    "and the preprint identifier are likewise taken from the paper on trust.  "
    "This program asserts nothing about any of that: not that the attribution "
    "is right, not that it is wrong.  Two things therefore have to be checked "
    "at the source and cannot be checked here.  First, that the quoted "
    "numbering is the source's numbering; note in particular that the cited "
    "work's title concerns colored Gaussian CYCLES while the conjecture "
    "refuted here is stated for colored PATHS.  Second, that the reading the "
    "paper adopts of the odd-order theorem -- taking the edge-index range in "
    "its condition (3) to be {1,...,m-1} -- is the intended one, since it is "
    "that reading, exactly as the paper states it, which this program "
    "implements and finds violated by the exhibited pair.  The algebra of the "
    "counterexample does not depend on any of this and is verified here in "
    "full; the identification of WHAT it refutes does depend on it and is not "
    "verified here at all.",
]

# ----------------------------------------------------------------- inputs ---
P_WORD = "bbcbcc"          # taken from the paper
Q_WORD = "bccbbc"          # taken from the paper
P_REVERSE_CLAIMED = "ccbcbb"
M_VERTICES = 7
PAPER_SIZE2 = ["13", "14", "15", "16", "24", "25", "26", "35", "36", "46"]
PAPER_SIZE3 = ["135", "136", "146", "246"]
PAPER_PARITY_TABLE = {
    "p": {"odd": ("b", "c", "c"), "even": ("b", "b", "c")},
    "q": {"odd": ("b", "c", "b"), "even": ("c", "b", "c")},
}
PAPER_CENSUS = {          # m -> (paths, equal pairs, reflections, constr, residual)
    3: (30, 14, 14, 0, 0),
    4: (470, 248, 232, 16, 0),
    5: (11562, 5780, 5772, 8, 0),
    6: (394092, 197496, 197016, 480, 0),
}


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)


# ------------------------------------------------- exact polynomial algebra --
# A polynomial is a dict {exponent tuple -> int}, with zero coefficients pruned.
# The exhibited object uses NV = 3 variables in the order (x, y_b, y_c).
NV = 3
X = (1, 0, 0)
YB = (0, 1, 0)
YC = (0, 0, 1)


def pzero():
    return {}


def pmono(exp, coeff=1):
    return {} if coeff == 0 else {tuple(exp): coeff}


def padd(a, b):
    out = dict(a)
    for k, v in b.items():
        nv = out.get(k, 0) + v
        if nv:
            out[k] = nv
        else:
            out.pop(k, None)
    return out


def pmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(ka[i] + kb[i] for i in range(len(ka)))
            nv = out.get(k, 0) + va * vb
            if nv:
                out[k] = nv
            else:
                out.pop(k, None)
    return out


def pscal(a, c):
    return {} if c == 0 else dict((k, v * c) for k, v in a.items())


def pstr(a, names=("x", "y_b", "y_c")):
    if not a:
        return "0"
    parts = []
    for k in sorted(a, reverse=True):
        fs = []
        if a[k] != 1 or all(e == 0 for e in k):
            fs.append(str(a[k]))
        for i, e in enumerate(k):
            if e == 1:
                fs.append(names[i])
            elif e > 1:
                fs.append("%s^%d" % (names[i], e))
        parts.append("*".join(fs))
    return " + ".join(parts)


def peval(a, vals):
    """Exact evaluation at rational values; vals is a tuple of Fractions."""
    tot = Fraction(0)
    for k, c in a.items():
        term = Fraction(c)
        for i, e in enumerate(k):
            term *= vals[i] ** e
        tot += term
    return tot


# --------------------------------------------------- paths and matchings -----
def matchings(n_edges):
    """All matchings of the path with edges 1..n_edges, as sorted index tuples.
    A matching is a set of edge indices no two of which are consecutive."""
    res = []

    def rec(i, cur):
        if i > n_edges:
            res.append(tuple(cur))
            return
        rec(i + 1, cur)                      # edge i unused
        if not cur or cur[-1] < i - 1:       # edge i usable
            cur.append(i)
            rec(i + 2, cur)
            cur.pop()

    rec(1, [])
    return sorted(res, key=lambda s: (len(s), s))


def kmatrix(vword, eword, varof):
    """The symmetric tridiagonal matrix of a colored path, entries as polys.
    varof maps a color letter to a variable exponent tuple."""
    m = len(vword)
    K = [[pzero() for _ in range(m)] for _ in range(m)]
    for i in range(m):
        K[i][i] = pmono(varof[vword[i]])
    for i in range(m - 1):
        e = pmono(varof[eword[i]])
        K[i][i + 1] = e
        K[i + 1][i] = e
    return K


def det_leibniz(K):
    """Determinant as the full signed sum over all m! permutations."""
    m = len(K)
    total = pzero()
    for perm in permutations(range(m)):
        term = None
        for i in range(m):
            e = K[i][perm[i]]
            if not e:
                term = None
                break
            term = e if term is None else pmul(term, e)
        if term is None:
            continue
        inv = 0
        for i in range(m):
            for j in range(i + 1, m):
                if perm[i] > perm[j]:
                    inv += 1
        total = padd(total, pscal(term, -1 if inv % 2 else 1))
    return total


def det_recurrence(K):
    """Determinant by the three-term recurrence for tridiagonal matrices."""
    m = len(K)
    prev2, prev1 = pmono((0,) * NV, 1), K[0][0]
    for i in range(1, m):
        cur = padd(pmul(K[i][i], prev1),
                   pscal(pmul(pmul(K[i][i - 1], K[i - 1][i]), prev2), -1))
        prev2, prev1 = prev1, cur
    return prev1


def det_matching_expansion(vword, eword, varof, Ms):
    """The paper's matching expansion, evaluated directly."""
    total = pzero()
    m = len(vword)
    for M in Ms:
        covered = set()
        term = pmono((0,) * NV, -1 if len(M) % 2 else 1)
        for i in M:
            term = pmul(term, pmul(pmono(varof[eword[i - 1]]),
                                   pmono(varof[eword[i - 1]])))
            covered.add(i)
            covered.add(i + 1)
        for j in range(1, m + 1):
            if j not in covered:
                term = pmul(term, pmono(varof[vword[j - 1]]))
        total = padd(total, term)
    return total


VAROF = {"a": X, "b": YB, "c": YC}


def vertex_word_of(eword):
    """The uniform vertex word the paper attaches to an edge word of a path:
    one vertex more than the number of edges, every vertex of color a."""
    return "a" * (len(eword) + 1)


def factorial(n):
    out = 1
    for i in range(2, n + 1):
        out *= i
    return out


def binom(n, k):
    """Exact binomial coefficient, zero outside the range."""
    if k < 0 or k > n or n < 0:
        return 0
    out = 1
    for i in range(k):
        out = out * (n - i) // (i + 1)
    return out


def check_wellformed(p, q):
    """Check 1: the exhibited object is what the paper says it is."""
    m = M_VERTICES
    facts = []
    ok = True
    for nm, w in (("p", p), ("q", q)):
        ok &= (len(w) == m - 1)
        ok &= set(w) == set("bc")
        facts.append("%s=%s len=%d colors=%s" %
                     (nm, w, len(w), "".join(sorted(set(w)))))
    ok &= (p != q)
    vword = "a" * m
    ok &= (len(vword) == m and len(set(vword)) == 1)
    ok &= set(vword).isdisjoint(set(p) | set(q))   # disjoint alphabets
    ck("exhibited_pair_wellformed", ok,
       "vertex word %s; %s; %s; distinct=%s; alphabets disjoint" %
       (vword, facts[0], facts[1], p != q))
    return ok


def check_multiplicities(p, q):
    """Check 2: the two necessary conditions the paper derives for equality."""
    mp = dict((d, p.count(d)) for d in sorted(set(p)))
    mq = dict((d, q.count(d)) for d in sorted(set(q)))
    same_edge = (mp == mq)
    vp, vq = vertex_word_of(p), vertex_word_of(q)
    nvp = dict((d, vp.count(d)) for d in sorted(set(vp)))
    nvq = dict((d, vq.count(d)) for d in sorted(set(vq)))
    same_vertex = (nvp == nvq)
    ck("color_multiplicities_match", same_edge and same_vertex,
       "edge multiplicities p=%s q=%s; vertex multiplicities %s and %s" %
       (mp, mq, nvp, nvq))
    return same_edge and same_vertex


def check_matching_lists():
    """Check 3: enumerate the matchings and compare with the printed lists.

    The size profile and the total are not merely pinned to the numbers a
    previous run happened to print: both are also DERIVED from closed forms
    computed here.  A matching of k edges in a path of n edges is a choice of k
    pairwise non-consecutive indices out of n, of which there are C(n-k+1, k),
    and the total over all k is the Fibonacci number F(n+2).  The enumeration
    must agree with both."""
    n_edges = M_VERTICES - 1
    Ms = matchings(n_edges)
    by_size = {}
    for M in Ms:
        by_size.setdefault(len(M), []).append("".join(str(i) for i in M))
    counts = dict((k, len(v)) for k, v in by_size.items())
    ok = (counts == {0: 1, 1: 6, 2: 10, 3: 4})
    ok &= (len(Ms) == 21)
    # the same two figures, derived instead of asserted
    derived = dict((k, binom(n_edges - k + 1, k)) for k in range(n_edges + 1)
                   if binom(n_edges - k + 1, k) > 0)
    fib = [1, 1]
    while len(fib) < n_edges + 2:
        fib.append(fib[-1] + fib[-2])       # fib[i] = F(i+1), so F(n+2)=fib[n+1]
    ok &= (counts == derived)
    ok &= (len(Ms) == sum(derived.values()) == fib[n_edges + 1])
    ok &= (sorted(by_size.get(2, [])) == sorted(PAPER_SIZE2))
    ok &= (sorted(by_size.get(3, [])) == sorted(PAPER_SIZE3))
    ck("matching_enumeration", ok,
       "sizes %s total %d; size2 %s; size3 %s; the profile and the total are "
       "also DERIVED, not just pinned: C(n-k+1,k) gives %s and the Fibonacci "
       "count F(%d)=%d agrees" %
       (counts, len(Ms), " ".join(sorted(by_size.get(2, []))),
        " ".join(sorted(by_size.get(3, []))), derived, n_edges + 2,
        fib[n_edges + 1]))
    return Ms


def paper_closed_form():
    """The polynomial the paper prints, assembled from the paper's coefficients:
       x^7 - 3x^5(B+C) + 2x^3(B^2+3BC+C^2) - 2xBC(B+C), B=y_b^2, C=y_c^2."""
    B = pmul(pmono(YB), pmono(YB))
    C = pmul(pmono(YC), pmono(YC))
    x = pmono(X)

    def xp(k):
        out = pmono((0,) * NV, 1)
        for _ in range(k):
            out = pmul(out, x)
        return out

    t1 = xp(7)
    t2 = pscal(pmul(xp(5), padd(B, C)), -3)
    t3 = pscal(pmul(xp(3), padd(padd(pmul(B, B), pscal(pmul(B, C), 3)),
                                pmul(C, C))), 2)
    t4 = pscal(pmul(x, pmul(pmul(B, C), padd(B, C))), -2)
    return padd(padd(t1, t2), padd(t3, t4))


def check_determinants(p, q, Ms):
    """Checks 4-8: three independent determinant algorithms, the closed form,
    and the equality that is the heart of the counterexample."""
    vword = "a" * M_VERTICES
    dets = {}
    for nm, w in (("p", p), ("q", q)):
        K = kmatrix(vword, w, VAROF)
        dets[nm] = {
            "leibniz": det_leibniz(K),
            "recur": det_recurrence(K),
            "match": det_matching_expansion(vword, w, VAROF, Ms),
        }
    agree = all(d["leibniz"] == d["recur"] == d["match"] for d in dets.values())
    ck("three_algorithms_agree", agree,
       "Leibniz over %d permutations, tridiagonal recurrence and matching "
       "expansion agree for both paths" % factorial(M_VERTICES))

    target = paper_closed_form()
    okp = dets["p"]["leibniz"] == target
    okq = dets["q"]["leibniz"] == target
    ck("determinant_of_p_matches_paper", okp, pstr(dets["p"]["leibniz"]))
    ck("determinant_of_q_matches_paper", okq, pstr(dets["q"]["leibniz"]))
    eq = dets["p"]["leibniz"] == dets["q"]["leibniz"]
    diff = padd(dets["p"]["leibniz"], pscal(dets["q"]["leibniz"], -1))
    ck("determinants_are_equal", eq and not diff,
       "computed difference D(a^7,p) - D(a^7,q) = %s" % pstr(diff))
    return dets["p"]["leibniz"]


def check_A_values(p, q, Ms):
    """Check 9: the auxiliary sums A_1, A_2, A_3 printed in the proof."""
    B = pmul(pmono(YB), pmono(YB))
    C = pmul(pmono(YC), pmono(YC))
    claimed = {
        1: padd(pscal(B, 3), pscal(C, 3)),
        2: padd(padd(pscal(pmul(B, B), 2), pscal(pmul(B, C), 6)),
                pscal(pmul(C, C), 2)),
        3: padd(pscal(pmul(pmul(B, B), C), 2), pscal(pmul(B, pmul(C, C)), 2)),
    }
    ok = True
    shown = []
    for nm, w in (("p", p), ("q", q)):
        for r in (1, 2, 3):
            acc = pzero()
            for M in Ms:
                if len(M) != r:
                    continue
                term = pmono((0,) * NV, 1)
                for i in M:
                    v = pmono(VAROF[w[i - 1]])
                    term = pmul(term, pmul(v, v))
                acc = padd(acc, term)
            ok &= (acc == claimed[r])
            if nm == "p":
                shown.append("A_%d=%s" % (r, pstr(acc)))
    ck("auxiliary_sums_A_r", ok, "; ".join(shown))
    return ok


def check_not_identity_or_reflection(p, q):
    """Check 10: the first two conjectured alternatives are computed to fail.
    Both vertex words are a^7, so a reflection of the pair is exactly a
    reversal of the edge word."""
    rp, rq = p[::-1], q[::-1]
    ok = (p != q) and (rp != q) and (rq != p)
    ok &= (rp == P_REVERSE_CLAIMED)      # the reverse word the paper prints
    ck("not_identity_and_not_reflection", ok,
       "p=%s q=%s reverse(p)=%s reverse(q)=%s : p!=q, rev(p)!=q, rev(q)!=p" %
       (p, q, rp, rq))
    return ok


def parity_classes(w):
    return tuple(w[0::2]), tuple(w[1::2])


def check_parity_table(p, q):
    """Check 11: the odd/even edge subsequences printed in the proof, and the
    computed fact that none of the four is monochromatic."""
    got = {}
    for nm, w in (("p", p), ("q", q)):
        o, e = parity_classes(w)
        got[nm] = {"odd": o, "even": e}
    ok = (got == PAPER_PARITY_TABLE)
    mono = [len(set(got[nm][k])) == 1 for nm in ("p", "q")
            for k in ("odd", "even")]
    ok &= (not any(mono))
    ck("parity_subsequence_table", ok,
       "p odd %s even %s | q odd %s even %s ; monochromatic flags %s" %
       (got["p"]["odd"], got["p"]["even"], got["q"]["odd"], got["q"]["even"],
        mono))
    return ok


def odd_construction_edge_condition(w1, w2):
    """The edge condition the paper quotes for the odd-order construction:
    each parity class of edge positions is monochromatic in each word, and the
    two colors are exchanged between the two words."""
    o1, e1 = parity_classes(w1)
    o2, e2 = parity_classes(w2)
    for cls in (o1, e1, o2, e2):
        if len(set(cls)) != 1:
            return False
    d1, d2 = o1[0], e1[0]
    if d1 == d2:
        return False
    return o2[0] == d2 and e2[0] == d1


def check_construction_fails(p, q):
    """Check 12 -- the load-bearing check.  Evaluate the quoted edge condition
    of the odd-order construction over every orientation of each path and both
    orderings of the pair.  All eight evaluations must come out false."""
    results = []
    for r1 in (False, True):
        for r2 in (False, True):
            w1 = p[::-1] if r1 else p
            w2 = q[::-1] if r2 else q
            for swap in (False, True):
                a, b = (w2, w1) if swap else (w1, w2)
                results.append(((r1, r2, swap),
                                odd_construction_edge_condition(a, b)))
    ok = not any(v for _, v in results)
    ck("odd_construction_condition_fails_in_all_orientations", ok,
       "%d of %d orientation/interchange choices satisfy the condition" %
       (sum(1 for _, v in results if v), len(results)))
    return ok


def leading_minors(vword, eword, xv, yv):
    """Exact leading principal minors of K, by the tridiagonal recurrence.
    xv maps a vertex color to a Fraction, yv an edge color to a Fraction."""
    d = [Fraction(1), xv[vword[0]]]
    for k in range(1, len(vword)):
        y = yv[eword[k - 1]]
        d.append(xv[vword[k]] * d[k] - y * y * d[k - 1])
    return d[1:]


def check_positive_definite_region(p, q, target):
    """Check 13: the remark's positive-definiteness claim, tested exactly at
    rational points of the stated region, with a control outside it."""
    F = Fraction
    inside = [(F(3), F(1), F(1)), (F(5), F(1), F(2)), (F(9, 4), F(1), F(1)),
              (F(41, 20), F(1), F(-1)), (F(7), F(-3), F(2)),
              (F(100, 3), F(7), F(-16))]
    in_region = all(xval > 2 * max(abs(b), abs(c)) for xval, b, c in inside)
    pd_everywhere = True
    matches_target = True
    for xval, b, c in inside:
        for w in (p, q):
            mins = leading_minors("a" * M_VERTICES, w, {"a": xval},
                                  {"b": b, "c": c})
            pd_everywhere &= all(mv > 0 for mv in mins)
            matches_target &= (mins[-1] == peval(target, (xval, b, c)))
        matches_target &= (peval(target, (xval, b, c)) ==
                           peval(target, (xval, c, b)))   # symmetry in b,c
    # control: a point outside the region where the test must report failure
    ctrl = leading_minors("a" * M_VERTICES, p, {"a": F(1)},
                          {"b": F(1), "c": F(1)})
    control_fires = not all(mv > 0 for mv in ctrl)
    ok = in_region and pd_everywhere and matches_target and control_fires
    ck("positive_definite_inside_stated_region", ok,
       "%d rational points, all inside the region=%s, all 7 leading minors "
       "positive for both paths=%s, last minor equals the determinant=%s; "
       "control x=1,y=1 detected as not positive definite=%s (minors %s)" %
       (len(inside), in_region, pd_everywhere, matches_target, control_fires,
        [str(v) for v in ctrl]))
    return ok


def canonical_words(n):
    """The paper's set of canonical words of length n: surjective words on
    {1..r} for some 1<=r<=n whose color multiplicities are weakly decreasing."""
    out = []
    for r in range(1, n + 1):
        for wd in product(range(1, r + 1), repeat=n):
            cnt = [0] * (r + 1)
            for ch in wd:
                cnt[ch] += 1
            if any(cnt[i] == 0 for i in range(1, r + 1)):
                continue
            if all(cnt[i] >= cnt[i + 1] for i in range(1, r)):
                out.append(wd)
    return out


def det_signature(u, w, Ms, m):
    """The determinant of a colored path as an exact coefficient dictionary on
    the full exponent vector (x_1..x_m, y_1..y_{m-1}), in canonical form."""
    d = {}
    for M in Ms:
        e = [0] * (2 * m - 1)
        cov = set()
        for i in M:
            e[m + w[i - 1] - 1] += 2
            cov.add(i)
            cov.add(i + 1)
        for j in range(1, m + 1):
            if j not in cov:
                e[u[j - 1] - 1] += 1
        k = tuple(e)
        d[k] = d.get(k, 0) + (-1 if len(M) % 2 else 1)
    return tuple(sorted((k, v) for k, v in d.items() if v))


def det_leibniz_numeric(u, w, m):
    """Leibniz determinant of a colored path given as numeric color words,
    returned in the same canonical form as det_signature."""
    nv = 2 * m - 1

    def unit(i):
        e = [0] * nv
        e[i] = 1
        return tuple(e)

    K = [[pzero() for _ in range(m)] for _ in range(m)]
    for i in range(m):
        K[i][i] = pmono(unit(u[i] - 1))
    for i in range(m - 1):
        e = pmono(unit(m + w[i] - 1))
        K[i][i + 1] = e
        K[i + 1][i] = e
    poly = det_leibniz(K)
    return tuple(sorted((k, v) for k, v in poly.items() if v))


def check_census_kernel():
    """Check: the signature the census groups on is the determinant.  Compared
    against the permutation expansion on every path of orders three and four
    and on a stride sample of order five."""
    ok = True
    tested = 0
    for m in (3, 4):
        Ms = matchings(m - 1)
        for u in canonical_words(m):
            for w in canonical_words(m - 1):
                tested += 1
                ok &= (det_signature(u, w, Ms, m) ==
                       det_leibniz_numeric(u, w, m))
    m = 5
    Ms = matchings(m - 1)
    pairs = [(u, w) for u in canonical_words(m) for w in canonical_words(m - 1)]
    for u, w in pairs[::57]:
        tested += 1
        ok &= (det_signature(u, w, Ms, m) == det_leibniz_numeric(u, w, m))
    ck("census_kernel_equals_permutation_expansion", ok,
       "%d paths compared coefficient by coefficient, signs included" % tested)
    return ok


def run_census(m, uwords=None, wwords=None):
    """Group all colored paths by exact determinant signature, then split the
    equal pairs into reflections and the rest."""
    U = canonical_words(m) if uwords is None else uwords
    W = canonical_words(m - 1) if wwords is None else wwords
    Ms = matchings(m - 1)
    groups = {}
    npaths = 0
    for u in U:
        for w in W:
            npaths += 1
            groups.setdefault(det_signature(u, w, Ms, m), []).append((u, w))
    eq = refl = 0
    rest = []
    for lst in groups.values():
        n = len(lst)
        for i in range(n):
            for j in range(i + 1, n):
                eq += 1
                a, b = lst[i], lst[j]
                if (a[0][::-1], a[1][::-1]) == b:
                    refl += 1
                else:
                    rest.append((a, b))
    return npaths, eq, refl, rest


def check_census_row(m):
    """Checks 14-17: regenerate one row of the minimality table."""
    npaths, eq, refl, rest = run_census(m)
    paths_c, eq_c, refl_c, con_c, res_c = PAPER_CENSUS[m]
    ok = (npaths == paths_c and eq == eq_c and refl == refl_c
          and len(rest) == con_c + res_c)
    ck("census_row_m%d" % m, ok,
       "paths %d, equal pairs %d, reflections %d, non-reflection pairs %d "
       "(paper: %d/%d/%d/%d+%d -- the constructions and residual columns are "
       "checked here only as their sum, not separated)" %
       (npaths, eq, refl, len(rest), paths_c, eq_c, refl_c, con_c, res_c))
    return rest


def parity_structure_holds(w1, w2):
    """Does the parity-monochromatic-exchange condition hold between the two
    words for some orientation of each and some ordering of the two?"""
    for r1 in (False, True):
        for r2 in (False, True):
            a = w1[::-1] if r1 else w1
            b = w2[::-1] if r2 else w2
            if odd_construction_edge_condition(a, b):
                return True
            if odd_construction_edge_condition(b, a):
                return True
    return False


def parity_admissible(w):
    """Can this word appear in ANY pair satisfying the condition?  The
    condition demands that both parity classes of both words be monochromatic
    in two distinct colors, and reversal cannot create that property (it either
    fixes the two classes or swaps them), so a word failing this cannot be
    rescued by orientation or by the ordering of the pair."""
    o, e = parity_classes(w)
    return len(set(o)) == 1 and len(set(e)) == 1 and o[0] != e[0]


def explained_by_parity_structure(pair, idx):
    """Does the parity-monochromatic-exchange condition hold on component idx
    of the pair -- idx 0 the vertex words, idx 1 the edge words -- for some
    orientation of each path and some ordering of the pair?"""
    return parity_structure_holds(pair[0][idx], pair[1][idx])


def explained_by_odd_construction(pair):
    """The quoted edge condition of the odd-order construction, tested on the
    edge words over all orientations and both orderings of the pair."""
    return explained_by_parity_structure(pair, 1)


def check_odd_rows_are_explained(rest_by_m):
    """Check 18: for the odd orders inside the census range, every pair the
    reflection test leaves over does satisfy the quoted edge condition -- the
    contrast that makes the seven-vertex pair a counterexample."""
    tot = good = 0
    for m in (3, 5):
        for pair in rest_by_m[m]:
            tot += 1
            if explained_by_odd_construction(pair):
                good += 1
    ck("odd_order_residual_pairs_below_seven_are_explained", good == tot,
       "%d of %d non-reflection equal pairs at m=3,5 satisfy the quoted edge "
       "condition" % (good, tot))
    return good == tot


def check_even_rows_have_parity_structure(rest_by_m):
    """The even rows of the minimality table, tested against a condition THIS
    PROGRAM INVENTED.  The paper cites Theorem 3.6 of the source for even m and
    never states its content, so there is no even-order condition in the paper
    to implement.  What is implemented instead: at even order the number of
    vertices is even, so it is the VERTEX positions that split into two parity
    classes of equal size, and the parity-monochromatic-exchange structure the
    paper does quote for the ODD-order construction is transplanted onto the
    vertex words.  That transplant is an inference of this program's, with no
    textual authority whatsoever, and it is not verified to follow from
    Theorem 3.6.  If it does not follow, the hit rates below say nothing about
    the even rows of the table.  What survives that caveat is a bare structural
    fact about the census output: every pair the reflection test leaves over at
    m = 4 and 6 has both vertex words alternating in two colors.  The same test
    is run on the odd rows m = 3, 5 as a discriminating control, where it must
    find nothing -- there the structure sits on the edge words instead."""
    ok = True
    parts = []
    for m in (4, 6):
        rest = rest_by_m[m]
        good = sum(1 for pr in rest if explained_by_parity_structure(pr, 0))
        ok &= (len(rest) > 0 and good == len(rest))
        parts.append("m=%d %d of %d" % (m, good, len(rest)))
    odd_tot = sum(len(rest_by_m[m]) for m in (3, 5))
    odd_hits = sum(1 for m in (3, 5) for pr in rest_by_m[m]
                   if explained_by_parity_structure(pr, 0))
    ok &= (odd_hits == 0 and odd_tot > 0)
    ck("even_order_residual_pairs_match_a_parity_analogue_this_program_invented",
       ok,
       "vertex-word parity structure %s; the same test on the odd rows finds "
       "%d of %d, as it must, since there the structure is on the edge words.  "
       "SCOPE: this vertex-word condition is NOT in the paper -- the paper "
       "cites Theorem 3.6 for even m without stating it, and this analogue is "
       "this program's own construction, not shown here to be implied by that "
       "theorem; if it is not implied, these hit rates corroborate nothing "
       "about the even rows and only the sum constructions+residual is "
       "reproduced there"
       % ("; ".join(parts), odd_hits, odd_tot))
    return ok


def check_even_analogue_selectivity(rest_by_m):
    """How much is the even-order analogue's perfect hit rate worth, taken
    purely as a statement about the census output?  A condition satisfied by
    almost every pair of vertex words would make 16 of 16 and 480 of 480
    meaningless.  Here the ambient rate is derived: over all unordered pairs
    (repetition allowed, since two paths may share a vertex word) of canonical
    vertex words of the given length, count those the analogue accepts.  At
    m = 4 the count is done twice -- once over every pair, once over the words
    parity_admissible keeps -- and the two must agree, which is what licenses
    the pruning used at m = 6.  Finally, the distinct vertex-word pairs actually
    carried by the left-over census pairs are exhibited and counted.  This check
    grants the invented analogue no authority it does not have: it measures the
    analogue's selectivity, it does not test Theorem 3.6."""
    ok = True
    parts = []
    for m in (4, 6):
        U = canonical_words(m)
        n = len(U)
        denom = n * (n + 1) // 2
        adm = [w for w in U if parity_admissible(w)]
        hits = sum(1 for i in range(len(adm)) for j in range(i, len(adm))
                   if parity_structure_holds(adm[i], adm[j]))
        if m == 4:                      # exhaustive audit of the pruning rule
            brute = sum(1 for i in range(n) for j in range(i, n)
                        if parity_structure_holds(U[i], U[j]))
            ok &= (brute == hits)
            parts.append("pruning audited exhaustively at m=4: %d accepted "
                         "pairs by brute force over all %d, %d by pruning" %
                         (brute, denom, hits))
        used = set()
        for a, b in rest_by_m[m]:
            used.add(tuple(sorted((a[0], b[0]))))
        inside = all(parity_structure_holds(x, y) for x, y in used)
        ok &= (0 < hits < denom) and inside and len(used) > 0
        parts.append("m=%d: %d of %d vertex-word pairs accepted (%d of %d words "
                     "can participate at all); the %d left-over census pairs "
                     "use %d distinct vertex-word pair(s) %s, all accepted=%s" %
                     (m, hits, denom, len(adm), n, len(rest_by_m[m]), len(used),
                      sorted("/".join("".join(str(t) for t in w) for w in pr)
                             for pr in used), inside))
    ck("even_order_analogue_is_selective_not_universal", ok, "  ".join(parts))
    return ok


def check_small_orders():
    """Check 19: the two closed forms the proof uses for m<=2, and that order
    two admits no unexplained equal pair."""
    v = {"1": (1, 0, 0, 0), "2": (0, 1, 0, 0), "e": (0, 0, 1, 0)}
    d1 = det_leibniz(kmatrix("1", "", v))
    ok = (d1 == {(1, 0, 0, 0): 1})
    d2 = det_leibniz(kmatrix("12", "e", v))
    want2 = {(1, 1, 0, 0): 1, (0, 0, 2, 0): -1}
    ok &= (d2 == want2)
    npaths, eq, refl, rest = run_census(2)
    ok &= (npaths == 3 and eq == 1 and refl == 1 and not rest)
    n2, n1 = len(canonical_words(2)), len(canonical_words(1))
    ok &= (npaths == n2 * n1)      # the path count, derived not pinned
    ck("orders_one_and_two", ok,
       "D(u1)=x_u1 and D(u1u2,w1)=x_u1*x_u2-y_w1^2 confirmed; order two: "
       "%d paths (= %d canonical vertex words x %d canonical edge words), "
       "%d equal pair, %d reflection, %d unexplained" %
       (npaths, n2, n1, eq, refl, len(rest)))
    return ok


def to_numeric(word):
    """Relabel a word by order of first appearance: b,c -> 1,2."""
    seen = {}
    out = []
    for ch in word:
        if ch not in seen:
            seen[ch] = len(seen) + 1
        out.append(seen[ch])
    return tuple(out)


def pair_present(rest, target):
    """Is the unordered pair `target` of (vertex word, edge word) in `rest`?"""
    t = set(target)
    for a, b in rest:
        if set((a, b)) == t:
            return True
    return False


def check_seven_vertex_subcensus(p, q):
    """Rerun the census machinery at seven vertices with a uniform vertex
    coloring.  It must rediscover the exhibited pair among the pairs the
    reflection test leaves over.

    A warning about the "none of the left-over pairs is explained" figure: at
    odd order the edge word has EVEN length, so a word whose two parity classes
    are monochromatic in distinct colors is alternating, and its color-exchanged
    partner is literally its reverse.  The uniform vertex word is a palindrome,
    so any pair meeting the quoted edge condition here is a strict reflection
    and has already been removed.  The zero is therefore structural, not
    evidence, and this function computes the structural fact rather than
    presenting the zero as a test: it exhibits the alternating canonical words
    and verifies that the exchange of each is its own reverse."""
    u = tuple([1] * M_VERTICES)
    W = canonical_words(M_VERTICES - 1)
    SCOPE_FIGURES["canonical_edge_words_6"] = len(W)
    pn, qn = to_numeric(p), to_numeric(q)
    ok = (pn in W) and (qn in W)
    npaths, eq, refl, rest = run_census(M_VERTICES, [u], W)
    ok &= (npaths == len(W))      # one vertex word times every canonical edge word
    found = pair_present(rest, ((u, pn), (u, qn)))
    expl = [pair for pair in rest if explained_by_odd_construction(pair)]
    alt = [w for w in W if len(set(w[0::2])) == 1 == len(set(w[1::2]))
           and w[0] != w[1]]
    exch_is_rev = all(
        tuple(w[1] if i % 2 == 0 else w[0] for i in range(len(w))) == w[::-1]
        for w in alt)
    ok &= found and len(rest) == 12
    ok &= (len(alt) == 2) and exch_is_rev and (not expl)
    # the "2 alternating words" figure, derived rather than pinned: a word whose
    # two parity classes are monochromatic in distinct colors IS an alternation
    # of two colors, and at even length both label orders have multiplicities
    # (n/2, n/2), so both are canonical -- exhibit them and demand equality.
    alt_expected = set([
        tuple(1 if i % 2 == 0 else 2 for i in range(M_VERTICES - 1)),
        tuple(2 if i % 2 == 0 else 1 for i in range(M_VERTICES - 1))])
    ok &= (set(alt) == alt_expected) and alt_expected.issubset(set(W))
    ck("seven_vertex_uniform_census_finds_the_pair", ok,
       "%d paths (= the %d canonical edge words, one vertex word), %d equal "
       "pairs, %d reflections, %d left over; exhibited pair "
       "%s present; %d alternating canonical edge words %s, derived "
       "independently as the two label orders of an alternation and each with "
       "exchange == reverse, so the %d left-over pairs are unexplained "
       "structurally" %
       (npaths, len(W), eq, refl, len(rest), "IS" if found else "IS NOT",
        len(alt),
        sorted("".join(str(t) for t in w) for w in alt),
        len(rest) - len(expl)))
    return ok


def check_seven_vertex_two_edge_colors(p, q):
    """Check 21 (beyond the paper): all seven-vertex paths whose edge word uses
    at most two colors, over every canonical vertex coloring."""
    U = canonical_words(M_VERTICES)
    W = [w for w in canonical_words(M_VERTICES - 1) if len(set(w)) <= 2]
    SCOPE_FIGURES["canonical_vertex_colorings_7"] = len(U)
    SCOPE_FIGURES["two_color_edge_words_6"] = len(W)
    npaths, eq, refl, rest = run_census(M_VERTICES, U, W)
    bad = [pair for pair in rest if not explained_by_odd_construction(pair)]
    u = tuple([1] * M_VERTICES)
    found = pair_present(bad, ((u, to_numeric(p)), (u, to_numeric(q))))
    ok = found and len(bad) == 48 and len(rest) == 84 and npaths == 482202
    ok &= (npaths == len(U) * len(W))    # the path count, derived not pinned
    ck("seven_vertex_two_edge_color_census", ok,
       "%d paths (= %d canonical vertex colorings x %d two-color edge words), "
       "%d equal pairs, %d reflections, %d left over, %d of those "
       "violate the quoted edge condition; exhibited pair among them: %s" %
       (npaths, len(U), len(W), eq, refl, len(rest), len(bad), found))
    return ok


def check_scope_note_figures():
    """The first NOT RE-RUN line is itself a quantitative claim -- it says the
    census that was NOT run is 11,481 x 1,602 = 18,392,562 paths, and that what
    WAS run covers 1,602 canonical edge words and the 42 edge words using at
    most two colors.  Those four counts must not be prose the program merely
    repeats: each is here the size of a set this program actually built while
    running the two seven-vertex sub-censuses, recorded at the moment it was
    built, and the product is multiplied out.  A disagreement between the note
    and the program is a FAIL, not an unnoticed inconsistency in the text."""
    parts = []
    ok = True
    for key in ("canonical_vertex_colorings_7", "canonical_edge_words_6",
                "two_color_edge_words_6"):
        got = SCOPE_FIGURES.get(key)
        want = NOTE_FIGURES[key]
        ok &= (got == want)
        parts.append("%s: derived %s, note %d" % (key, got, want))
    a = SCOPE_FIGURES.get("canonical_vertex_colorings_7")
    b = SCOPE_FIGURES.get("canonical_edge_words_6")
    prod = a * b if (isinstance(a, int) and isinstance(b, int)) else None
    ok &= (prod == NOTE_FIGURES["un_run_product"])
    parts.append("un-run full census %s x %s = %s, note %d" %
                 (a, b, prod, NOTE_FIGURES["un_run_product"]))
    ck("scope_note_figures_are_derived_not_quoted", ok, "; ".join(parts))
    return ok


def check_hypotheses(p, q, target):
    """Check 22: every hypothesis of the clause being refuted holds for the
    exhibited pair, so the clause does apply to it."""
    m = M_VERTICES
    hyp = []
    hyp.append(("order is odd", m % 2 == 1))
    hyp.append(("both paths have %d vertices" % m,
                len(p) + 1 == m and len(q) + 1 == m))
    hyp.append(("vertex and edge alphabets disjoint",
                set(vertex_word_of(p)).isdisjoint(set(p) | set(q))))
    hyp.append(("vertex colorings identical",
                vertex_word_of(p) == vertex_word_of(q)))
    hyp.append(("edge-color multiplicities identical",
                sorted(p) == sorted(q)))
    K1 = kmatrix("a" * m, p, VAROF)
    K2 = kmatrix("a" * m, q, VAROF)
    hyp.append(("determinant polynomials equal",
                det_leibniz(K1) == det_leibniz(K2) == target))
    hyp.append(("the two paths are distinct objects", p != q))
    ok = all(v for _, v in hyp)
    ck("hypotheses_of_the_refuted_clause", ok,
       "; ".join("%s=%s" % (n, v) for n, v in hyp))
    return ok


def run_all():
    p, q = P_WORD, Q_WORD
    print("# exhibited object, decoded from the paper")
    print("#   vertices %d, all of color a; edge words p=%s q=%s"
          % (M_VERTICES, p, q))

    check_wellformed(p, q)
    check_multiplicities(p, q)
    Ms = check_matching_lists()
    target = check_determinants(p, q, Ms)
    check_A_values(p, q, Ms)
    check_hypotheses(p, q, target)
    check_not_identity_or_reflection(p, q)
    check_parity_table(p, q)
    check_construction_fails(p, q)
    check_positive_definite_region(p, q, target)
    check_small_orders()
    check_census_kernel()
    rest_by_m = {}
    for m in (3, 4, 5, 6):
        rest_by_m[m] = check_census_row(m)
    check_odd_rows_are_explained(rest_by_m)
    check_even_rows_have_parity_structure(rest_by_m)
    check_even_analogue_selectivity(rest_by_m)
    check_seven_vertex_subcensus(p, q)
    check_seven_vertex_two_edge_colors(p, q)
    check_scope_note_figures()          # must run last: it audits the two above


def main():
    try:
        run_all()
    except Exception as exc:                       # a corrupted input must not
        ck("completed_without_unhandled_error", False,   # break the contract
           "%s: %s" % (type(exc).__name__, exc))
    for line in NOTES:
        print("NOT RE-RUN: " + line)
    nfail = sum(1 for _, ok in CHECKS if not ok)
    if nfail:
        print("VERDICT: %d OF %d CHECKS FAILED" % (nfail, len(CHECKS)))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % len(CHECKS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
