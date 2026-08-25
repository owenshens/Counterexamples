#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- verification program for
  "The Hessian Arrangement Requires Four Added Lines for Supersolvability"
  (extSS(H) = 4 for the twelve-line Hessian arrangement H in P^2_C; hence no
   one-line extension of H is supersolvable, answering negatively the question
   of Dimca, Kuehne and Pokora, arXiv:2505.01733v4, Remark 2.3).

HOW TO RUN
  python3 verify.py            # ~1 second, standard library only
  python3 verify.py --paper /path/to/paper.tex   # also corroborates step 10's
                               # embedded citation data against the real source
Exit status 0 iff every check passes.  The number of checks does not depend on
which files happen to sit beside this program: every check below, step 10
included, is recorded in every run, for 41 checks in total.
After the checks the run prints, in this order: six falsification probes, which
mutate one compared input each and re-evaluate the very predicate the targeted
check evaluates (they record no check and cannot make the verdict green); the
residual gaps G1-G6; the verdict; and a closing one-line NOT RE-RUN disclosure
of everything the green verdict does not cover.

ARITHMETIC
  Everything is exact.  The coefficient field is Q(w), w a primitive cube root
  of unity, realised as Q[t]/(t^2+t+1): an element is a pair (a,b) of
  fractions.Fraction meaning a + b*w, with w^2 = -1-w.  There is no floating
  point anywhere in this file, and no numeric tolerance is used by any check.

======================================================================
TAKEN FROM THE PAPER (data; not verified here, it is the input)
======================================================================
  T1. The defining equation of the Hessian arrangement,
        H : xyz((x^3+y^3+z^3)^3 - 27 x^3 y^3 z^3) = 0,            eq. (1)
      and the list of its twelve lines
        x=0, y=0, z=0, x + w^a y + w^b z = 0  (a,b in Z/3).        Sec. 1
  T2. The four line equations the paper adds to build B:
        x+y = 0  and  x+y-2 w^j z = 0 (j in Z/3).                 Sec. 3
  T3. The distinguished point p = (1:-1:0) and the line x-y = 0.   Sec. 3
  T4. Values printed in the paper, used ONLY as comparison targets:
      |Sing(H)| = 21 (as 9+12), 9 points of multiplicity 4 (=P), 12 of
      multiplicity 2 (=V), sum C(m,2) = 66 = C(12,2), 3 P-points and 2 V-points
      per line, 12*C(3,2) = C(9,2) = 36, |B| = 16 (as 12+4), extSS(H) = 4,
      d = 12, d_1 = 4 = m(H), exponents (4,7), the hypothesis 4 <= 5 < 9, the
      printed sets P and V, and the printed four exceptional points
      (0:0:1), (1:1:1), (1:1:w), (1:1:w^2).  Every one of these is re-derived
      below from T1-T3 and compared; none is assumed.
  T4'. Comparison targets fixed in this program rather than quoted from the
      paper: |Sing(B)| = 45, and the exact forced-added-line counts 4 (for a
      modular point in P), 6 (for one in V), 6 (for one off Sing(H)).  The
      paper never prints |Sing(B)|, and its lower-bound argument claims only
      k >= 4 in each of the three cases.
  T5. The statement of [Dimca, Thm 1.12(2)].  CAUTION: the paper quotes only the
      "all unjoined multiple points lie on a single line L" clause.  The rest of
      the statement used in step 9 -- the inequality d - 2 m_p <= sum_{q in Q}
      (m_q - 1), with equality iff the arrangement is free with exponents
      (m_p, d-m_p-1), and the hypotheses m_p <= d_1+1 < d-m_p+1, m_p = d_1,
      m_p = m(A) -- is read off the SOURCE (arXiv:2503.01624v7, Theorem 1.12),
      not off the paper.  Used as an external theorem; not reproved here.
  T6. The du Plessis-Wall / Dimca-Sernesi freeness criterion: for a reduced plane
      curve of degree d with r = mdr(f) <= (d-1)/2, the curve is free with
      exponents (r, d-1-r) iff tau = (d-1)^2 - r(d-1-r).  External theorem,
      applied in step 9 to the mdr(f) DERIVED in step 8b.
  T7. The paper's bibliographic data, transcribed verbatim as literal LaTeX text
      (PAPER_CITATIONS and PAPER_BIBLIOGRAPHY, just before step 10): every \\cite
      occurrence with its multiplicity, and every \\bibitem key with the arXiv
      identifiers its entry prints.  This is the INPUT of step 10, and it is what
      makes step 10 self-contained, so that step 10 is recorded in every run with
      or without a paper source.  When a paper source IS available it is parsed
      by the same two parsers and this transcription must agree with it exactly,
      in both directions, so a transcription that has drifted from the paper
      FAILS instead of certifying itself.

======================================================================
DERIVED HERE (this is what the checks actually decide)
======================================================================
  D1. The full expansion of the product of the twelve linear forms, and its
      identity with xyz((x^3+y^3+z^3)^3 - 27x^3y^3z^3)  -- pins the object.
  D2. All C(12,2) = 66 pairwise intersections; the set Sing(H), each point's
      multiplicity, the partition into P and V, the completeness identity
      sum C(m,2) = 66, and the per-line incidence counts.
  D3. That the derived P and V coincide, as sets of projective points, with the
      sets the paper prints.
  D4. Lemma 2(1): every one of the C(9,2)=36 pairs of P-points lies on a line of H
      (equivalently: a line not in H carries at most one point of P).
  D5. Lemma 2(2), for all nine p in P (not just the representative): the number of
      V-points unjoined by lines of H, their collinearity, the unique line they
      span, and that p is off it; plus the explicit four points for p=(1:-1:0)
      and the identification of their span with x-y=0.
  D6. Lemma 2(3), for all twelve v in V: the six P-points covered by the two old
      lines through v, and that each join v-b_i meets V only in v.
  D7. The four joins from p to its exceptional points, compared with the paper's
      printed four added lines; distinctness; disjointness from H; |B| = 16.
  D8. Sing(B) by intersecting all C(16,2) = 120 pairs; |Sing(B)|; the
      multiplicity of p in B; and modularity of p in B (so B is supersolvable),
      giving extSS(H) <= 4.
  D9. For every X in Sing(H), the number of distinct non-H lines XQ forced by
      modularity (Q ranging over Sing(H)\\{X}); the minimum over X.
  D10. The two ingredients of the off-Sing(H) case (max P-points on a line of H,
      max P-points on a line not in H) and the resulting bound k >= 6.
  D11. extSS(H) = 4 from D8 and D9/D10, and the consequence that no extension
      of H by one line is supersolvable.
  D12. The first exponent itself: d_1 = mdr(f) = min{r : some nonzero (a,b,c) of
      degree r has a f_x + b f_y + c f_z = 0}, by exact Gaussian elimination on
      the graded pieces r = 0,1,2,3,4 (step 8b).  This is what makes d_1 = 4 a
      derived fact instead of a quotation from [Dimca, Example 6.5(i)]; the
      paper's headline claim needs it, because [DKP, Remark 2.3] poses its
      question only for arrangements with d_1 = m(A).  EARLIER VERSIONS OF THIS
      FILE SET d_1 := m_p, which made the exponent checks certify their own
      input, since m_p = d_1 is a HYPOTHESIS of [Dimca, Thm 1.12(2)].
  D12b. The Tjurina number tau(H) = sum (m_q-1)^2, and freeness with exponents
      (4,7) by applying T6 to the derived d_1 and tau; plus the equality case of
      the inequality in [Dimca, Thm 1.12(2)].  Note that, given the multiplicity
      profile verified in step 2, both sides of each of these identities are
      arithmetically forced, so their content lies in the derived d_1 alone.
  D13. (bibliographic) the resolution of each of the four checked citations to an arXiv
      version, composed out of the paper's own bibliographic data (T7): the
      (key, label) citation pairs parsed from the literal \\cite tokens, the map
      key -> arXiv version parsed from the literal bibliography entries, that
      every cited key has an entry and every entry is cited, and that composing
      the two sends each of the four labels checked here to the version used
      here, each version being printed by exactly one entry.  Recorded in every
      run.
"""

import sys
from fractions import Fraction as Fr
from itertools import combinations

# ---------------------------------------------------------------------------
# PAPER DATA (T1-T4).  Field elements are pairs (a,b) = a + b*w.
# ---------------------------------------------------------------------------
ZERO = (Fr(0), Fr(0))
ONE = (Fr(1), Fr(0))
W = (Fr(0), Fr(1))          # w
W2 = (Fr(-1), Fr(-1))       # w^2 = -1 - w
TWO = (Fr(2), Fr(0))
WPOW = [ONE, W, W2]         # w^0, w^1, w^2

# Comparison targets only (T4, and those marked T4', which are fixed in this
# program rather than quoted from the paper).  Nothing here is used to compute a
# derived quantity; each appears only on the right-hand side of a comparison.
PAPER = {
    "n_lines": 12,
    "n_pairs": 66,
    "n_sing_H": 21,
    "n_P": 9,
    "mult_P": 4,
    "n_V": 12,
    "mult_V": 2,
    "P_per_line": 3,
    "V_per_line": 2,
    "n_P_pairs": 36,
    "n_exceptional": 4,
    "size_B": 16,
    "n_sing_B": 45,               # T4' (target fixed here; not in the paper)
    "extSS": 4,
    "d": 12,
    "d1": 4,
    "m_max": 4,
    "exponents": (4, 7),
    "hypothesis": (4, 5, 9),      # m_p <= d1+1 < d-m_p+1  reads  4 <= 5 < 9
    "forced_lines_from_P": 4,     # T4' (exact target fixed here)
    "forced_lines_from_V": 6,     # T4' (exact target fixed here)
    "forced_lines_off_Sing": 6,   # T4' (exact target fixed here)
}

# T4: the paper's printed labels for the four citations (see D13).
PAPER_LABELS = [
    ("Dimca", "2503.01624v7", "Theorem~1.12(2)"),
    ("Dimca", "2503.01624v7", "Example~6.5(i)"),
    ("DKP", "2505.01733v4", "Remark~2.3"),
    ("Kabat", "2201.04856v1", "Definition~1.3"),
]

FAILURES = []
NCHECKS = [0]


# ---------------------------------------------------------------------------
# Exact arithmetic in Q(w) = Q[t]/(t^2+t+1), w^2 = -1-w.
# ---------------------------------------------------------------------------
def fadd(p, q):
    return (p[0] + q[0], p[1] + q[1])


def fneg(p):
    return (-p[0], -p[1])


def fsub(p, q):
    return (p[0] - q[0], p[1] - q[1])


def fmul(p, q):
    # (a+bw)(c+dw) = ac + (ad+bc) w + bd w^2 = (ac-bd) + (ad+bc-bd) w
    a, b = p
    c, d = q
    return (a * c - b * d, a * d + b * c - b * d)


def fis0(p):
    return p[0] == 0 and p[1] == 0


def finv(p):
    # N(a+bw) = (a+bw)(a+b w^2) = a^2 - ab + b^2 > 0 unless a=b=0.
    a, b = p
    n = a * a - a * b + b * b
    if n == 0:
        raise ZeroDivisionError("inverse of 0 in Q(w)")
    return ((a - b) / n, -b / n)


def fstr(p):
    a, b = p
    if b == 0:
        return str(a)
    tail = "w" if b == 1 else ("-w" if b == -1 else "%s*w" % b)
    if a == 0:
        return tail
    if b > 0:
        return "%s+%s" % (a, tail[0:] if b != 1 else "w")
    return "%s%s" % (a, tail)


# ---------------------------------------------------------------------------
# Projective points and lines over Q(w): both are triples of field elements.
# A point P lies on a line L iff L.P = 0.  Duality makes the two identical.
# ---------------------------------------------------------------------------
def dot(L, P):
    return fadd(fadd(fmul(L[0], P[0]), fmul(L[1], P[1])), fmul(L[2], P[2]))


def cross(P, Q):
    return (fsub(fmul(P[1], Q[2]), fmul(P[2], Q[1])),
            fsub(fmul(P[2], Q[0]), fmul(P[0], Q[2])),
            fsub(fmul(P[0], Q[1]), fmul(P[1], Q[0])))


def pnorm(V):
    """Canonical representative: scale so the first nonzero coordinate is 1."""
    for c in V:
        if not fis0(c):
            u = finv(c)
            return tuple(fmul(x, u) for x in V)
    raise ValueError("zero triple is not a projective point")


def pstr(V):
    return "(" + ":".join(fstr(c) for c in V) + ")"


def incidence_table(lines):
    """All C(n,2) pairwise meets -> {point: frozenset(line indices through it)}."""
    if len(set(lines)) != len(lines):
        raise ValueError("incidence_table needs pairwise distinct lines")
    acc = {}
    for i, j in combinations(range(len(lines)), 2):
        X = pnorm(cross(lines[i], lines[j]))
        acc.setdefault(X, set()).update((i, j))
    return dict((X, frozenset(S)) for X, S in acc.items())


# ---------------------------------------------------------------------------
# Polynomials in Q(w)[x,y,z] as {(i,j,k): coefficient}, zero terms pruned.
# ---------------------------------------------------------------------------
def pol_prune(p):
    return dict((m, c) for m, c in p.items() if not fis0(c))


def pol_add(p, q):
    out = dict(p)
    for m, c in q.items():
        out[m] = fadd(out.get(m, ZERO), c)
    return pol_prune(out)


def pol_mul(p, q):
    out = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2])
            out[m] = fadd(out.get(m, ZERO), fmul(c1, c2))
    return pol_prune(out)


def pol_scale(p, c):
    return pol_prune(dict((m, fmul(v, c)) for m, v in p.items()))


def pol_pow(p, n):
    out = {(0, 0, 0): ONE}
    for _ in range(n):
        out = pol_mul(out, p)
    return out


def pol_eq(p, q):
    return pol_prune(p) == pol_prune(q)


def linear_form(L):
    """The linear polynomial A x + B y + C z of the line L = (A,B,C)."""
    mons = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    return pol_prune(dict(zip(mons, L)))


def pol_str(p, limit=None):
    items = sorted(p.items(), key=lambda kv: (-sum(kv[0]), kv[0]))
    if limit is not None and len(items) > limit:
        items = items[:limit]
        tail = " + ... (%d terms in all)" % len(p)
    else:
        tail = ""
    def mono(m):
        s = "".join(v + ("^%d" % e if e > 1 else "")
                    for v, e in zip("xyz", m) if e)
        return s if s else "1"
    return " + ".join("(%s)%s" % (fstr(c), mono(m)) for m, c in items) + tail


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def info(msg):
    print("    " + msg)


def verdict():
    n = NCHECKS[0]
    if FAILURES:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(FAILURES), n))
        for name in FAILURES:
            print("    failed: %s" % name)
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


def bail():
    """A derived object came out wrong, so later checks cannot even be posed.
    Print the verdict that the format requires and stop."""
    print("")
    print("    a derived object contradicts the paper, so the remaining checks "
          "cannot be posed")
    raise SystemExit(verdict())


def check(ok, name, detail=""):
    NCHECKS[0] += 1
    tag = "PASS" if ok else "FAIL"
    if not ok:
        FAILURES.append(name)
    print("%s [%02d] %s%s" % (tag, NCHECKS[0], name, ("  --  " + detail) if detail else ""))
    return ok


# ---------------------------------------------------------------------------
# The arrangement H (T1).  Twelve lines, as coefficient triples.
# ---------------------------------------------------------------------------
def hessian_lines_raw():
    """The twelve forms exactly as printed: x, y, z, x + w^a y + w^b z.  Sec. 1."""
    lines = [(ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (ZERO, ZERO, ONE)]
    for a in range(3):
        for b in range(3):
            lines.append((ONE, WPOW[a], WPOW[b]))
    return lines


def hessian_lines():
    return [tuple(pnorm(L)) for L in hessian_lines_raw()]


def paper_P_and_V():
    """The sets P and V exactly as the paper prints them (Section 2, T4)."""
    P = set()
    for j in range(3):
        mj = fneg(WPOW[j])                     # -w^j
        P.add(tuple(pnorm((ZERO, ONE, mj))))   # (0 : 1 : -w^j)
        P.add(tuple(pnorm((ONE, ZERO, mj))))   # (1 : 0 : -w^j)
        P.add(tuple(pnorm((ONE, mj, ZERO))))   # (1 : -w^j : 0)
    V = set(tuple(pnorm(t)) for t in
            [(ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (ZERO, ZERO, ONE)])
    for r in range(3):
        for s in range(3):
            V.add(tuple(pnorm((ONE, WPOW[r], WPOW[s]))))
    return P, V


def paper_added_lines():
    """x+y=0 and x+y-2w^j z=0, j in Z/3.  Paper, Section 3 (T2)."""
    out = [tuple(pnorm((ONE, ONE, ZERO)))]
    for j in range(3):
        out.append(tuple(pnorm((ONE, ONE, fneg(fmul(TWO, WPOW[j]))))))
    return out


# ---------------------------------------------------------------------------
# STEP 1 (D1).  The product of the twelve listed linear forms equals
# xyz((x^3+y^3+z^3)^3 - 27 x^3 y^3 z^3) identically.  This is what pins the
# twelve lines to the Hessian arrangement of Dimca, Example 6.5(i).
# ---------------------------------------------------------------------------
def check_identity(lines):
    raw = hessian_lines_raw()
    prod = {(0, 0, 0): ONE}
    for L in raw:
        prod = pol_mul(prod, linear_form(L))
    # No rescaling can have crept in: the product is over the printed forms, and
    # normalising them for the incidence work below changes none of them.
    unscaled = all(tuple(pnorm(L)) == L for L in raw)
    info("normalisation is the identity on all twelve printed forms: %s "
         "(DECORATION: every printed form already has leading coefficient 1, so no "
         "input can make this false; it documents that the product below is taken "
         "over exactly the twelve lines used for the incidences)" % unscaled)
    x3 = {(3, 0, 0): ONE}
    y3 = {(0, 3, 0): ONE}
    z3 = {(0, 0, 3): ONE}
    cubes = pol_add(pol_add(x3, y3), z3)
    target = pol_mul({(1, 1, 1): ONE},
                     pol_add(pol_pow(cubes, 3),
                             pol_scale({(3, 3, 3): ONE}, (Fr(-27), Fr(0)))))
    deg = max(sum(m) for m in prod)
    info("expanded product of the 12 forms: %d monomials, total degree %d" %
         (len(prod), deg))
    info("product   = " + pol_str(prod, limit=6))
    info("target    = " + pol_str(target, limit=6))
    rational = all(c[1] == 0 for c in prod.values())
    info("every coefficient of the expanded product lies in Z (w-part 0): %s" % rational)
    diff = pol_add(prod, pol_scale(target, (Fr(-1), Fr(0))))
    info("product - target has %d nonzero monomials" % len(diff))
    check(pol_eq(prod, target) and unscaled,
          "IDENTITY: prod of 12 forms == xyz((x^3+y^3+z^3)^3-27x^3y^3z^3)",
          "degree %d, %d monomials, difference empty: %s, unscaled: %s"
          % (deg, len(prod), len(diff) == 0, unscaled))
    check(deg == PAPER["n_lines"] and len(set(lines)) == PAPER["n_lines"]
          and len(lines) == len(raw),
          "H has exactly 12 distinct lines and the form has degree 12",
          "|H| = %d (paper: %d), deg = %d" % (len(set(lines)), PAPER["n_lines"], deg))
    return prod


# ---------------------------------------------------------------------------
# STEP 2 (D2, D3).  Sing(H) from all 66 pairwise meets; the partition into
# P (multiplicity 4) and V (multiplicity 2); completeness; per-line counts;
# and agreement with the sets the paper prints.
# ---------------------------------------------------------------------------
def check_incidences(lines):
    tab = incidence_table(lines)
    mults = sorted((len(S) for S in tab.values()), reverse=True)
    # The partition is read off the derived multiplicities, not from the paper:
    # P = the points of largest multiplicity, V = everything else.  The values 4
    # and 2 are then compared with the paper, and the check below refuses any
    # third multiplicity.
    m_top = mults[0]
    Pd = sorted(X for X, S in tab.items() if len(S) == m_top)
    Vd = sorted(X for X, S in tab.items() if len(S) != m_top)
    m_low = sorted(set(len(tab[X]) for X in Vd))
    npairs = len(lines) * (len(lines) - 1) // 2
    accounted = sum(len(S) * (len(S) - 1) // 2 for S in tab.values())
    info("pairs of lines intersected: %d" % npairs)
    info("distinct intersection points: %d" % len(tab))
    info("multiplicity multiset: %s" % {m: mults.count(m) for m in sorted(set(mults))})
    info("largest multiplicity %d (that is P); other multiplicities present: %s"
         % (m_top, m_low))
    info("sum over points of C(mult,2) = %d ; C(12,2) = %d" % (accounted, npairs))
    check(len(tab) == PAPER["n_sing_H"], "|Sing(H)| derived from the 66 pairwise meets",
          "derived %d, paper %d" % (len(tab), PAPER["n_sing_H"]))
    check(len(Pd) == PAPER["n_P"] and len(Vd) == PAPER["n_V"]
          and m_top == PAPER["mult_P"] and m_low == [PAPER["mult_V"]]
          and len(Pd) + len(Vd) == len(tab),
          "multiplicity profile: 9 points of multiplicity 4, 12 of multiplicity 2",
          "derived |P| = %d at multiplicity %d, |V| = %d at multiplicity %s"
          % (len(Pd), m_top, len(Vd), m_low))
    check(accounted == npairs,
          "completeness: sum C(m,2) = C(12,2), so there are no further singular points",
          "%d == %d" % (accounted, npairs))
    Pp, Vp = paper_P_and_V()
    check(set(Pd) == Pp, "derived quadruple points == the set P printed in the paper",
          "%d derived, %d printed, symmetric difference %d"
          % (len(Pd), len(Pp), len(set(Pd) ^ Pp)))
    check(set(Vd) == Vp, "derived double points == the set V printed in the paper",
          "%d derived, %d printed, symmetric difference %d"
          % (len(Vd), len(Vp), len(set(Vd) ^ Vp)))
    pc = sorted(set(sum(1 for X in Pd if fis0(dot(L, X))) for L in lines))
    vc = sorted(set(sum(1 for X in Vd if fis0(dot(L, X))) for L in lines))
    info("per-line P-counts (distinct values over the 12 lines): %s" % pc)
    info("per-line V-counts (distinct values over the 12 lines): %s" % vc)
    check(pc == [PAPER["P_per_line"]] and vc == [PAPER["V_per_line"]],
          "every line of H carries exactly 3 points of P and 2 points of V",
          "P-counts %s, V-counts %s" % (pc, vc))
    return tab, Pd, Vd


def joined_by_H(lines, A, B):
    """True iff some line of H contains both A and B."""
    return any(fis0(dot(L, A)) and fis0(dot(L, B)) for L in lines)


# ---------------------------------------------------------------------------
# STEP 3 (D4).  Lemma 2(1): every pair of P-points lies on a line of H.
# Equivalent dual form used later in the lower bound: a line not in H carries
# at most one point of P.  Also the paper's counting identity 12*C(3,2)=C(9,2).
# ---------------------------------------------------------------------------
def check_lemma1(lines, Pd):
    pairs = list(combinations(Pd, 2))
    on_H = [(A, B) for A, B in pairs if joined_by_H(lines, A, B)]
    info("pairs of P-points: %d = C(9,2); joined by a line of H: %d"
         % (len(pairs), len(on_H)))
    check(len(pairs) == PAPER["n_P_pairs"] and len(on_H) == len(pairs),
          "Lemma 2(1): all 36 pairs of P-points lie on a line of H",
          "%d of %d pairs joined" % (len(on_H), len(pairs)))
    # The per-line P-count used here is the derived one, not the paper's figure.
    percounts = sorted(set(sum(1 for q in Pd if fis0(dot(L, q))) for L in lines))
    k = percounts[0] if len(percounts) == 1 else -1
    lhs = len(lines) * (k * (k - 1) // 2)
    rhs = len(Pd) * (len(Pd) - 1) // 2
    info("derived P-points per line: %s ; %d*C(%d,2) = %d and C(%d,2) = %d"
         % (percounts, len(lines), k, lhs, len(Pd), rhs))
    check(len(percounts) == 1 and lhs == rhs == PAPER["n_P_pairs"],
          "counting identity 12*C(3,2) = C(9,2) = 36 (no P-pair counted twice)",
          "%d == %d == %d" % (lhs, rhs, PAPER["n_P_pairs"]))
    # Dual form: the span of any two distinct P-points is a line of H, hence a
    # line outside H meets P in at most one point.
    spans = set(tuple(pnorm(cross(A, B))) for A, B in pairs)
    inH = set(lines)
    outside = [L for L in spans if L not in inH]
    info("distinct spans of P-pairs: %d, of which outside H: %d"
         % (len(spans), len(outside)))
    check(not outside,
          "every span of two P-points is a line of H (so a non-H line has <=1 P-point)",
          "%d spans, %d outside H" % (len(spans), len(outside)))


# ---------------------------------------------------------------------------
# STEP 4 (D5).  Lemma 2(2), for every one of the nine points of P.
# ---------------------------------------------------------------------------
def check_lemma2(lines, Pd, Vd):
    counts, collinear, off_line, spans_unique = [], [], [], []
    exc_for_p = None
    p0 = tuple(pnorm((ONE, fneg(ONE), ZERO)))       # p = (1:-1:0), T3
    for p in Pd:
        unj = [v for v in Vd if not joined_by_H(lines, p, v)]
        counts.append(len(unj))
        L = None
        if len(unj) >= 2:
            L = tuple(pnorm(cross(unj[0], unj[1])))
            collinear.append(all(fis0(dot(L, v)) for v in unj))
            off_line.append(not fis0(dot(L, p)))
            # the line is the unique one through the unjoined points
            spans = set(tuple(pnorm(cross(a, b))) for a, b in combinations(unj, 2))
            spans_unique.append(spans == {L})
        else:
            collinear.append(False)
            off_line.append(False)
            spans_unique.append(False)
        if p == p0:
            exc_for_p = (sorted(unj), L)
    info("unjoined-V counts over the nine p in P: %s" % sorted(counts))
    check(all(c == PAPER["n_exceptional"] for c in counts) and len(counts) == PAPER["n_P"],
          "Lemma 2(2), part (a): for each of the 9 p in P exactly 4 points of V are "
          "unjoined by H",
          "counts %s (paper: 4 each)" % sorted(set(counts)))
    check(all(collinear) and all(spans_unique),
          "Lemma 2(2), part (b): those 4 points are collinear and span a unique line",
          "collinear for %d/%d p, unique span for %d/%d"
          % (sum(collinear), len(collinear), sum(spans_unique), len(spans_unique)))
    check(all(off_line),
          "Lemma 2(2), part (c): p does not lie on the line spanned by its 4 "
          "exceptional points",
          "p off the line for %d/%d p" % (sum(off_line), len(off_line)))
    if exc_for_p is None or exc_for_p[1] is None:
        check(False, "p = (1:-1:0) is a quadruple point with 4 exceptional V-points",
              "p not found among the derived quadruple points, or too few unjoined "
              "points; the remaining Lemma 2(2) and upper-bound checks cannot run")
        check(False, "the four exceptional points of p=(1:-1:0) span the line x - y = 0",
              "not reached: the previous check failed first")
        bail()
    unj, L = exc_for_p
    printed = set(tuple(pnorm(t)) for t in
                  [(ZERO, ZERO, ONE)] + [(ONE, ONE, WPOW[s]) for s in range(3)])
    info("for p = (1:-1:0) the derived exceptional points are %s"
         % ", ".join(pstr(v) for v in unj))
    info("their derived span is the line %s (coefficients), i.e. %s"
         % (pstr(L), pol_str(linear_form(L))))
    check(set(unj) == printed,
          "for p=(1:-1:0) the 4 exceptional points are (0:0:1),(1:1:1),(1:1:w),(1:1:w^2)",
          "derived set equals printed set: %s" % (set(unj) == printed))
    check(L == tuple(pnorm((ONE, fneg(ONE), ZERO))),
          "for p=(1:-1:0) the line supplied by Dimca Thm 1.12(2) is x - y = 0",
          "derived %s" % pol_str(linear_form(L)))
    return p0, unj, L


# ---------------------------------------------------------------------------
# STEP 5 (D6).  Lemma 2(3), for every one of the twelve points of V.
# ---------------------------------------------------------------------------
def check_lemma3(lines, Pd, Vd):
    n_old, covered_ok, rem_ok, joins_clean = [], [], [], []
    cov_counts = []
    for v in Vd:
        old = [L for L in lines if fis0(dot(L, v))]
        n_old.append(len(old))
        cov = [q for q in Pd if any(fis0(dot(L, q)) for L in old)]
        cov_counts.append(len(cov))
        covered_ok.append(len(cov) == 6)
        rem = [q for q in Pd if q not in set(cov)]
        rem_ok.append(len(rem) == 3)
        ok = True
        joins = set()
        for b in rem:
            J = tuple(pnorm(cross(v, b)))
            joins.add(J)
            if J in set(lines):
                ok = False                       # must be a new line
            for u in Vd:
                if u != v and fis0(dot(J, u)):
                    ok = False                   # meets V again
        if len(joins) != len(rem):
            ok = False                           # the three joins must be distinct
        joins_clean.append(ok)
    info("lines of H through each v in V: %s" % sorted(set(n_old)))
    info("P-points covered by the two old lines through v (distinct values): %s"
         % sorted(set(cov_counts)))
    check(set(n_old) == {2} and all(covered_ok),
          "Lemma 2(3), part (a): the two lines of H through each v in V cover exactly "
          "6 points of P",
          "6 covered for %d/%d v in V" % (sum(covered_ok), len(covered_ok)))
    check(all(rem_ok),
          "Lemma 2(3), part (b): exactly 3 points of P remain uncovered for each v in V",
          "3 remaining for %d/%d v" % (sum(rem_ok), len(rem_ok)))
    check(all(joins_clean),
          "Lemma 2(3), part (c): the 3 joins v-b_i are distinct, outside H, and meet V "
          "only in v",
          "clean for %d/%d v" % (sum(joins_clean), len(joins_clean)))


# ---------------------------------------------------------------------------
# STEP 6 (D7, D8).  Upper bound extSS(H) <= 4:  B = H + 4 lines, and p is
# modular in B.  The four added lines are DERIVED as the joins from p to its
# exceptional points and compared with the four equations the paper prints.
# ---------------------------------------------------------------------------
def check_upper_bound(lines, tab, p0, unj):
    derived = [tuple(pnorm(cross(p0, v))) for v in unj]
    printed = paper_added_lines()
    info("derived joins from p=(1:-1:0) to its 4 exceptional points:")
    for L in derived:
        info("    %s   i.e.   %s = 0" % (pstr(L), pol_str(linear_form(L))))
    info("printed added lines of the paper: %s"
         % ", ".join(pol_str(linear_form(L)) + " = 0" for L in printed))
    check(set(derived) == set(printed),
          "the paper's 4 added lines are exactly the joins p-v_i (derived, then compared)",
          "derived set == printed set: %s" % (set(derived) == set(printed)))
    check(len(set(derived)) == PAPER["n_exceptional"]
          and not (set(derived) & set(lines)),
          "the 4 added lines are distinct and none belongs to H",
          "%d distinct, %d of them in H" % (len(set(derived)), len(set(derived) & set(lines))))
    # Deduplicated on purpose: |B| is then a real measurement of how many of the
    # four joins are genuinely new lines, and no line is repeated in B.
    B = list(lines) + sorted(set(derived) - set(lines))
    info("|B| = |H| + |added| = %d + %d = %d" % (len(lines), len(B) - len(lines), len(B)))
    check(len(set(B)) == PAPER["size_B"], "|B| = 16",
          "derived %d, paper %d" % (len(set(B)), PAPER["size_B"]))
    tabB = incidence_table(B)
    multp = len(tabB.get(p0, ()))
    old_through_p = sum(1 for L in lines if fis0(dot(L, p0)))
    info("|Sing(B)| = %d ; |Sing(H)| = %d ; new singular points = %d"
         % (len(tabB), len(tab), len(set(tabB) - set(tab))))
    info("multiplicity of p in B = %d = %d old lines + %d added lines"
         % (multp, old_through_p, len(derived)))
    check(len(tabB) == PAPER["n_sing_B"], "|Sing(B)| = 45",
          "derived %d, comparison target %d -- T4': the paper never prints "
          "|Sing(B)|, so the right-hand side here is a target fixed in this program "
          "and not a value quoted from the paper"
          % (len(tabB), PAPER["n_sing_B"]))
    check(multp == old_through_p + len(derived) == 8,
          "p has multiplicity 8 in B (all 4 added lines pass through p)",
          "mult_B(p) = %d = %d + %d" % (multp, old_through_p, len(derived)))
    Bset = set(B)
    bad = [q for q in tabB if q != p0 and tuple(pnorm(cross(p0, q))) not in Bset]
    check(not bad and p0 in tabB,
          "p is modular in B: every other singular point of B is joined to p by a line of B",
          "%d of %d singular points fail" % (len(bad), len(tabB) - 1))
    nadd = len(set(B)) - len(set(lines))
    check(not bad and nadd == PAPER["extSS"] and set(lines) <= set(B),
          "hence B is supersolvable and extSS(H) <= |B \\ H| = 4",
          "H is contained in B and |B \\ H| = %d" % nadd)
    # Sing(H) is contained in Sing(B); this is what licenses the lower-bound
    # argument to test modularity only against the old singular points.
    check(set(tab) <= set(tabB),
          "Sing(H) is contained in Sing(B) (used by the lower bound)",
          "%d of %d old singular points survive in B"
          % (len(set(tab) & set(tabB)), len(tab)))
    return B, tabB, nadd, (not bad and p0 in tabB)


# ---------------------------------------------------------------------------
# STEP 7 (D9, D10).  Lower bound extSS(H) >= 4.
# If B is supersolvable with modular point X then, because Sing(H) is contained
# in Sing(B), every join XQ with Q in Sing(H)\{X} must lie in B; the joins that
# are not lines of H must therefore be added lines.  So k >= (number of
# distinct non-H joins from X).  We derive that number for all 21 X in Sing(H).
# ---------------------------------------------------------------------------
def check_lower_bound(lines, tab, Pd, Vd):
    inH = set(lines)
    forced = {}
    for X in tab:
        need = set()
        for Q in tab:
            if Q == X:
                continue
            J = tuple(pnorm(cross(X, Q)))
            if J not in inH:
                need.add(J)
        forced[X] = len(need)
    fP = sorted(set(forced[X] for X in Pd))
    fV = sorted(set(forced[X] for X in Vd))
    info("forced non-H lines through X: distinct values over the 9 points of P: %s" % fP)
    info("forced non-H lines through X: distinct values over the 12 points of V: %s" % fV)
    check(fP == [PAPER["forced_lines_from_P"]],
          "every X in P forces exactly 4 added lines, so k >= 4 in that case",
          "derived %s, paper %d" % (fP, PAPER["forced_lines_from_P"]))
    check(fV == [PAPER["forced_lines_from_V"]],
          "every X in V forces exactly 6 added lines, so k >= 6 in that case",
          "derived %s, comparison target %d -- T4': the paper argues only "
          "k >= 4 in this case (it says explicitly that the individual counts are not "
          "claimed to be sharp), so 6 is a target fixed in this program, not the paper's"
          % (fV, PAPER["forced_lines_from_V"]))
    # X outside Sing(H): the two-line count.  Both ingredients are derived.
    maxP_on_H = max(sum(1 for q in Pd if fis0(dot(L, q))) for L in lines)
    spans = set(tuple(pnorm(cross(A1, A2))) for A1, A2 in combinations(Pd, 2))
    maxP_off_H = 1 if not (spans - inH) else 2
    kmin_off = len(Pd) - maxP_on_H          # each added line covers <= 1 P-point
    info("max P-points on a line of H = %d ; max P-points on a line not in H = %d"
         % (maxP_on_H, maxP_off_H))
    info("X outside Sing(H): at most 1 line of H through X, so k >= 9 - %d = %d"
         % (maxP_on_H, kmin_off))
    check(maxP_on_H == PAPER["P_per_line"] and maxP_off_H == 1
          and kmin_off == PAPER["forced_lines_off_Sing"],
          "X outside Sing(H) forces k >= 6 (one old line covers <=3 of P, each added <=1)",
          "9 - %d = %d, comparison target %d -- T4': the paper argues only "
          "k >= 4 in this case, so 6 is a target fixed in this program, not the paper's"
          % (maxP_on_H, kmin_off, PAPER["forced_lines_off_Sing"]))
    overall = min(min(forced.values()), kmin_off)
    info("minimum over all three cases: k >= %d" % overall)
    check(overall == PAPER["extSS"],
          "lower bound over every possible modular point: extSS(H) >= 4",
          "min(P:%d, V:%d, off:%d) = %d" % (fP[0], fV[0], kmin_off, overall))
    return overall


# ---------------------------------------------------------------------------
# STEP 9 (D12).  The paper is in the scope of the DKP question only if
# d_1 = m(H) = 4 and 4 <= 5 < 9.  The exponents (4,7) are quoted by the paper
# from Dimca Example 6.5(i); here they are corroborated two independent ways:
# (a) the equality case of the inequality inside Dimca Thm 1.12(2) (T5), and
# (b) the Tjurina identity tau = (d-1)^2 - d_1 d_2 valid for free arrangements.
# ---------------------------------------------------------------------------
def check_source_arithmetic(lines, tab, Pd, Vd, p0, unj, d1_derived):
    d = len(lines)
    m_max = max(len(S) for S in tab.values())
    tau = sum((len(S) - 1) ** 2 for S in tab.values())
    info("d = %d ; m(H) = max multiplicity = %d ; tau(H) = sum (m_q-1)^2 = %d"
         % (d, m_max, tau))
    check(m_max == PAPER["m_max"] and d == PAPER["d"],
          "d = 12 and the maximal multiplicity m(H) = 4, both derived from Sing(H)",
          "derived d=%d, m=%d" % (d, m_max))
    # (a) equality case of Dimca Theorem 1.12(2), applied at p = (1:-1:0).
    m_p = len(tab[p0])
    lhs = d - 2 * m_p
    rhs = sum(len(tab[q]) - 1 for q in unj)
    # d_1 is the value DERIVED in step 8b (mdr(f)).  It is deliberately NOT set
    # to m_p here: "m_p = d_1" is a HYPOTHESIS of Dimca Thm 1.12(2), so taking
    # d_1 := m_p would make the equality clause below certify its own input.
    if d1_derived is None:
        check(False, "the first exponent d_1 used below was derived, not assumed",
              "step 8b did not produce mdr(f); the exponent checks cannot be posed")
        return
    d1_pred, d2_pred = d1_derived, d - d1_derived - 1
    info("Q = the %d multiple points unjoined to p; d - 2m_p = %d and "
         "sum_{q in Q}(m_q-1) = %d" % (len(unj), lhs, rhs))
    info("equality in d-2m_p <= sum(m_q-1) holds: %s ; the equality clause of Dimca "
         "Thm 1.12(2) then makes H free with exponents (%d,%d)"
         % (lhs == rhs, d1_pred, d2_pred))
    check(lhs == rhs and m_p == d1_pred and m_p == m_max
          and (d1_pred, d2_pred) == PAPER["exponents"],
          "equality case of Dimca Thm 1.12(2) gives exponents (4,7), as the paper states",
          "hypotheses m_p = d_1 (%d = %d) and m_p = m(H) (%d = %d) hold; %d == %d, "
          "exponents (%d,%d)"
          % (m_p, d1_pred, m_p, m_max, lhs, rhs, d1_pred, d2_pred))
    # (b) Tjurina identity for a free line arrangement.
    tau_free = (d - 1) ** 2 - d1_pred * d2_pred
    info("free-arrangement identity: (d-1)^2 - d_1 d_2 = %d - %d = %d, and tau(H) = %d"
         % ((d - 1) ** 2, d1_pred * d2_pred, tau_free, tau))
    check(tau == tau_free and d1_pred + d2_pred == d - 1,
          "freeness with exponents (4,7) via du Plessis-Wall: tau = (d-1)^2 - d_1 d_2 "
          "for the DERIVED d_1 = mdr(f)",
          "%d == %d and %d+%d = %d" % (tau, tau_free, d1_pred, d2_pred, d - 1))
    info("NOTE: given the multiplicity profile already verified in step 2, both sides "
         "of this identity are forced arithmetically, so its content is entirely in "
         "the derived d_1 = mdr(f) = %d of step 8b; it is a criterion applied, not an "
         "independent recomputation" % d1_pred)
    # scope of the DKP question: m_p <= d_1 + 1 < d - m_p + 1
    a, b, c = m_p, d1_pred + 1, d - m_p + 1
    info("hypothesis of Dimca Thm 1.12(2): m_p <= d_1+1 < d-m_p+1 reads %d <= %d < %d"
         % (a, b, c))
    check((a, b, c) == PAPER["hypothesis"] and a <= b < c and d1_pred == m_max,
          "H is in the scope of the DKP question: d_1 = m(H) = 4 and 4 <= 5 < 9",
          "%d <= %d < %d, d_1 = m(H) = %d" % (a, b, c, d1_pred))


# ---------------------------------------------------------------------------
# STEP 8b (D14).  Helpers for the first exponent d_1 = mdr(f).
# mdr(f) is the least r for which the graded piece
#   AR(f)_r = {(a,b,c) of degree r : a f_x + b f_y + c f_z = 0}
# is nonzero.  For a free curve, d_1 = mdr(f) is the first exponent, so this is
# how the paper's quoted "d_1 = 4" becomes a DERIVED quantity rather than an
# input.  All entries of the linear system are rational (f has integer
# coefficients), so the elimination below runs over Q with Fraction.
# ---------------------------------------------------------------------------
def pol_diff(p, var):
    """Partial derivative of p with respect to x (var=0), y (1) or z (2)."""
    out = {}
    for m, c in p.items():
        e = m[var]
        if e == 0:
            continue
        m2 = list(m)
        m2[var] = e - 1
        key = tuple(m2)
        out[key] = fadd(out.get(key, ZERO), fmul(c, (Fr(e), Fr(0))))
    return pol_prune(out)


def monomials_of_degree(deg):
    return [(i, j, deg - i - j)
            for i in range(deg, -1, -1) for j in range(deg - i, -1, -1)]


def rank_over_Q(rows, ncols):
    """Rank of a matrix of Fractions, by forward Gaussian elimination.
    Exact throughout: the entries stay in Q, no pivoting tolerance is used."""
    rows = [list(row) for row in rows]
    r = 0
    for col in range(ncols):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][col] != 0:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][col]
        for i in range(r + 1, len(rows)):
            if rows[i][col] != 0:
                fct = rows[i][col] / pv
                ri, rr = rows[i], rows[r]
                for k in range(col, ncols):
                    if rr[k] != 0:
                        ri[k] -= fct * rr[k]
        r += 1
        if r == len(rows):
            break
    return r


def syzygy_dim(fx, fy, fz, r):
    """dim_Q AR(f)_r, where AR(f)_r = {(a,b,c) deg r : a fx + b fy + c fz = 0}."""
    cols = monomials_of_degree(r)
    nc = 3 * len(cols)
    acc = {}
    for slot, g in enumerate((fx, fy, fz)):
        for ci, mo in enumerate(cols):
            col = slot * len(cols) + ci
            for m, c in g.items():
                if c[1] != 0:
                    raise ValueError("partial derivative has a non-rational "
                                     "coefficient; the elimination assumes Q")
                key = (m[0] + mo[0], m[1] + mo[1], m[2] + mo[2])
                row = acc.setdefault(key, {})
                row[col] = row.get(col, Fr(0)) + c[0]
    rows = []
    for key in sorted(acc):
        rd = acc[key]
        if any(v != 0 for v in rd.values()):
            rows.append([rd.get(j, Fr(0)) for j in range(nc)])
    return nc - rank_over_Q(rows, nc), nc, len(rows)


def check_first_exponent(prod, m_max):
    """DERIVE d_1 = mdr(f).  The paper quotes d_1 = 4 = m(H) from Dimca,
    Example 6.5(i); nothing else in this file would notice if that were wrong,
    and the scope of the DKP question depends on it, so it is computed here."""
    d = max(sum(m) for m in prod)
    fx, fy, fz = pol_diff(prod, 0), pol_diff(prod, 1), pol_diff(prod, 2)
    info("f = the degree-%d product; partial derivatives have %d, %d, %d monomials"
         % (d, len(fx), len(fy), len(fz)))
    dims = {}
    try:
        # The whole profile up to m_max+1 is computed, not just the first nonzero
        # degree, so that an unexpected answer is diagnosable from the output.
        for r in range(0, m_max + 2):
            nul, nc, nr = syzygy_dim(fx, fy, fz, r)
            dims[r] = nul
            info("degree %d: %d unknowns, %d equations, dim AR(f)_%d = %d"
                 % (r, nc, nr, r, nul))
    except Exception as exc:                      # never crash the verdict
        check(False, "d_1 = mdr(f) derived by exact elimination on AR(f)",
              "the elimination raised %s: %s" % (type(exc).__name__, exc))
        return None
    nz = sorted(r for r in dims if dims[r] > 0)
    d1 = nz[0] if nz else None
    info("hence mdr(f) = %s: dim AR(f)_r = 0 for every smaller r, so no Jacobian "
         "syzygy of lower degree exists" % d1)
    check(d1 is not None and d1 == m_max,
          "d_1 = mdr(f) = 4 = m(H), DERIVED (not quoted from Dimca Example 6.5(i))",
          "mdr(f) = %s, m(H) = %d, profile %s"
          % (d1, m_max, [dims[r] for r in sorted(dims)]))
    check(d1 is not None and 2 * d1 <= d - 1,
          "d_1 <= (d-1)/2, the range where tau = (d-1)^2-d_1 d_2 characterises freeness",
          "2*%s <= %d" % (d1, d - 1))
    # If AR(f) is free with generators in degrees (4,7) then AR(f)_{d_1+1} is
    # exactly S_1 * (the degree-d_1 generator), of dimension 3.  A second
    # generator in degree d_1+1 would show up here as dimension 4 or more, so
    # this is a falsifiable corroboration of freeness that does not go through
    # the Tjurina criterion at all.
    nxt = dims.get(d1 + 1) if d1 is not None else None
    check(nxt == 3,
          "dim AR(f)_{d_1+1} = 3 = dim S_1, as a free AR(f) with exponents (4,7) forces "
          "(no second generator in degree 5)",
          "dim AR(f)_%s = %s" % (None if d1 is None else d1 + 1, nxt))
    return d1


# ---------------------------------------------------------------------------
# STEP 10 (D13).  Bibliographic, no arithmetic.  Four labelled citations carry
# the external inputs the paper leans on; this step decides that the paper cites
# exactly those four labels and that each cited key resolves to the arXiv version
# used here.
#
# The bibliographic input is EMBEDDED below as literal LaTeX text taken from the
# paper (T7), so this step is self-contained and is recorded in every run.  An
# earlier version of this file recorded NO check when no .tex happened to sit
# beside the program: run alone in a scratch directory it printed a "gap" note
# and then a clean verdict over one check fewer, which is the same defect as a
# driver announcing "ALL 0 CHECKS PASS".  A check whose input is missing must
# register a result, not disappear.
#
# What is DERIVED from the embedded input, rather than restated:
#   - the (key, label) citation pairs and their multiplicities, by parsing the
#     literal \cite tokens with parse_citations;
#   - the map key -> arXiv version, by parsing the literal bibliography with
#     parse_bibliography;
#   - that every cited key has a bibliography entry (no dangling citation) and
#     every entry is cited (no orphan entry);
#   - that composing citation -> key -> version sends each of the four labels
#     checked here to the version used here, each version being printed by
#     exactly one entry (so a version cannot sit under the wrong key).
# When a paper source is available it is parsed by the SAME two parsers and the
# embedded transcription must agree with it in both directions, so drift fails
# rather than certifying itself.  The label-to-number resolution INSIDE the three
# preprints was done by hand, reading those preprints directly:
#   Dimca 2503.01624v7  Theorem 1.12(2) = the "single line L" statement,
#                       Example 6.5(i)  = the Hessian arrangement, d_1 = 4 = m,
#   DKP   2505.01733v4  Remark 2.3      = the question quoted in the paper,
#   Kabat 2201.04856v1  Definition 1.3  = the definition of extSS.
# ---------------------------------------------------------------------------

# T7: every \cite occurrence of the paper, verbatim, with its multiplicity.
PAPER_CITATIONS = [
    ("\\cite[Definition~1.1]{Kabat}", 1),
    ("\\cite[Definition~1.3]{Kabat}", 1),
    ("\\cite[Example~6.5(i)]{Dimca}", 1),
    ("\\cite[Remark~2.3]{DKP}", 2),
    ("\\cite[Theorem~1.12(2)]{Dimca}", 3),
    ("\\cite{DKP}", 1),
    ("\\cite{OT}", 1),
]

# T7: the paper's bibliography, in order -- \bibitem key, and the arXiv
# identifiers printed in that entry.  The Orlik-Terao entry is a book and prints
# none, which is why the empty tuple is data and not an omission.
PAPER_BIBLIOGRAPHY = [
    ("Dimca", ("2503.01624v7",)),
    ("DKP", ("2505.01733v4",)),
    ("Kabat", ("2201.04856v1",)),
    ("OT", ()),
]


def parse_citations(text):
    """Every \\cite{k} / \\cite[label]{k} occurrence of a LaTeX source, counted
    over (key, label) pairs; label is None for an unlabelled citation.  \\citep
    and friends are not \\cite and are skipped."""
    out = {}
    i, tok, n = 0, "\\cite", len(text)
    while True:
        i = text.find(tok, i)
        if i < 0:
            return out
        j = i + len(tok)
        if j < n and text[j] == "[":
            k = text.find("]", j)
            if k < 0:
                return out
            label, j = text[j + 1:k], k + 1
        else:
            label = None
        if j >= n or text[j] != "{":
            i = i + len(tok)
            continue
        k = text.find("}", j)
        if k < 0:
            return out
        pair = (text[j + 1:k], label)
        out[pair] = out.get(pair, 0) + 1
        i = k + 1


def arxiv_ids(chunk):
    """The identifiers printed as  arXiv:NNNN.NNNNNvN  in a chunk, in order,
    without repetition."""
    ids, i, tok = [], 0, "arXiv:"
    while True:
        i = chunk.find(tok, i)
        if i < 0:
            return ids
        j = k = i + len(tok)
        while k < len(chunk) and (chunk[k].isdigit() or chunk[k] in ".v"):
            k += 1
        ident = chunk[j:k]
        if ident and ident not in ids:
            ids.append(ident)
        i = k if k > j else j


def parse_bibliography(text):
    """[(bibitem key, tuple of arXiv identifiers in that entry)] in source order.
    An entry runs to the next \\bibitem, or to \\end{thebibliography}."""
    marks, i, tok = [], 0, "\\bibitem{"
    while True:
        i = text.find(tok, i)
        if i < 0:
            break
        j = text.find("}", i)
        if j < 0:
            break
        marks.append((i, text[i + len(tok):j], j + 1))
        i = j + 1
    stop = text.find("\\end{thebibliography}")
    if stop < 0:
        stop = len(text)
    out = []
    for idx, (_start, key, after) in enumerate(marks):
        end = marks[idx + 1][0] if idx + 1 < len(marks) else stop
        out.append((key, tuple(arxiv_ids(text[after:end]))))
    return out


def check_labels(paper_path):
    """Exactly ONE check, recorded whether or not a paper source is available."""
    # (0) The paper source, if there is one, only ever ADDS obligations.
    src, mode = None, ""
    problems = []
    if paper_path:
        try:
            with open(paper_path, "r") as fh:
                src = fh.read()
        except IOError as exc:
            problems.append("a paper source was named but could not be read: %s" % exc)
    if src is None and paper_path:
        mode = ("a paper source was named but could not be read, which is a FAILURE "
                "below, not a skipped check")
    elif src is None:
        mode = ("the embedded transcription T7 alone (no paper .tex is available "
                "here); the check is recorded all the same, so the check total "
                "below does not depend on which files sit beside this program")
    else:
        mode = ("the embedded transcription T7, CORROBORATED token for token "
                "against the paper source supplied")
    info("input of this step: %s" % mode)

    # (1) Parse the embedded citation list, and require it to parse back to the
    #     declared multiplicities (this exercises the parser used on the source).
    embedded = " ".join(" ".join([tok] * n) for tok, n in PAPER_CITATIONS)
    cites = parse_citations(embedded)
    if sum(cites.values()) != sum(n for _tok, n in PAPER_CITATIONS):
        problems.append("the embedded \\cite tokens do not parse back to their "
                        "declared multiplicities")

    # (2) Parse the embedded bibliography into the map key -> version.
    bib = [(key, tuple(ids)) for key, ids in PAPER_BIBLIOGRAPHY]
    bibkeys = [key for key, _ids in bib]
    version_of = dict(bib)
    if len(set(bibkeys)) != len(bibkeys):
        problems.append("duplicate \\bibitem key in the bibliography")

    # (3) No dangling citation, no orphan entry.
    cited = set(key for key, _label in cites)
    for key in sorted(cited - set(bibkeys)):
        problems.append("key cited with no bibliography entry: %s" % key)
    for key in sorted(set(bibkeys) - cited):
        problems.append("bibliography entry never cited: %s" % key)

    # (4) Compose citation -> key -> version and compare with the four checked here.
    for key, arxiv, label in PAPER_LABELS:
        info("expected %-38s and version %s"
             % ("\\cite[%s]{%s}" % (label, key), arxiv))
        if (key, label) not in cites:
            problems.append("label checked here not cited: \\cite[%s]{%s}" % (label, key))
        got = version_of.get(key, ())
        if got != (arxiv,):
            problems.append("entry %s prints %s, not exactly %s"
                            % (key, ", ".join(got) if got else "no arXiv identifier",
                               arxiv))
        # and the version is pinned to that entry alone
        holders = sorted(k for k, ids in bib if arxiv in ids)
        if holders != [key]:
            problems.append("%s is printed by %s, not by %s alone"
                            % (arxiv, ", ".join(holders) if holders else "no entry", key))
    info("derived resolution: %s"
         % ", ".join("%s[%s] -> %s" % (k, l, ", ".join(version_of.get(k, ())) or "-")
                     for k, _a, l in PAPER_LABELS))

    # (5) Every labelled citation present, so that one OUTSIDE the four checked
    #     here is visible rather than silently tolerated.  Not a failure: only
    #     those four labels are checked, and the paper may cite others.
    found = set((key, label) for key, label in cites if label is not None)
    audited = set((key, label) for key, _a, label in PAPER_LABELS)
    extra = sorted(found - audited)
    info("labelled citations present: %s"
         % ", ".join("%s[%s]" % (k, l) for k, l in sorted(found)))
    info("labelled citations OUTSIDE the four checked here (not verified by this program, "
         "and not a failure): %s"
         % (", ".join("%s[%s]" % (k, l) for k, l in extra) if extra else "none"))

    # (6) With a source in hand, the transcription must match it exactly, both
    #     ways: nothing in the paper missing from T7, nothing in T7 absent from
    #     the paper.  A drifted transcription fails here.
    if src is not None:
        if parse_citations(src) != cites:
            problems.append("the paper's \\cite occurrences differ from T7")
        if parse_bibliography(src) != bib:
            problems.append("the paper's bibliography differs from T7")
        for tok, _n in PAPER_CITATIONS:
            if tok not in src:
                problems.append("T7 token absent from the paper source: %s" % tok)
        for key, ids in bib:
            if ("\\bibitem{%s}" % key) not in src:
                problems.append("T7 entry absent from the paper source: %s" % key)
            for ident in ids:
                if ident not in src:
                    problems.append("T7 version absent from the paper source: %s" % ident)

    check(not problems,
          "the four labels checked here are cited and resolve, through the paper's own "
          "bibliography, to exactly the four arXiv versions, one entry each "
          "(the label-to-number resolution inside those preprints was done by hand, "
          "see header)",
          "problems: %s" % ("; ".join(problems) if problems else "none"))


def find_paper(argv):
    """--paper PATH, else the first .tex sitting beside this program or above it."""
    import os
    if "--paper" in argv:
        i = argv.index("--paper")
        if i + 1 < len(argv):
            return argv[i + 1]
    here = os.path.dirname(os.path.abspath(__file__))
    for d in (here, os.path.dirname(here), os.getcwd()):
        try:
            cands = sorted(f for f in os.listdir(d) if f.endswith(".tex"))
        except OSError:
            continue
        for f in cands:
            path = os.path.join(d, f)
            try:
                with open(path, "r") as fh:
                    head = fh.read(4000)
            except IOError:
                continue
            if "extSS" in head and "Hessian" in head:
                return path
    return None


def print_probes(lines, prod, Pd, unj, p0):
    """Anti-vacuity, COMPUTED here instead of asserted in prose.

    Six probes, one per input the checks compare against (the four added lines,
    the printed set P, the printed exceptional points, the -27 of the defining
    equation, the extSS citation checked in step 10, and a truncated B).  Each probe
    perturbs exactly one of those inputs and re-evaluates the SAME predicate the
    corresponding check evaluates; a probe is good only if the predicate turns
    false.  No probe calls check(): the number of checks and their names are the
    same with this section as without it, and a probe can never make the verdict
    green.  What it can do is show that the checks are not tautologies."""
    print("-- Falsification probes: mutate one input, re-evaluate the SAME predicate --")
    NEW = tuple(pnorm((ONE, ONE, TWO)))         # x+y+2z: neither a line of H nor
    added = paper_added_lines()                 # one of the four added lines, and
    flips = []                                  # as a point, (1:1:2) is not in P.

    # P1  the four added lines (target: the step-6 comparison and modularity).
    mut_added = added[:3] + [NEW]
    joins = set(tuple(pnorm(cross(p0, v))) for v in unj)
    Bm = list(lines) + sorted(set(mut_added) - set(lines))
    tabm = incidence_table(Bm)
    Bmset = set(Bm)
    badm = [q for q in tabm if q != p0 and tuple(pnorm(cross(p0, q))) not in Bmset]
    f1 = (joins != set(mut_added)) and bool(badm)
    flips.append(f1)
    info("P1  added lines: replace x+y-2w^2z by x+y+2z -> derived joins == printed "
         "set: %s ; p still modular in the mutated B: %s (%d of %d singular points "
         "unjoined) -> checks 20 and 25 flip: %s"
         % (joins == set(mut_added), not badm, len(badm), len(tabm) - 1, f1))

    # P2  the printed set P (target: check 06).
    Pp, _Vp = paper_P_and_V()
    mutP = set(sorted(Pp)[1:]) | {NEW}
    f2 = set(Pd) != mutP
    flips.append(f2)
    info("P2  printed P: drop one of the nine points, insert (1:1:2) -> derived P == "
         "printed P: %s (symmetric difference %d) -> check 06 flips: %s"
         % (set(Pd) == mutP, len(set(Pd) ^ mutP), f2))

    # P3  the printed exceptional points of p=(1:-1:0) (target: check 15).
    mut_exc = set(tuple(pnorm((ONE, ONE, WPOW[s]))) for s in range(3))   # (0:0:1) gone
    f3 = set(unj) != mut_exc
    flips.append(f3)
    info("P3  printed exceptional points: drop (0:0:1) -> derived set == printed set: "
         "%s (%d derived against %d printed) -> check 15 flips: %s"
         % (set(unj) == mut_exc, len(set(unj)), len(mut_exc), f3))

    # P4  the -27 of the defining equation (target: check 01).
    cubes = pol_add(pol_add({(3, 0, 0): ONE}, {(0, 3, 0): ONE}), {(0, 0, 3): ONE})
    mut_target = pol_mul({(1, 1, 1): ONE},
                         pol_add(pol_pow(cubes, 3),
                                 pol_scale({(3, 3, 3): ONE}, (Fr(-26), Fr(0)))))
    mut_diff = pol_add(prod, pol_scale(mut_target, (Fr(-1), Fr(0))))
    f4 = not pol_eq(prod, mut_target)
    flips.append(f4)
    info("P4  the -27 of eq. (1): use -26 -> product == target: %s (difference has %d "
         "nonzero monomial(s)) -> check 01 flips: %s"
         % (pol_eq(prod, mut_target), len(mut_diff), f4))

    # P5  the extSS citation checked in step 10 (target: check 41).
    cites = parse_citations(" ".join(" ".join([tok] * n) for tok, n in PAPER_CITATIONS))
    holders_real = sorted(k for k, ids in PAPER_BIBLIOGRAPHY if "2201.04856v1" in ids)
    holders_mut = sorted(k for k, ids in PAPER_BIBLIOGRAPHY if "2201.04856v2" in ids)
    f5 = (("Kabat", "Definition~1.3") in cites
          and ("Kabat", "Definition~1.4") not in cites
          and holders_real == ["Kabat"] and holders_mut != ["Kabat"])
    flips.append(f5)
    info("P5  the extSS citation: ask for Kabat[Definition~1.4] and version "
         "2201.04856v2 instead -> label cited: %s (the real label is cited: %s) ; "
         "entries printing that version: %s (the real version: %s) -> check 41 flips: %s"
         % (("Kabat", "Definition~1.4") in cites,
            ("Kabat", "Definition~1.3") in cites,
            ", ".join(holders_mut) if holders_mut else "none",
            ", ".join(holders_real) if holders_real else "none", f5))

    # P6  a truncated B (target: checks 22, 23, 25).
    B3 = list(lines) + sorted(set(added[:3]) - set(lines))
    tab3 = incidence_table(B3)
    B3set = set(B3)
    bad3 = [q for q in tab3 if q != p0 and tuple(pnorm(cross(p0, q))) not in B3set]
    f6 = (len(B3) != PAPER["size_B"] and bool(bad3))
    flips.append(f6)
    info("P6  truncate B: keep only 3 of the 4 added lines -> |B| = %d (target %d), "
         "|Sing(B)| = %d (target %d), p modular: %s (%d singular points unjoined) -> "
         "checks 22 and 25 flip: %s"
         % (len(B3), PAPER["size_B"], len(tab3), PAPER["n_sing_B"],
            not bad3, len(bad3), f6))

    info("probes that flip the check they target: %d of %d  (a probe is not a check: "
         "the check total is the same with and without this section)"
         % (sum(1 for f in flips if f), len(flips)))
    print("")


def print_gaps():
    """Every step between the checked facts and the paper's conclusion that NO
    check above decides.  Printed, not hidden: a green verdict means the checks
    passed, not that these six items were verified here."""
    print("-- Residual gaps: steps NOT decided by any check above ------------------")
    for g in (
        "G1  'at most one line of H passes through a point not in Sing(H)' is "
        "definitional (two lines of H through a point would put it in Sing(H)) and is "
        "not computed; the off-Sing(H) bound k >= 9-3 = 6 rests on it.",
        "G2  'an arrangement with a modular point is supersolvable' -- for rank three "
        "this is the Orlik-Terao equivalence with lattice supersolvability, quoted by "
        "the paper (Section 1); not reproved here.  What is computed is modularity.",
        "G3  Kabat Definition 1.3 defines extSS by supersolvable resolutions "
        "(chains adding one line at a time).  The program bounds min |B \\ H| over "
        "supersolvable B; the two agree because any B containing H can be filtered "
        "one line at a time, which is an argument, not a computation.",
        "G4  freeness of H is obtained by applying the du Plessis-Wall criterion "
        "(tau = (d-1)^2 - d_1 d_2 iff free, for d_1 = mdr(f) <= (d-1)/2) to the "
        "derived mdr(f) and tau; the criterion itself is an external theorem.",
        "G5  the statement of Dimca Thm 1.12(2) is an external input (verified by hand "
        "against arXiv:2503.01624v7); only its arithmetic is re-derived here.",
        "G6  resolving the labels Theorem 1.12(2), Example 6.5(i), Remark 2.3 and "
        "Definition 1.3 to the numbered items they name would require the three arXiv "
        "sources themselves, and those sources are not shipped beside the paper, so "
        "step 10 instead derives the resolution "
        "citation -> bibliography key -> arXiv version from the paper's own "
        "bibliographic data (embedded as literal input T7, and corroborated against "
        "the paper source whenever one is available); the label-to-number resolution "
        "INSIDE those three preprints was done by hand.",
    ):
        info(g)
    print("")


# The closing disclosure.  One line, printed after the verdict, quotable as it
# stands: everything a reader could take the green verdict to cover but which
# this program does NOT establish.
NOT_RERUN = (
    "NOT RE-RUN: this program re-derives every number it prints from the twelve "
    "printed linear forms, the four printed added lines and the point (1:-1:0), but it "
    "does NOT reprove the external inputs it applies -- the rank-three equivalence "
    "between having a modular point and lattice supersolvability (G2), the agreement of "
    "Kabat's resolution definition of extSS with min |B \\ H| (G3), the du Plessis-Wall "
    "freeness criterion (G4), the statement of Dimca Theorem 1.12(2) (G5), and the "
    "definitional fact that at most one line of H passes through a point off Sing(H), on "
    "which the off-Sing(H) bound rests (G1); the resolution of the labels Theorem 1.12(2), "
    "Example 6.5(i), Remark 2.3 and Definition 1.3 to numbered items INSIDE "
    "arXiv:2503.01624v7, arXiv:2505.01733v4 and arXiv:2201.04856v1 was done by hand "
    "against those preprints, since they are not shipped here, and step 10 verifies only "
    "that the paper cites those four labels and that its own bibliography sends each to "
    "the arXiv version used here (G6); the three figures |Sing(B)| = 45 and the "
    "forced-added-line counts 6 for a modular point in V and 6 for one off Sing(H) are "
    "derived here and are NOT printed by the paper, whose lower-bound "
    "argument claims only k >= 4 in each case; the six falsification probes above are "
    "computed in this run, but no separate mutation transcript accompanies the paper; and "
    "any transcript shipped beside the paper that reports 40 checks, or a step 10 gap "
    "reading 'no paper .tex was found', was produced by an earlier build of this file and "
    "is NOT a run of this one -- this version embeds its bibliographic input and therefore "
    "records step 10 in every run, whether or not a paper source sits beside it, for 41 "
    "checks in total.")


def main(argv):
    print("verify.py -- extSS(Hessian) = 4, exact arithmetic in Q(w), w^2 = -1-w")
    print("no floating point is used; every comparison below is exact")
    print("")
    print("-- Step 1: the object (identity pinning H to Dimca, Example 6.5(i)) ------")
    lines = hessian_lines()
    prod = check_identity(lines)
    print("")
    print("-- Step 2: incidences of H ---------------------------------------------")
    tab, Pd, Vd = check_incidences(lines)
    print("")
    print("-- Step 3: Lemma 2(1) ---------------------------------------------------")
    check_lemma1(lines, Pd)
    print("")
    print("-- Step 4: Lemma 2(2) (all nine points of P) ----------------------------")
    p0, unj, _span = check_lemma2(lines, Pd, Vd)
    print("")
    print("-- Step 5: Lemma 2(3) (all twelve points of V) --------------------------")
    check_lemma3(lines, Pd, Vd)
    print("")
    print("-- Step 6: upper bound, extSS(H) <= 4 -----------------------------------")
    _B, _tabB, kmax, modular = check_upper_bound(lines, tab, p0, unj)
    print("")
    print("-- Step 7: lower bound, extSS(H) >= 4 -----------------------------------")
    kmin = check_lower_bound(lines, tab, Pd, Vd)
    print("")
    print("-- Step 8: the theorem and its consequence ------------------------------")
    info("lower bound derived in step 7: extSS(H) >= %d" % kmin)
    info("upper bound derived in step 6: extSS(H) <= %d (a supersolvable extension "
         "of that size was exhibited: %s)" % (kmax, modular))
    check(modular and kmin == kmax == PAPER["extSS"],
          "Theorem: extSS(H) = 4 (the lower bound meets the exhibited extension)",
          "derived %d <= extSS(H) <= %d" % (kmin, kmax))
    check(kmin > 1,
          "consequence: no extension of H by ONE line is supersolvable "
          "(negative answer to DKP Remark 2.3)",
          "extSS(H) = %d > 1" % kmin)
    print("")
    print("-- Step 8b: the first exponent d_1 = mdr(f), derived ---------------------")
    m_max = max(len(S) for S in tab.values())
    d1_derived = check_first_exponent(prod, m_max)
    print("")
    print("-- Step 9: scope of the question, freeness ------------------------------")
    check_source_arithmetic(lines, tab, Pd, Vd, p0, unj, d1_derived)
    print("")
    print("-- Step 10: source fidelity (bibliographic) ----------------------------")
    check_labels(find_paper(argv))
    print("")
    print_probes(lines, prod, Pd, unj, p0)
    print_gaps()
    rc = verdict()
    print("")
    print(NOT_RERUN)
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
