#!/usr/bin/env python3
"""Verification program for the refutation of the deep-hole completeness
conjecture for Reed-Solomon-type codes on elliptic-curve point sets.

Scope: the DECISIVE upper endpoint k = N-4 over the small field named in the
paper.  Everything else the paper claims is named in the closing disclosure.

Standard library only.  Exact integer / Fraction arithmetic throughout.
"""

from fractions import Fraction
from itertools import product

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         ("  -- " + detail) if detail else ""))


def finish():
    bad = [n for n, ok in CHECKS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), len(CHECKS)))
    else:
        print("VERDICT: ALL %d CHECKS PASS" % len(CHECKS))
    return 1 if bad else 0


# ----------------------------------------------------------------------
# exact arithmetic over F_p (p prime), integers only
# ----------------------------------------------------------------------

def inv(a, p):
    a %= p
    if a == 0:
        raise ZeroDivisionError("no inverse of 0")
    return pow(a, p - 2, p)


def rref(rows, p):
    """Row-reduce a list of vectors over F_p.  Returns (reduced, pivots)."""
    mat = [list(r) for r in rows]
    ncol = len(mat[0]) if mat else 0
    pivots = []
    r = 0
    for c in range(ncol):
        pr = None
        for i in range(r, len(mat)):
            if mat[i][c] % p:
                pr = i
                break
        if pr is None:
            continue
        mat[r], mat[pr] = mat[pr], mat[r]
        iv = inv(mat[r][c], p)
        mat[r] = [(v * iv) % p for v in mat[r]]
        for i in range(len(mat)):
            if i != r and mat[i][c] % p:
                f = mat[i][c] % p
                mat[i] = [(mat[i][j] - f * mat[r][j]) % p for j in range(ncol)]
        pivots.append(c)
        r += 1
        if r == len(mat):
            break
    return [tuple(row) for row in mat[:r]], pivots


def canon(v, basis, pivots, p):
    """Canonical representative of the coset v + span(basis), for basis in rref."""
    w = [x % p for x in v]
    for row, c in zip(basis, pivots):
        f = w[c]
        if f:
            for j in range(len(w)):
                w[j] = (w[j] - f * row[j]) % p
    return tuple(w)


def wt(v):
    return sum(1 for x in v if x)


def span(basis, p):
    """All linear combinations of basis vectors over F_p."""
    out = [tuple([0] * len(basis[0]))] if basis else [()]
    for row in basis:
        new = []
        for c in range(p):
            for v in out:
                new.append(tuple((v[j] + c * row[j]) % p for j in range(len(row))))
        out = new
    return sorted(set(out))


# ----------------------------------------------------------------------
# the curve, its points, and the Riemann-Roch spaces L(mO)
# ----------------------------------------------------------------------

def affine_points(a, b, p):
    """Points of y^2 = x^3 + a x + b over F_p, in lexicographic (x, y) order."""
    pts = []
    for x in range(p):
        rhs = (x * x * x + a * x + b) % p
        for y in range(p):
            if (y * y) % p == rhs:
                pts.append((x, y))
    return pts


def singular_projective_points(a, b, p):
    """Singular points of F = Y^2 Z - X^3 - a X Z^2 - b Z^3 in PG(2,p)."""
    bad = []
    for X, Y, Z in product(range(p), repeat=3):
        if (X, Y, Z) == (0, 0, 0):
            continue
        F = (Y * Y * Z - X ** 3 - a * X * Z * Z - b * Z ** 3) % p
        FX = (-3 * X * X - a * Z * Z) % p
        FY = (2 * Y * Z) % p
        FZ = (Y * Y - 2 * a * X * Z - 3 * b * Z * Z) % p
        if F == 0 == FX and FY == 0 == FZ:
            bad.append((X, Y, Z))
    return bad


def lm_monomials(m):
    """Basis monomials x^i y^j (j <= 1) of L(mO): 2i + 3j <= m."""
    mons = []
    for j in (0, 1):
        i = 0
        while 2 * i + 3 * j <= m:
            mons.append((i, j))
            i += 1
    return sorted(mons, key=lambda t: (2 * t[0] + 3 * t[1], t[1]))


def eval_basis(m, pts, p):
    """Evaluation vectors of the L(mO) basis on the ordered point list."""
    return [tuple((pow(x, i, p) * pow(y, j, p)) % p for (x, y) in pts)
            for (i, j) in lm_monomials(m)]


# ======================================================================
# INPUTS TAKEN FROM THE PAPER.  Everything else below is derived.
# ======================================================================
# the field, the curve y^2 = x^3 + A x + B, the ordering of D = E \ {O},
# and the exhibited word.  No count, dimension, radius or coset total is
# taken as input.
Q = 3
A, B = 2, 1
POINT_ORDER = [(0, 1), (0, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
EXHIBITED_WORD = (0, 1, 0, 1, 0, 2)

print("inputs: q=%d, curve y^2 = x^3 + %dx + %d, point order %s, word %s"
      % (Q, A, B, POINT_ORDER, EXHIBITED_WORD))
print("a verification program does accompany this paper: this one. From the field, the curve, "
      "the point order and the word printed in the paper it re-derives the point set, the codes "
      "C_L(D,kO), the covering radius, the deep-hole coset count, the proposed layer and the "
      "k=2 coset census, and compares each derived number with the value the paper prints.")
print("")

# --- check 1: the exhibited data are well formed -----------------------
ok = (len(POINT_ORDER) == 6
      and len(set(POINT_ORDER)) == len(POINT_ORDER)
      and all(0 <= c < Q for pt in POINT_ORDER for c in pt)
      and len(EXHIBITED_WORD) == len(POINT_ORDER)
      and all(0 <= c < Q for c in EXHIBITED_WORD))
ck("exhibited data parse: %d distinct affine points over F_%d, word of that length"
   % (len(POINT_ORDER), Q), ok)

# --- check 2: the curve is nonsingular (derived from A, B) -------------
sing = singular_projective_points(A, B, Q)
ck("curve is nonsingular over F_%d (Jacobian criterion in PG(2,%d))" % (Q, Q),
   not sing, "singular points found: %s" % (sing,) if sing else "none")

# --- check 3: derive E(F_q) and N; the paper's order is exactly D ------
aff = affine_points(A, B, Q)
N = len(aff) + 1                      # + the point at infinity O
n = N - 1
ok3 = (sorted(aff) == sorted(POINT_ORDER))
ck("derived affine point set of E equals the paper's ordered list of D",
   ok3, "derived %s" % (aff,))
ck("derived N = #E(F_%d) = %d satisfies the hypothesis N >= q+4 = %d"
   % (Q, N, Q + 4), N >= Q + 4, "N=%d, n=|D|=%d" % (N, n))

# --- check 4: Hasse bound, tested with integers only ------------------
ck("derived N obeys Hasse: (N-q-1)^2 = %d <= 4q = %d"
   % ((N - Q - 1) ** 2, 4 * Q), (N - Q - 1) ** 2 <= 4 * Q)


def coset_table(basis, pivots, p, length):
    """Brute-force coset census: canonical rep -> minimum weight in the coset."""
    tab = {}
    for v in product(range(p), repeat=length):
        r = canon(v, basis, pivots, p)
        w = wt(v)
        if r not in tab or w < tab[r]:
            tab[r] = w
    return tab


def code_of(m, pts, p):
    """rref basis and pivots of C_L(D, m O) for the given point order."""
    rows, piv = rref(eval_basis(m, pts, p), p)
    return rows, piv


# --- check 5: Riemann-Roch dimensions and injectivity of evaluation ----
dims_ok = True
detail = []
for m in range(1, n + 1):
    if len(lm_monomials(m)) != m:
        dims_ok = False
    detail.append("l(%dO)=%d" % (m, len(lm_monomials(m))))
ck("derived dim L(mO) = m for 1 <= m <= n (genus one, Riemann-Roch)",
   dims_ok, ", ".join(detail))

K = N - 4                              # the decisive upper endpoint
C_rows, C_piv = code_of(K, POINT_ORDER, Q)
Cp_rows, Cp_piv = code_of(K + 1, POINT_ORDER, Q)
ck("evaluation on D is injective at m=%d and m=%d: derived dim C = %d and dim C' = %d"
   % (K, K + 1, len(C_rows), len(Cp_rows)),
   len(C_rows) == K and len(Cp_rows) == K + 1,
   "expected %d and %d from l(mO)" % (K, K + 1))

# --- check 6: the endpoint really is codimension three, with d >= 3 ----
words_C = span(C_rows, Q)
dmin = min(wt(v) for v in words_C if any(v))
ck("upper endpoint k = N-4 = %d is codimension three (n-k = %d) and derived d(C) = %d >= n-k"
   % (K, n - K, dmin), n - K == 3 and dmin >= n - K,
   "|C| = %d words" % len(words_C))

# --- check 7: the Lemma's hypothesis, and the derived covering radius ---
ck("Lemma hypothesis holds on derived data: q+1 = %d < n = %d < q^2+q+1 = %d"
   % (Q + 1, n, Q * Q + Q + 1), Q + 1 < n < Q * Q + Q + 1)

tab = coset_table(C_rows, C_piv, Q, n)
rho = max(tab.values())
ck("derived covering radius of C = C_L(D,%dO) is rho = %d (Lemma predicts 2)"
   % (K, rho), rho == 2, "%d cosets censused" % len(tab))

# --- check 8: the deep-hole coset count, derived vs the Lemma formula --
by_w = {}
for r, w in tab.items():
    by_w[w] = by_w.get(w, 0) + 1
deep = by_w.get(rho, 0)
formula = (Q - 1) * (Q * Q + Q + 1 - n)
ck("derived number of deep-hole cosets = %d equals (q-1)(q^2+q+1-n) = %d"
   % (deep, formula), deep == formula,
   "coset weight profile %s" % (sorted(by_w.items()),))
ck("coset partition is exact: 1 + n(q-1) + deep = %d + %d + %d = q^3 = %d"
   % (by_w.get(0, 0), by_w.get(1, 0), deep, Q ** 3),
   by_w.get(0, 0) == 1 and by_w.get(1, 0) == n * (Q - 1)
   and 1 + n * (Q - 1) + deep == Q ** 3 == len(tab))
ck("deep-hole coset fraction is exactly %s (Fraction arithmetic, no floats)"
   % Fraction(deep, len(tab)),
   Fraction(deep, len(tab)) == Fraction((Q - 1) * (Q * Q + Q + 1 - n), Q ** 3))

# --- check 9: the proposed layer C'/C -- size, and that it IS deep holes -
words_Cp = span(Cp_rows, Q)
layer_reps = set(canon(v, C_rows, C_piv, Q) for v in words_Cp)
zero = tuple([0] * n)
layer_nonzero = sorted(layer_reps - {zero})
ck("derived layer C_L(D,%dO)/C_L(D,%dO) has exactly q-1 = %d nonzero cosets"
   % (K + 1, K, Q - 1), len(layer_nonzero) == Q - 1,
   "%d cosets including the zero coset" % len(layer_reps))
ck("every coset of the derived layer has coset weight exactly rho = %d, so the layer does consist of deep holes"
   % rho, all(tab[r] == rho for r in layer_nonzero),
   "layer coset weights %s" % [tab[r] for r in layer_nonzero])

# --- check 10: DECISIVE -- a derived deep hole outside the layer --------
witness = None
for v in product(range(Q), repeat=n):
    r = canon(v, C_rows, C_piv, Q)
    if tab[r] == rho and r not in layer_reps and wt(v) == rho:
        witness = v
        break
outside = deep - len(layer_nonzero)
if witness is None:
    ck("a deep hole of C outside the proposed layer exists", False,
       "no witness found")
else:
    in_Cp = witness in set(words_Cp)
    dw = min(wt(tuple((witness[i] - c[i]) % Q for i in range(n))) for c in words_C)
    print("derived endpoint witness w* = %s" % (witness,))
    ck("derived witness w* has d(w*, C) = %d = rho, so w* is a deep hole of C_L(D,%dO)"
       % (dw, K), dw == rho)
    ck("derived witness w* lies outside C_L(D,%dO), hence outside the proposed layer"
       % (K + 1), (not in_Cp) and canon(witness, C_rows, C_piv, Q) not in layer_reps)

ck("conjecture fails at k = N-4 = %d: derived deep-hole cosets %d > layer cosets %d, leaving %d outside"
   % (K, deep, len(layer_nonzero), outside),
   deep > len(layer_nonzero) and outside == formula - (Q - 1) and outside > 0)

# ----------------------------------------------------------------------
# The paper's exhibited word, and the finite census it reports, at k = 2.
# ----------------------------------------------------------------------
K2 = 2
C2_rows, C2_piv = code_of(K2, POINT_ORDER, Q)
words_C2 = span(C2_rows, Q)
stated = set()
for u0, u1 in product(range(Q), repeat=2):
    u2 = (-u0 - u1) % Q
    stated.add((u0, u0, u1, u1, u2, u2))
ck("derived C_L(D,2O) equals the paper's stated set {(u0,u0,u1,u1,u2,u2): u0+u1+u2=0}",
   set(words_C2) == stated, "%d derived words, %d stated" % (len(words_C2), len(stated)))

tab2 = coset_table(C2_rows, C2_piv, Q, n)
rho2 = max(tab2.values())
ck("derived covering radius of C_L(D,2O) is rho = %d" % rho2, rho2 == 3,
   "%d cosets censused" % len(tab2))

d_w = min(wt(tuple((EXHIBITED_WORD[i] - c[i]) % Q for i in range(n)))
          for c in words_C2)
ck("the paper's exhibited word has d(w, C_L(D,2O)) = %d = rho, so w is a deep hole at k=2"
   % d_w, d_w == rho2)

pair_diffs = [(EXHIBITED_WORD[0] - EXHIBITED_WORD[1]) % Q,
              (EXHIBITED_WORD[2] - EXHIBITED_WORD[3]) % Q,
              (EXHIBITED_WORD[4] - EXHIBITED_WORD[5]) % Q]
const_in_code = all(len({(c[0] - c[1]) % Q, (c[2] - c[3]) % Q,
                         (c[4] - c[5]) % Q}) == 1 for c in words_C)
ck("every word of C_L(D,3O) has constant pair difference while w has pair differences %s, "
   "so w is outside the k=2 layer" % (pair_diffs,),
   const_in_code and len(set(pair_diffs)) > 1
   and EXHIBITED_WORD not in set(words_C))

layer2 = set(canon(v, C2_rows, C2_piv, Q) for v in words_C)
layer2_nonzero = sorted(layer2 - {zero})
census = sum(1 for w in tab2.values() if w == rho2)
ck("census recomputed: %d of the %d cosets of C_L(D,2O) have coset weight %d, "
   "and the layer accounts for %d of them, leaving %d outside"
   % (census, len(tab2), rho2, len(layer2_nonzero), census - len(layer2_nonzero)),
   census == 24 and len(tab2) == Q ** (n - K2) and len(layer2_nonzero) == Q - 1
   and census - len(layer2_nonzero) == 22
   and all(tab2[r] == rho2 for r in layer2_nonzero))

print("note: the 81-coset census just recomputed is the one finite computation the paper quotes, "
      "and the source of its numbers 24 and 22. It needs no program of its own: the paper's "
      "decisive steps -- rho(C_2) = 3, d(w,C_2) = 3, w outside C_L(D,3O), and the endpoint count "
      "of the covering-radius lemma -- are hand-sized and are carried out in the paper itself, "
      "and the paper states that nothing in it depends on the census, which is reported only to "
      "quantify how far the layer falls short at k = 2. It is recomputed here regardless, and "
      "the recomputation agrees with the printed 24 and 22. So an absent program would not have "
      "meant an absent computation here; it would have meant a computation the reader can redo "
      "by hand from the printed data.")

print("note: at the endpoint k = %d the paper argues by count, not by a word; "
      "the derived d(w, C_L(D,%dO)) = %d shows the k=2 word is not itself an endpoint "
      "witness, which is why w* above was derived instead."
      % (K, K,
         min(wt(tuple((EXHIBITED_WORD[i] - c[i]) % Q for i in range(n)))
             for c in words_C)))

# ----------------------------------------------------------------------
# A bounded sweep: the same endpoint conclusion on other small odd PRIME
# fields.  This does not prove the general theorem; it tests it on a
# finite family, and each instance can fail on its own.
# ----------------------------------------------------------------------

def endpoint_instance(a, b, p):
    """Derive (n, rho<=2?, deep-hole count, layer count) at k = N-4 for one curve."""
    pts = affine_points(a, b, p)
    nn = len(pts)                       # n = N - 1
    kk = nn - 3
    rows, piv = rref(eval_basis(kk, pts, p), p)
    if len(rows) != kk:
        return None
    zed = tuple([0] * nn)
    s1 = {canon(zed, rows, piv, p)}
    for i in range(nn):
        for c in range(1, p):
            v = list(zed)
            v[i] = c
            s1.add(canon(v, rows, piv, p))
    s2 = set(s1)
    for i in range(nn):
        for j in range(i + 1, nn):
            for c1 in range(1, p):
                for c2 in range(1, p):
                    v = list(zed)
                    v[i] = c1
                    v[j] = c2
                    s2.add(canon(v, rows, piv, p))
    new = [m for m in lm_monomials(kk + 1) if m not in lm_monomials(kk)]
    if len(new) != 1:
        return None
    i0, j0 = new[0]
    g = tuple((pow(x, i0, p) * pow(y, j0, p)) % p for (x, y) in pts)
    layer = set()
    for c in range(1, p):
        layer.add(canon([(c * t) % p for t in g], rows, piv, p))
    ok_inst = (len(s1) == 1 + nn * (p - 1)          # d(C) >= 3
               and len(s2) == p ** 3                # rho <= 2
               and len(layer) == p - 1              # layer is q-1 cosets
               and all(r not in s1 for r in layer)  # layer cosets are deep
               and p ** 3 - len(s1) == (p - 1) * (p * p + p + 1 - nn))
    return (nn, ok_inst, p ** 3 - len(s1), p - 1)


sweep_total = 0
sweep_bad = []
sweep_primes = set()
for p in (3, 5, 7, 11):
    budget = 2 if p >= 11 else 10 ** 6
    used = 0
    for a, b in product(range(p), repeat=2):
        if used >= budget:
            break
        if singular_projective_points(a, b, p):
            continue
        nn = len(affine_points(a, b, p))
        if nn + 1 < p + 4:
            continue
        used += 1
        res = endpoint_instance(a, b, p)
        sweep_total += 1
        sweep_primes.add(p)
        if res is None or not res[1] or not res[2] > res[3]:
            sweep_bad.append((p, a, b, res))
ck("bounded sweep over small odd PRIME fields %s: %d curve instances with N >= q+4 all give "
   "rho = 2, deep-hole count (q-1)(q^2+q+1-n), and a layer of only q-1 cosets that is a "
   "proper subset of the deep holes"
   % (sorted(sweep_primes), sweep_total),
   sweep_total >= 4 and len(sweep_primes) >= 3 and not sweep_bad,
   "failures: %s" % (sweep_bad[:3],) if sweep_bad else "no failures")

rc = finish()
print("NOT RE-RUN HERE: the endpoint theorem in the generality claimed -- it is asserted "
      "for every odd prime power q, every elliptic curve E/F_q with #E(F_q) >= q+4 and every "
      "choice of O, whereas this program derives it in full only for the single curve "
      "y^2 = x^3 + 2x + 1 over F_3 with O the point at infinity, plus a bounded sweep of "
      "curves over the PRIME fields F_3, F_5, F_7, F_11 (capped at two curves for q=11) with "
      "O at infinity: no non-prime q (F_9, F_25, ...) is touched at all, which is exactly the "
      "case the paper says the source's side conditions miss, and no origin O other than the "
      "point at infinity is used; that the projective model has exactly one point at infinity, "
      "namely O = [0:1:0], is assumed rather than derived -- N is formed as (number of affine "
      "solutions) + 1 -- and the paper's own particular singularity computations (F_X = Z^2 over "
      "F_3, and F_Z(O) = 1) are not reproduced step by step, being replaced by an exhaustive "
      "Jacobian scan of PG(2,q); the codimension-three lemma as a general statement about "
      "arbitrary [n,n-3] codes over F_q with d >= 3, and its PG(2,q) line-counting proof, are "
      "only instantiated here, never proved; the general theorems the paper invokes -- "
      "Riemann-Roch giving l(mO) = m in genus one, injectivity of evaluation when deg(kO) < n, "
      "the designed-distance bound d >= n-k, Hasse's bound, and the cited textbook -- are used "
      "as the monomial-basis construction and confirmed only by rank and weight computations on "
      "the finitely many codes built here; every statement about the source paper is unverified "
      "because no external document was fetched or parsed, including that its Conjecture 1.4, "
      "Theorem 1.2(ii)-(iv), Remark 1.3 and Corollary 4.11 say what is quoted, that their count "
      "(#E(F_q)-n)(q-1)q^k specialises to the layer, that their side conditions are n >= q+k or "
      "q prime or k <= sqrt(q) and that at k = n-3 the first forces q <= 3 while the third "
      "fails, that Corollary 4.11 concerns residue codes and does not conflict, and the remark "
      "about preprint versus printed numbering; the claim that at k = n-2 = N-3 the covering "
      "radius drops to 1 so completeness fails trivially is not computed; at k = 2 the paper's "
      "hand arguments -- that every word of F_3^6 agrees with a codeword in at least three "
      "coordinates, and the step A_0 + A_1 = F_3 -- are not re-derived, only their conclusions "
      "rho = 3, d(w,C) = 3, w outside C_L(D,3O) and the 24/22 census are recomputed, and only "
      "for this one F_3 curve; nothing is checked about the interior of the range, where the "
      "paper itself decides nothing for q > 3 and any k < N-4; on every curve of the bounded "
      "sweep only the single value k = N-4 is examined, so no k < N-4 -- in particular not k = 2 "
      "-- is examined for any curve other than y^2 = x^3 + 2x + 1 over F_3, and the sweep "
      "enumerates only models in short Weierstrass form y^2 = x^3 + a x + b, one instance per "
      "(a,b) pair with no reduction to isomorphism classes, so it covers neither every elliptic "
      "curve over those four fields nor every model of the curves it does reach; the prior-art "
      "and novelty "
      "statements -- no resolution of the conjecture in the literature, no earlier appearance of "
      "the exhibited word, no novelty claimed for the lemma -- involve no computation and no "
      "literature search, of the citing literature or otherwise, was performed; and every "
      "provenance and record-keeping assertion the paper makes about material outside itself is "
      "uncorroborated here: that a first implementation of the census exists, that a second one "
      "was written independently from the definitions of the note, that its 81 coset-weight "
      "determinations and the counts 24 and 22 agreed throughout with the printed values, and "
      "that both programs and the transcript of that run are retained in an auxiliary archive "
      "available on request, are statements this program cannot test -- no such archive, program "
      "or transcript was read by it, none is distributed beside it, and neither earlier program "
      "is re-run here; this program is a further, independent implementation written from the "
      "printed data alone, and its agreement with the printed numbers is the only evidence for "
      "the census offered in this bundle.")
raise SystemExit(rc)
