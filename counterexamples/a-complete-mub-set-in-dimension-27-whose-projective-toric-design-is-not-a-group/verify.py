#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify.py -- re-derives every quantity claimed in paper.tex / paper.pdf,
   "A complete set of mutually unbiased bases in dimension 27 whose projective
    toric design is not a group".

WHAT IT TAKES AS INPUT.  Only what the paper prints:

  * the presentation  F = GF(27) = F_3[t]/(t^3 - t - 1),  element index
    e = a_0 + 3 a_1 + 9 a_2  <->  a_0 + a_1 t + a_2 t^2   (paper, Section 2);
  * the 27 x 27 digit table of Table 1, transcribed verbatim below into
    PRINTED_FORMS (row = one quadratic form Q, column w = the value Q(w) in
    F_3, coordinates ordered by the index e above);
  * the three generators of Lambda printed after Table 1, transcribed into
    PRINTED_LAMBDA_GENS.

Everything else -- the field arithmetic, the trace, the spread maps, the design
X, the character sums, the closure sweeps -- is rebuilt here from that data and
from the formulas quoted in the paper, and compared against it.

CONTRACT.  Python 3.9+, standard library only.  All arithmetic is exact integer
arithmetic in Z, F_3 = Z/3 or the ring Z[zeta_8]; no floating point value ever
decides a check.  Character and angle sums over cube roots of unity are decided
by the integer identity

    n_0 + n_1 w + n_2 w^2 = 0   in C, w = exp(2 pi i / 3)
    <=>   n_0 = n_1 = n_2                       (n_i integers),
    |n_0 + n_1 w + n_2 w^2|^2 = n_0^2 + n_1^2 + n_2^2 - n_0 n_1 - n_1 n_2 - n_2 n_0.

One line `PASS <name> [detail]` is printed per check, then the closing verdict.
Exit status is 0 if and only if every check passed.  Runs in a few seconds.
"""

import itertools
import sys

# ---------------------------------------------------------------------------
# 0.  THE OBJECT, TRANSCRIBED FROM THE PAPER
# ---------------------------------------------------------------------------
# Table 1 of the paper: the 27 quadratic forms Q of the Ball-Bamberg-Lavrauw-
# Penttila symplectic spread over GF(27) (the Hering plane of order 27), each as
# the digit string (Q(w))_{w} in F_3, w running over GF(27) in the index order
# e = a_0 + 3 a_1 + 9 a_2.  Labels are the paper's: Q0; Qm<m> for the 13
# nonsquares m; Qs<s> for the 13 chosen representatives s of F* / {+-1}.
PRINTED_FORMS = """
Q0    000000000000000000000000000
Qm2   000111111210021021201012012
Qm6   000210201201111102210120111
Qm7   000102120021120111012111102
Qm8   000021012111102120111102120
Qm10  022112121112010100121100001
Qm14  022121112001211010010001211
Qm17  022211211100100211100211100
Qm18  011110101011002212011221020
Qm19  011002020101011221110212011
Qm21  011200200110221110101101212
Qm23  011011011020212101002110221
Qm24  011020002212110011221011101
Qm26  011101110122101002122020110
Qs1   000222222012201201021210210
Qs3   022022022220001112202121010
Qs4   022202220010112022001022121
Qs5   022220202121022001112010022
Qs9   000201210222012210222201021
Qs10  011212221221200020212002200
Qs11  022001010211220202211220202
Qs12  022010001022121220022202112
Qs13  000012021102222012120021222
Qs14  011122122200122200200200122
Qs15  022100100202202121220112220
Qs16  000120102120210222102222201
Qs17  011221212002020122020122002
"""

# The three printed generators of Lambda, the group of trace functionals:
# lambda_v(w) = tr(v w) for v = 1, t, t^2, i.e. for element indices 1, 3, 9.
PRINTED_LAMBDA_GENS = """
L1    000000000222222222111111111
Lt    000222111000222111000222111
Lt2   021021021210210210102102102
"""

# The non-closure witness printed in the paper (Section 4).
PRINTED_WITNESS_PHI = "000111111210021021201012012"          # the row Qm2
PRINTED_WITNESS_SUM = "000222222120012012102021021"          # phi + phi mod 3

# The d = 2 remark of Section 5: the four points of P(T^2), second angle in
# units of pi/4.  X_2 = {pi/4, 3pi/4, 5pi/4, 7pi/4}.
PRINTED_D2_ANGLES = (1, 3, 5, 7)

D = 27

# ---------------------------------------------------------------------------
# 1.  CHECK PLUMBING
# ---------------------------------------------------------------------------
_n_pass = 0
_n_fail = 0


def check(name, ok, detail=""):
    global _n_pass, _n_fail
    if ok:
        _n_pass += 1
        print("PASS %s%s" % (name, (" " + detail) if detail else ""))
    else:
        _n_fail += 1
        print("FAIL %s%s" % (name, (" " + detail) if detail else ""))


def popcount(n):
    try:
        return n.bit_count()
    except AttributeError:
        return bin(n).count("1")


def parse_table(block):
    labels, rows = [], []
    for line in block.strip().splitlines():
        parts = line.split()
        if len(parts) != 2:
            raise SystemExit("malformed transcribed row: %r" % line)
        labels.append(parts[0])
        rows.append(tuple(int(c) for c in parts[1]))
    return labels, rows


# ---------------------------------------------------------------------------
# 2.  GF(27) = F_3[t]/(t^3 - t - 1), REBUILT FROM THE PRINTED PRESENTATION
# ---------------------------------------------------------------------------
def digits(e):
    return (e % 3, (e // 3) % 3, (e // 9) % 3)


def index(a):
    return (a[0] % 3) + 3 * (a[1] % 3) + 9 * (a[2] % 3)


ADD = [[index(tuple(x + y for x, y in zip(digits(a), digits(b)))) for b in range(27)]
       for a in range(27)]
NEG = [index(tuple(-x for x in digits(a))) for a in range(27)]


def _raw_mul(a, b):
    A, B = digits(a), digits(b)
    c = [0] * 5
    for i in range(3):
        for j in range(3):
            c[i + j] = (c[i + j] + A[i] * B[j]) % 3
    # reduce with t^3 = t + 1, hence t^4 = t^2 + t
    c0, c1, c2, c3, c4 = c
    return index(((c0 + c3) % 3, (c1 + c3 + c4) % 3, (c2 + c4) % 3))


MUL = [[_raw_mul(a, b) for b in range(27)] for a in range(27)]


def add(a, b):
    return ADD[a][b]


def mul(a, b):
    return MUL[a][b]


def power(a, k):
    r = 1
    for _ in range(k):
        r = mul(r, a)
    return r


FROB = [power(x, 3) for x in range(27)]
NONZERO = [x for x in range(27) if x != 0]

# --- check 1: the presentation really is a field ---------------------------
roots = [x for x in (0, 1, 2) if (x ** 3 - x - 1) % 3 == 0]
check("field-irreducible", roots == [],
      "[t^3-t-1 sends 0,1,2 in F_3 to %s, so it has no root and, being cubic, is irreducible]"
      % ",".join(str((x ** 3 - x - 1) % 3) for x in (0, 1, 2)))

_units = all(any(mul(a, b) == 1 for b in NONZERO) for a in NONZERO)
_comm = all(mul(a, b) == mul(b, a) for a in range(27) for b in range(27))
_assoc = all(mul(mul(a, b), c) == mul(a, mul(b, c))
             for a in range(27) for b in range(27) for c in range(27))
_distr = all(mul(a, add(b, c)) == add(mul(a, b), mul(a, c))
             for a in range(27) for b in range(27) for c in range(27))
_cyclic = False
for g in NONZERO:
    seen, x = set(), 1
    for _ in range(26):
        x = mul(x, g)
        seen.add(x)
    if len(seen) == 26:
        _cyclic = True
        break
check("field-axioms", _units and _comm and _assoc and _distr and _cyclic,
      "[27 elements; multiplication commutative, associative, distributive over the "
      "componentwise mod-3 addition; all 26 nonzero elements invertible; F* cyclic of order 26]")

# --- check 2: Frobenius ----------------------------------------------------
_fr_add = all(FROB[add(a, b)] == add(FROB[a], FROB[b]) for a in range(27) for b in range(27))
_fr_mul = all(FROB[mul(a, b)] == mul(FROB[a], FROB[b]) for a in range(27) for b in range(27))
_fr_ord = all(FROB[FROB[FROB[x]]] == x for x in range(27)) and any(FROB[x] != x for x in range(27))
_fr_fix = sorted(x for x in range(27) if FROB[x] == x) == [0, 1, 2]
check("frobenius", _fr_add and _fr_mul and _fr_ord and _fr_fix,
      "[x -> x^3 is an additive, multiplicative automorphism of order 3 fixing exactly F_3 = {0,1,2}, "
      "hence F_3-linear]")

# --- check 3: the trace ----------------------------------------------------
TR = [add(add(x, FROB[x]), FROB[FROB[x]]) for x in range(27)]
_tr_vals = all(TR[x] in (0, 1, 2) for x in range(27))
_tr_lin = all(TR[add(a, b)] == (TR[a] + TR[b]) % 3 for a in range(27) for b in range(27))
_tr_bal = [TR.count(0), TR.count(1), TR.count(2)] == [9, 9, 9]
check("trace", _tr_vals and _tr_lin and _tr_bal,
      "[tr(x)=x+x^3+x^9 lands in F_3, is F_3-linear, and is balanced: fibres of size 9,9,9]")

FUNCTIONALS = [tuple(TR[mul(v, w)] for w in range(27)) for v in range(27)]
check("trace-nondegenerate", len(set(FUNCTIONALS)) == 27 and FUNCTIONALS[0] == tuple([0] * 27),
      "[v -> tr(v .) is injective, so the 27 maps w -> tr(vw) are exactly all 27 = |Hom_{F_3}(F,F_3)| "
      "F_3-functionals on F]")

# --- check 4: FACT 1, -1 is a nonsquare ------------------------------------
SQUARES = sorted({mul(x, x) for x in NONZERO})
NONSQUARES = sorted(set(NONZERO) - set(SQUARES))
_sq_group = all(mul(a, b) in SQUARES for a in SQUARES for b in SQUARES)
_order_of_minus_one = next(k for k in range(1, 27) if power(NEG[1], k) == 1)
check("fact1-minus-one-nonsquare",
      len(SQUARES) == 13 and len(NONSQUARES) == 13 and _sq_group
      and _order_of_minus_one == 2 and NEG[1] in NONSQUARES,
      "[13 squares forming a subgroup of the ODD order 13; -1 = index %d has multiplicative order 2, "
      "and 2 does not divide 13, so -1 is a nonsquare]" % NEG[1])

# --- check 5: beta inverts alpha ------------------------------------------
# beta(v) = (1/2)(-v + v^3 + v^9) = 2(-v + v^3 + v^9) since 1/2 = 2 in F_3;
# alpha(u) = u^3 + u^9.
BETA = [mul(2, add(add(NEG[v], FROB[v]), FROB[FROB[v]])) for v in range(27)]
_beta_bad = [v for v in range(27) if add(FROB[BETA[v]], FROB[FROB[BETA[v]]]) != v]
_beta_lin = all(BETA[add(a, b)] == add(BETA[a], BETA[b]) for a in range(27) for b in range(27))
check("beta-inverts-alpha", _beta_bad == [] and _beta_lin,
      "[beta(v)=2(-v+v^3+v^9) is F_3-linear and alpha(beta(v)) = beta(v)^3 + beta(v)^9 = v for all 27 v, "
      "0 failures]")

# ---------------------------------------------------------------------------
# 3.  THE PRINTED TABLE IS WHAT THE PAPER'S FORMULAS SAY
# ---------------------------------------------------------------------------
LABELS, FORMS = parse_table(PRINTED_FORMS)
GLABELS, LGENS = parse_table(PRINTED_LAMBDA_GENS)

check("table-shape",
      len(FORMS) == 27 and all(len(r) == 27 for r in FORMS)
      and all(c in (0, 1, 2) for r in FORMS for c in r) and len(set(FORMS)) == 27,
      "[Table 1 transcribes as 27 rows of 27 digits in {0,1,2}, all 27 rows distinct]")

check("table-zero-row", FORMS[0] == tuple([0] * 27) and LABELS[0] == "Q0",
      "[the row labelled Q0 is identically zero]")

_m_labels = [int(l[2:]) for l in LABELS if l.startswith("Qm")]
_m_ok = sorted(_m_labels) == NONSQUARES
for m in _m_labels:
    if FORMS[LABELS.index("Qm%d" % m)] != tuple(TR[mul(power(m, 3), power(w, 4))] for w in range(27)):
        _m_ok = False
check("table-nonsquare-rows", _m_ok,
      "[the 13 labels Qm* are exactly the 13 nonsquares %s, and each printed row equals "
      "w -> tr(m^3 w^4), 13 of 13]" % ",".join(map(str, NONSQUARES)))

_s_labels = [int(l[2:]) for l in LABELS if l.startswith("Qs")]
_reps, _seen = [], set()
for x in NONZERO:
    if x in _seen:
        continue
    _reps.append(x)
    _seen.update((x, NEG[x]))
_s_ok = len(_s_labels) == 13 and len({frozenset((s, NEG[s])) for s in _s_labels}) == 13
for s in _s_labels:
    if FORMS[LABELS.index("Qs%d" % s)] != tuple(
            TR[mul(2, mul(mul(w, s), BETA[mul(w, s)]))] for w in range(27)):
        _s_ok = False
check("table-class-rows", _s_ok,
      "[the 13 labels Qs* are 13 pairwise distinct classes of F*/{+-1}, and each printed row equals "
      "w -> tr(2 w s beta(w s)), 13 of 13]")

_lam_ok = (len(LGENS) == 3
           and LGENS[0] == FUNCTIONALS[1] and LGENS[1] == FUNCTIONALS[3] and LGENS[2] == FUNCTIONALS[9])
LAMBDA = sorted({tuple((c0 * LGENS[0][w] + c1 * LGENS[1][w] + c2 * LGENS[2][w]) % 3 for w in range(27))
                 for c0 in range(3) for c1 in range(3) for c2 in range(3)})
check("lambda-generators", _lam_ok and len(LAMBDA) == 27 and set(LAMBDA) == set(FUNCTIONALS),
      "[the 3 printed generators are w -> tr(vw) for v = 1, t, t^2; their F_3-span has 27 elements and "
      "equals the full set of trace functionals]")

# The 27 F_3-linear spread maps h, and Q_h(w) = tr((1/2) w h(w)) = tr(2 w h(w)).
HMAPS, HLABELS = [], []
HMAPS.append(tuple(0 for _ in range(27)))
HLABELS.append("h0")
for m in NONSQUARES:
    HMAPS.append(tuple(add(mul(m, FROB[FROB[x]]), mul(power(m, 3), FROB[x])) for x in range(27)))
    HLABELS.append("hm%d" % m)
for s in _s_labels:
    HMAPS.append(tuple(mul(s, BETA[mul(x, s)]) for x in range(27)))
    HLABELS.append("hs%d" % s)

_order = ["Q0"] + ["Qm%d" % m for m in NONSQUARES] + ["Qs%d" % s for s in _s_labels]
_perm = [LABELS.index(nm) for nm in _order]
_qh_ok = (len(HMAPS) == 27 and len(set(HMAPS)) == 27
          and all(tuple(TR[mul(2, mul(w, HMAPS[i][w]))] for w in range(27)) == FORMS[_perm[i]]
                  for i in range(27)))
check("forms-from-spread-maps", _qh_ok,
      "[the 27 maps h_0 = 0, h_m(x) = m x^9 + m^3 x^3, h_s(x) = s beta(xs) are distinct, and "
      "Q_h(w) = tr(2 w h(w)) reproduces the printed row of Table 1 for every one of them]")

# ---------------------------------------------------------------------------
# 4.  EACH ROW IS A QUADRATIC FORM; THE 27 OF THEM ARE A SYMPLECTIC SPREAD
# ---------------------------------------------------------------------------
check("rows-are-even", all(Q[0] == 0 and all(Q[NEG[w]] == Q[w] for w in range(27)) for Q in FORMS),
      "[every printed row satisfies Q(0) = 0 and Q(-w) = Q(w), 27 of 27]")


def polar(Q, w, wp):
    return (Q[add(w, wp)] - Q[w] - Q[wp]) % 3


_bil = True
for Q in FORMS:
    for w in range(27):
        for wp in range(27):
            if polar(Q, w, wp) != polar(Q, wp, w):
                _bil = False
    for w in range(27):
        for a in range(27):
            for b in range(27):
                if polar(Q, w, add(a, b)) != (polar(Q, w, a) + polar(Q, w, b)) % 3:
                    _bil = False
                    break
            if not _bil:
                break
        if not _bil:
            break
check("polar-form-bilinear", _bil,
      "[B_Q(w,w') = Q(w+w') - Q(w) - Q(w') is symmetric and additive in each argument for all 27 rows, "
      "so each row is a genuine quadratic form on the F_3-space F]")

_pol_tr = all(polar(FORMS[_perm[i]], w, wp) == TR[mul(w, HMAPS[i][wp])]
              for i in range(27) for w in range(27) for wp in range(27))
check("polar-form-is-trace-pairing", _pol_tr,
      "[B_{Q_h}(w,w') = tr(w h(w')) for all 27 maps h and all 729 pairs (w,w'); this step uses "
      "char F = 3 odd, via 1/2 in F_3]")

BASIS = (1, 3, 9)                     # 1, t, t^2
MATS = [tuple(tuple(polar(Q, BASIS[i], BASIS[j]) for j in range(3)) for i in range(3)) for Q in FORMS]
check("spread-matrices-distinct", len(set(MATS)) == 27
      and all(M[i][j] == M[j][i] for M in MATS for i in range(3) for j in range(3)),
      "[the 27 polar forms give 27 distinct symmetric 3x3 matrices over F_3]")


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) % 3


_sing, _pairs = 0, 0
for i in range(27):
    for j in range(i + 1, 27):
        Dm = tuple(tuple((MATS[i][a][b] - MATS[j][a][b]) % 3 for b in range(3)) for a in range(3))
        _pairs += 1
        if det3(Dm) == 0:
            _sing += 1
check("spread-differences-nonsingular", _pairs == 351 and _sing == 0,
      "[all %d = C(27,2) pairwise differences of the 27 symmetric matrices are nonsingular over F_3, "
      "0 singular; with the component at infinity this is a symplectic spread of 28 components of "
      "GF(27) x GF(27)]" % _pairs)

# ---------------------------------------------------------------------------
# 5.  THE DESIGN X, AND THE TWO HYPOTHESES OF THE CONJECTURE
# ---------------------------------------------------------------------------
XLIST = sorted({tuple((Q[w] + L[w]) % 3 for w in range(27)) for Q in FORMS for L in LAMBDA})
XSET = set(XLIST)
check("design-size", len(XLIST) == 729 == D * D,
      "[X = {Q + lambda : Q in Table 1, lambda in Lambda} has 729 = 27^2 distinct points, so "
      "(Q,lambda) -> Q + lambda is a bijection onto X]")

check("design-normalised", all(x[0] == 0 for x in XLIST),
      "[the w = 0 coordinate of every one of the 729 points is 0, the normalisation of P(T^d) = T^d/T]")

ZERO = tuple([0] * 27)
check("design-contains-identity", ZERO in XSET,
      "[the identity of P(T^27), the all-zero point, lies in X (it is Q0 + 0)]")

# --- the toric angle condition -------------------------------------------
MASK = []
for x in XLIST:
    m = [0, 0, 0]
    for w, v in enumerate(x):
        m[v] |= 1 << w
    MASK.append(tuple(m))

_ang_viol, _ang_pairs, _ang_vals, _ang_shapes = 0, 0, set(), {}
for i in range(729):
    a0, a1, a2 = MASK[i]
    for j in range(i + 1, 729):
        b0, b1, b2 = MASK[j]
        n0 = popcount(a0 & b0) + popcount(a1 & b1) + popcount(a2 & b2)
        n1 = popcount(b0 & a1) + popcount(b1 & a2) + popcount(b2 & a0)
        n2 = 27 - n0 - n1
        val = n0 * n0 + n1 * n1 + n2 * n2 - n0 * n1 - n1 * n2 - n2 * n0
        _ang_vals.add(val)
        key = (tuple(sorted((n0, n1, n2))), val)
        _ang_shapes[key] = _ang_shapes.get(key, 0) + 1
        _ang_pairs += 1
        if val not in (0, 27):
            _ang_viol += 1
_want_shapes = {((9, 9, 9), 0): 9477, ((6, 9, 12), 27): 255879}
check("toric-angle",
      _ang_pairs == 265356 and _ang_viol == 0 and sorted(_ang_vals) == [0, 27]
      and _ang_shapes == _want_shapes,
      "[over all %d = C(729,2) unordered pairs phi != theta of X, |sum_w exp(i(phi_w-theta_w))|^2 lies in "
      "{0,27} = {0,d}: 0 violations. The residue counts (n_0,n_1,n_2) of the difference take exactly two "
      "shapes: (9,9,9) on 9477 pairs, giving 0, and permutations of (6,9,12) on 255879 pairs, giving 27]"
      % _ang_pairs)

# --- the 2-design condition ----------------------------------------------
# k = e_{a1} + e_{a2} - e_{b1} - e_{b2} with {a1,a2}, {b1,b2} multisets from the
# 27 coordinates.  There are 378 such multisets, hence 378^2 = 142884 ordered
# pairs; k = 0 exactly when the two multisets coincide.
MULTISETS = [(a, b) for a in range(27) for b in range(a, 27)]
check("characters-enumerated", len(MULTISETS) == 378 and len(MULTISETS) ** 2 == 142884,
      "[the 27 coordinates give %d size-2 multisets and %d ordered pairs of them, i.e. %d characters "
      "chi_k of Definition 2.6]" % (len(MULTISETS), len(MULTISETS) ** 2, len(MULTISETS) ** 2))

_kzero_int, _kzero_mod3 = 0, 0
for (a1, a2) in MULTISETS:
    for (b1, b2) in MULTISETS:
        kv = [0] * 27
        kv[a1] += 1
        kv[a2] += 1
        kv[b1] -= 1
        kv[b2] -= 1
        if all(c == 0 for c in kv):
            _kzero_int += 1
        if all(c % 3 == 0 for c in kv):
            _kzero_mod3 += 1
check("characters-nonzero-mod-3", _kzero_int == 378 and _kzero_mod3 == 378,
      "[exactly %d of the %d characters have k = 0 as an integer vector, and exactly the same %d have "
      "k = 0 mod 3; so 'k nonzero' is unambiguous and %d characters must be tested]"
      % (_kzero_int, len(MULTISETS) ** 2, _kzero_mod3, len(MULTISETS) ** 2 - _kzero_int))

UM = []
for (a, b) in MULTISETS:
    m = [0, 0, 0]
    for idx, x in enumerate(XLIST):
        m[(x[a] + x[b]) % 3] |= 1 << idx
    UM.append(tuple(m))

_des_bad, _des_tested = 0, 0
for i in range(378):
    u0, u1, u2 = UM[i]
    for j in range(378):
        if i == j:
            continue
        v0, v1, v2 = UM[j]
        n0 = popcount(u0 & v0) + popcount(u1 & v1) + popcount(u2 & v2)
        n1 = popcount(u1 & v0) + popcount(u2 & v1) + popcount(u0 & v2)
        _des_tested += 1
        if not (n0 == 243 and n1 == 243):
            _des_bad += 1
check("two-design", _des_tested == 142506 and _des_bad == 0,
      "[all %d nonzero characters chi_k sum to zero over X: every one has residue counts "
      "(n_0,n_1,n_2) = (243,243,243), so sum_{x in X} exp(i k.x) = 243(1 + w + w^2) = 0; 0 failures]"
      % _des_tested)

# ---------------------------------------------------------------------------
# 6.  X IS NEITHER A SUBGROUP NOR A COSET OF ONE
# ---------------------------------------------------------------------------
_phi = tuple(int(c) for c in PRINTED_WITNESS_PHI)
_sum = tuple(int(c) for c in PRINTED_WITNESS_SUM)
check("witness-transcription",
      _phi == FORMS[LABELS.index("Qm2")] and _sum == tuple((2 * c) % 3 for c in _phi)
      and _phi in XSET and _sum not in XSET,
      "[the printed witness phi is the row Qm2, it lies in X, the printed phi+phi is its digitwise "
      "double mod 3, and that point is not in X]")

_out = 0
for x in XLIST:
    for y in XLIST:
        if tuple((x[w] + y[w]) % 3 for w in range(27)) not in XSET:
            _out += 1
check("not-a-subgroup", _out == 435942,
      "[%d of the %d ordered sums x + y with x,y in X leave X; X is closed under neither, so it is not "
      "a subgroup of P(T^27)]" % (_out, 729 * 729))

_coset_closed, _worst = 0, 0
for x0 in XLIST:
    found, seen = False, 0
    for x in XLIST:
        for y in XLIST:
            seen += 1
            if tuple((x[w] + y[w] - x0[w]) % 3 for w in range(27)) not in XSET:
                found = True
                break
        if found:
            break
    _worst = max(_worst, seen)
    if not found:
        _coset_closed += 1
check("not-a-coset", _coset_closed == 0,
      "[for every one of the 729 base points x0 in X, the translate X - x0 fails additive closure: "
      "0 of 729 base points give a subgroup (worst case %d pairs examined before a failure). Since "
      "X = x0 + G with G a subgroup and 0 in X forces G = X, this is also what the single check "
      "'X is not a subgroup' already gives]" % _worst)

FSET = set(FORMS)
_fout = [(i, j) for i in range(27) for j in range(27)
         if tuple((FORMS[i][w] + FORMS[j][w]) % 3 for w in range(27)) not in FSET]
_first = _fout[0] if _fout else None
check("form-set-not-closed",
      len(_fout) == 598 and _first is not None
      and LABELS[_first[0]] == "Qm2" and LABELS[_first[1]] == "Qm2",
      "[%d of the 729 sums of two rows of Table 1 are not rows of Table 1; in the table's own order the "
      "first failure is %s + %s, i.e. m = 2 = -1, the nonsquare the argument names]"
      % (len(_fout), LABELS[_first[0]], LABELS[_first[1]]))

# --- FACT 2 ---------------------------------------------------------------
def h_m(m):
    return tuple(add(mul(m, FROB[FROB[x]]), mul(power(m, 3), FROB[x])) for x in range(27))


check("fact2-doubling", all(tuple(add(h_m(m)[x], h_m(m)[x]) for x in range(27)) == h_m(NEG[m])
                            for m in range(27)),
      "[in characteristic 3, h_m + h_m = 2m x^9 + 2 m^3 x^3 = (2m) x^9 + (2m)^3 x^3 = h_{-m}, "
      "verified as maps for all 27 values of m]")

# --- FACT 3 ---------------------------------------------------------------
LIN = {}
for a0 in range(27):
    for a1 in range(27):
        for a2 in range(27):
            LIN[tuple(add(add(mul(a0, x), mul(a1, FROB[x])), mul(a2, FROB[FROB[x]]))
                      for x in range(27))] = (a0, a1, a2)
check("fact3-linearised-basis", len(LIN) == 19683 == 3 ** 9,
      "[the 27^3 = 3^9 maps a_0 x + a_1 x^3 + a_2 x^9 are pairwise distinct and all F_3-linear, so they "
      "exhaust End_{F_3}(GF(27)), whose cardinality is 3^9 = 19683: the linearised representation is unique]")

_coef = [LIN[h] for h in HMAPS]
_m_a0 = all(_coef[1 + i][0] == 0 for i in range(13))
_s_a0 = all(_coef[14 + i][0] == mul(_s_labels[i], _s_labels[i]) and _coef[14 + i][0] != 0
            for i in range(13))
check("fact3-families-disjoint", _m_a0 and _s_a0,
      "[every h_m has linearised coefficient a_0 = 0 (and a_1 = m^3, a_2 = m); every h_s has "
      "a_0 = s^2 != 0. The two families are separated by that one coefficient]")

HSET = set(HMAPS)
_wit = [m for m in NONSQUARES if NEG[m] in SQUARES and h_m(NEG[m]) not in HSET]
check("non-closure-mechanism", len(_wit) == 13 == len(NONSQUARES),
      "[for each of the 13 nonsquares m: -m is a square (FACT 1), h_m + h_m = h_{-m} (FACT 2), and "
      "h_{-m} is neither 0 nor any h_{m'} nor any h_s (FACT 3) -- 13 of 13. So the 27-element map set "
      "H is not closed under addition]")

# ---------------------------------------------------------------------------
# 7.  THE CONJECTURE, ASSEMBLED
# ---------------------------------------------------------------------------
_hyp = (len(XLIST) == D * D and _des_bad == 0 and _ang_viol == 0 and _out > 0 and _coset_closed == 0)
check("conjecture-4.9-false", _hyp,
      "[X is a P(T^27) 2-design (0 of 142506 characters fail) of size |X| = 729 = 27^2 that is not a "
      "subgroup of P(T^27) (indeed not a coset of one) and yet satisfies eq:toric-angle (0 of 265356 "
      "pairs fail). Conjecture 4.9 of arXiv:2311.13479, main.tex L786-788, is therefore false at d = 27]")

# ---------------------------------------------------------------------------
# 8.  THE d = 2 REMARK: THE LITERAL STATEMENT ALSO DIES TRIVIALLY
# ---------------------------------------------------------------------------
# Exact arithmetic in Z[zeta_8] = Z[z]/(z^4 + 1); an element is a 4-tuple of
# integer coefficients of 1, z, z^2, z^3.  zeta_8 = exp(i pi / 4).
def z8(m):
    m %= 8
    v = [0, 0, 0, 0]
    if m < 4:
        v[m] = 1
    else:
        v[m - 4] = -1
    return tuple(v)


def z8_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def z8_mul(a, b):
    c = [0] * 7
    for i in range(4):
        for j in range(4):
            c[i + j] += a[i] * b[j]
    return (c[0] - c[4], c[1] - c[5], c[2] - c[6], c[3])


def z8_conj(a):
    r = (0, 0, 0, 0)
    for i in range(4):
        r = z8_add(r, z8_mul((a[i], 0, 0, 0), z8(-i)))
    return r


def z8_norm(a):
    n = z8_mul(a, z8_conj(a))
    if n[1] != 0 or n[2] != 0 or n[3] != 0:
        return None                     # not a rational integer
    return n[0]


X2 = [(0, t) for t in PRINTED_D2_ANGLES]        # points of P(T^2), angles in units of pi/4
_a_ok, _a_vals = True, set()
for p, q in itertools.combinations(X2, 2):
    s = (0, 0, 0, 0)
    for j in range(2):
        s = z8_add(s, z8(p[j] - q[j]))
    n = z8_norm(s)
    if n is None:
        _a_ok = False
    else:
        _a_vals.add(n)
        if n not in (0, 2):
            _a_ok = False
_d2_chars, _d2_bad = 0, 0
for (a1, a2) in [(0, 0), (0, 1), (1, 1)]:
    for (b1, b2) in [(0, 0), (0, 1), (1, 1)]:
        kv = [0, 0]
        kv[a1] += 1
        kv[a2] += 1
        kv[b1] -= 1
        kv[b2] -= 1
        if kv == [0, 0]:
            continue
        _d2_chars += 1
        s = (0, 0, 0, 0)
        for x in X2:
            s = z8_add(s, z8(kv[0] * x[0] + kv[1] * x[1]))
        if s != (0, 0, 0, 0):
            _d2_bad += 1
check("d2-literal-statement-fails",
      _a_ok and sorted(_a_vals) == [0, 2] and _d2_chars == 6 and _d2_bad == 0
      and (0, 0) not in X2,
      "[X_2 = {pi/4, 3pi/4, 5pi/4, 7pi/4} in P(T^2) has |X_2| = 4 = 2^2, all %d of its nonzero "
      "characters vanish, all C(4,2) = 6 pairs give |sum|^2 in {0,2} = {0,d}, and it does not contain "
      "the identity, so it is not a subgroup: the LITERAL statement already fails at d = 2. But X_2 IS "
      "a coset of the subgroup {0, pi/2, pi, 3pi/2}, which is why the d = 27 object above, which "
      "contains the identity, is the substantive witness]" % _d2_chars)

# ---------------------------------------------------------------------------
# 9.  SCOPE AND VERDICT
# ---------------------------------------------------------------------------
print("")
print("NOT RE-RUN: nothing here concerns d = 6 or d = 8, or any characteristic-2 dimension. "
      "The polar-form step 'polar-form-is-trace-pairing' needs char F odd, and no even-characteristic "
      "(Galois-ring / Kerdock) construction was built or tested.")
print("NOT RE-RUN: no MUB matrices are constructed. The program verifies the spread and quadratic-form "
      "data from which Theorem 4.4 of arXiv:2311.13479 and Abdukhalikov's construction produce a "
      "complete set of 28 MUBs in C^27; it does not multiply out the 28 bases.")
print("NOT RE-RUN: MUB-equivalence is untouched. The program does not decide whether some complete set "
      "of MUBs MUB-equivalent to this one has a subgroup (or coset) projective toric design; the "
      "MUB-equivalence group is a continuum and was not searched. Nor does it recoordinatise the spread "
      "by sending each of the 28 components to infinity in turn.")
print("NOT RE-RUN: no minimality is claimed or tested. Nothing here says d = 27 is the least dimension, "
      "or the least prime-power dimension, admitting a non-group projective toric 2-design of size d^2 "
      "that satisfies the angle condition; no census of other dimensions was run.")
print("NOT RE-RUN: the bridging dictionary between group toric designs and semifield spreads is used "
      "only in the direction needed here, and only for this X. The program verifies steps (i)-(iii) of "
      "Lemma 3.1 computationally on the printed object, not as a general theorem.")
print("")
if _n_fail:
    print("VERDICT: %d CHECKS FAILED" % _n_fail)
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % _n_pass)
sys.exit(0)
