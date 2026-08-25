#!/usr/bin/env python3
"""Verification of a length-100 binary periodic Golay pair whose 50-compression
refutes a conjectured classification of compressions of periodic Golay pairs.

TAKEN FROM THE PAPER (inputs, transcribed literally into the block below):
  - the four length-10 seed sequences U, V, C, D as +/- strings;
  - the autocorrelation vector c, the base sequences p, q and their row sums;
  - the definitions of N, P, of r, e and of the (v/2)-compression;
  - the definition of Turyn multiplication T and the recursion (A_1,B_1)=(C,D);
  - the two conjectured forms of the compression, as functions of the row sums;
  - the five equivalence generators (swap, cyclic shift or reversal of the
    first sequence, simultaneous decimation, simultaneous alternating
    negation); the operations they induce at length two are written down as
    the paper states them but are then re-derived below at length 100;
  - the numbers the paper reports: the compressions ([12,2],[4,-6]) and
    ([12,2],[-4,6]), row sums (14,2), the forms ([0,14],[0,2]) and
    ([8,6],[8,-6]), the maxima 12/14/8, z_1, w_1, z_2, w_2, t = 1+3i and the
    ratio 4/5 + 3/5 i.

DERIVED HERE (computed from the seeds; no number below is assumed):
  - the seeds' autocorrelations, and that (U,V) and (C,D) are Golay pairs;
  - p, q recomputed from U, V, with p_j q_j = 0, p_j^2 + q_j^2 = 1 and the
    identity that makes the Turyn factor identically 10;
  - the length-100 pair itself, built by Turyn multiplication, printed back as
    +/- strings, counted, and shown binary;
  - that it is a Golay pair AND a periodic Golay pair, from the periodic
    definition directly, at all 99 nonzero shifts (every hypothesis of the
    conjecture);
  - its 50-compression from the definition, its row sums, and the fact that the
    compression is itself a length-two periodic Golay pair of energy 2v;
  - the two conjectured forms, instantiated from the computed row sums;
  - the load-bearing conclusion: the full equivalence class of the compression
    under the induced operations is enumerated by breadth-first search and
    contains neither form; independently, the unit-condition obstruction is
    evaluated in exact Gaussian-integer arithmetic and fails for the
    compression while holding for both forms;
  - the action of the five generators on the compression: each generator is
    applied to the actual length-100 pair (all 100 cyclic shifts, the
    reversal composed with every shift, all admissible decimations, the swap
    and simultaneous alternating negation), the image is recompressed, and the
    result is required to equal the transcribed operation and to lie in the
    enumerated class -- and every image is checked to be a periodic Golay pair
    still, so these really are equivalence operations;
  - the tower: lengths 10^1..10^6 built, compressed and refuted, with the
    autocorrelation conditions recomputed in full up to the ceiling the
    printed NOT RE-RUN line states (no exponent is written down twice);
  - the number-theoretic step: the unit group of Z[i], the factorisation of t,
    non-associateness of 2+i and 2-i and of the powers t^k, conj(t)^k, and the
    ratio (-i z_1)/conj(z_1) as an exact rational pair.

Exact integer and rational arithmetic throughout; no floating point anywhere.
"""
import sys
import time
from fractions import Fraction
from operator import mul

CHECKS = []

# ---------------------------------------------------------------- paper inputs
U_STR = "+--++++++-"
V_STR = "+--+-+---+"
C_STR = "-++-+-+++-"
D_STR = "-++++++--+"
C_AUTO = (3, 0, 1, 0, -1, -2, 1, 2, -1)          # claimed (N_U(s))_{s=1..9}
P_CLAIM = (1, -1, -1, 1, 0, 1, 0, 0, 0, 0)        # claimed p
Q_CLAIM = (0, 0, 0, 0, 1, 0, 1, 1, 1, -1)         # claimed q
P_SUM_CLAIM, Q_SUM_CLAIM = 1, 3                   # eq. (sums)
H2_CLAIM = ((12, 2), (4, -6))                     # compression of (A_2, B_2)
H2T_CLAIM = ((12, 2), (-4, 6))                    # after negating B_2
AB_CLAIM = (14, 2)                                # row sums a, b of H2T
FORM1_CLAIM = ((0, 14), (0, 2))
FORM2_CLAIM = ((8, 6), (8, -6))
MAXABS_CLAIM = (12, 14, 8)                        # of H2T, form 1, form 2
Z1_CLAIM, W1_CLAIM = (2, 4), (4, -2)              # z_1, w_1 as (re, im)
Z2_CLAIM, W2_CLAIM = (14, -2), (10, 10)           # z_2, w_2
RATIO_CLAIM = (Fraction(4, 5), Fraction(3, 5))    # (-i z_1) / conj(z_1)
T_CLAIM = (1, 3)                                  # t = 1 + 3i
UNITS = ((1, 0), (-1, 0), (0, 1), (0, -1))        # the group U


def decode(s):
    """'+'/'-' string -> tuple of +-1 integers."""
    return tuple(1 if ch == "+" else -1 for ch in s)


def N(X, s):
    """Aperiodic autocorrelation sum_{j=0}^{v-1-s} x_j x_{j+s}."""
    return sum(map(mul, X, X[s:]))


def P(X, s):
    """Periodic autocorrelation sum_{j=0}^{v-1} x_j x_{j+s mod v}."""
    return sum(map(mul, X, X[s:] + X[:s]))


def rowsum(X):
    return sum(X)


def altsum(X):
    return sum(X[j] if j % 2 == 0 else -X[j] for j in range(len(X)))


def compress(X):
    """(v/2)-compression [(r+e)/2, (r-e)/2]; exact integer division."""
    r, e = rowsum(X), altsum(X)
    if (r + e) % 2 or (r - e) % 2:
        raise ValueError("compression is not integral")
    return ((r + e) // 2, (r - e) // 2)


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)
    return bool(ok)


def finish():
    n = len(CHECKS)
    bad = [c for c in CHECKS if not c[1]]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


def check_seed_pairs():
    """Check 1: the four length-10 seeds and their aperiodic autocorrelations."""
    U, V, Cs, D = decode(U_STR), decode(V_STR), decode(C_STR), decode(D_STR)
    wf = all(len(X) == 10 and set(X) <= {1, -1} for X in (U, V, Cs, D))
    ck("seed_sequences_well_formed", wf,
       "U=%s V=%s C=%s D=%s (length 10, entries +-1)" % (U_STR, V_STR, C_STR, D_STR))
    nU = tuple(N(U, s) for s in range(1, 10))
    nV = tuple(N(V, s) for s in range(1, 10))
    nC = tuple(N(Cs, s) for s in range(1, 10))
    nD = tuple(N(D, s) for s in range(1, 10))
    neg = tuple(-x for x in C_AUTO)
    ck("seed_autocorrelations_equal_paper_vector",
       nU == C_AUTO and nD == C_AUTO and nV == neg and nC == neg,
       "N_U=%s N_D=%s N_V=%s N_C=%s" % (nU, nD, nV, nC))
    ck("seed_pairs_are_golay",
       all(nU[i] + nV[i] == 0 for i in range(9))
       and all(nC[i] + nD[i] == 0 for i in range(9)),
       "N_U+N_V=0 and N_C+N_D=0 for s=1..9")
    return U, V, Cs, D


def check_pq(U, V):
    """Check 2: the base sequences p, q and the identity making T binary."""
    p = tuple((U[j] + V[j]) // 2 for j in range(10))
    q = tuple((U[j] - V[j]) // 2 for j in range(10))
    ck("pq_from_seeds_match_paper", p == P_CLAIM and q == Q_CLAIM,
       "p=%s q=%s" % (p, q))
    ck("pq_disjoint_support_unit_norm",
       all(p[j] * q[j] == 0 and p[j] ** 2 + q[j] ** 2 == 1 for j in range(10)),
       "p_j q_j = 0 and p_j^2 + q_j^2 = 1 for all j")
    ck("pq_row_sums", sum(p) == P_SUM_CLAIM and sum(q) == Q_SUM_CLAIM,
       "sum p = %d, sum q = %d" % (sum(p), sum(q)))
    energy = sum(x * x for x in p) + sum(x * x for x in q)
    off = [N(p, s) + N(q, s) for s in range(1, 10)]
    ck("pq_turyn_factor_is_ten", energy == 10 and all(x == 0 for x in off),
       "P(z)P(1/z)+Q(z)Q(1/z) = 10 identically; energy=%d, off-shifts=%s"
       % (energy, off))
    return p, q


def turyn(X, Y, p, q):
    """T(X,Y) of the paper: concatenate p_j X + q_j Y^R and -q_j X^R + p_j Y."""
    XR, YR = X[::-1], Y[::-1]
    A, B = [], []
    for j in range(10):
        A.extend(p[j] * X[i] + q[j] * YR[i] for i in range(len(X)))
        B.extend(-q[j] * XR[i] + p[j] * Y[i] for i in range(len(Y)))
    return tuple(A), tuple(B)


def check_exhibited_object(Cs, D, p, q):
    """Check 3: the length-100 pair is well formed, decoded and counted."""
    A2, B2 = turyn(Cs, D, p, q)
    ck("exhibited_pair_is_binary_length_100",
       len(A2) == 100 and len(B2) == 100
       and set(A2) == {1, -1} and set(B2) == {1, -1},
       "|A|=%d |B|=%d, entries +-1, +1 count A=%d B=%d"
       % (len(A2), len(B2), A2.count(1), B2.count(1)))
    enc = lambda X: "".join("+" if x == 1 else "-" for x in X)
    print("  A = " + enc(A2))
    print("  B = " + enc(B2))
    ck("exhibited_pair_decodes_reversibly",
       decode(enc(A2)) == A2 and decode(enc(B2)) == B2,
       "printed +/- strings re-decode to the same sequences")
    return A2, B2


def check_golay(A2, B2):
    """Check 4: hypotheses -- (A,B) is a Golay pair, hence a periodic one."""
    v = len(A2)
    off = [N(A2, s) + N(B2, s) for s in range(1, v)]
    ck("exhibited_pair_is_a_golay_pair",
       len(A2) == len(B2) and all(x == 0 for x in off),
       "N_A(s)+N_B(s)=0 for all %d shifts; max |value| = %d"
       % (len(off), max(abs(x) for x in off)))
    poff = [P(A2, s) + P(B2, s) for s in range(1, v)]
    ck("exhibited_pair_is_a_periodic_golay_pair", all(x == 0 for x in poff),
       "P_A(s)+P_B(s)=0 for all %d nonzero shifts, computed from the "
       "periodic definition; max |value| = %d"
       % (len(poff), max(abs(x) for x in poff)))
    nB2 = tuple(-x for x in B2)
    poff2 = [P(A2, s) + P(nB2, s) for s in range(1, v)]
    ck("negated_partner_is_still_a_periodic_golay_pair",
       all(x == 0 for x in poff2)
       and all(N(A2, s) + N(nB2, s) == 0 for s in range(1, v)),
       "(A, -B) also satisfies both Golay conditions at all %d nonzero shifts"
       % len(poff2))
    return nB2


def check_compression(A2, B2, nB2):
    """Check 5: the 50-compression equals the pair the paper exhibits."""
    H2 = (compress(A2), compress(B2))
    H2t = (compress(A2), compress(nB2))
    ck("compression_of_exhibited_pair_matches_paper",
       H2 == H2_CLAIM and H2t == H2T_CLAIM,
       "H_2 = (%s,%s), after negating the partner (%s,%s)"
       % (list(H2[0]), list(H2[1]), list(H2t[0]), list(H2t[1])))
    a, b = sum(H2t[0]), sum(H2t[1])
    ck("row_sums_of_compression", (a, b) == AB_CLAIM, "a = %d, b = %d" % (a, b))
    # the compression is itself a periodic Golay pair of length two, with
    # total energy 2v -- both are hypotheses the conjecture's objects satisfy.
    paf1 = 2 * (H2t[0][0] * H2t[0][1] + H2t[1][0] * H2t[1][1])
    energy = sum(x * x for x in H2t[0]) + sum(x * x for x in H2t[1])
    ck("compression_is_a_length_two_periodic_golay_pair",
       paf1 == 0 and energy == 2 * len(A2),
       "P(1) sum = %d, energy = %d = 2v with v = %d"
       % (paf1, energy, len(A2)))
    return H2, H2t, a, b, energy


def forms(a, b):
    """The two pairs of eq. (forms), built from the row sums, not transcribed."""
    if (a + b) % 2 or (a - b) % 2:
        raise ValueError("row sums of opposite parity")
    f1 = ((0, a), (0, b))
    f2 = (((a + b) // 2, (a - b) // 2), ((a + b) // 2, (b - a) // 2))
    return f1, f2


def check_forms(a, b, energy):
    """Check 6: the two conjectured representatives for these row sums.

    `energy` is the energy computed from the exhibited compression, not a
    constant: the forms must match the object they are being compared with."""
    f1, f2 = forms(a, b)
    ck("conjectured_forms_instantiated",
       f1 == FORM1_CLAIM and f2 == FORM2_CLAIM,
       "([0,a],[0,b]) = (%s,%s) and the second form = (%s,%s)"
       % (list(f1[0]), list(f1[1]), list(f2[0]), list(f2[1])))
    ok = True
    for nm, F in (("form1", f1), ("form2", f2)):
        if sum(F[0]) != a or sum(F[1]) != b:
            ok = False
        if 2 * (F[0][0] * F[0][1] + F[1][0] * F[1][1]) != 0:
            ok = False
        if sum(x * x for x in F[0]) + sum(x * x for x in F[1]) != energy:
            ok = False
    ck("conjectured_forms_are_admissible_candidates", ok,
       "both forms have row sums (%d,%d), vanishing shift-1 sum and the same "
       "energy %d as the exhibited compression, so the refutation is not an "
       "arithmetic mismatch" % (a, b, energy))
    return f1, f2


# --------------------------------------- exact Gaussian-integer arithmetic
def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def gconj(u):
    return (u[0], -u[1])


def gpow(u, k):
    r = (1, 0)
    for _ in range(k):
        r = gmul(r, u)
    return r


def gdiv(a, b):
    """Exact quotient a/b in Z[i], or None if b does not divide a."""
    nb = b[0] ** 2 + b[1] ** 2
    num = gmul(a, gconj(b))
    if num[0] % nb or num[1] % nb:
        return None
    return (num[0] // nb, num[1] // nb)


def zw(H):
    (x0, x1), (y0, y1) = H
    return (x0 + x1, y0 + y1), (x0 - x1, y0 - y1)


def unit_condition(H):
    """True iff w(H) lies in U z(H) union U conj(z(H)); exact integers."""
    z, w = zw(H)
    orb = [gmul(u, z) for u in UNITS] + [gmul(u, gconj(z)) for u in UNITS]
    return w in orb


def check_obstruction(H2, H2t, f1, f2):
    """Check 7: the obstruction of the lemma, computed on all four objects."""
    z, w = zw(H2t)
    ck("z_and_w_of_compression",
       zw(H2) == (Z2_CLAIM, W2_CLAIM)
       and (z, w) == (gconj(Z2_CLAIM), gconj(W2_CLAIM)),
       "for H_2: z = %s, w = %s; for the negated partner z = %s, w = %s"
       % (Z2_CLAIM, W2_CLAIM, z, w))
    ck("both_conjectured_forms_satisfy_the_obstruction",
       unit_condition(f1) and unit_condition(f2),
       "w = -z for form 1 and w = i conj(z) for form 2, so any pair "
       "equivalent to a form must satisfy the unit condition")
    ck("compression_violates_the_obstruction",
       (not unit_condition(H2t)) and (not unit_condition(H2)),
       "w(H) is in neither U z(H) nor U conj(z(H)): w = %s, U z = %s, "
       "U conj z = %s" % (w, [gmul(u, z) for u in UNITS],
                          [gmul(u, gconj(z)) for u in UNITS]))


# ------------------- the induced equivalence operations on length-two pairs
# These five maps are written down as the paper states them, but they are NOT
# taken on trust: check_induced_action below applies the corresponding
# generators to the actual length-100 pair, recompresses, and requires the
# result to agree with the map used here (and to stay a periodic Golay pair).
def op_swap(H):
    return (H[1], H[0])


def op_shift(H):            # cyclic shift of the first sequence
    return ((H[0][1], H[0][0]), H[1])


def op_reverse(H):          # reversal of the first sequence
    return ((H[0][1], H[0][0]), H[1])


def op_decimate(H):         # admissible decimations are odd: trivial here
    return H


def op_altneg(H):           # simultaneous alternating negation
    return ((H[0][0], -H[0][1]), (H[1][0], -H[1][1]))


def op_negfirst(H):
    return ((-H[0][0], -H[0][1]), H[1])


def op_negsecond(H):
    return (H[0], (-H[1][0], -H[1][1]))


CORE_OPS = (op_swap, op_shift, op_reverse, op_decimate, op_altneg)
EXT_OPS = CORE_OPS + (op_negfirst, op_negsecond)


def orbit(H, ops):
    """Full equivalence class of a length-two pair under the given operations."""
    seen, frontier = {H}, [H]
    while frontier:
        nxt = []
        for X in frontier:
            for g in ops:
                Y = g(X)
                if Y not in seen:
                    seen.add(Y)
                    nxt.append(Y)
        frontier = nxt
    return seen


def check_inequivalence(H2t, H2, f1, f2, a, b):
    """Check 8 (load bearing): the exhibited compression is equivalent to
    neither conjectured form -- decided by exhausting the equivalence class."""
    core = orbit(H2t, CORE_OPS)
    ext = orbit(H2t, EXT_OPS)
    ck("equivalence_class_closed_and_negation_generated",
       core == ext and op_negsecond(H2t) in core and op_negfirst(H2t) in core,
       "class size %d under the five stated operations, and the same set "
       "under the larger system that also negates one sequence outright"
       % len(core))
    ck("compression_is_equivalent_to_neither_conjectured_form",
       f1 not in ext and f2 not in ext,
       "([0,%d],[0,%d]) and (%s,%s) are both absent from the %d-element "
       "equivalence class of (%s,%s)"
       % (a, b, list(f2[0]), list(f2[1]), len(ext),
          list(H2t[0]), list(H2t[1])))
    g1, g2 = forms(sum(H2[0]), sum(H2[1]))
    orb2 = orbit(H2, EXT_OPS)
    ck("unnegated_compression_also_refutes_the_conjecture",
       g1 not in orb2 and g2 not in orb2,
       "for H_2 the row sums are (%d,%d), forms (%s,%s) and (%s,%s), class "
       "size %d, neither form present"
       % (sum(H2[0]), sum(H2[1]), list(g1[0]), list(g1[1]),
          list(g2[0]), list(g2[1]), len(orb2)))
    ck("obstruction_is_an_equivalence_invariant",
       all(not unit_condition(X) for X in ext)
       and all(unit_condition(X) for X in orbit(f1, EXT_OPS))
       and all(unit_condition(X) for X in orbit(f2, EXT_OPS)),
       "the unit condition is constant on each class: false throughout the "
       "class of the exhibited pair, true throughout both form classes")
    inv = {max(abs(t) for s in X for t in s) for X in ext}
    trip = (max(abs(t) for s in H2t for t in s),
            max(abs(t) for s in f1 for t in s),
            max(abs(t) for s in f2 for t in s))
    ck("maximum_absolute_entry_separates_the_three_pairs",
       len(inv) == 1 and trip == MAXABS_CLAIM,
       "invariant on the class of the exhibited pair (value %s); the three "
       "maxima are %s, so no two of the pairs can be equivalent"
       % (sorted(inv), trip))
    return ext


# ------------- the same five generators, acting on the actual sequences
def seq_shift(X, k):
    """Cyclic shift of a single sequence."""
    return X[k:] + X[:k]


def seq_reverse(X):
    return X[::-1]


def seq_decimate(X, d):
    """Decimation x_j -> x_{dj mod v}; admissible means gcd(d, v) = 1."""
    v = len(X)
    return tuple(X[(d * j) % v] for j in range(v))


def seq_altneg(X):
    """Alternating negation x_j -> (-1)^j x_j."""
    return tuple(X[j] if j % 2 == 0 else -X[j] for j in range(len(X)))


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def check_induced_action(A2, B2, H2t, cls):
    """Check 8b (load bearing, DERIVED not transcribed): the five equivalence
    generators, applied to the actual length-100 pair, induce on its
    50-compression exactly the operations op_* used to enumerate the class."""
    v = len(A2)
    dec = [d for d in range(1, v) if gcd(d, v) == 1]
    named = [("swap_the_sequences", (B2, A2), op_swap),
             ("cyclic_shift_of_the_first", (seq_shift(A2, 1), B2), op_shift),
             ("reversal_of_the_first", (seq_reverse(A2), B2), op_reverse),
             ("simultaneous_alternating_negation",
              (seq_altneg(A2), seq_altneg(B2)), op_altneg)]
    named += [("simultaneous_decimation_by_%d" % d,
               (seq_decimate(A2, d), seq_decimate(B2, d)), op_decimate)
              for d in dec]
    hom = []
    for nm, (X, Y), op in named:
        got = (compress(X), compress(Y))
        if got != op(H2t):
            hom.append((nm, got, op(H2t)))
    ck("induced_action_on_the_compression_is_derived_not_assumed",
       not hom and dec and all(d % 2 for d in dec),
       "for the swap, a cyclic shift and a reversal of the first sequence, "
       "simultaneous alternating negation and all %d admissible decimations "
       "(every d with gcd(d,%d)=1, each of them odd), the %d-compression of "
       "the transformed length-%d pair equals the transcribed operation "
       "applied to the compression; %d mismatches %s"
       % (len(dec), v, v // 2, v, len(hom), hom[:4]))
    images = [(seq_shift(A2, k), B2) for k in range(v)]
    images += [(seq_reverse(seq_shift(A2, k)), B2) for k in range(v)]
    images += [pr for _, pr, _ in named]
    imgs = [(compress(X), compress(Y)) for X, Y in images]
    outside = [h for h in imgs if h not in cls]
    distinct = set(imgs)
    # len(distinct) > 1 is the non-vacuity condition: if every generator fixed
    # the compression, membership in the class would be trivially satisfied.
    ck("every_generated_compression_lands_in_the_enumerated_class",
       not outside and len(distinct) > 1,
       "%d images of the length-%d pair (all %d cyclic shifts of the first "
       "sequence, its reversal composed with every shift, all %d admissible "
       "decimations, the swap and alternating negation) recompress into the "
       "%d-element class, realising %d distinct compressions so the generators "
       "do move it; %d outside %s"
       % (len(images), v, v, len(dec), len(cls), len(distinct),
          len(outside), outside[:4]))
    bad = [nm for nm, (X, Y), _ in named
           if any(P(X, s) + P(Y, s) for s in range(1, v))
           or set(X) | set(Y) != {1, -1}]
    pg = ["image_%d" % i for i, (X, Y) in enumerate(images)
          if any(P(X, s) + P(Y, s) for s in range(1, v))
          or set(X) | set(Y) != {1, -1} or len(X) != v or len(Y) != v]
    ck("the_generators_preserve_the_periodic_golay_property",
       not bad and not pg and len(images) == 2 * v + len(named) and images,
       "every one of the %d images is binary of length %d with "
       "P_A(s)+P_B(s)=0 at all %d nonzero shifts, recomputed from the "
       "periodic definition, so these are equivalence operations on the "
       "objects of the conjecture; failures: %s"
       % (len(images), v, v - 1, (bad + pg)[:4]))


TOWER_TOP = 6          # sequences of length 10^6 are built and compressed
GOLAY_TOP = 4          # autocorrelations re-run in full up to length 10^4
ASSOC_TOP = 40         # exponents tested for the Gaussian-associate step


def check_tower(Cs, D, p, q):
    """Check 9: the recurrences, closed forms and refutation for n = 1..TOP."""
    A, B = Cs, D
    zs, ws, rec_ok, golay_ok, refute_ok, energy_ok = [], [], True, True, True, True
    for n in range(1, TOWER_TOP + 1):
        t0 = time.perf_counter()
        z, w = (rowsum(A), rowsum(B)), (altsum(A), altsum(B))
        zs.append(z)
        ws.append(w)
        H = (compress(A), compress(B))
        if len(A) != 10 ** n or len(B) != 10 ** n or set(A) | set(B) != {1, -1}:
            golay_ok = False
        if n <= GOLAY_TOP:
            print("  tower level n=%d: recomputing all %d aperiodic and "
                  "periodic autocorrelations ..." % (n, 2 * (10 ** n - 1)))
            sys.stdout.flush()
            if any(N(A, s) + N(B, s) for s in range(1, len(A))):
                golay_ok = False
            if any(P(A, s) + P(B, s) for s in range(1, len(A))):
                golay_ok = False
        if zw(H) != (z, w):
            rec_ok = False
        e2 = sum(x * x for x in H[0]) + sum(x * x for x in H[1])
        p1 = H[0][0] * H[0][1] + H[1][0] * H[1][1]
        if e2 != 2 * 10 ** n or p1 != 0:
            energy_ok = False
        if n >= 2:
            f1, f2 = forms(sum(H[0]), sum(H[1]))
            cls = orbit(H, EXT_OPS)
            if unit_condition(H) or f1 in cls or f2 in cls:
                refute_ok = False
        print("  tower level n=%d: length %d built and compressed to (%s,%s)"
              "%s [%.1f s]"
              % (n, 10 ** n, list(H[0]), list(H[1]),
                 "; autocorrelations recomputed in full"
                 if n <= GOLAY_TOP else
                 "; autocorrelations not recomputed at this length",
                 time.perf_counter() - t0))
        sys.stdout.flush()
        if n < TOWER_TOP:
            A, B = turyn(A, B, p, q)
    ck("tower_is_binary_and_golay_where_recomputed",
       golay_ok and min(GOLAY_TOP, TOWER_TOP) >= 3 and TOWER_TOP >= 2,
       "lengths 10^1..10^%d are binary of the right length; Golay and "
       "periodic Golay conditions recomputed in full up to length 10^%d"
       % (TOWER_TOP, GOLAY_TOP))
    ck("compression_matches_row_and_alternating_sums", rec_ok
       and zs[0] == Z1_CLAIM and ws[0] == W1_CLAIM and zs[1] == Z2_CLAIM
       and ws[1] == W2_CLAIM,
       "z_n = r(A_n)+i r(B_n) and w_n = e(A_n)+i e(B_n) reproduce z(H_n), "
       "w(H_n); z_1 = %s, w_1 = %s, z_2 = %s, w_2 = %s"
       % (zs[0], ws[0], zs[1], ws[1]))
    tb = gconj(T_CLAIM)
    closed = all(zs[n] == gmul(gpow(tb, n), zs[0])
                 and ws[n] == gmul(gmul((0, -1), gpow(T_CLAIM, n)), zs[0])
                 for n in range(len(zs)))
    rec = all(zs[n + 1] == gmul((1, -3), zs[n])
              and ws[n + 1] == gmul((1, 3), ws[n]) for n in range(len(zs) - 1))
    ck("recurrences_and_closed_forms_hold", closed and rec and ws[0] == gmul((0, -1), zs[0]),
       "z_{n+1} = (1-3i)z_n, w_{n+1} = (1+3i)w_n, z_n = conj(t)^{n-1} z_1, "
       "w_n = -i t^{n-1} z_1, and w_1 = -i z_1, for n up to %d" % TOWER_TOP)
    ck("energy_and_shift_sum_of_every_compression", energy_ok,
       "each compression is a length-two periodic Golay pair of energy 2v")
    ck("refutation_holds_at_every_length_recomputed", refute_ok,
       "for n = 2..%d the compression violates the unit condition and its "
       "equivalence class contains neither conjectured form" % TOWER_TOP)


def check_number_theory():
    """Check 10: the two arithmetic facts the infinitude argument rests on."""
    box = [(a, b) for a in range(-4, 5) for b in range(-4, 5)
           if a * a + b * b == 1]
    closed = all(gmul(u, w) in UNITS for u in UNITS for w in UNITS)
    inv = all(any(gmul(u, w) == (1, 0) for w in UNITS) for u in UNITS)
    ck("unit_group_of_the_gaussian_integers", set(box) == set(UNITS)
       and len(UNITS) == 4 and closed and inv
       and all(gconj(u) in UNITS for u in UNITS),
       "U = %s is exactly the set of Gaussian integers of norm 1, and is "
       "closed under multiplication, inversion and conjugation" % (UNITS,))
    t = T_CLAIM
    g = gdiv(t, (1, 1))                 # the odd prime factor of t, derived
    ck("factorisation_of_t", g is not None and gmul((1, 1), g) == t
       and gmul((1, -1), gconj(g)) == gconj(t)
       and g[0] ** 2 + g[1] ** 2 == 5,
       "t/(1+i) = %s is a Gaussian integer of norm %s; conj(t) = (1-i) times "
       "its conjugate" % (g, None if g is None else g[0] ** 2 + g[1] ** 2))
    nonassoc = g is not None and all(gmul(u, gconj(g)) != g for u in UNITS)
    ck("odd_prime_factor_and_its_conjugate_are_nonassociate", nonassoc,
       "no unit u satisfies u conj(g) = g for g = %s, and norm 5 is prime so "
       "g is a Gaussian prime" % (g,))
    powers = all(all(gmul(u, gpow(gconj(t), k)) != gpow(t, k) for u in UNITS)
                 for k in range(1, ASSOC_TOP + 1))
    ck("t_and_its_conjugate_have_nonassociate_powers", powers,
       "t^k is not a unit multiple of conj(t)^k for k = 1..%d, so the case "
       "w_n = u z_n is impossible in those lengths" % ASSOC_TOP)
    z1 = Z1_CLAIM
    num = gmul((0, -1), z1)
    den = gconj(z1)
    nn = den[0] ** 2 + den[1] ** 2
    q = (Fraction(gmul(num, gconj(den))[0], nn),
         Fraction(gmul(num, gconj(den))[1], nn))
    ck("ratio_minus_i_z1_over_conj_z1", q == RATIO_CLAIM
       and all(q != (Fraction(u[0]), Fraction(u[1])) for u in UNITS),
       "(-i z_1)/conj(z_1) = %s + %s i, exact rationals, not a unit, so the "
       "case w_n = u conj(z_n) is impossible for every n" % (q[0], q[1]))


def main():
    U, V, Cs, D = check_seed_pairs()
    p, q = check_pq(U, V)
    A2, B2 = check_exhibited_object(Cs, D, p, q)
    nB2 = check_golay(A2, B2)
    H2, H2t, a, b, energy = check_compression(A2, B2, nB2)
    f1, f2 = check_forms(a, b, energy)
    check_obstruction(H2, H2t, f1, f2)
    cls = check_inequivalence(H2t, H2, f1, f2, a, b)
    check_induced_action(A2, nB2, H2t, cls)
    check_tower(Cs, D, p, q)
    check_number_theory()
    print("NOT RE-RUN: the claim that the conjecture holds for every length "
          "v < 100, which is what makes length 100 least possible, is a "
          "published exhaustive classification and is not recomputed here; "
          "inside the tower, Golay and periodic-Golay autocorrelations are "
          "recomputed in full only up to length 10^%d (the refuting length "
          "100 = 10^2 is itself verified in full above, while lengths 10^%d to "
          "10^%d are built, compressed and refuted but their autocorrelations "
          "are not recomputed -- the per-level progress lines above report how "
          "many autocorrelations each level recomputed and how long it took); "
          "the Gaussian-associate step "
          "is tested for exponents up to %d rather than all n; and the action "
          "of the five equivalence generators on a length-two compression, "
          "although derived here rather than assumed, is derived only for the "
          "exhibited length-100 pair (each generator applied once, plus every "
          "cyclic shift, every shift-then-reversal and every admissible "
          "decimation) -- the paper's general statement of that action at "
          "arbitrary even length, and its cited claim that negation of one "
          "sequence is generated by the five operations, are not proved here, "
          "though the latter is confirmed as a set identity on the class of "
          "the exhibited compression."
          % (GOLAY_TOP, GOLAY_TOP + 1, TOWER_TOP, ASSOC_TOP))


if __name__ == "__main__":
    try:
        main()
        ck("program_ran_to_completion", True, "no unhandled exception")
    except Exception as exc:                     # a crash must not read as a pass
        ck("program_ran_to_completion", False, "unhandled exception: %r" % (exc,))
    finish()
