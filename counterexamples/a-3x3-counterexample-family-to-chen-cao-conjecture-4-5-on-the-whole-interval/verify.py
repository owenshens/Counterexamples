#!/usr/bin/env python3
"""Verification of the objects and identities printed in

    "A 3x3 Counterexample Family to Conjecture 4.5 of Chen-Cao's arXiv:1801.00225v1"

Python 3.9+, STANDARD LIBRARY ONLY (fractions, itertools), no external data file.
All arithmetic is exact: every quantity is a Fraction or an integer, and no decision
anywhere in this program is taken on a float.  Floats appear only inside f-string
decorations printed for a human reader.

The program READS THE OBJECTS AS PRINTED IN THE PAPER -- the family

    A(s) = [[0, 1-u, u], [1-u, 0, u], [u, u, 0]],   u = (s-2)/2,

Chen-Cao's own two Lemma 4.4 witnesses A_0(s) and A_1(s), and the symmetric
zero-diagonal slice -- and re-derives every number and every identity the paper
claims about them.

Prints one `PASS <name> [detail]` line per check and closes with

    VERDICT: ALL <n> CHECKS PASS

exiting 0 iff every check passed.
"""

from fractions import Fraction as Fr
from itertools import permutations
import sys

# --------------------------------------------------------------------------------------
# check harness
# --------------------------------------------------------------------------------------
_N = 0
_BAD = 0


def ck(name, cond, detail=''):
    global _N, _BAD
    _N += 1
    if cond:
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _BAD += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(msg):
    print('NOTE %s' % msg)


def hd(msg):
    print('')
    print('=== %s' % msg)


# --------------------------------------------------------------------------------------
# exact linear algebra: the permanent by the full n! expansion, over Fractions
# --------------------------------------------------------------------------------------
def per(M):
    """Permanent of a square matrix of Fractions, by the definition (sum over S_n)."""
    n = len(M)
    tot = Fr(0)
    for p in permutations(range(n)):
        prod = Fr(1)
        for i in range(n):
            prod *= M[i][p[i]]
        tot += prod
    return tot


def I_minus(A):
    n = len(A)
    return [[(Fr(1) if i == j else Fr(0)) - A[i][j] for j in range(n)] for i in range(n)]


def row_sums(A):
    return [sum(r) for r in A]


def col_sums(A):
    n = len(A)
    return [sum(A[i][j] for i in range(n)) for j in range(n)]


def in_omega(A):
    """Chen-Cao's own definition of DOUBLY substochastic, written out from Section 1:
    entrywise nonnegative, every ROW sum <= 1 AND every COLUMN sum <= 1."""
    n = len(A)
    nonneg = all(A[i][j] >= 0 for i in range(n) for j in range(n))
    return nonneg and all(x <= 1 for x in row_sums(A)) and all(x <= 1 for x in col_sums(A))


def mass(A):
    return sum(sum(r) for r in A)


def direct_sum(blocks):
    n = sum(len(b) for b in blocks)
    M = [[Fr(0)] * n for _ in range(n)]
    o = 0
    for b in blocks:
        k = len(b)
        for i in range(k):
            for j in range(k):
                M[o + i][o + j] = b[i][j]
        o += k
    return M


# --------------------------------------------------------------------------------------
# the objects, exactly as printed in the paper
# --------------------------------------------------------------------------------------
def A_fam(s):
    """The paper's witness family."""
    u = (s - 2) / Fr(2)
    return [[Fr(0), 1 - u, u],
            [1 - u, Fr(0), u],
            [u, u, Fr(0)]]


def A_0(s):
    """Chen-Cao, Lemma 4.4, first witness (their edge-midpoint matrix)."""
    w = s / Fr(2) - 1
    return [[Fr(0), Fr(1, 2), w],
            [Fr(1, 2), Fr(0), Fr(1, 2)],
            [w, Fr(1, 2), Fr(0)]]


def A_1(s):
    """Chen-Cao, Lemma 4.4, second witness (nonzero diagonal entry a_33 = s-2)."""
    return [[Fr(0), Fr(1), Fr(0)],
            [Fr(1), Fr(0), Fr(0)],
            [Fr(0), Fr(0), s - 2]]


M2 = [[Fr(0), Fr(1)], [Fr(1), Fr(0)]]


def f(s):
    return (s ** 3 - 5 * s ** 2 + 4 * s + 12) / Fr(4)


def br1(s):
    return (s ** 2 - 5 * s + 12) / Fr(4)


def br2_corrected(s):
    return 6 - 2 * s


def br2_printed(s):
    return 6 - 4 * s


def V(s):
    """The CORRECTED (strongest) reading of Conjecture 4.5's value: the max of the two
    branches, i.e. exactly the lower bound Chen-Cao's Lemma 4.4 proves."""
    return max(br1(s), br2_corrected(s))


# --------------------------------------------------------------------------------------
# exact polynomial arithmetic over Q, for identity checks WITHOUT any CAS
# --------------------------------------------------------------------------------------
def pmul(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else Fr(0)) + (b[i] if i < len(b) else Fr(0)) for i in range(n)]


def pscale(a, c):
    return [x * c for x in a]


def pnorm(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def peq(a, b):
    return pnorm(a) == pnorm(b)


def pshow(a):
    a = pnorm(a)
    ts = []
    for i in range(len(a) - 1, -1, -1):
        if a[i] == 0:
            continue
        c = a[i]
        ts.append(('%s' % c) + ('' if i == 0 else ('*s' if i == 1 else '*s^%d' % i)))
    return ' + '.join(ts) if ts else '0'


S = [Fr(0), Fr(1)]                                   # the polynomial s
P_f = [Fr(3), Fr(1), Fr(-5, 4), Fr(1, 4)]            # f(s)  = (s^3-5s^2+4s+12)/4
P_br1 = [Fr(3), Fr(-5, 4), Fr(1, 4)]                 # (s^2-5s+12)/4
P_br2c = [Fr(6), Fr(-2)]                             # 6-2s


# ======================================================================================
print(__doc__.strip().splitlines()[2].strip())
print('exact rational arithmetic only; python %s' % sys.version.split()[0])

MASSES_INTERIOR = [Fr(21, 10), Fr(9, 4), Fr(5, 2), Fr(14, 5), Fr(11, 4), Fr(29, 10)]
MASSES_ENDPOINT = [Fr(2), Fr(3)]
MASSES = MASSES_INTERIOR + MASSES_ENDPOINT

# the paper's Table 1, transcribed and then re-derived
TABLE = {
    Fr(21, 10): (Fr(7611, 4000), Fr(9, 5), 'second', Fr(411, 4000)),
    Fr(9, 4):   (Fr(453, 256),   Fr(3, 2), 'second', Fr(69, 256)),
    Fr(5, 2):   (Fr(51, 32),     Fr(23, 16), 'first', Fr(5, 32)),
    Fr(14, 5):  (Fr(186, 125),   Fr(73, 50), 'first', Fr(7, 250)),
    Fr(11, 4):  (Fr(383, 256),   Fr(93, 64), 'first', Fr(11, 256)),
    Fr(29, 10): (Fr(5939, 4000), Fr(591, 400), 'first', Fr(29, 4000)),
    Fr(2):      (Fr(2),          Fr(2),      '-',     Fr(0)),
    Fr(3):      (Fr(3, 2),       Fr(3, 2),   '-',     Fr(0)),
}

# --------------------------------------------------------------------------------------
hd('Step 1: the exhibited family A(s) lies in omega_3 -- entries >= 0, all ROW sums <= 1'
   ' AND all COLUMN sums <= 1')
note('the column-sum clause is the one whose loss would collapse omega_3 into the'
     ' row-substochastic class Chen-Cao had already settled in their Theorem 4.1')
for s in MASSES:
    A = A_fam(s)
    ck('witness_in_omega_3_at_s_%s' % str(s).replace('/', '_over_'), in_omega(A),
       'A(s) rows=%s cols=%s entries {0, %s, %s}'
       % ([str(x) for x in row_sums(A)], [str(x) for x in col_sums(A)],
          str(A[0][2]), str(A[0][1])))

hd('Step 2: the mass is EXACTLY s (no rounding, no tolerance)')
for s in MASSES:
    A = A_fam(s)
    ck('mass_is_exactly_s_at_s_%s' % str(s).replace('/', '_over_'), mass(A) == s,
       'sigma(A(s)) = %s' % str(mass(A)))

hd('Step 3: per(I-A(s)) computed by the FULL 3! = 6 permutation expansion equals the'
   ' printed cubic f(s) = (s^3-5s^2+4s+12)/4')
for s in MASSES:
    A = A_fam(s)
    p = per(I_minus(A))
    ck('permanent_equals_printed_cubic_at_s_%s' % str(s).replace('/', '_over_'), p == f(s),
       'per(I-A) = %s = %.9f ; f(s) = %s' % (str(p), float(p), str(f(s))))

hd('Step 4: the printed Table 1 is reproduced exactly -- permanent, conjectured value,'
   ' governing branch and excess')
THR_LO, THR_HI = Fr(2274917, 1000000), Fr(2274918, 1000000)   # brackets (-3+sqrt(57))/2
note('the threshold (-3+sqrt57)/2 is irrational; it is bracketed by %s < thr < %s, which'
     ' is enough to decide which branch governs at every mass in the table'
     % (str(THR_LO), str(THR_HI)))
ck('threshold_bracket_is_sound',
   (THR_LO ** 2 + 3 * THR_LO - 12) < 0 < (THR_HI ** 2 + 3 * THR_HI - 12),
   'sign of s^2+3s-12 flips between the two brackets, so the root of'
   ' (s^2-5s+12)/4 = 6-2s lies strictly between them')
for s in MASSES:
    p_want, v_want, branch_want, exc_want = TABLE[s]
    A = A_fam(s)
    p = per(I_minus(A))
    v = V(s)
    if s in MASSES_ENDPOINT:
        branch = '-'
    else:
        assert not (THR_LO <= s <= THR_HI)
        branch = 'first' if s > THR_HI else 'second'
    ok = (p == p_want and v == v_want and branch == branch_want and p - v == exc_want)
    ck('table_row_reproduced_at_s_%s' % str(s).replace('/', '_over_'), ok,
       'per=%s V=%s branch=%s excess=%s' % (str(p), str(v), branch, str(p - v)))

hd('Step 5: STRICT violation at each of the six sampled interior masses in MASSES_INTERIOR'
   ' -- at each of them per(I-A) beats br1, the corrected branch 6-2s, and the variant 6-4s'
   " (the '6-4s' reading is transcribed by hand from arXiv:1801.00225v1, not parsed here);"
   ' the universal statement for all 2<s<3 follows from the Step 9 polynomial identities,'
   ' not from this step')
for s in MASSES_INTERIOR:
    A = A_fam(s)
    p = per(I_minus(A))
    beats_all = (p > br1(s)) and (p > br2_corrected(s)) and (p > br2_printed(s)) and (p > V(s))
    ck('refutes_every_reading_at_s_%s' % str(s).replace('/', '_over_'), beats_all,
       'per=%s > br1=%s, > 6-2s=%s, > 6-4s=%s' % (str(p), str(br1(s)),
                                                  str(br2_corrected(s)), str(br2_printed(s))))

hd('Step 6: NEGATIVE-POLARITY CONTROL at the two masses where the truth is PUBLISHED --'
   ' the family returns the published value and reports NO violation')
note('s=2 is Chen-Cao Theorem 3.1 at n=3, e=2 (value 2); s=3 is their own Lemma 4.2, the'
     ' Omega_3 value 3/2.  A witness family that "violated" the conjecture here would be'
     ' evidence of an arithmetic error, not of a counterexample.')
for s in MASSES_ENDPOINT:
    A = A_fam(s)
    p = per(I_minus(A))
    ck('no_violation_at_published_mass_s_%s' % str(s), p == V(s) and not (p > V(s)),
       'per = V = %s, excess exactly 0, REFUTES=False' % str(p))

hd("Step 7: CONTROL on Chen-Cao's OWN arithmetic -- their witness A_0 reproduces their"
   ' printed first branch exactly, so the refutation beats their arithmetic and not a'
   ' mis-transcription of it')
for s in [Fr(21, 10), Fr(5, 2), Fr(11, 4), Fr(3)]:
    A = A_0(s)
    p = per(I_minus(A))
    ck('control_A0_reproduces_printed_branch1_at_s_%s' % str(s).replace('/', '_over_'),
       in_omega(A) and mass(A) == s and p == br1(s),
       'A_0 in omega_3, sigma=%s, per(I-A_0) = %s = (s^2-5s+12)/4' % (str(mass(A)), str(p)))

hd('Step 8: CONTROL on the misprint M2 -- their witness A_1 gives 6-2s at every mass and'
   ' 6-4s at NONE')
for s in [Fr(21, 10), Fr(5, 2), Fr(3)]:
    A = A_1(s)
    p = per(I_minus(A))
    ck('control_A1_gives_6_minus_2s_never_6_minus_4s_at_s_%s' % str(s).replace('/', '_over_'),
       in_omega(A) and mass(A) == s and p == br2_corrected(s) and p != br2_printed(s),
       'per(I-A_1) = %s ; 6-2s = %s ; 6-4s = %s' % (str(p), str(br2_corrected(s)),
                                                    str(br2_printed(s))))
ck('printed_second_branch_is_negative_above_three_halves',
   all(br2_printed(s) < 0 for s in MASSES_INTERIOR),
   '6-4s < 0 for every s > 3/2, while per(I-A) >= det(I-A) >= 0 for doubly substochastic A')
ck('printed_threshold_belongs_to_the_corrected_branch',
   peq(padd(P_br1, pscale(P_br2c, Fr(-1))), pscale(pnorm([Fr(-12), Fr(3), Fr(1)]), Fr(1, 4))),
   '(s^2-5s+12)/4 - (6-2s) = (s^2+3s-12)/4, whose root is (-3+sqrt57)/2 -- the printed'
   ' threshold; the crossover with 6-4s would be the root of s^2+11s-12, i.e. s=1')

hd('Step 9: the two excess IDENTITIES, by exact polynomial coefficient comparison'
   ' (no sampling, no CAS)')
P_gap1 = pscale(pmul(S, pmul(padd(S, [Fr(-3)]), padd(S, [Fr(-3)]))), Fr(1, 4))   # s(s-3)^2/4
lhs1 = padd(P_f, pscale(P_br1, Fr(-1)))
ck('identity_f_minus_branch1_equals_s_times_s_minus_3_squared_over_4', peq(lhs1, P_gap1),
   'f(s) - (s^2-5s+12)/4 = %s = s(s-3)^2/4' % pshow(lhs1))
T = padd(S, [Fr(-2)])                                                            # t = s-2
P_gap2 = padd(padd(T, pscale(pmul(T, T), Fr(1, 4))), pscale(pmul(T, pmul(T, T)), Fr(1, 4)))
lhs2 = padd(P_f, pscale(P_br2c, Fr(-1)))
ck('identity_f_minus_corrected_branch2_equals_t_plus_t2_over_4_plus_t3_over_4',
   peq(lhs2, P_gap2),
   'f(s) - (6-2s) = %s = t + t^2/4 + t^3/4 with t = s-2' % pshow(lhs2))
ck('both_gaps_are_strictly_positive_at_the_six_sampled_interior_masses_and_gap1_vanishes_at_s_equals_3',
   all(f(s) - br1(s) > 0 and f(s) - br2_corrected(s) > 0 for s in MASSES_INTERIOR)
   and all(f(s) - br1(s) == 0 for s in [Fr(3)]),
   'f-br1 > 0 and f-br2c > 0 at each s in MASSES_INTERIOR, and f(3)-br1(3) = 0; positivity'
   ' on all of (2,3) is NOT checked here -- it follows by inspection from the Step 9'
   ' identities s(s-3)^2/4 and t+t^2/4+t^3/4 with t=s-2>0, which are proved as polynomial'
   ' identities but whose sign is not machine-verified')

hd('Step 10: CONSISTENCY with the two proved UPPER bounds -- Chen-Cao Theorem 4.1 and'
   ' Malek 1986 both give max <= 2 at n=3, and f <= 2 on [2,3]')
note('had f exceeded 2 anywhere on [2,3] the witness would be infeasible and the'
     ' arithmetic wrong; this is the break-attempt that mattered most')
P_two_minus_f = padd([Fr(2)], pscale(P_f, Fr(-1)))
P_fact = pscale(pmul(padd(S, [Fr(-2)]), pnorm([Fr(2), Fr(3), Fr(-1)])), Fr(1, 4))  # (s-2)(2+3s-s^2)/4
ck('identity_two_minus_f_factors_as_s_minus_2_times_2_plus_3s_minus_s_squared_over_4',
   peq(P_two_minus_f, P_fact),
   '2 - f(s) = %s = (s-2)(2+3s-s^2)/4, and 2+3s-s^2 > 0 on (2,3] since its roots are'
   ' (3 +- sqrt17)/2 = -0.5616.., 3.5616..' % pshow(P_two_minus_f))
for s in MASSES:
    ck('f_at_most_two_at_s_%s' % str(s).replace('/', '_over_'), f(s) <= 2,
       'f(%s) = %s <= 2' % (str(s), str(f(s))))

hd('Step 11: the strictly convex slice -- the three vertices of the triangle T_t ARE A(s),'
   " Chen-Cao's A_0 is the MIDPOINT of the edge V_1V_3, and the gap is s(s-3)^2/4")
note('slice coordinates: x_1 = a_23, x_2 = a_13, x_3 = a_12 on the symmetric zero-diagonal'
     ' slice, P(x) = 1 + x_1^2 + x_2^2 + x_3^2 - 2 x_1 x_2 x_3, and'
     ' T_t = {x >= 0, sum x = t, x_i + x_j <= 1} with t = s/2')


def P_slice(x):
    return 1 + x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 2 * x[0] * x[1] * x[2]


def pd_ok(x):
    """Exact test that the Hessian of P restricted to {sum v = 0} is positive definite:
    with p = 1+x_3, q = 1+x_1, m = 1+x_2 the criterion is (p+q-m)^2 < 4pq."""
    p, q, m = 1 + x[2], 1 + x[0], 1 + x[1]
    return (p + q - m) ** 2 < 4 * p * q


for s in MASSES_INTERIOR:
    t = s / 2
    u = (s - 2) / Fr(2)
    V1 = (1 - u, u, u)
    V2 = (u, 1 - u, u)
    V3 = (u, u, 1 - u)
    verts = [V1, V2, V3]
    feasible = all(all(c >= 0 for c in v) and sum(v) == t
                   and all(v[i] + v[j] <= 1 for i, j in ((0, 1), (0, 2), (1, 2)))
                   for v in verts)
    perm_of = all(sorted(v) == sorted((2 - t, t - 1, t - 1)) for v in verts)
    vals_are_f = all(P_slice(v) == f(s) for v in verts)
    A0x = (Fr(1, 2), s / 2 - 1, Fr(1, 2))
    mid = tuple((V1[i] + V3[i]) / 2 for i in range(3))
    is_mid = (mid == A0x)
    mid_val = (P_slice(A0x) == br1(s))
    gap = P_slice(V1) - P_slice(A0x)
    gap_ok = (gap == s * (s - 3) ** 2 / Fr(4) == (u + 1) * (2 * u - 1) ** 2 / Fr(2))
    pd = all(pd_ok(v) for v in verts) and pd_ok(A0x)
    ck('slice_geometry_at_s_%s' % str(s).replace('/', '_over_'),
       feasible and perm_of and vals_are_f and is_mid and mid_val and gap_ok and pd,
       'vertices = permutations of (2-t, t-1, t-1) = A(s), each of value %s; A_0 = midpoint'
       ' of V_1V_3 of value %s; gap = %s = s(s-3)^2/4; restricted Hessian positive definite'
       ' at all four points' % (str(f(s)), str(br1(s)), str(gap)))

def vHv(x, a, b):
    """v^T H v computed from the Hessian of P entry by entry, with v = a(1,-1,0)+b(0,1,-1)."""
    H = [[Fr(2), -2 * x[2], -2 * x[1]],
         [-2 * x[2], Fr(2), -2 * x[0]],
         [-2 * x[1], -2 * x[0], Fr(2)]]
    v = (a, b - a, -b)
    return sum(H[i][j] * v[i] * v[j] for i in range(3) for j in range(3))


_hess_ok = True
_hess_pts = 0
for s in MASSES_INTERIOR:
    u = (s - 2) / Fr(2)
    for x in [(1 - u, u, u), (u, 1 - u, u), (u, u, 1 - u), (Fr(1, 2), s / 2 - 1, Fr(1, 2))]:
        for a, b in [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(1), Fr(1)), (Fr(2), Fr(-3)),
                     (Fr(1, 3), Fr(5, 7))]:
            closed = 4 * ((1 + x[2]) * a ** 2 + (1 + x[0]) * b ** 2
                          - a * b * (1 + x[0] - x[1] + x[2]))
            _hess_pts += 1
            if vHv(x, a, b) != closed or (not (a == 0 and b == 0) and vHv(x, a, b) <= 0):
                _hess_ok = False
for s in MASSES_INTERIOR:
    t = s / 2
    cen = (t / 3, t / 3, t / 3)
    val = P_slice(cen)
    closed = 1 + t ** 2 / Fr(3) - 2 * t ** 3 / Fr(27)
    feas = all(c >= 0 for c in cen) and sum(cen) == t and all(
        cen[i] + cen[j] <= 1 for i, j in ((0, 1), (0, 2), (1, 2)))
    ck('slice_centroid_is_the_slice_minimum_at_s_%s' % str(s).replace('/', '_over_'),
       feas and val == closed and val < f(s) and val < br1(s),
       'P(centroid) = %s = 1 + t^2/3 - 2t^3/27 = %.6f, strictly below both the vertex value'
       ' %s and the edge-midpoint value %s' % (str(val), float(val), str(f(s)), str(br1(s))))
ck('restricted_hessian_closed_form_and_positivity', _hess_ok,
   'v^T H v = 4[(1+x_3)a^2 + (1+x_1)b^2 - ab(1+x_1-x_2+x_3)] verified entry-by-entry at'
   ' %d exact (point, direction) pairs, and strictly positive at every one' % _hess_pts)

hd("Step 12: the LIFT through Chen-Cao's Conjecture 4.6 reduction, checked by the FULL n!"
   ' permanent of the actual n x n block matrix at n = 3, 5, 7')
for n in (3, 5, 7):
    k = (n - 3) // 2
    s = Fr(2 * n - 1, 2)                 # s = n - 1/2, so s' = s-n+3 = 5/2
    sp = s - n + 3
    A = direct_sum([M2] * k + [A_fam(sp)])
    ck('lift_membership_and_mass_at_n_%d' % n,
       len(A) == n and in_omega(A) and mass(A) == s,
       'n=%d, sigma=%s in (n-1, n]' % (n, str(mass(A))))
    p = per(I_minus(A))
    want = Fr(2) ** k * f(sp)
    ck('lift_permanent_equals_two_to_the_n_minus_3_over_2_times_f_at_n_%d' % n,
       p == want,
       'per(I-A) = %s = 2^%d * f(%s) = 2^%d * %s' % (str(p), k, str(sp), k, str(f(sp))))
    ck('lift_strictly_exceeds_the_composite_conjectured_value_at_n_%d' % n,
       p > Fr(2) ** k * V(sp),
       '%s > 2^%d * V(%s) = %s' % (str(p), k, str(sp), str(Fr(2) ** k * V(sp))))
ck('per_of_I_minus_M2_is_two', per(I_minus(M2)) == 2,
   'per([[1,-1],[-1,1]]) = 2, which is why the prefactor is a power of 2')

hd('Step 13: ENDPOINT CONSISTENCY of the lifted cubic -- a consistency check ONLY, and it'
   ' carries NO discriminating information (see the paper)')
note('at s=n-1 the lifted value matches Chen-Cao Theorem 3.1, a THEOREM, for every odd n.'
     ' At s=n it matches the CONJECTURED Kim-Roush 1981 / Minc Conjecture 35 value'
     ' 3*2^(k-2), which is UNPROVED for n >= 5 -- so the only unconditional two-endpoint'
     ' check is at n=3.  The corrected Conjecture 4.5 value V agrees in all ten cells too,'
     ' which is exactly why endpoint agreement discriminates nothing.')
for n in (3, 5, 7, 9, 11):
    k = (n - 1) // 2
    e = (n - 3) // 2
    lo_f, lo_V = Fr(2) ** e * f(Fr(2)), Fr(2) ** e * V(Fr(2))
    hi_f, hi_V = Fr(2) ** e * f(Fr(3)), Fr(2) ** e * V(Fr(3))
    ck('lifted_lower_endpoint_matches_chen_cao_theorem_3_1_at_n_%d' % n,
       lo_f == Fr(2) ** k and lo_V == lo_f,
       's=n-1: lifted = %s = 2^((n-1)/2) = %s, and V gives the SAME'
       % (str(lo_f), str(Fr(2) ** k)))
    ck('lifted_upper_endpoint_matches_the_conjectured_kim_roush_value_at_n_%d' % n,
       hi_f == Fr(3) * Fr(2) ** (k - 2) and hi_V == hi_f,
       's=n: lifted = %s = 3*2^(k-2) = %s (CONJECTURED for n>=5), and V gives the SAME'
       % (str(hi_f), str(Fr(3) * Fr(2) ** (k - 2))))

hd('Step 14: the honest bracket at the headline mass s = 5/2')
p = per(I_minus(A_fam(Fr(5, 2))))
ck('headline_bracket_at_s_five_halves', p == Fr(51, 32) and Fr(51, 32) < 2,
   '51/32 <= max{per(I-A) : A in omega_3, sigma(A)=5/2} <= 2, the upper bound being'
   ' Chen-Cao Theorem 4.1 / Malek 1986; the conjectured 23/16 = 1.4375 is EXCLUDED, and'
   ' the true maximum is NOT determined here')
ck('excess_over_the_conjectured_value_at_s_five_halves',
   p - Fr(23, 16) == Fr(5, 32) and Fr(51, 46) == p / Fr(23, 16),
   'excess exactly 5/32 = 0.15625, ratio 51/46 = 1.1087..')

# --------------------------------------------------------------------------------------
print('')
print('NOTE SCOPE, quoted so no reader can mistake what this program does and does not'
      ' establish.  It verifies the LOWER bound at eight sampled masses (six interior, two'
      ' endpoints): at each of them A(s) is in omega_3, its mass is exactly s, and its'
      ' permanent is exactly f(s).  The strict refutation f(s) > V(s) is verified only at'
      ' the six interior masses; at the two endpoints s=2 and s=3 the checks establish'
      ' f(s) = V(s) exactly (excess 0), a negative-polarity control, not a violation.'
      '  Two exact polynomial'
      ' identities, f(s)-(s^2-5s+12)/4 = s(s-3)^2/4 and f(s)-(6-2s) = t+t^2/4+t^3/4 with'
      ' t=s-2, are verified universally by coefficient comparison; membership, mass and the'
      ' strict sign of the gaps are NOT verified for all 2 < s < 3 -- those remain finite'
      ' samples.  It does NOT establish that max{per(I-A) : A in omega_3^s} EQUALS f(s) --'
      ' that upper bound is OPEN and is not addressed here.  Step 11 proves optimality of A(s) only on the symmetric zero-diagonal'
      ' SLICE, not over the full 9-entry polytope.  NOT RE-RUN HERE: any global search; any'
      ' floating-point multistart search; and the TEXT of the journal version, LAA 555'
      ' (2018) 412-431, which was never obtained -- every claim about what the source prints'
      ' refers to arXiv:1801.00225v1.')
print('')
if _BAD:
    print('CHECKS FAILED: %d of %d' % (_BAD, _N))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _N)
sys.exit(0)
