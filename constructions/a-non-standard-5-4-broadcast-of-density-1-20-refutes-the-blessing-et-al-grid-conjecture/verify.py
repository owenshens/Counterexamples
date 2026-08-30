#!/usr/bin/env python3
"""Verification of every computational claim in

    "A Non-Standard (5,4) Broadcast of the Infinite Grid of Density 1/20"

Python 3.9+, STANDARD LIBRARY ONLY (fractions, math, collections).  No third-party package, no external
data file, no floating point in any decision: every comparison is between integers or Fraction objects.

WHAT IT READS.  The object is not recomputed from a search.  Everything below the header block
`PAPER` is derived from the two lines of data that the paper prints -- the integer basis
[[10,0],[4,2]] and the membership predicate `y % 2 == 0 and (x - 2*y) % 10 == 0` -- and is then
compared against the numbers the paper displays, which are also transcribed literally into `PAPER`.

DEFINITIONS USED, as restated in Section 1 of the paper from Blessing-Insko-Johnson-Mauretour
(arXiv:1401.2499, source lines 335-350) and Drews-Harris-Randolph (arXiv:1712.00150, lines 109/132/138):

    d(u,v)       = graph distance on Z x Z = the L1 distance |x_u - x_v| + |y_u - y_v|
    sig_t(u,v)   = max(0, t - d(u,v))                         (a tower at u gives ITSELF t)
    rec_S(u)     = sum of sig_t(u,v) over all v in S           (the tower at u is included)
    S is (t,r) broadcast dominating   iff   rec_S(u) >= r for EVERY u in Z^2, towers included
    density      = lim_n |S cap [-n,n]^2| / (2n+1)^2 , the box centred at the origin
    T(d,e)       = span{(d,0),(e,1)} = {(x,y) : x = e*y mod d}, the STANDARD broadcasts, density 1/d

Run:  python3 verify.py     (a few seconds)
"""

import math
from collections import Counter
from fractions import Fraction

# ---------------------------------------------------------------------------
# PAPER -- every number this program is allowed to compare against, transcribed
#          from the printed paper and from nowhere else.
# ---------------------------------------------------------------------------
PAPER = {
    'basis':            ((10, 0), (4, 2)),
    'hermite':          (10, 4, 2),                 # (a, b, c) with rows (a,0), (b,c)
    'index':            20,
    't':                5,
    'r':                4,
    'coset_row_y0':     [5, 4, 5, 4, 5, 4, 5, 4, 5, 4],
    'coset_row_y1':     [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    'row_sums':         (45, 40),
    'cone_weight_W5':   85,
    'cone_weight_W4':   44,
    'min_reception':    4,
    'tight_classes':    15,
    'density':          Fraction(1, 20),
    'towers_in_window': 337,                        # |S cap [-40,40]^2|
    'window_histogram': {4: 4880, 5: 1681},
    'min_nonzero_norm': 6,
    'smith':            (2, 10),
    'delta_4_2':        Fraction(1, 18),            # Drews-Harris-Randolph: 1/(2(t-1)^2) at t=4
    'index20_hits':     [(10, 4, 2), (10, 6, 2)],
    'sigma':            {19: 20, 20: 42, 21: 32, 22: 36},
    'D_max_standard':   21,                         # forced by the checksum lemma: 4D <= 85
    'best_standard_d':  18,
    'best_standard_es': [4, 5, 7, 11, 13, 14],
    'C_5_4':            84,                         # Shlomi's usable-signal constant
    'shlomi_bound':     Fraction(1, 21),            # 4/84
    'T18_4_receptions': [5, 5, 6, 5, 5, 4, 5, 4, 4, 4, 4, 4, 5, 4, 5, 5, 6, 5],
    'T18_5_receptions': [4, 3, 2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 3],
    'family':           {2: (5, 20, 18), 3: (13, 52, 50), 4: (25, 100, 98),
                         5: (41, 164, 162), 6: (61, 244, 242)},   # t' -> (D, 4D, 2(2t'-1)^2)
    'hvh_t0':           12,
    'shlomi_table_row': {(5, 4): 18, (7, 4): 50, (9, 4): 98},
    'shlomi_table_ub':  {(5, 4): 20, (7, 4): 52, (9, 4): 100},
}

# ---------------------------------------------------------------------------
# the check harness
# ---------------------------------------------------------------------------
_state = {'n': 0, 'bad': 0}


def PASS(name, detail=''):
    _state['n'] += 1
    print('PASS %-34s %s' % (name, detail))


def FAIL(name, detail=''):
    _state['bad'] += 1
    print('FAIL %-34s %s' % (name, detail))


def check(name, condition, detail=''):
    if condition:
        PASS(name, detail)
    else:
        FAIL(name, detail)


# ---------------------------------------------------------------------------
# the primitives
# ---------------------------------------------------------------------------
def ball(t):
    """The cells that a tower of strength t reaches: L1 norm <= t-1."""
    return [(dx, dy) for dx in range(-(t - 1), t) for dy in range(-(t - 1), t)
            if abs(dx) + abs(dy) <= t - 1]


def reception(in_s, u, t, B=None):
    """rec_S(u) with S given by its membership predicate; the tower at u, if any, contributes t."""
    B = ball(t) if B is None else B
    x0, y0 = u
    tot = 0
    for dx, dy in B:
        if in_s(x0 + dx, y0 + dy):
            tot += t - abs(dx) - abs(dy)
    return tot


def cone_weight(t):
    """W_t = sum over ALL z in Z^2 of max(0, t - |z|_1)."""
    return sum(t - abs(dx) - abs(dy) for dx, dy in ball(t))


def cone_weight_formula(t):
    """The same, counted by shells: 1 cell of norm 0 and 4m cells of norm m."""
    return t + sum(4 * m * (t - m) for m in range(1, t))


def lattice_pred(a, b, c):
    """Membership in the sublattice with rows (a,0) and (b,c), a,c > 0."""
    return lambda x, y: (y % c == 0) and ((x - (y // c) * b) % a == 0)


def standard(d, e):
    """T(d,e) = {(x,y) : x = e*y mod d}."""
    return lambda x, y: (x - e * y) % d == 0


def lattice_min_reception(a, b, c, t):
    """min over the a*c classes of Z^2 / <(a,0),(b,c)>; the reps are (x,y), 0<=x<a, 0<=y<c."""
    B = ball(t)
    p = lattice_pred(a, b, c)
    return min(reception(p, (x, y), t, B) for y in range(c) for x in range(a))


def sublattices(D):
    """Hermite normal form: rows (a,0),(b,c) with a*c = D and 0 <= b < a.  There are sigma(D) of them."""
    out = []
    for c in range(1, D + 1):
        if D % c:
            continue
        a = D // c
        for b in range(a):
            out.append((a, b, c))
    return out


def sigma(D):
    return sum(k for k in range(1, D + 1) if D % k == 0)


def smith_2x2(m):
    """Invariant factors (d1, d2) of a 2x2 integer matrix: d1 = gcd of the entries, d1*d2 = |det|."""
    (p, q), (r, s) = m
    d1 = math.gcd(math.gcd(abs(p), abs(q)), math.gcd(abs(r), abs(s)))
    det = abs(p * s - q * r)
    return d1, det // d1


def count_congruent(a, n):
    """#{x in [-n,n] : x = a mod 10}, in closed form so no O(n^2) loop is needed."""
    return (n - a) // 10 - (-n - 1 - a) // 10


# ---------------------------------------------------------------------------
# THE WITNESS, read off the paper
# ---------------------------------------------------------------------------
A, Bh, C = PAPER['hermite']
T, R = PAPER['t'], PAPER['r']
B5 = ball(T)


def in_S(x, y):
    """The paper's printed predicate, verbatim: y even and x = 2y mod 10."""
    return (y % 2 == 0) and ((x - 2 * y) % 10 == 0)


def in_S_span(x, y):
    """The paper's printed basis, as an integer span: (x,y) = m*(10,0) + n*(4,2)."""
    (a1, a2), (b1, b2) = PAPER['basis']
    if (y - 0) % b2 != 0:
        return False
    n = y // b2
    return (x - n * b1) % a1 == 0


print('=' * 100)
print('VERIFICATION OF THE PRINTED OBJECT   S = {(x,y) in Z^2 : y = 0 mod 2, x = 2y mod 10}'
      ' = span{(10,0),(4,2)}')
print('=' * 100)

# --- 1. the object is the lattice the paper says it is ----------------------
W = 40
pred_pts = set((x, y) for x in range(-W, W + 1) for y in range(-W, W + 1) if in_S(x, y))
span_pts = set((x, y) for x in range(-W, W + 1) for y in range(-W, W + 1) if in_S_span(x, y))
gen_pts = set()
for m in range(-30, 31):
    for n in range(-30, 31):
        x, y = 10 * m + 4 * n, 2 * n
        if -W <= x <= W and -W <= y <= W:
            gen_pts.add((x, y))
check('predicate-equals-basis-span', pred_pts == span_pts == gen_pts,
      'three independent constructions of S agree over [-40,40]^2, %d points each' % len(pred_pts))
check('tower-count-in-window', len(pred_pts) == PAPER['towers_in_window'],
      '|S cap [-40,40]^2| = %d, paper prints %d' % (len(pred_pts), PAPER['towers_in_window']))

det = PAPER['basis'][0][0] * PAPER['basis'][1][1] - PAPER['basis'][0][1] * PAPER['basis'][1][0]
check('index-is-the-determinant', det == PAPER['index'],
      'det [[10,0],[4,2]] = %d = index of S in Z^2' % det)

# --- 2. the 20 printed coset representatives really are a full set ----------
reps = [(x, y) for y in range(C) for x in range(A)]
check('reps-count', len(reps) == PAPER['index'], '%d representatives (x in 0..9, y in 0..1)' % len(reps))
pairwise = all(not in_S(p[0] - q[0], p[1] - q[1]) for i, p in enumerate(reps) for q in reps[i + 1:])
covering = True
for x in range(-25, 26):
    for y in range(-25, 26):
        hits = [u for u in reps if in_S(x - u[0], y - u[1])]
        if len(hits) != 1:
            covering = False
check('reps-are-a-transversal', pairwise and covering,
      'the 20 reps are pairwise inequivalent mod S, and every cell of [-25,25]^2 is congruent to '
      'exactly one of them')

# --- 3. reception is S-periodic, so the 20 classes decide everything --------
periodic = all(reception(in_S, u, T, B5)
               == reception(in_S, (u[0] + 10, u[1]), T, B5)
               == reception(in_S, (u[0] + 4, u[1] + 2), T, B5) for u in reps)
check('reception-is-S-periodic', periodic,
      'rec(u) = rec(u+(10,0)) = rec(u+(4,2)) at all 20 representatives')

# --- 4. the printed coset table --------------------------------------------
row0 = [reception(in_S, (x, 0), T, B5) for x in range(A)]
row1 = [reception(in_S, (x, 1), T, B5) for x in range(A)]
check('coset-table-row-y0', row0 == PAPER['coset_row_y0'], 'recomputed %s' % row0)
check('coset-table-row-y1', row1 == PAPER['coset_row_y1'], 'recomputed %s' % row1)
check('coset-table-row-sums', (sum(row0), sum(row1)) == PAPER['row_sums'],
      'row sums %d and %d' % (sum(row0), sum(row1)))

# --- 5. the checksum lemma ------------------------------------------------
w5a, w5b = cone_weight(T), cone_weight_formula(T)
check('cone-weight-W5-two-ways', w5a == w5b == PAPER['cone_weight_W5'],
      'W_5 = sum_z max(0,5-|z|_1) = %d by direct summation and = 5 + sum_m 4m(5-m) = %d' % (w5a, w5b))
total = sum(row0) + sum(row1)
check('checksum-lemma-instance', total == w5a,
      'the 20 coset receptions sum to %d = W_5; one slip anywhere in the table would break this' % total)

# --- 6. the decisive inequality on the table ------------------------------
allvals = row0 + row1
check('min-reception-equals-r', min(allvals) == PAPER['min_reception'] == R,
      'min over the 20 classes = %d = r, so S IS a (5,4) broadcast dominating set of Z^2' % min(allvals))
check('tight-class-count', allvals.count(R) == PAPER['tight_classes'],
      '%d of the 20 classes receive exactly 4; the other %d receive 5'
      % (allvals.count(R), 20 - allvals.count(R)))

# --- 7. an independent recheck over a window, no lattice arithmetic -------
hist = Counter()
wmin = None
for x in range(-W, W + 1):
    for y in range(-W, W + 1):
        v = reception(in_S, (x, y), T, B5)
        hist[v] += 1
        wmin = v if wmin is None else min(wmin, v)
check('window-recheck-81x81', wmin == R and dict(hist) == PAPER['window_histogram'],
      'min = %d, histogram = %s over the 6561 cells of [-40,40]^2' % (wmin, dict(sorted(hist.items()))))
even_even = sum(1 for x in range(-W, W + 1) for y in range(-W, W + 1) if x % 2 == 0 and y % 2 == 0)
check('window-histogram-accounting', sum(hist.values()) == 6561 and hist[5] == even_even == 41 * 41,
      '6561 = %d + %d, and the %d cells receiving 5 are exactly the 41^2 cells with both coordinates even'
      % (hist[4], hist[5], hist[5]))

# --- 8. the density, exactly ---------------------------------------------
brute = sum(1 for x in range(-50, 51) if (x - 0) % 10 == 0)
check('closed-form-count-agrees', brute == count_congruent(0, 50),
      'the O(1) congruence count is validated against a brute-force count at n = 50 (%d)' % brute)
dens_ok = True
report = []
for n in (50, 200, 500, 2000, 10 ** 6):
    cnt = sum(count_congruent((2 * y) % 10, n) for y in range(-n, n + 1) if y % 2 == 0)
    frac = Fraction(cnt, (2 * n + 1) ** 2)
    if abs(frac - PAPER['density']) > Fraction(1, n):
        dens_ok = False
    report.append('n=%d:%s' % (n, frac.limit_denominator(10 ** 7)))
check('density-is-one-twentieth', dens_ok,
      'exact rational counts |S cap [-n,n]^2|/(2n+1)^2 within 1/n of 1/20: ' + ', '.join(report))

# --- 9. structure: minimum distance, doubling, non-standardness ----------
norms = [abs(x) + abs(y) for (x, y) in pred_pts if (x, y) != (0, 0)]
check('min-nonzero-L1-norm', min(norms) == PAPER['min_nonzero_norm'],
      'the shortest nonzero vector of S has L1 norm %d > 4, so a tower receives only its own 5' % min(norms))

lee1 = standard(5, 2)
perfect = all(sum(1 for dx in range(-1, 2) for dy in range(-1, 2)
                  if abs(dx) + abs(dy) <= 1 and lee1(X + dx, Y + dy)) == 1
              for X in range(-30, 31) for Y in range(-30, 31))
doubled = set((2 * X, 2 * Y) for X in range(-30, 31) for Y in range(-30, 31)
              if lee1(X, Y) and -W <= 2 * X <= W and -W <= 2 * Y <= W)
check('S-is-twice-the-Lee-code', perfect and doubled == pred_pts,
      'T(5,2) is a perfect Lee code of radius 1 (every cell within distance 1 of exactly one codeword) '
      'and 2*T(5,2) = S over [-40,40]^2')

on_line_y1 = sum(1 for x in range(-10 ** 4, 10 ** 4 + 1) if in_S(x, 1))
std_meets = all(any(standard(d, e)(x, 1) for x in range(d)) for d in range(1, 26) for e in range(d))
check('S-is-not-standard-by-a-line', on_line_y1 == 0 and std_meets,
      'S has no tower on the line y = 1, while every T(d,e), d <= 25, meets every horizontal line')

sm = smith_2x2(PAPER['basis'])
std_smith = set(smith_2x2(((d, 0), (e, 1))) for d in range(1, 26) for e in range(d))
check('quotient-is-not-cyclic', sm == PAPER['smith'] and all(f[0] == 1 for f in std_smith),
      'Z^2/S has invariant factors (%d,%d), i.e. Z/2 + Z/10, not cyclic; every Z^2/T(d,e) has invariant '
      'factors (1,d) and is cyclic. Non-cyclicity is GL_2(Z)-invariant, so S is no relabelling of a T(d,e).'
      % sm)

# --- 10. the refutation --------------------------------------------------
delta42 = Fraction(1, 2 * (4 - 1) ** 2)
check('delta-4-2-from-DHR-formula', delta42 == PAPER['delta_4_2'],
      'the Drews-Harris-Randolph value 1/(2(t-1)^2) at t = 4 is %s' % delta42)
check('strict-inequality-1-20-lt-1-18', PAPER['density'] < delta42,
      'delta_{5,4} <= 1/20 < 1/18 = delta_{4,2} as exact rationals, so the equality '
      'delta_{t,r} = delta_{t+1,r+2} FAILS at (t,r) = (4,2)')

# --- 11. the index censuses ---------------------------------------------
for D in (19, 20, 21, 22):
    lats = sublattices(D)
    hits = [(a, b, c) for (a, b, c) in lats if lattice_min_reception(a, b, c, T) >= R]
    if D == 20:
        check('index-20-census', len(lats) == sigma(D) == PAPER['sigma'][D]
              and hits == PAPER['index20_hits'],
              'all %d sublattices of index 20 examined (= sigma(20)); exactly two are (5,4) broadcasts: '
              '%s -- our witness and its mirror image' % (len(lats), hits))
    else:
        check('index-%d-census-empty' % D, len(lats) == sigma(D) == PAPER['sigma'][D] and hits == [],
              'all %d sublattices of index %d examined (= sigma(%d)); none is a (5,4) broadcast'
              % (len(lats), D, D))

# --- 12. the elementary density bound the checksum lemma gives -----------
check('checksum-bound-D-at-most-21', 4 * 21 <= w5a < 4 * 22 and 4 * PAPER['D_max_standard'] <= w5a,
      'a sublattice with one tower per period and rec >= 4 everywhere obeys 4D <= W_5 = 85, hence '
      'D <= 21; this is what makes the standard-broadcast search below finite')

capped = sum(min(T - abs(dx) - abs(dy), R) for dx, dy in B5)
shlomi_C = Fraction(2, 3) * R ** 3 - 2 * R ** 2 * T + 2 * R * T ** 2 + Fraction(R, 3)
check('shlomi-constant-C-5-4', capped == shlomi_C == PAPER['C_5_4'],
      'sum over the ball of min(sig,r) = 4+16+24+24+16 = %d, and Shlomi\'s closed form '
      '2r^3/3 - 2r^2 t + 2rt^2 + r/3 = %s -- equal' % (capped, shlomi_C))
check('delta-5-4-not-determined',
      PAPER['shlomi_bound'] == Fraction(R, PAPER['C_5_4']) < PAPER['density'],
      'Shlomi\'s bound gives delta_{5,4} >= 4/84 = 1/21, and our witness gives <= 1/20; '
      '1/21 < 1/20 strictly, so the exact value of delta_{5,4} is NOT settled here')

# --- 13. controls, both polarities --------------------------------------
t184 = [reception(standard(18, 4), (x, 0), T, B5) for x in range(18)]
check('control-positive-T(18,4)', min(t184) >= R and t184 == PAPER['T18_4_receptions']
      and sum(t184) == w5a,
      'T(18,4) IS a (5,4) broadcast: receptions %s, min %d, sum %d = W_5. So 1/18 is achieved at (5,4) '
      'and the comparison 1/20 < 1/18 is against a realisable density.' % (t184, min(t184), sum(t184)))

best_d, best_es = 0, []
for d in range(1, 26):
    es = [e for e in range(d) if min(reception(standard(d, e), (x, 0), T, B5) for x in range(d)) >= R]
    if es:
        best_d, best_es = d, es
check('control-best-standard-is-18', best_d == PAPER['best_standard_d'] and best_es == PAPER['best_standard_es'],
      'over d = 1..25 the largest d admitting a standard (5,4) broadcast is %d, at e = %s -- reproducing '
      'the published T(18,4) including its lowest-e tie-break' % (best_d, best_es))

bad = [(d, e) for d in range(19, 23) for e in range(d)
       if min(reception(standard(d, e), (x, 0), T, B5) for x in range(d)) >= R]
check('control-negative-d-19-to-22', bad == [],
      'no T(d,e) with d in 19..22 is a (5,4) broadcast (must be empty, and is: %s)' % bad)

t185 = [reception(standard(18, 5), (x, 0), 4, ball(4)) for x in range(18)]
check('control-third-party-T(18,5)', t185 == PAPER['T18_5_receptions'] and min(t185) >= 2
      and sum(t185) == cone_weight(4) == PAPER['cone_weight_W4'],
      'our (4,2) reception vector for T(18,5) is %s, identical entry by entry to the published one; '
      'sum %d = W_4. This pins our orientation, offset and weight conventions against a third party.'
      % (t185, sum(t185)))

nonstd23 = lambda x, y: ((x + 2 * y) % 4) in (0, 1)
v23 = [reception(nonstd23, (x, y), 2, ball(2)) for y in range(1) for x in range(4)]
check('control-DHR-non-standard-(2,3)', min(v23) >= 3 and v23 == [3, 3, 3, 3],
      'the known non-standard optimum at (t,r) = (2,3), {(x,y): x+2y mod 4 in {0,1}} of density 1/2, '
      'validates with receptions %s -- so this checker can recognise a non-standard optimum' % v23)

# --- 14. by-product: no standard (5,4) broadcast is optimal --------------
check('no-standard-5-4-is-optimal',
      Fraction(1, best_d) > PAPER['density'] and best_d <= PAPER['D_max_standard'],
      'the least density of a standard (5,4) broadcast is 1/%d, and 1/%d > 1/20; the search is complete '
      'because 4D <= 85 forces D <= 21. So (5,4) is a pair for which no standard broadcast is optimal, '
      'answering the Drews-Harris-Randolph Further Question affirmatively.' % (best_d, best_d))

# --- 15. by-product: the corpus correction ------------------------------
tab_ok = all(Fraction(1, PAPER['shlomi_table_row'][k]) > Fraction(1, PAPER['shlomi_table_ub'][k])
             for k in PAPER['shlomi_table_row'])
check('table-entries-exceed-our-bounds', tab_ok,
      'the printed reciprocal-density entries 18, 50, 98 at (5,4),(7,4),(9,4) correspond to densities '
      '1/18, 1/50, 1/98, each STRICTLY GREATER than our 1/20, 1/52, 1/100 -- so as captioned '
      '("the reciprocal of the minimum density") those three entries are wrong')

# --- 16. the doubling family, t' = 2..6 --------------------------------
hvh = []
for tp in sorted(PAPER['family']):
    rho = tp - 1
    D = 2 * rho * rho + 2 * rho + 1
    Dref, idx_ref, cmp_ref = PAPER['family'][tp]
    lee = standard(D, (-(2 * rho + 1)) % D)
    perf = all(sum(1 for dx in range(-rho, rho + 1) for dy in range(-rho, rho + 1)
                   if abs(dx) + abs(dy) <= rho and lee(X + dx, Y + dy)) == 1
               for X in range(-2 * D, 2 * D + 1) for Y in range(-rho - 2, rho + 3))
    tt = 2 * tp + 1
    Btt = ball(tt)
    fam = (lambda D=D, lee=lee: (lambda x, y: x % 2 == 0 and y % 2 == 0 and lee(x // 2, y // 2)))()
    vals = [reception(fam, (x, y), tt, Btt) for y in range(2) for x in range(2 * D)]
    Wtt = cone_weight_formula(tt)
    ours = Fraction(1, 4 * D)
    theirs = Fraction(1, 2 * (2 * tp - 1) ** 2)
    good = (perf and D == Dref and 4 * D == idx_ref and 2 * (2 * tp - 1) ** 2 == cmp_ref
            and len(vals) == 4 * D and min(vals) == R and sum(vals) == Wtt and ours < theirs)
    if good:
        hvh.append(2 * tp)
    check('family-t-prime-%d' % tp, good,
          "L is a perfect Lee code of radius %d (index %d); S = 2L has index %d, min reception %d over all "
          "%d classes, checksum %d = W_%d; density 1/%d < 1/%d = delta_{%d,2}. So the conjecture also fails "
          "at (t,r) = (%d,2)." % (rho, D, 4 * D, min(vals), len(vals), sum(vals), tt, 4 * D,
                                  2 * (2 * tp - 1) ** 2, 2 * tp, 2 * tp))

check('herrman-van-hintum-range', hvh == [4, 6, 8, 10, 12] and max(hvh) == PAPER['hvh_t0'],
      'counterexamples to delta_{t+1,4} = delta_{t,2} are exhibited at t = %s, so the Herrman-van Hintum '
      'weakening "there exists t_0 such that for all t >= t_0 ..." is refuted at r = 2 for EVERY t_0 <= %d '
      'and for no larger t_0' % (hvh, max(hvh)))

# ---------------------------------------------------------------------------
# scope: what this program does NOT cover
# ---------------------------------------------------------------------------
print('-' * 100)
print('NOT RE-RUN: the lower bound delta_{t,2} >= 1/(2(t-1)^2) of Drews-Harris-Randolph is QUOTED, not')
print('            reproved. Only its value at t = 4 is arithmetic-checked. The strict inequality')
print('            1/20 < delta_{4,2} rests entirely on that published theorem.')
print('NOT RE-RUN: Shlomi\'s lower bound delta_{t,r} >= r / C_{t,r} is QUOTED, not reproved; only the')
print('            arithmetic C_{5,4} = 84 and 4/84 = 1/21 is checked here.')
print('NOT RE-RUN: no search above index 22. Nothing in this program excludes a (5,4) broadcast of')
print('            density below 1/20, so delta_{5,4} is bounded, not determined: 1/21 <= d <= 1/20.')
print('NOT RE-RUN: the censuses at index 19..22 enumerate SUBLATTICES ONLY (one tower per period).')
print('            Periodic sets with two or more towers per period are not enumerated at any index.')
print('NOT RE-RUN: the doubling family is enumerated only for 2 <= t\' <= 6. The general lemma')
print('            "S_{t\'} = 2 L_{t\'} is a (2t\'+1,4) broadcast for every t\' >= 2" is NOT proved,')
print('            so a Herrman-van Hintum threshold t_0 >= 13 is untouched.')
print('NOT RE-RUN: no literature search. In particular Harris-Insko-Johnson, "Projects in (t,r) Broadcast')
print('            Domination" (Springer 2020, doi 10.1007/978-3-030-37853-0_8), is paywalled and unread.')
print('-' * 100)

if _state['bad']:
    print('VERDICT: %d CHECK(S) FAILED' % _state['bad'])
    raise SystemExit(1)
print('VERDICT: ALL %d CHECKS PASS' % _state['n'])
raise SystemExit(0)
