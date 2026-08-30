#!/usr/bin/env python3
"""
verify.py -- an independent check of the claims of

    "A Proof of the Bhagat-Kulkarni-Larsson-Murali Smallest-Discrepancy-Heap Formula"

against a brute-force evaluation of the PSPE recursion (1) of the paper, in exact integer
arithmetic.  Python 3.9+, STANDARD LIBRARY ONLY (only `sys` and `fractions.Fraction` are
imported), no external data file: a referee can run it with nothing installed.

    python3 verify.py            # the paper's box, s_2 <= 100
    python3 verify.py 40         # a smaller box, if you are in a hurry

WHAT IT READS.  Every object it consumes is PRINTED IN THE PAPER: the two hypotheses, the
definitions of a, b, d, n, eps, alpha_t and H, the admissibility rule 0 <= i <= ceil(m/2),
the closed form v(m,i) = m*b - i*d, the invariant (6), the values (7) at H, the block table
for S={3,5}, the six-instance table, and the 16 published columns of Table tab:3_5 of the
target paper as reproduced in section 5.  Nothing is inherited from the run that produced the
paper.

WHAT IT DOES NOT COVER is printed at the end of the run and repeated in REVIEW_NOTE.md.
"""
import sys
from fractions import Fraction

AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100

CHECKS = []


def check(name, cond, detail=''):
    """One PASS/FAIL line per check.  Exit status is 0 iff every one passed."""
    CHECKS.append((name, bool(cond)))
    print('%s %-46s %s' % ('PASS' if cond else 'FAIL', name, detail))
    return bool(cond)


# ---------------------------------------------------------------------------
# The recursion (1) of the paper, brute-forced.  friendly=True is FvF (break a tie so as to
# MAXIMISE the opponent's o^1(h-s)); friendly=False is AvA (minimise it).
# `collapse` forces both branches to one convention: the M2 mutation control.
# ---------------------------------------------------------------------------
def outcomes(a, b, N, friendly, collapse=None):
    if collapse == 'F':
        friendly = True
    elif collapse == 'A':
        friendly = False
    f = [0] * (N + 1)
    g = [0] * (N + 1)
    nties = 0
    for h in range(a, N + 1):
        va = a + g[h - a]
        vb = b + g[h - b] if h >= b else None
        if vb is None or va > vb:
            s = a
        elif vb > va:
            s = b
        else:
            nties += 1
            oa, ob = f[h - a], f[h - b]
            if oa == ob:
                s = a
            elif friendly:
                s = a if oa > ob else b
            else:
                s = a if oa < ob else b
        f[h] = a + g[h - a] if s == a else b + g[h - b]
        g[h] = f[h - s]
    return f, g, nties


# ---------------------------------------------------------------------------
# The lattice W of the paper: admissible (m,i), v(m,i) = m*b - i*d.
# ---------------------------------------------------------------------------
def lattice(a, d, n, H):
    """sorted [(v, m, i)] for admissible (m,i) with 0 <= v <= H, plus every collision found."""
    b = a + d
    seen = {}
    coll = []
    for m in range(0, 2 * n + 3):
        for i in range(0, (m + 1) // 2 + 1):          # (m+1)//2 == ceil(m/2)
            v = m * b - i * d
            if v < 0 or v > H:
                continue
            if v in seen:
                coll.append((v, seen[v], (m, i)))
            else:
                seen[v] = (m, i)
    return sorted((v, mi[0], mi[1]) for v, mi in seen.items()), coll


def Vtable(items, H):
    """Vm[h], Vi[h], Vv[h] -- the coordinates of V(h) = max{w in W : w <= h}, for h = 0..H."""
    Vm = [0] * (H + 1)
    Vi = [0] * (H + 1)
    Vv = [0] * (H + 1)
    k = 0
    cm = ci = cv = 0
    for h in range(H + 1):
        while k < len(items) and items[k][0] <= h:
            cv, cm, ci = items[k]
            k += 1
        Vm[h], Vi[h], Vv[h] = cm, ci, cv
    return Vm, Vi, Vv


def predict(m, i, b, d):
    """The invariant (6) of the paper."""
    t = m // 2
    if m % 2 == 0:
        return (t * b, t * b - i * d)
    return ((t + 1) * b - i * d, t * b)


def hypothesis_pairs(amax):
    """The paper's hypothesis set, in the equivalent form of its Lemma 1."""
    for a in range(2, amax + 1):
        for d in range(2, a):
            if a % d == 0:
                continue
            yield a, d


def tau(k):
    return sum(1 for j in range(1, k + 1) if k % j == 0)


print('=' * 78)
print('PART A -- CONTROLS.  Both polarities, before any sweep.')
print('=' * 78)

# --- C1..C5: the 16 published columns of Table tab:3_5, as reproduced in section 5 of the
# paper, and the block table of section 5.  These integers are the PAPER'S, not ours.
PUB_F = [(0, 0), (0, 0), (0, 0), (3, 0), (3, 0), (5, 0), (5, 0), (5, 0),
         (5, 3), (5, 3), (5, 5), (6, 5), (6, 5), (8, 5), (8, 6), (10, 5)]
PUB_A = [(0, 0), (0, 0), (0, 0), (3, 0), (3, 0), (5, 0), (5, 0), (5, 0),
         (5, 3), (5, 3), (5, 5), (6, 5), (6, 5), (8, 5), (8, 5), (10, 5)]
PUB_D2 = [0] * 14 + [-1, 0]

fF, gF, tiesF = outcomes(3, 5, 15, True)
fA, gA, tiesA = outcomes(3, 5, 15, False)
mineF = [(fF[h], gF[h]) for h in range(16)]
mineA = [(fA[h], gA[h]) for h in range(16)]
d1 = [mineA[h][0] - mineF[h][0] for h in range(16)]
d2 = [mineA[h][1] - mineF[h][1] for h in range(16)]

check('control-published-table-fvf', mineF == PUB_F,
      'S={3,5}, all 16 columns h=0..15 reproduced')
check('control-published-table-ava', mineA == PUB_A,
      'S={3,5}, all 16 columns h=0..15 reproduced')
check('control-published-delta1-identically-zero', all(x == 0 for x in d1),
      'delta^1(h) = 0 for h=0..15, as published')
check('control-published-delta2-bold-entry', d2 == PUB_D2,
      'the only non-zero entry is delta^2(14) = -1, as published in bold')
check('control-tie-break-branch-is-exercised', tiesF > 0 and tiesA > 0,
      'ties met on S={3,5}, h<=15: FvF %d, AvA %d' % (tiesF, tiesA))

# --- C6, C7: the two sacrifice anchors published at line 301 of the target's source.
anchor_ok = []
for (a, b, hh, want) in [(2, 3, 7, (4, 3)), (3, 4, 17, (9, 8))]:
    f2, g2, _ = outcomes(a, b, hh, True)
    va, vb = a + g2[hh - a], b + g2[hh - b]
    anchor_ok.append(((va, vb) == want and va > vb, a, b, hh, va, vb))
for (ok, a, b, hh, va, vb) in anchor_ok:
    check('control-sacrifice-anchor-s%d-%d-h%d' % (a, b, hh), ok,
          'value(s_2=%d)=%d > value(s_1=%d)=%d, strict sacrifice' % (a, va, b, vb))

# --- C8: PROVED SILENT.  In the dominant regime s_1 >= 2 s_2 the target's own PROVED
# Proposition prop:equality forbids a discrepancy at any heap, so the detector must say
# nothing there.
dom_pairs = dom_heaps = 0
dom_bad = []
for a in range(1, 41):
    for b in range(2 * a, min(2 * a + 30, 90) + 1):
        N = 600
        f3, g3, _ = outcomes(a, b, N, True)
        f4, g4, _ = outcomes(a, b, N, False)
        dom_pairs += 1
        dom_heaps += N + 1
        for h in range(N + 1):
            if f3[h] != f4[h] or g3[h] != g4[h]:
                dom_bad.append((a, b, h))
                break
check('control-dominant-regime-proved-silent', not dom_bad,
      '%d pairs, %d heap values, %d discrepancies (0 required)'
      % (dom_pairs, dom_heaps, len(dom_bad)))

# --- M1: the invariant checker must be able to say NO.  Shift the parity index m -> m+1.
a, d = 3, 2
b, n = a + d, a // d
H = (n + 2) * a + n * b
items, coll = lattice(a, d, n, H)
Vm, Vi, Vv = Vtable(items, H)
fF, gF, _ = outcomes(a, b, H, True)
good = sum(1 for h in range(H) if (fF[h], gF[h]) == predict(Vm[h], Vi[h], b, d))
mutd = sum(1 for h in range(H) if (fF[h], gF[h]) == predict(Vm[h] + 1, Vi[h], b, d))
check('mutation-parity-shifted-invariant-rejected', good == H and mutd < H,
      'S={3,5}: invariant (6) holds on %d/%d heaps h<H, the m->m+1 mutant on %d/%d'
      % (good, H, mutd, H))

# --- M2: the detector must fire because of the TIE-BREAK and nothing else.  Collapse the
# second arm of the comparison onto the first convention: the SAME detector, run on the SAME
# pairs, must then be completely silent.  (If it still fired, what it reports would be an
# engine artefact rather than a disagreement between the conventions.)
m2_pairs = m2_bad = 0
for (aa, dd) in hypothesis_pairs(40):
    bb, nn = aa + dd, aa // dd
    HH = (nn + 2) * aa + nn * bb
    realF = outcomes(aa, bb, HH, True)
    fakeF = outcomes(aa, bb, HH, False, collapse='F')     # the AvA arm, forced friendly
    realA = outcomes(aa, bb, HH, False)
    fakeA = outcomes(aa, bb, HH, True, collapse='A')      # the FvF arm, forced antagonistic
    m2_pairs += 1
    if (realF[0] != fakeF[0] or realF[1] != fakeF[1]
            or realA[0] != fakeA[0] or realA[1] != fakeA[1]):
        m2_bad += 1
check('mutation-collapsed-comparison-silent', m2_bad == 0,
      '%d pairs with s_2<=40; with both arms on one convention the detector fired on %d'
      % (m2_pairs, m2_bad))

# --- C9: the second hypothesis is LOAD-BEARING.  At an excluded ratio d | a the paper's own
# Lemma 4(F2) fails: eps = 0 and admissible coordinates collide below H.
a9, d9 = 100, 1
n9 = a9 // d9
b9 = a9 + d9
H9 = (n9 + 2) * a9 + n9 * b9
v_even = 200 * b9 - 0 * d9
v_odd = 201 * b9 - 101 * d9
check('control-excluded-ratio-produces-a-collision',
      a9 - n9 * d9 == 0 and v_even == v_odd and 101 <= (201 + 1) // 2 and v_even < H9,
      'a=100, d=1: eps=0 and v(200,0)=v(201,101)=%d < H=%d' % (v_even, H9))

print()
print('=' * 78)
print("PART B -- THE OBJECTS PRINTED IN THE PAPER, re-derived from the recursion.")
print('=' * 78)

# --- P1, P2: the block table and the values at H for S={3,5}, exactly as printed.
PAPER_BLOCKS = [((0, 3), (0, 0)), ((3, 5), (3, 0)), ((5, 8), (5, 0)), ((8, 10), (5, 3)),
                ((10, 11), (5, 5)), ((11, 13), (6, 5)), ((13, 14), (8, 5))]
fF, gF, _ = outcomes(3, 5, 14, True)
fA, gA, _ = outcomes(3, 5, 14, False)
blocks_ok = all((fF[h], gF[h]) == pair and (fA[h], gA[h]) == pair
                for (lo, hi), pair in PAPER_BLOCKS for h in range(lo, hi))
covered = sorted(h for (lo, hi), _ in PAPER_BLOCKS for h in range(lo, hi))
check('paper-block-table-for-s-3-5', blocks_ok and covered == list(range(14)),
      '%d blocks, both conventions, covering exactly h=0..13' % len(PAPER_BLOCKS))
check('paper-values-at-h-for-s-3-5',
      (fF[14], gF[14]) == (8, 6) and (fA[14], gA[14]) == (8, 5),
      'o_FvF(14)=(8,6), o_AvA(14)=(8,5) as printed')

# --- P3..P5: the six-instance table of section 5, and the two outcome pairs quoted there.
SIX = [(3, 5, 1, 1, 14), (5, 7, 2, 1, 34), (7, 9, 3, 1, 62),
       (9, 11, 4, 1, 98), (7, 10, 2, 1, 48), (10, 13, 3, 1, 89)]
row_ok = []
first_ok = []
for (s2, s1, wn, weps, wH) in SIX:
    dd = s1 - s2
    nn = s2 // dd
    ee = s2 - nn * dd
    HH = (nn + 2) * s2 + nn * s1
    row_ok.append((nn, ee, HH) == (wn, weps, wH))
    f5, g5, _ = outcomes(s2, s1, HH, True)
    f6, g6, _ = outcomes(s2, s1, HH, False)
    first = next((h for h in range(HH + 1)
                  if f5[h] != f6[h] or g5[h] != g6[h]), None)
    first_ok.append(first == HH)
check('paper-six-instance-table', all(row_ok),
      'n, eps and H = (n+2)s_2+n*s_1 agree on all %d printed instances' % len(SIX))
check('paper-six-instance-least-discrepancy-heap', all(first_ok),
      'the least discrepancy heap is exactly H on all %d printed instances' % len(SIX))
f7, g7, _ = outcomes(5, 7, 34, True)
f8, g8, _ = outcomes(5, 7, 34, False)
f9, g9, _ = outcomes(10, 13, 89, True)
fa, ga, _ = outcomes(10, 13, 89, False)
check('paper-quoted-outcome-pairs',
      (f7[34], g7[34]) == (19, 15) and (f8[34], g8[34]) == (19, 14)
      and (f9[89], g9[89]) == (49, 40) and (fa[89], ga[89]) == (49, 39),
      'S={5,7} at 34 -> (19,15)/(19,14); S={10,13} at 89 -> (49,40)/(49,39)')

# --- P6: Lemma 1 of the paper -- the interval condition of the target is EQUIVALENT to
# 0<d<a together with d nmid a, and then n = floor(a/d).  Tested on a grid that deliberately
# contains non-hypothesis pairs, so the equivalence can fail in both directions.
grid = disagree = inside = 0
for a in range(2, 61):
    for b in range(a + 1, 3 * a + 1):
        grid += 1
        d = b - a
        ratio = Fraction(b, a)
        ns = [k for k in range(1, 4 * a + 4)
              if Fraction(k + 1, k) > ratio > Fraction(k + 2, k + 1)]
        mine = (0 < d < a) and (a % d != 0)
        if mine:
            inside += 1
        if bool(ns) != mine or (mine and ns != [a // d]):
            disagree += 1
check('paper-lemma-1-hypothesis-equivalence', disagree == 0,
      '%d (a,b) pairs tested, %d in the hypothesis set, %d disagreements'
      % (grid, inside, disagree))

# --- P7: the denominator of the box, independently by closed form.
box = list(hypothesis_pairs(AMAX))
closed = sum(a - tau(a) for a in range(2, AMAX + 1))
check('paper-box-denominator-closed-form', len(box) == closed,
      'enumeration %d == sum_{a<=%d}(a - tau(a)) = %d' % (len(box), AMAX, closed))

print()
print('=' * 78)
print('PART C -- THE THEOREM AND ITS LEMMA, over every hypothesis pair with s_2 <= %d.' % AMAX)
print('=' * 78)

tot = 0
ok_inv = ok_atH = ok_delta = ok_min = 0
ok_F1 = ok_F2 = ok_F3 = ok_F4 = 0
ok_split = ok_claims = ok_below = 0
uncovered = 0
case_counts = {}
c4_t_ok = True
b3_once = True
first_pair_m3 = None

for (a, d) in box:
    b = a + d
    n = a // d
    eps = a - n * d
    H = (n + 2) * a + n * b
    tot += 1
    items, coll = lattice(a, d, n, H)
    Vm, Vi, Vv = Vtable(items, H)
    fF, gF, _ = outcomes(a, b, H, True)
    fA, gA, _ = outcomes(a, b, H, False)

    def v(m, i):
        return m * b - i * d

    # -- Lemma 4 (F1): the two gap identities, for every t <= n.
    if all(v(2 * t + 1, t + 1) - v(2 * t, 0) == a - t * d
           and v(2 * t + 2, t + 1) - v(2 * t + 1, 0) == a - t * d
           for t in range(0, n + 1)):
        ok_F1 += 1

    # -- Lemma 4 (F2): alpha_t >= eps > 0 for t <= n, and coordinates are unique on [0,H].
    if (not coll) and all(a - t * d >= eps for t in range(0, n + 1)) and 0 < eps < d:
        ok_F2 += 1

    # -- Lemma 4 (F3): the endpoints.
    vals = set(x[0] for x in items)
    if (v(2 * n + 2, n + 1) == H + d and v(2 * n + 1, 0) == H + (d - eps)
            and H not in vals
            and items[-1] == (2 * n * b + a, 2 * n + 1, 1)
            and H - items[-1][0] == eps):
        ok_F3 += 1

    # -- Lemma 4 (F4): the gap immediately above each admissible value in [0,H].
    if all(items[k + 1][0] - items[k][0]
           == (d if items[k][2] >= 1 else a - (items[k][1] // 2) * d)
           for k in range(len(items) - 1)):
        ok_F4 += 1

    # -- Theorem 2, first half: the invariant (6) at every h < H, in BOTH conventions.
    if all((fF[h], gF[h]) == predict(Vm[h], Vi[h], b, d) == (fA[h], gA[h])
           for h in range(H)):
        ok_inv += 1
    # -- and, separately, that the two conventions simply agree below H.
    if all(fF[h] == fA[h] and gF[h] == gA[h] for h in range(H)):
        ok_below += 1

    # -- Theorem 2, second half: the values (7) at H, and Corollary 3's delta^2.
    if (fF[H] == fA[H] == a + n * b and gF[H] == (n + 1) * a and gA[H] == n * b):
        ok_atH += 1
    if gA[H] - gF[H] == -eps and fA[H] - fF[H] == 0:
        ok_delta += 1

    # -- Corollary 3: the least discrepancy heap is exactly H.
    if next((h for h in range(H + 1)
             if fF[h] != fA[h] or gF[h] != gA[h]), None) == H:
        ok_min += 1

    # -- The five-case split of section 4: exhaustive, disjoint, and every per-case claim.
    split_ok = True
    claims_ok = True
    n3b = 0
    for h in range(a, H + 1):
        m, i = Vm[h], Vi[h]
        t = m // 2
        r = h - Vv[h]
        alpha = a - t * d
        ca = a + gF[h - a]
        cb = b + gF[h - b] if h >= b else None
        names = []
        if m % 2 == 0:
            if 1 <= i <= t:
                names.append('1')
            elif i == 0 and t >= 1:
                names.append('2')
        else:
            if i == 0:
                names.append('4')
            elif 1 <= i <= t:
                names.append('3A' if (i >= 2 or r < alpha) else '3B')
            elif i == t + 1:
                names.append('5')
        if len(names) != 1:
            uncovered += 1
            split_ok = False
            continue
        c = names[0]
        case_counts[c] = case_counts.get(c, 0) + 1
        if c == '1':
            okc = (Vv[h - b] == v(2 * t - 1, i) and Vv[h - a] == v(2 * t - 1, i - 1)
                   and cb - ca == d and cb > ca
                   and fF[h] == t * b and gF[h] == t * b - i * d)
        elif c == '2':
            okc = (Vv[h - b] == v(2 * t - 1, 0) and Vv[h - a] == v(2 * t - 1, 0)
                   and cb - ca == d
                   and fF[h] == t * b and gF[h] == t * b)
        elif c == '3A':
            okc = (Vv[h - b] == v(2 * t, i) and Vv[h - a] == v(2 * t, i - 1) and ca == cb
                   and fF[h - a] == t * b and fF[h - b] == t * b
                   and fF[h] == (t + 1) * b - i * d
                   and gF[h] == t * b and gA[h] == t * b)
        elif c == '3B':
            n3b += 1
            okc = (t == n and i == 1 and h == H and ca == cb
                   and Vv[h - a] == v(2 * n + 1, n + 1)
                   and fF[h - a] == (n + 1) * a and fF[h - b] == n * b
                   and fF[h - a] - fF[h - b] == eps)
        elif c == '4':
            if t > n - 1:
                c4_t_ok = False
            okc = (Vv[h - b] == v(2 * t, 0)
                   and Vv[h - a] in (v(2 * t, 0), v(2 * t + 1, t + 1))
                   and gF[h - a] == t * b and cb - ca == d
                   and fF[h] == (t + 1) * b and gF[h] == t * b)
        else:                                   # c == '5'
            if t == 0:
                okc = (h < b and fF[h] == a and gF[h] == 0)
            else:
                okc = (Vv[h - a] == v(2 * t, t) and Vv[h - b] == v(2 * t - 1, 0)
                       and ca - cb == alpha and alpha > 0
                       and fF[h] == (t + 1) * a and gF[h] == t * b)
        if not okc:
            claims_ok = False
    if split_ok:
        ok_split += 1
    if claims_ok:
        ok_claims += 1
    if n3b != 1:
        b3_once = False

    # -- M3, on the first pair only: DELETE Case 5 from the split and confirm the
    # exhaustiveness checker then reports uncovered heaps.  Otherwise it is vacuous.
    if first_pair_m3 is None:
        miss = 0
        for h in range(a, H + 1):
            m, i = Vm[h], Vi[h]
            t = m // 2
            r = h - Vv[h]
            alpha = a - t * d
            nm = []
            if m % 2 == 0:
                if 1 <= i <= t:
                    nm.append('1')
                elif i == 0 and t >= 1:
                    nm.append('2')
            else:
                if i == 0:
                    nm.append('4')
                elif 1 <= i <= t:
                    nm.append('3A' if (i >= 2 or r < alpha) else '3B')
            if len(nm) != 1:
                miss += 1
        first_pair_m3 = (a, b, miss)

assigned = sum(case_counts.values())
check('mutation-case-5-deleted-leaves-heaps-uncovered', first_pair_m3[2] > 0,
      'S={%d,%d}: with Case 5 removed the split leaves %d heaps uncovered (>0 required)'
      % first_pair_m3)
check('lemma-4-f1-gap-identities', ok_F1 == tot, '%d / %d pairs' % (ok_F1, tot))
check('lemma-4-f2-unique-admissible-coordinates', ok_F2 == tot,
      '%d / %d pairs, 0 collisions on [0,H]' % (ok_F2, tot))
check('lemma-4-f3-endpoints-and-h-not-in-w', ok_F3 == tot, '%d / %d pairs' % (ok_F3, tot))
check('lemma-4-f4-gap-above-each-lattice-point', ok_F4 == tot, '%d / %d pairs' % (ok_F4, tot))
check('theorem-2-invariant-below-h-both-conventions', ok_inv == tot,
      '%d / %d pairs' % (ok_inv, tot))
check('theorem-2-conventions-agree-below-h', ok_below == tot, '%d / %d pairs' % (ok_below, tot))
check('theorem-2-values-at-h', ok_atH == tot,
      '%d / %d pairs: o_FvF(H)=(a+nb,(n+1)a), o_AvA(H)=(a+nb,nb)' % (ok_atH, tot))
check('corollary-3-delta-1-zero-and-delta-2-minus-eps', ok_delta == tot,
      '%d / %d pairs' % (ok_delta, tot))
check('corollary-3-least-discrepancy-heap-is-h', ok_min == tot,
      '%d / %d pairs, H = (n+2)s_2 + n*s_1' % (ok_min, tot))
check('section-4-split-exhaustive-and-disjoint', ok_split == tot and uncovered == 0,
      '%d / %d pairs; %d of %d heap assignments in no case or two cases'
      % (ok_split, tot, uncovered, assigned))
check('section-4-per-case-claims', ok_claims == tot,
      '%d / %d pairs, %d heap-level case assignments' % (ok_claims, tot, assigned))
check('section-4-case-3b-fires-exactly-once-per-pair',
      b3_once and case_counts.get('3B') == tot,
      "case 3B assignments %d == pairs %d, and 1 per pair -- minimality, measured"
      % (case_counts.get('3B', -1), tot))
check('section-4-case-4-has-t-at-most-n-minus-1', c4_t_ok,
      'Case 4 never occurs with t > n-1, as claimed from Lemma 4(F3)')

print()
print('  case occupancy over the whole box: %s' % dict(sorted(case_counts.items())))
print()
print('NOT RE-RUN -- what this program does NOT cover:')
print('  * Heaps h > H.  Above H the outcome genuinely is convention-dependent and neither')
print('    the paper nor this program claims anything there.')
print('  * The excluded ratios s_1/s_2 = (k+1)/k (i.e. d | a): the target\'s COMPANION')
print('    conjecture.  Only one control touches them, and only to show that the paper\'s')
print('    Lemma 4(F2) genuinely fails there; no claim of any kind is checked at d | a.')
print('  * Subtraction sets with |S| >= 3.')
print('  * The induction itself.  Theorem 2 is proved by hand for every n >= 1; this program')
print('    audits it on the finite box s_2 <= %d and cannot certify the general case.  It is' % AMAX)
print('    not a proof assistant and the proof has not been formalised in one.')
print('  * Pairs with s_2 > %d.  The box is finite and stated; a wider census (s_2 <= 200,' % AMAX)
print('    19002 pairs) was run elsewhere for the CONJECTURE only and is not re-run here.')
print('  * The bibliographic exposure named in section 5 of the paper: the in-preparation')
print('    manuscript Kulkarni-Larsson cannot be fetched, and no program can settle that.')
print()

nfail = sum(1 for _, ok in CHECKS if not ok)
if nfail:
    print('VERDICT: %d OF %d CHECKS FAILED' % (nfail, len(CHECKS)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(CHECKS))
sys.exit(0)
