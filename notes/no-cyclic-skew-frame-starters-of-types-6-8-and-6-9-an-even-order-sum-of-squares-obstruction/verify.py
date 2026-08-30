#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- independent re-derivation of every quantity claimed in

    "No Cyclic Skew Frame Starters of Types 6^8 and 6^9:
     An Even-Order Sum-of-Squares Obstruction"   (paper.tex / paper.pdf)

The input of this program is exactly what the paper prints and nothing else:

  * the 35 (type, verdict) rows of Table 1 of arXiv:2211.12367v1 as transcribed
    into Table 1 of the paper, plus the off-table published nonexistence 5^5;
  * the nine rows whose authority column in the source reads "exhaustive search"
    (the dagger-marked rows of the paper's table);
  * the three cyclic skew frame starters printed in the paper's subsection
    "Objects for the negative polarity to be checked against", and the Z_21
    patterned anti-control described there;
  * the parameters of the 4^11 exact-cover model stated in the paper's
    Section "Reproduction and verification".

Nothing is imported beyond the Python standard library and no external data file
is read.  All arithmetic is exact integer arithmetic on Python ints; there is no
floating point anywhere and no randomness, so the run is deterministic.

Python 3.9 or later.

WHAT THIS PROGRAM DOES NOT DO is printed by the program itself, at the end, as a
`NOT RE-RUN:` line; the same statement appears in REVIEW_NOTE.md under `## Scope`.
"""

import itertools
import sys

# ---------------------------------------------------------------------------
# 0.  the check driver:  the VERDICT number is the number of PASS lines, by
#     construction, because both come from this one counter.
# ---------------------------------------------------------------------------

_N = 0
_BAD = 0


def check(name, cond, detail=''):
    """Print one PASS line, or one FAIL line and remember it."""
    global _N, _BAD
    _N += 1
    if cond:
        print('PASS %s%s' % (name, (' ' + detail) if detail else ''))
    else:
        _BAD += 1
        print('FAIL %s%s' % (name, (' ' + detail) if detail else ''))


# ---------------------------------------------------------------------------
# 1.  THE PAPER'S DATA, TRANSCRIBED.  (h, u) pairs; the type is h^u, g = h*u.
# ---------------------------------------------------------------------------

# Table 1 of arXiv:2211.12367v1, verdict "yes" (a skew frame starter is exhibited)
T1_YES = [(2, 5), (2, 13), (2, 17), (2, 25), (2, 29), (3, 13), (3, 19),
          (4, 5), (4, 13), (5, 7), (5, 11), (8, 5)]

# Table 1, verdict "no"
T1_NO = [(2, 8), (2, 9), (2, 12), (2, 16), (2, 20), (2, 21), (2, 24), (2, 28),
         (3, 7), (3, 9), (3, 11), (3, 15), (3, 17),
         (4, 7), (4, 8), (4, 9), (4, 10), (4, 12), (5, 9), (6, 5)]

# Table 1, verdict "?" -- the three open cells
T1_OPEN = [(4, 11), (6, 8), (6, 9)]

# published nonexistence omitted from Table 1 by its tablenote 2 (paper: 5^5 row)
OFF_TABLE_NO = [(5, 5)]

# the paper's dagger-marked rows: authority column reads "exhaustive search"
EXHAUSTIVE = [(2, 8), (2, 9), (2, 16), (3, 7), (4, 7), (4, 8), (4, 9), (4, 10), (6, 5)]

# The three starters printed in the paper, verbatim, as unordered pairs of
# residues.  Keyed by (h, u); g = h*u and H = <g/h> in every case.
PRINTED_STARTERS = {
    (2, 5): [(1, 2), (6, 8), (4, 7), (3, 9)],
    (4, 5): [(1, 2), (8, 19), (7, 9), (3, 11), (13, 16), (14, 18), (6, 12), (4, 17)],
    (5, 7): [(1, 2), (3, 9), (6, 19), (5, 32), (10, 30), (8, 33), (12, 24), (26, 29),
             (11, 15), (17, 34), (13, 18), (22, 31), (23, 25), (16, 27), (4, 20)],
}

# the paper's numbers for the three open cells and for 5^5:  (h, u) -> (N, Q mod g)
PAPER_CELL_NUMBERS = {
    (4, 11): (3510, 0),
    (6, 8): (4025, 40),
    (6, 9): (5176, 36),
    (5, 5): (996, 0),
}

# the H sets the paper prints for the two decided cells
PAPER_H = {
    (6, 8): [0, 8, 16, 24, 32, 40],
    (6, 9): [0, 9, 18, 27, 36, 45],
}

# the paper's 4^11 exact-cover model
MODEL_411 = dict(g=44, h=4, u=11, H=[0, 11, 22, 33], columns=80, rows=640,
                 pairs_total=780, diff_in_H=60, sum_in_H=80, degree=32, m=20,
                 level1_shards=32, subshards=869,
                 subshard_split=[(7, 26), (13, 27), (12, 28)])

ALL_CELLS = T1_YES + T1_NO + T1_OPEN + OFF_TABLE_NO


# ---------------------------------------------------------------------------
# 2.  the objects and invariants of the paper, recomputed from the definitions
# ---------------------------------------------------------------------------

def subgroup(g, h):
    """The unique subgroup of order h of Z_g:  H = <u> with u = g/h."""
    assert g % h == 0
    u = g // h
    return sorted({(u * k) % g for k in range(h)})


def upper_half(g, H):
    """The set called \\mathcal H in the paper (Lemma 'parity')."""
    Hs = set(H)
    return [j for j in range(1, (g - 1) // 2 + 1) if j not in Hs]


def Q_direct(g, H):
    """Q = sum of z^2 over z in Z_g \\ H, as a plain integer (not reduced)."""
    Hs = set(H)
    return sum(z * z for z in range(g) if z not in Hs)


def N_of(h, u):
    """N = (g-1)(2g-1) - u(h-1)(2h-1)."""
    g = h * u
    return (g - 1) * (2 * g - 1) - u * (h - 1) * (2 * h - 1)


def Q_closed(h, u):
    """Q = g*N/6, as an exact integer (the division is asserted to be exact)."""
    g = h * u
    N = N_of(h, u)
    num = g * N
    assert num % 6 == 0, ('g*N not divisible by 6', h, u, N)
    return num // 6


def clauses(h, u):
    """(i) h even => u odd;   (ii) 3|g => 3|h and u = 1 mod 3."""
    g = h * u
    c1 = (u % 2 == 1) if h % 2 == 0 else True
    c2 = (h % 3 == 0 and u % 3 == 1) if g % 3 == 0 else True
    return c1, c2


def fires(h, u, weak=False):
    """The criterion FIRES (forces nonexistence) when Q is nonzero mod g."""
    g = h * u
    q = Q_direct(g, subgroup(g, h)) % g
    return ((2 * q) % g != 0) if weak else (q != 0)


def starter_report(g, h, pairs):
    """Check (P1),(P2),(P3) of the paper on a candidate set of pairs.

    Returns a dict with a boolean per condition, plus the internal quantities
    Theta, sum d_i^2, sum s_i^2 and sum (x_i^2 + y_i^2) used in the proof of the
    paper's Theorem 'crit'.  Everything is an exact integer.
    """
    H = subgroup(g, h)
    T = sorted(set(range(g)) - set(H))
    m = (g - h) // 2
    ent, dif, sm = [], [], []
    sq_entries = 0
    sum_d2 = 0
    sum_s2 = 0
    for x, y in pairs:
        x, y = x % g, y % g
        ent += [x, y]
        d = (x - y) % g
        s = (x + y) % g
        dif += [d, (-d) % g]
        sm += [s, (-s) % g]
        sq_entries += x * x + y * y
        sum_d2 += d * d
        sum_s2 += s * s
    return {
        'count_ok': len(pairs) == m,
        'm': m,
        'P1': sorted(ent) == T,
        'P2': sorted(dif) == T,
        'P3': sorted(sm) == T,
        'Theta': sum(j * j for j in upper_half(g, H)),
        'sum_d2': sum_d2,
        'sum_s2': sum_s2,
        'sq_entries': sq_entries,
        'Q': Q_direct(g, H),
    }


# ---------------------------------------------------------------------------
# 3.  the checks
# ---------------------------------------------------------------------------

def main():
    print('verify.py -- "No Cyclic Skew Frame Starters of Types 6^8 and 6^9"')
    print('python: %d.%d.%d, standard library only, exact integer arithmetic'
          % sys.version_info[:3])
    print()

    # --- 3.1 the paper's transcription of Table 1 ---------------------------
    print('--- Table 1 of arXiv:2211.12367v1 as transcribed by the paper ---')
    check('table1-row-count',
          len(T1_YES) + len(T1_NO) + len(T1_OPEN) == 35,
          '[35 rows = %d yes + %d no + %d open]' % (len(T1_YES), len(T1_NO), len(T1_OPEN)))
    check('table1-partition-12-20-3',
          (len(T1_YES), len(T1_NO), len(T1_OPEN)) == (12, 20, 3),
          '[12 / 20 / 3, as the paper states]')
    check('table1-rows-distinct',
          len(set(T1_YES + T1_NO + T1_OPEN)) == 35,
          '[no type appears in two verdict classes]')
    check('table1-open-cells-are-4pow11-6pow8-6pow9',
          sorted(T1_OPEN) == [(4, 11), (6, 8), (6, 9)],
          '[g = 44, 48, 54]')
    check('exhaustive-search-rows-are-table1-no-rows',
          len(EXHAUSTIVE) == 9 and all(c in T1_NO for c in EXHAUSTIVE),
          '[9 dagger-marked rows, all with verdict no]')
    check('off-table-cell-not-a-table1-row',
          all(c not in T1_YES + T1_NO + T1_OPEN for c in OFF_TABLE_NO),
          '[5^5 is omitted from Table 1 by its tablenote 2]')
    print()

    # --- 3.2 Lemma "parity" ------------------------------------------------
    print('--- Lemma (parity) on every cell ---')
    bad = []
    for h, u in ALL_CELLS:
        g = h * u
        H = subgroup(g, h)
        if len(H) != h or (u % g) not in H:
            bad.append((h, u))
    check('subgroup-order-h', not bad, '[H = <g/h> has order h on all %d cells]' % len(ALL_CELLS))

    bad = [(h, u) for h, u in ALL_CELLS if (h * u) % 2 != h % 2]
    check('parity-g-equals-parity-h', not bad,
          '[g = h (mod 2) on all %d cells; equivalently g-h is even]' % len(ALL_CELLS))

    even = [(h, u) for h, u in ALL_CELLS if (h * u) % 2 == 0]
    bad = [(h, u) for h, u in even if (h * u // 2) not in subgroup(h * u, h)]
    check('even-order-involution-lies-in-H', not bad,
          '[g/2 in H on all %d even-order cells, so Z_g \\ H has no self-inverse element]' % len(even))

    bad = []
    for h, u in ALL_CELLS:
        g = h * u
        H = subgroup(g, h)
        half = upper_half(g, H)
        if len(half) != (g - h) // 2:
            bad.append((h, u))
        elif sorted(half + [(-j) % g for j in half]) != sorted(set(range(g)) - set(H)):
            bad.append((h, u))
    check('half-decomposition-of-T', not bad,
          '[T = half U (-half) with |half| = m = (g-h)/2 on all %d cells]' % len(ALL_CELLS))

    for h, u in T1_OPEN:
        g = h * u
        H = subgroup(g, h)
        check('parity-step-%d^%d' % (h, u), (g // 2) in H,
              '[g=%d, g/2=%d in H=%s]' % (g, g // 2, H))
    print()

    # --- 3.3 Q two ways, and Lemma "elementary form" ------------------------
    print('--- Q computed two independent ways, and its elementary form ---')
    bad = []
    for h, u in ALL_CELLS:
        g = h * u
        if Q_direct(g, subgroup(g, h)) != Q_closed(h, u):
            bad.append((h, u))
    check('Q-direct-equals-g-N-over-6', not bad,
          '[sum of squares over Z_g \\ H equals g*N/6 on all %d cells]' % len(ALL_CELLS))

    bad = []
    for h, u in ALL_CELLS:
        g = h * u
        H = subgroup(g, h)
        Theta = sum(j * j for j in upper_half(g, H))
        if (Q_direct(g, H) - 2 * Theta) % g != 0:
            bad.append((h, u))
    check('Q-congruent-2-Theta', not bad,
          '[the starter-free identity Q = 2*Theta (mod g) on all %d cells]' % len(ALL_CELLS))

    bad = []
    for h, u in ALL_CELLS:
        g = h * u
        q0 = Q_direct(g, subgroup(g, h)) % g == 0
        if q0 != (N_of(h, u) % 6 == 0):
            bad.append((h, u))
    check('Q-zero-iff-six-divides-N', not bad,
          '[Q = 0 (mod g)  <=>  6 | N  on all %d cells]' % len(ALL_CELLS))

    bad = []
    for h, u in ALL_CELLS:
        g = h * u
        q0 = Q_direct(g, subgroup(g, h)) % g == 0
        c1, c2 = clauses(h, u)
        if q0 != (c1 and c2):
            bad.append((h, u))
    check('clause-form-agrees', not bad,
          '[clauses (i) and (ii) hold  <=>  Q = 0 (mod g)  on all %d cells]' % len(ALL_CELLS))

    # the exhaustive check of the mod-3 half of Lemma "elementary form", over a
    # range far wider than the paper's cells
    bad = []
    for h in range(2, 41):
        for u in range(2, 41):
            g = h * u
            if (g - h) % 2:
                continue
            c1, c2 = clauses(h, u)
            if (N_of(h, u) % 6 == 0) != (c1 and c2):
                bad.append((h, u))
    check('clause-form-agrees-on-a-wide-range', not bad,
          '[6 | N  <=>  clauses (i) and (ii), over every admissible 2<=h,u<=40]')
    print()

    # --- 3.4 both polarities on Table 1 ------------------------------------
    print('--- both polarities of the criterion on Table 1 (paper items (a)-(e)) ---')
    yes_fire = [c for c in T1_YES if fires(*c)]
    check('a-criterion-never-fires-on-a-published-yes-row',
          len(yes_fire) == 0,
          '[0 of 12 -- a necessary condition must be silent where a starter is exhibited]')

    no_fire = [c for c in T1_NO if fires(*c)]
    check('b-criterion-fires-on-18-of-20-no-rows',
          len(no_fire) == 18, '[%d of 20]' % len(no_fire))
    silent_no = sorted('%d^%d' % c for c in T1_NO if not fires(*c))
    check('b-silent-exactly-on-3pow7-and-4pow7',
          silent_no == ['3^7', '4^7'], '[silent on %s]' % ' '.join(silent_no))
    check('b-off-table-5pow5-silent',
          not fires(5, 5) and N_of(5, 5) == 996 and 996 == 6 * 166,
          '[N = 996 = 6*166, so Q = 0 mod 25 and the criterion is silent]')

    open_fire = sorted('%d^%d' % c for c in T1_OPEN if fires(*c))
    check('c-criterion-fires-on-exactly-two-open-cells',
          open_fire == ['6^8', '6^9'], '[fires on %s, silent on 4^11]' % ' '.join(open_fire))

    ex_fire = [c for c in EXHAUSTIVE if fires(*c)]
    check('d-supersedes-seven-of-nine-exhaustive-searches',
          len(ex_fire) == 7 and sorted('%d^%d' % c for c in ex_fire) ==
          sorted(['2^8', '2^9', '2^16', '4^8', '4^9', '4^10', '6^5']),
          '[%s]' % ' '.join('%d^%d' % c for c in ex_fire))
    check('d-all-seven-have-even-order',
          all((h * u) % 2 == 0 for h, u in ex_fire),
          '[g = %s]' % ' '.join(str(h * u) for h, u in ex_fire))
    check('d-the-two-not-superseded-are-the-silent-no-rows',
          sorted('%d^%d' % c for c in EXHAUSTIVE if not fires(*c)) == ['3^7', '4^7'],
          '[3^7 and 4^7 remain the source\'s own searches; the criterion adds nothing there]')

    weak_no = [c for c in T1_NO if fires(*c, weak=True)]
    check('e-weak-form-fires-on-only-12-of-20-no-rows',
          len(weak_no) == 12, '[%d of 20, versus 18 for the strong form]' % len(weak_no))
    weak_ex = [c for c in EXHAUSTIVE if fires(*c, weak=True)]
    check('e-weak-form-supersedes-only-three-exhaustive-searches',
          sorted('%d^%d' % c for c in weak_ex) == sorted(['2^9', '4^9', '6^5']),
          '[%s]' % ' '.join('%d^%d' % c for c in weak_ex))
    check('e-weak-form-still-kills-6pow8-and-6pow9',
          fires(6, 8, weak=True) and fires(6, 9, weak=True),
          '[2Q = 32 mod 48 and 2Q = 18 mod 54]')
    check('e-weak-form-implied-by-strong-form',
          all((not fires(*c, weak=True)) or fires(*c) for c in ALL_CELLS),
          '[every cell the weak form kills, the strong form kills too]')
    print()

    # --- 3.5 the two decided cells, number by number ------------------------
    print('--- Theorem: the cells 6^8 and 6^9 ---')
    for h, u in [(6, 9), (6, 8)]:
        g = h * u
        H = subgroup(g, h)
        N = N_of(h, u)
        Qi = Q_direct(g, H)
        pN, pQ = PAPER_CELL_NUMBERS[(h, u)]
        check('cell-%d^%d-subgroup-as-printed' % (h, u), H == PAPER_H[(h, u)],
              '[H = %s]' % H)
        check('cell-%d^%d-m' % (h, u), (g - h) // 2 == (24 if u == 9 else 21),
              '[m = (g-h)/2 = %d pairs]' % ((g - h) // 2))
        check('cell-%d^%d-N' % (h, u), N == pN,
              '[N = %d, paper says %d; 6 does not divide it (N mod 6 = %d)]'
              % (N, pN, N % 6))
        check('cell-%d^%d-Q' % (h, u), Qi == g * N // 6 and Qi % g == pQ,
              '[Q = %d = %d*%d + %d, so Q = %d mod %d, nonzero]'
              % (Qi, g, Qi // g, Qi % g, Qi % g, g))
        c1, c2 = clauses(h, u)
        check('cell-%d^%d-clause-that-fails' % (h, u), not (c1 and c2),
              '[clause (i) %s, clause (ii) %s]'
              % ('holds' if c1 else 'FAILS', 'holds' if c2 else 'FAILS'))
        check('cell-%d^%d-nonexistence' % (h, u), fires(h, u),
              '[NO CYCLIC SKEW FRAME STARTER OF TYPE %d^%d]' % (h, u))
    check('cell-6pow9-arithmetic-as-printed',
          53 * 107 - 9 * 5 * 11 == 5176 and 54 * 5176 // 6 == 46584
          and 46584 == 54 * 862 + 36,
          '[53*107 - 9*5*11 = 5176; 54*5176/6 = 46584 = 54*862 + 36]')
    check('cell-6pow8-arithmetic-as-printed',
          47 * 95 - 8 * 5 * 11 == 4025 and 4025 % 2 == 1
          and 48 * 4025 // 6 == 32200 and 32200 == 48 * 670 + 40,
          '[47*95 - 8*5*11 = 4025 is odd; 48*4025/6 = 32200 = 48*670 + 40]')
    print()

    # --- 3.6 the cell the paper does NOT settle ----------------------------
    print('--- Remark: the criterion is provably silent on 4^11 ---')
    N411 = N_of(4, 11)
    check('cell-4pow11-N-divisible-by-six',
          N411 == 3510 and N411 == 6 * 585,
          '[N = 43*87 - 11*3*7 = 3741 - 231 = 3510 = 6*585]')
    check('cell-4pow11-criterion-silent',
          Q_direct(44, subgroup(44, 4)) % 44 == 0 and not fires(4, 11),
          '[Q = 0 mod 44: both clauses hold, so the criterion decides nothing]')
    print()

    # --- 3.7 the printed objects -------------------------------------------
    print('--- the three starters printed in the paper, against (P1),(P2),(P3) ---')
    for (h, u) in sorted(PRINTED_STARTERS):
        g = h * u
        r = starter_report(g, h, PRINTED_STARTERS[(h, u)])
        check('printed-starter-%d^%d-is-a-cyclic-skew-frame-starter' % (h, u),
              r['count_ok'] and r['P1'] and r['P2'] and r['P3'],
              '[Z_%d, m=%d pairs, P1 %s, P2 %s, P3 %s]'
              % (g, r['m'], r['P1'], r['P2'], r['P3']))
        check('printed-starter-%d^%d-proof-quantities' % (h, u),
              (r['sq_entries'] - r['Q']) % g == 0
              and (r['sum_d2'] - r['Theta']) % g == 0
              and (r['sum_s2'] - r['Theta']) % g == 0
              and (r['sum_s2'] + r['sum_d2'] - 2 * r['Q']) % g == 0,
              '[sum(x^2+y^2) = Q, sum d^2 = sum s^2 = Theta = %d, and '
              'sum s^2 + sum d^2 = 2Q, all mod %d]' % (r['Theta'] % g, g))
        check('printed-starter-%d^%d-satisfies-the-theorem' % (h, u),
              r['Q'] % g == 0,
              '[Q = 0 mod %d, as Theorem (crit) requires of any type carrying a starter]' % g)

    print('--- the Z_21 anti-control: a frame starter that is not skew ---')
    anti = [(x, (-x) % 21) for x in range(1, 11) if x % 7]
    r = starter_report(21, 3, anti)
    check('anti-control-has-the-right-size', r['count_ok'] and len(anti) == 9,
          '[m = (21-3)/2 = 9 pairs {x,-x}, 1 <= x <= 10, 7 does not divide x]')
    check('anti-control-satisfies-P1-and-P2', r['P1'] and r['P2'],
          '[it really is a frame starter of type 3^7]')
    check('anti-control-is-rejected-on-P3-alone', not r['P3'],
          '[every sum x + (-x) is 0, which lies in H, so (P3) fails]')
    check('anti-control-cell-is-a-silent-row', not fires(3, 7),
          '[3^7 has Q = 0 mod 21: the obstruction there is not Q, consistently with item (b)]')

    check('source-Z10-example',
          Q_direct(10, [0, 5]) == 260 and 260 == 10 * 26,
          '[sum of squares over Z_10 \\ {0,5} is 260 = 10*26, so Q = 0 mod 10]')
    print()

    # --- 3.8 the 4^11 exact-cover model arithmetic --------------------------
    print('--- the 4^11 exact-cover model stated in the paper (arithmetic only) ---')
    M = MODEL_411
    g, h = M['g'], M['h']
    H = subgroup(g, h)
    check('model-subgroup', H == M['H'], '[H = <11> = %s]' % H)
    T = [z for z in range(g) if z not in set(H)]
    check('model-T-size', len(T) == 40, '[|Z_44 \\ H| = %d]' % len(T))
    check('model-no-self-inverse-in-T', 22 in H,
          '[22 in H, so every class {z,-z} in T has two elements]')
    nclass = len({frozenset((z, (-z) % g)) for z in T})
    check('model-column-count',
          len(T) + nclass + nclass == M['columns'] and nclass == 20,
          '[%d element + %d difference + %d sum classes = %d columns]'
          % (len(T), nclass, nclass, M['columns']))

    allpairs = list(itertools.combinations(T, 2))
    bad_d = [p for p in allpairs if (p[0] - p[1]) % g in set(H)]
    bad_s = [p for p in allpairs if (p[0] + p[1]) % g in set(H)]
    both = [p for p in allpairs if p in set(bad_d) and p in set(bad_s)]
    good = [p for p in allpairs
            if (p[0] - p[1]) % g not in set(H) and (p[0] + p[1]) % g not in set(H)]
    check('model-pair-total', len(allpairs) == M['pairs_total'],
          '[C(40,2) = %d]' % len(allpairs))
    check('model-excluded-counts',
          len(bad_d) == M['diff_in_H'] and len(bad_s) == M['sum_in_H'] and len(both) == 0,
          '[%d with difference in H, %d with sum in H, %d with both]'
          % (len(bad_d), len(bad_s), len(both)))
    check('model-row-count-by-inclusion-exclusion',
          len(good) == M['rows']
          and M['pairs_total'] - M['diff_in_H'] - M['sum_in_H'] == M['rows'],
          '[780 - 60 - 80 = %d admissible rows]' % len(good))
    degs = {x: sum(1 for p in good if x in p) for x in T}
    check('model-row-count-by-degree',
          set(degs.values()) == {M['degree']}
          and len(T) * M['degree'] // 2 == M['rows'],
          '[every element has exactly %d admissible partners, and 40*32/2 = %d]'
          % (M['degree'], M['rows']))
    check('model-solution-size', (g - h) // 2 == M['m']
          and M['m'] * 4 == M['columns'],
          '[a solution is %d rows, each covering 2 element + 1 difference + 1 sum '
          'column, so 20*4 = %d columns]' % (M['m'], M['columns']))
    check('model-shard-arithmetic',
          sum(a for a, _ in M['subshard_split']) == M['level1_shards']
          and sum(a * b for a, b in M['subshard_split']) == M['subshards'],
          '[7 + 13 + 12 = %d level-one shards; 7*26 + 13*27 + 12*28 = %d sub-shards]'
          % (M['level1_shards'], M['subshards']))
    print()

    # --- 3.9 verdict --------------------------------------------------------
    print('NOT RE-RUN: the exhaustive census of the 4^11 exact-cover instance.  The paper')
    print('  does not claim the 4^11 cell and this program does not decide it.  Only the MODEL')
    print('  arithmetic above (80 columns, 640 rows, 869 sub-shards) is checked; no search is')
    print('  performed here, and the engine that ran that census is not part of this folder.')
    print('NOT RE-RUN: nothing about NONCYCLIC groups of order 44, 48 or 54, about skew Room')
    print('  FRAMES as opposed to starters, or about STRONG (non-skew) frame starters.  Those')
    print('  boundaries are asserted in the paper on the cited literature, not recomputed here.')
    print('NOT RE-RUN: the transcription itself.  The 35 verdicts, the nine "exhaustive search"')
    print('  authorities and the three printed starters are read off arXiv:2211.12367v1 by hand;')
    print('  this program checks that they are internally consistent with the criterion, which')
    print('  is a strong test of the transcription but not a substitute for reading the source.')
    print()

    if _BAD:
        print('VERDICT: %d of %d CHECKS FAILED' % (_BAD, _N))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _N)
    return 0


if __name__ == '__main__':
    sys.exit(main())
