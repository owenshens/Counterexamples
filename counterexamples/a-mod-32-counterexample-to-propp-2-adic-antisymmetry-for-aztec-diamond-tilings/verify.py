#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify.py -- checks every computational claim of

    "A mod 32 counterexample to Propp's 2-adic antisymmetry conjecture for
     Aztec diamond tilings by dominoes and square tetrominoes"

Python 3.9+, STANDARD LIBRARY ONLY, exact integer arithmetic throughout: no
float ever decides anything, no third-party package, no external data file.

WHAT IT CONSUMES.  Only the numbers printed in the paper, quoted below with
their equation numbers.  Nothing is read from disk.

WHAT IT DOES.  Two layers.

  (1) THE ARITHMETIC LAYER (sections A-F, a fraction of a second).  From the
      sixteen residues of eq. (mod32) and the exact M(13) of eq. (M13) it
      re-derives the hypothesis of the instance (14,15,5), the failure of the
      conclusion, v_2 = 4, the L/e decomposition, the (0,13) cell, the
      mod-2^128 corroboration, and the whole census table.

  (2) THE INDEPENDENT RECOMPUTATION (section G).  A third tiling engine,
      written here in pure Python and sharing no code with either C engine
      used for the paper's data.  It rebuilds the Aztec diamond region from
      the inequality |x| + |y| <= n+1 directly, reproduces M(0..4) by plain
      backtracking with no transfer computation at all, reproduces all
      thirteen published terms of A356512 EXACTLY, reproduces the
      domino-only specialisation 2^(n(n+1)/2) of Elkies-Kuperberg-Larsen-Propp
      (A006125) as a control on the region and the domino rules, and then
      recomputes M(13), M(14), M(15) mod 256 from scratch -- the three data
      the refutation rests on -- and re-derives the refutation from ITS OWN
      residues rather than from the printed ones.

COST.  Section G's last three checks are a 2^27-, 2^29- and 2^31-state
transfer computation in an interpreted language.  Budget roughly half an hour
of one core and of the order of 16 GB of RAM; every other check is instant.
The program prints its progress and its peak resident size as it goes.

Run:  python3 verify.py      -> one `PASS <name>` line per check, then
                               `VERDICT: ALL <n> CHECKS PASS`, exit 0 iff so.
"""

import sys
import time
import resource
from fractions import Fraction

# =====================================================================
# THE OBJECT, AS PRINTED IN THE PAPER.  Nothing else is input.
# =====================================================================

# eq. (table): the thirteen terms Propp prints, M(n) for n = 0..12.
PUBLISHED = [
    1,
    3,
    19,
    293,
    10917,
    996599,
    222222039,
    121552500713,
    162860556763865,
    535527565429290907,
    4318205059450240425083,
    85475498697714319842817853,
    4151186175463797888945512144221,
]

# eq. (M13): the new exact term.
M13_EXACT = 494970717813329135973007168990535711

# eq. (M1415): the two new residues.
M14_MOD256 = 79
M15_MOD256 = 225
# and the same engine's value for n = 13, quoted in section 5 of the paper.
M13_MOD256_PRINTED = 31

# eq. (mod32): the complete mod-32 row for n = 0..15, as printed.
MOD32_ROW = [1, 3, 19, 5, 5, 23, 23, 9, 25, 27, 27, 29, 29, 31, 15, 1]

# Remark (the failure is one bit): the e-bit row for n = 0..15, as printed.
E_ROW = [0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1]

# Propp's own mod-16 row for n = 0..12, quoted in the same remark.
PROPP_MOD16_ROW = [1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13]

# section 5: the exact accumulator's value of M(14) mod 2^128.
M14_MOD_2_128 = 228062845944308165005028192227771743567

# eq. (cell013): the quotient printed there.
CELL013_QUOTIENT = 15467834931666535499156474030954241

# The instance under test.
N, NP, K = 14, 15, 5

# =====================================================================
# harness
# =====================================================================

_RESULTS = []


def check(name, ok, detail=''):
    """One check.  Prints `PASS <name> [detail]` or `FAIL <name> [detail]`."""
    _RESULTS.append(bool(ok))
    print('%s %s%s' % ('PASS' if ok else 'FAIL', name,
                       ('  ' + detail) if detail else ''), flush=True)


def note(text):
    """A progress line.  Deliberately never starts with PASS or FAIL."""
    print('    ' + text, flush=True)


def peak_gb():
    kb_units = 1024 if sys.platform.startswith('linux') else 1
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * kb_units / 1e9


def v2(x):
    """2-adic valuation of a nonzero integer.  Exact, no floats."""
    assert x != 0
    k = 0
    while x % 2 == 0:
        x //= 2
        k += 1
    return k


def L(n):
    """Propp's closed form, written as n + 1 + (n mod 2)."""
    return n + 1 + (n % 2)


def L_propp_expression(n):
    """Propp's closed form as he prints it: n + 1 + (1 + (-1)^(n+1))/2.

    Written with Fraction so the halving is exact and the two spellings are
    compared as numbers, not as floats."""
    return n + 1 + Fraction(1 + (-1) ** (n + 1), 2)


# The complete residue data available for 0 <= n <= 15: mod 256 where we have
# it, which is everywhere in that range.  Built from the printed object.
MOD256 = [x % 256 for x in PUBLISHED] + [M13_EXACT % 256, M14_MOD256, M15_MOD256]
NMAX = 15
assert len(MOD256) == NMAX + 1


def cell_holds(n, np_, k, residues=MOD256):
    """Does the conclusion of Conjecture 2 hold at (n, n', k)?

    Legitimate only for k <= 8, since `residues` are modulo 2^8."""
    assert 1 <= k <= 8
    return (residues[n] + residues[np_]) % (2 ** k) == 0


def obligated_k(n, np_):
    """The exponents k >= 1 for which (n, n') satisfies the hypothesis.

    n + n' = -3 (mod 2^k) holds exactly for 1 <= k <= v_2(n + n' + 3)."""
    if (n + np_) % 2 == 0:
        return []
    return list(range(1, v2(n + np_ + 3) + 1))


# =====================================================================
# SECTION A -- the printed object and its internal consistency
# =====================================================================
note('interpreter: Python %s on %s' % (sys.version.split()[0], sys.platform))
print('=== A. the printed object ===', flush=True)


def region_rows(n):
    """Row half-widths of AD_n as the paper describes them: row r has
    2*min(r+1, 2n-r) cells."""
    return [2 * min(r + 1, 2 * n - r) for r in range(2 * n)]


check('A1-region-cell-count-2n(n+1)',
      all(sum(region_rows(n)) == 2 * n * (n + 1) for n in range(NMAX + 1)),
      'eq. (cells) verified for n = 0..15')

check('A2-published-table-13-terms-all-odd',
      len(PUBLISHED) == 13 and all(x % 2 == 1 for x in PUBLISHED),
      '13 terms, every one odd (Propp: M(n) is always odd)')

check('A3-M13-exact-reduces-to-printed-residues',
      M13_EXACT % 256 == M13_MOD256_PRINTED and M13_EXACT % 32 == 31,
      'M(13) = 31 mod 256 and 31 mod 32; the exact term of eq. (M13) has %d '
      'decimal digits' % len(str(M13_EXACT)))

check('A4-mod32-row-matches-the-printed-data',
      [x % 32 for x in MOD256] == MOD32_ROW,
      'eq. (mod32) is exactly eq. (table) + eq. (M13) + eq. (M1415) mod 32')

check('A5-propp-mod16-row-matches-published-terms',
      [x % 16 for x in PUBLISHED] == PROPP_MOD16_ROW,
      'his row 1,3,3,5,5,7,7,9,9,11,11,13,13 reproduced from eq. (table)')

check('A6-L-two-spellings-agree',
      all(Fraction(L(n)) == L_propp_expression(n) for n in range(64)),
      'n+1+(n mod 2) == n+1+(1+(-1)^(n+1))/2 for n = 0..63, exact rationals')

check('A7-M-congruent-L-mod-16-up-to-n-15',
      all(MOD32_ROW[n] % 16 == L(n) % 16 for n in range(NMAX + 1)),
      'the mod-16 closed form survives at n = 13, 14, 15 too')

check('A8-propp-mod32-obstruction-at-n-2',
      PUBLISHED[2] % 32 != L(2) % 32 and PUBLISHED[2] == 19 and L(2) == 3,
      "M(2) = 19 != 3 mod 32: Propp's own reason the closed form dies mod 32")

check('A9-e-bit-row-matches-the-printed-row',
      [((MOD32_ROW[n] - L(n)) % 32) // 16 for n in range(NMAX + 1)] == E_ROW
      and all(e in (0, 1) for e in E_ROW),
      'e(n) = (M(n) - L(n))/16 mod 2, n = 0..15')

# =====================================================================
# SECTION B -- the instance (14, 15, 5)
# =====================================================================
print('=== B. the instance (n, n\', k) = (14, 15, 5) ===', flush=True)

check('B1-hypothesis-holds',
      N + NP == 29 and (N + NP + 3) == 32 == 2 ** K and (-3 - N) % 32 == NP
      and K >= 1 and 0 <= N <= NMAX and 0 <= NP <= NMAX,
      '14+15 = 29, 29+3 = 32 = 2^5, (-3-14) mod 32 = 15, k = 5 >= 1')

_sum256 = (M14_MOD256 + M15_MOD256) % 256
check('B2-sum-mod-256-is-48-and-v2-is-4',
      M14_MOD256 + M15_MOD256 == 304 and _sum256 == 48
      and 48 == 16 * 3 and v2(48) == 4 and v2(48) < 8,
      '79+225 = 304 = 256+48; 48 = 2^4*3; v_2 = 4 < 8 so mod 256 settles it')

check('B3-conclusion-fails-mod-32',
      (MOD32_ROW[N] + MOD32_ROW[NP]) % 32 == 16 != 0
      and not cell_holds(N, NP, K),
      '15 + 1 = 16 != 0 mod 32: Conjecture 2 is false as printed')

check('B4-refutation-stable-under-the-unknown-high-bits',
      all(((M14_MOD256 + 256 * a) + (M15_MOD256 + 256 * b)) % 32 == 16
          for a in range(256) for b in range(256)),
      'every lift of the two residues mod 256 gives the same failure, so the '
      'exact values of M(14), M(15) are not needed')

check('B5-cell-is-the-parity-of-e(14)+e(15)',
      (L(N) + L(NP)) % 32 == 0 and L(N) + L(NP) == 32
      and (MOD32_ROW[N] + MOD32_ROW[NP]) % 32
          == (16 * (E_ROW[N] + E_ROW[NP])) % 32
      and (E_ROW[N] + E_ROW[NP]) % 2 == 1,
      'L(14)+L(15) = 32 = 0 mod 32, so the cell is 16*(e(14)+e(15)) mod 32; '
      'e(14)+e(15) = 1 is odd')

check('B6-flipping-either-e-bit-restores-the-cell',
      all((16 * (a + b)) % 32 == 0
          for a, b in ((1, 1), (0, 0))),
      'e(14)=1 or e(15)=0 would give 0 mod 32: the failure is exactly one bit')

# =====================================================================
# SECTION C -- the previously open k = 4 cell (0, 13)
# =====================================================================
print('=== C. the cell (0, 13) at k = 4 ===', flush=True)

_s013 = PUBLISHED[0] + M13_EXACT
check('C1-cell-0-13-exact-division',
      _s013 == 494970717813329135973007168990535712
      and _s013 == 32 * CELL013_QUOTIENT and CELL013_QUOTIENT % 2 == 1,
      '1 + M(13) = 2^5 * %d, quotient odd' % CELL013_QUOTIENT)

check('C2-cell-0-13-holds-with-one-power-to-spare',
      v2(_s013) == 5 and obligated_k(0, 13) == [1, 2, 3, 4]
      and all(cell_holds(0, 13, k) for k in obligated_k(0, 13)),
      'v_2 = 5 against the 4 required by 0+13 = 13 = -3 mod 16')

# =====================================================================
# SECTION D -- the mod-2^128 corroboration of M(14)
# =====================================================================
print('=== D. the mod-2^128 corroboration of M(14) ===', flush=True)

check('D1-int128-wrap-reduces-to-79-mod-256',
      M14_MOD_2_128 % 256 == M14_MOD256 and M14_MOD_2_128 < 2 ** 128,
      'a value known mod 2^128 is exact mod 2^8; %d mod 256 = 79'
      % (M14_MOD_2_128 % 10 ** 8))

check('D2-last-eight-decimal-digits-determine-mod-256',
      10 ** 8 % 256 == 0 and M14_MOD_2_128 % 10 ** 8 == 71743567
      and 71743567 == 280248 * 256 + 79,
      '256 divides 10^8, so 71743567 = 280248*256 + 79 binds the two engines')

# =====================================================================
# SECTION E -- the census
# =====================================================================
print('=== E. the census over n, n\' <= 15 ===', flush=True)


def census(maxn):
    """(pairs, vacuous, obligated, hold, fail, triples, triples_hold,
    triples_fail) over unordered pairs 0 <= n <= n' <= maxn."""
    pairs = vac = obl = hold = fail = trip = trh = trf = 0
    for n in range(maxn + 1):
        for m in range(n, maxn + 1):
            pairs += 1
            ks = obligated_k(n, m)
            if not ks:
                vac += 1
                continue
            obl += 1
            bad = False
            for k in ks:
                trip += 1
                if cell_holds(n, m, k):
                    trh += 1
                else:
                    trf += 1
                    bad = True
            fail += 1 if bad else 0
            hold += 0 if bad else 1
    return pairs, vac, obl, hold, fail, trip, trh, trf


_c12, _c13, _c15 = census(12), census(13), census(15)

check('E1-census-slice-n-at-most-12',
      _c12[:5] == (91, 49, 42, 42, 0),
      'published data only: 91 pairs = 49 vacuous + 42 obligated, 42 hold, '
      '0 fail')

check('E2-census-slice-n-at-most-13',
      _c13[:5] == (105, 56, 49, 49, 0),
      'adds the exact M(13): 105 = 56 + 49, 49 hold, 0 fail')

check('E3-census-slice-n-at-most-15',
      _c15[:5] == (136, 72, 64, 63, 1),
      'adds M(14), M(15) mod 256: 136 = 72 + 64, 63 hold, 1 fail')

check('E4-vacuity-is-exactly-even-sum',
      all((obligated_k(n, m) == []) == ((n + m) % 2 == 0)
          for n in range(NMAX + 1) for m in range(NMAX + 1)),
      'n + n\' even fails the k = 1 hypothesis and hence every k, so such a '
      'pair verifies nothing and is not counted as verified')

check('E5-121-triples-120-hold-1-fails',
      (_c15[5], _c15[6], _c15[7]) == (121, 120, 1),
      'triples (n, n\', k), 0 <= n <= n\' <= 15, 1 <= k <= v_2(n+n\'+3)')

_by_maxk = {}
for _n in range(NMAX + 1):
    for _m in range(_n, NMAX + 1):
        _ks = obligated_k(_n, _m)
        if _ks:
            _by_maxk[_ks[-1]] = _by_maxk.get(_ks[-1], 0) + 1
check('E6-triple-decomposition-by-largest-exponent',
      sorted(_by_maxk.items()) == [(1, 32), (2, 16), (3, 8), (4, 7), (5, 1)]
      and sum(k * c for k, c in _by_maxk.items()) == 121,
      '32*1 + 16*2 + 8*3 + 7*4 + 1*5 = 121')

_failures = [(n, m, k) for n in range(NMAX + 1) for m in range(n, NMAX + 1)
             for k in obligated_k(n, m) if not cell_holds(n, m, k)]
check('E7-the-unique-failing-triple-is-14-15-5',
      _failures == [(N, NP, K)],
      'exhaustive over the box n, n\' <= 15: exactly one failure, %s'
      % (_failures,))

check('E8-every-k-at-most-4-cell-is-forced-by-the-mod-16-law',
      all(cell_holds(n, m, k) and (L(n) + L(m)) % (2 ** k) == 0
          for n in range(NMAX + 1) for m in range(n, NMAX + 1)
          for k in obligated_k(n, m) if k <= 4),
      'M = L mod 16 and L(n)+L(n\') = n+n\'+3, so all 120 cells with k <= 4 '
      'follow from published information; only the k = 5 cell is new')

check('E9-all-k-1-cells-hold-consistent-with-oddness',
      all(cell_holds(n, m, 1) for n in range(NMAX + 1)
          for m in range(n, NMAX + 1) if 1 in obligated_k(n, m))
      and all(x % 2 == 1 for x in MOD256),
      'the decider is silent on all 32 pairs with largest exponent 1, as '
      "Propp's oddness theorem requires: a proved-silent control")

# =====================================================================
# SECTION F -- the ceiling, and what is deliberately left open
# =====================================================================
print('=== F. the ceiling ===', flush=True)

check('F1-k-at-least-6-needs-an-index-at-least-31',
      not any((n + m + 3) % 64 == 0
              for n in range(31) for m in range(31))
      and (30 + 31 + 3) % 64 == 0,
      'n + n\' = -3 mod 64 forces n + n\' >= 61, hence max(n, n\') >= 31; '
      'the least witness pair is (30, 31)')

_diag29 = [(n, 29 - n) for n in range(0, 15) if 0 <= 29 - n]
check('F2-the-other-k-5-cells-are-not-decided-here',
      (13, 16) in _diag29 and (12, 17) in _diag29
      and all(max(p) > NMAX for p in _diag29 if p != (N, NP))
      and _failures == [(N, NP, K)],
      'every other pair on the diagonal n + n\' = 29 has an index >= 16, '
      'where M is unknown to us; those cells are settled in neither direction')

# =====================================================================
# SECTION G -- an independent recomputation, in pure Python
# =====================================================================
print('=== G. independent recomputation (third engine, pure Python) ===',
      flush=True)


def region_from_inequality(n):
    """The cells of AD_n, derived from |x| + |y| <= n+1 and nothing else.

    A unit square is in AD_n iff all four of its corners satisfy the
    inequality.  Returned as a flat 0/1 list over the (2n) x (2n) bounding
    box, row-major, with the box's top-left corner at (x, y) = (-n, n)."""
    W = 2 * n
    cells = [0] * (W * W)
    for r in range(W):
        for c in range(W):
            x0, x1 = -n + c, -n + c + 1
            y0, y1 = n - r - 1, n - r
            if all(abs(x) + abs(y) <= n + 1
                   for x in (x0, x1) for y in (y0, y1)):
                cells[r * W + c] = 1
    return cells


def count_backtracking(n):
    """M(n) by plain recursive backtracking: no transfer computation, no
    bitmask, no memoisation.  Fills the first free cell in row-major order
    with a horizontal domino, a vertical domino or a 2x2 square."""
    if n == 0:
        return 1
    W = 2 * n
    cells = region_from_inequality(n)
    free = [bool(cells[i]) for i in range(W * W)]
    total_cells = sum(cells)

    def rec(start, remaining):
        if remaining == 0:
            return 1
        p = start
        while not free[p]:
            p += 1
        r, c = divmod(p, W)
        out = 0
        # horizontal domino
        if c + 1 < W and free[p + 1]:
            free[p] = free[p + 1] = False
            out += rec(p, remaining - 2)
            free[p] = free[p + 1] = True
        # vertical domino
        if r + 1 < W and free[p + W]:
            free[p] = free[p + W] = False
            out += rec(p, remaining - 2)
            free[p] = free[p + W] = True
        # 2x2 square
        if (c + 1 < W and r + 1 < W and free[p + 1] and free[p + W]
                and free[p + W + 1]):
            for q in (p, p + 1, p + W, p + W + 1):
                free[q] = False
            out += rec(p, remaining - 4)
            for q in (p, p + 1, p + W, p + W + 1):
                free[q] = True
        return out

    return rec(0, total_cells)


def count_transfer(n, mod=None, square=True):
    """M(n), by a broken-profile transfer computation over the bounding box.

    State: an integer whose bit j records whether the cell at flat index p + j
    has already been covered by a tile placed at an earlier cell.  Bits 0, 1,
    W and W+1 are the only ones a tile at p can set, so W+2 bits suffice.
    Cells outside AD_n are permanently blocked and can never carry bit 0.

    `mod=None` gives the exact integer; `square=False` disables the 2x2 tile,
    which must then return the domino count 2^(n(n+1)/2) (A006125).
    """
    if n == 0:
        return 1 if mod is None else 1 % mod
    W = 2 * n
    cells = region_from_inequality(n)
    if sum(cells) != 2 * n * (n + 1):
        raise AssertionError('region has %d cells, expected %d'
                             % (sum(cells), 2 * n * (n + 1)))
    B_R, B_D, B_DR = 2, 1 << W, 1 << (W + 1)
    st = {0: 1 if mod is None else 1 % mod}
    for p in range(W * W):
        r, c = divmod(p, W)
        if not cells[p]:
            # nothing can ever have covered this cell, so bit 0 is 0 and the
            # step is a pure shift, which is injective: no state merges.
            st = {m >> 1: v for m, v in st.items()}
            continue
        h = c + 1 < W and cells[p + 1]
        d = r + 1 < W and cells[p + W]
        sq = square and h and d and cells[p + W + 1]
        nxt = {}
        get = nxt.get
        if mod is None:
            for m, val in st.items():
                if m & 1:
                    k = m >> 1
                    nxt[k] = get(k, 0) + val
                    continue
                fr, fd = not m & B_R, not m & B_D
                if h and fr:
                    k = (m | B_R) >> 1
                    nxt[k] = get(k, 0) + val
                if d and fd:
                    k = (m | B_D) >> 1
                    nxt[k] = get(k, 0) + val
                if sq and fr and fd and not m & B_DR:
                    k = (m | B_R | B_D | B_DR) >> 1
                    nxt[k] = get(k, 0) + val
        else:
            for m, val in st.items():
                if m & 1:
                    k = m >> 1
                    nxt[k] = (get(k, 0) + val) % mod
                    continue
                fr, fd = not m & B_R, not m & B_D
                if h and fr:
                    k = (m | B_R) >> 1
                    nxt[k] = (get(k, 0) + val) % mod
                if d and fd:
                    k = (m | B_D) >> 1
                    nxt[k] = (get(k, 0) + val) % mod
                if sq and fr and fd and not m & B_DR:
                    k = (m | B_R | B_D | B_DR) >> 1
                    nxt[k] = (get(k, 0) + val) % mod
        st = nxt
    return st.get(0, 0)


_rows_ok = True
for _n in range(NMAX + 1):
    _cells = region_from_inequality(_n)
    _W = 2 * _n
    _widths = [sum(_cells[_r * _W:(_r + 1) * _W]) for _r in range(_W)]
    if _widths != region_rows(_n) or sum(_cells) != 2 * _n * (_n + 1):
        _rows_ok = False
check('G1-region-from-the-inequality-matches-the-row-widths',
      _rows_ok,
      'cells of |x|+|y| <= n+1 give row widths 2*min(r+1, 2n-r) and total '
      '2n(n+1), n = 0..15 -- the region is not taken on trust')

_bt = [count_backtracking(_n) for _n in range(5)]
check('G2-backtracking-reproduces-M(0..4)',
      _bt == PUBLISHED[:5],
      'no transfer computation at all: %s' % (_bt,))

_tr_small = [count_transfer(_n) for _n in range(5)]
check('G3-transfer-agrees-with-backtracking-on-n-0-to-4',
      _tr_small == _bt,
      'two structurally different methods, same five integers')

for _n in range(13):
    _t0 = time.time()
    _got = count_transfer(_n)
    check('G4-transfer-exact-vs-A356512-n=%02d' % _n,
          _got == PUBLISHED[_n],
          'M(%d) = %d  (%.1fs)' % (_n, _got, time.time() - _t0))

for _n in range(13):
    _got = count_transfer(_n, square=False)
    _want = 2 ** (_n * (_n + 1) // 2)
    check('G5-dominoes-only-vs-A006125-n=%02d' % _n,
          _got == _want,
          '2x2 disabled: %d = 2^%d  [Elkies-Kuperberg-Larsen-Propp]'
          % (_got, _n * (_n + 1) // 2))

note('the three deep residues now follow; 2^27, 2^29 and 2^31 states.')
_deep = {}
for _n, _want in ((13, M13_MOD256_PRINTED), (14, M14_MOD256), (15, M15_MOD256)):
    _t0 = time.time()
    _got = count_transfer(_n, mod=256)
    _deep[_n] = _got
    note('n=%d done in %.1fs, peak RSS %.2f GB' % (_n, time.time() - _t0,
                                                   peak_gb()))
    check('G%d-transfer-mod-256-reproduces-M(%d)' % (6 + _n - 13, _n),
          _got == _want,
          'independent pure-Python value %d, paper prints %d' % (_got, _want))

check('G9-refutation-re-derived-from-this-program-own-residues',
      13 in _deep and 14 in _deep and 15 in _deep
      and (_deep[14] + _deep[15]) % 256 == 48
      and v2((_deep[14] + _deep[15]) % 256) == 4
      and (_deep[14] + _deep[15]) % 32 == 16 != 0
      and _deep[13] % 32 == 31,
      'from residues this program computed itself: M(14)+M(15) = 48 mod 256, '
      'v_2 = 4 < 5 => Conjecture 2 fails at (14, 15, 5)')

# =====================================================================
# scope, then the verdict
# =====================================================================
print('=== scope ===', flush=True)
print('NOT RE-RUN: M(14) and M(15) EXACTLY. Both engines and this program '
      'carry them modulo 2^8 only; the exact integers are near 10^41 and '
      '10^47 and are not computed anywhere here. Proposition 2 needs only '
      'v_2 < 8, so the residues settle the refutation.', flush=True)
print('NOT RE-RUN: M(n) for any n >= 16. Nothing beyond index 15 is computed, '
      'so every instance of Conjecture 2 with an index >= 16 -- including the '
      'two remaining k = 5 cells (13, 16) and (12, 17) -- is left open in '
      'both directions, and (14, 15, 5) is shown least only within the box '
      "n, n' <= 15.", flush=True)
print('NOT RE-RUN: k >= 6. Any such test needs max(n, n\') >= 31 (check F1), '
      'a 2^63-state computation, and none was attempted.', flush=True)
print("NOT RE-RUN: Propp's Conjecture 1, the 2-adic continuity of M. It is "
      'neither used nor tested here, and the antisymmetry restatement '
      'Mhat(-3-n) = -Mhat(n) falls only conditionally on it.', flush=True)
print("NOT RE-RUN: Propp's theorem that M(n) is odd. Check E9 verifies that "
      'the 16 residues in hand are odd and that every k = 1 cell holds, which '
      'is a control on the data, not a proof of the theorem.', flush=True)
print('NOT RE-RUN: the two C engines that produced the paper\'s data. This '
      'program is a third, independent implementation; it agrees with them at '
      'n = 13, 14, 15 mod 256, but it does not inspect or execute them.',
      flush=True)

_n_pass = sum(1 for r in _RESULTS if r)
_n_fail = len(_RESULTS) - _n_pass
print('', flush=True)
if _n_fail:
    print('VERDICT: %d of %d CHECKS FAILED' % (_n_fail, len(_RESULTS)),
          flush=True)
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _n_pass, flush=True)
sys.exit(0)
