#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Referee checker for `Two-Circulant Butson Hadamard Matrices of Orders 46 and 74`.

Python 3.9+, STANDARD LIBRARY ONLY, no external data file: every object it consumes is
transcribed below from the paper itself.  All arithmetic is exact integer arithmetic in
Z[zeta_6] = Z[z]/(z^2 - z + 1); there is no floating point anywhere and no decision is
taken on an inexact quantity.

Prints one `PASS <name> [detail]` line per check and exits 0 iff every check passed.
"""

import sys
from itertools import product

# ---------------------------------------------------------------------------
# The objects, transcribed from Section 3 of paper.tex.
# ---------------------------------------------------------------------------
W1 = {
    'name': 'W1_BH46_split_3_43',
    'claim': 'BH(46,6)', 'p': 23, 'k': 11, 'g': 5, 'G0': [1, 22],
    'a': [1, 5, 5, 2, 4, 4, 3, 2, 2, 1, 0],
    'b': [2, 1, 4, 3, 2, 1, 5, 3, 3, 0, 2],
    'e_x': [0, 1, 5, 4, 4, 5, 2, 2, 3, 0, 2, 1, 1, 2, 0, 3, 2, 2, 5, 4, 4, 5, 1],
    'e_y': [0, 2, 4, 1, 2, 1, 3, 3, 5, 2, 3, 0, 0, 3, 2, 5, 3, 3, 1, 2, 1, 4, 2],
    'nx': 3, 'ny': 43,
}
W2 = {
    'name': 'W2_BH46_split_9_37',
    'claim': 'BH(46,6)', 'p': 23, 'k': 11, 'g': None, 'G0': [1, 22],
    'a': None, 'b': None,
    'e_x': [0, 1, 2, 3, 4, 1, 4, 3, 4, 1, 4, 0, 0, 4, 1, 4, 3, 4, 1, 4, 3, 2, 1],
    'e_y': [0, 5, 3, 3, 5, 4, 1, 0, 1, 0, 4, 0, 0, 4, 0, 1, 0, 1, 4, 5, 3, 3, 5],
    'nx': 9, 'ny': 37,
}
W3 = {
    'name': 'W3_BH74',
    'claim': 'BH(74,6)', 'p': 37, 'k': 12, 'g': 2, 'G0': [1, 10, 26],
    'a': [1, 4, 2, 0, 4, 3, 4, 3, 2, 5, 5, 0],
    'b': [1, 5, 4, 0, 0, 2, 1, 2, 4, 2, 0, 3],
    'e_x': [0, 1, 4, 2, 2, 0, 0, 2, 0, 4, 1, 4, 4, 0, 5, 4, 4, 3, 3, 0, 4, 5, 3,
            0, 3, 5, 1, 4, 5, 5, 2, 5, 3, 2, 2, 3, 4],
    'e_y': [0, 1, 5, 4, 4, 3, 0, 4, 0, 0, 1, 1, 0, 3, 2, 5, 0, 2, 2, 3, 5, 0, 2,
            0, 2, 0, 1, 1, 0, 2, 4, 2, 2, 4, 4, 2, 1],
    'nx': 31, 'ny': 43,
}
WITNESSES = [W1, W2, W3]

# other printed quantities
W2_TWO_GAMMA_X = [-4, 2, -4, 6, -8, -2, 0, -16, 2, 6, 4]      # 2*gamma_s(x), s = 1..11
W3_GAMMA_1 = (-4, 1)                                          # -4 + zeta_6
W3_GAMMA_14 = (6, -7)                                         # 6 - 7 zeta_6
ODD_LOESCHIAN_BELOW_46 = [1, 3, 7, 9, 13, 19, 21, 25, 27, 31, 37, 39, 43]
SPLITS_46 = [(3, 43), (7, 39), (9, 37), (19, 27), (21, 25)]
INDEX1_NORMS_23 = [529, 507, 463, 441, 463, 507]
INDEX1_NORMS_47 = [2209, 2163, 2071, 2025, 2071, 2163]

# ---------------------------------------------------------------------------
# Exact Z[zeta_6] arithmetic.  An element is a pair (A, B) meaning A + B*z, z^2 = z - 1.
# ---------------------------------------------------------------------------
ZPOW = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]    # z^m, m = 0..5

def add(u, v):
    return (u[0] + v[0], u[1] + v[1])

def neg(u):
    return (-u[0], -u[1])

def mul(u, v):
    a1, b1 = u
    a2, b2 = v
    return (a1 * a2 - b1 * b2, a1 * b2 + a2 * b1 + b1 * b2)

def conj(u):
    # conj(z) = z^5 = 1 - z, so conj(A + Bz) = (A + B) - Bz
    return (u[0] + u[1], -u[1])

def norm(u):
    return u[0] * u[0] + u[0] * u[1] + u[1] * u[1]

def reduce_counts(cnt):
    """sum_m cnt[m] * z^m, via the paper's A = n0-n2-n3+n5, B = n1+n2-n4-n5."""
    return (cnt[0] - cnt[2] - cnt[3] + cnt[5], cnt[1] + cnt[2] - cnt[4] - cnt[5])

def is_loeschian(n):
    """n = A^2 + AB + B^2 for some integers A, B, by the classical criterion: every prime
    congruent to 2 mod 3 occurs to an even power.  Trial division, exact."""
    if n < 0:
        return False
    if n == 0:
        return True
    m, d = n, 2
    while d * d <= m:
        if m % d == 0:
            e = 0
            while m % d == 0:
                m //= d
                e += 1
            if d % 3 == 2 and e % 2 == 1:
                return False
        d += 1
    if m > 1 and m % 3 == 2:
        return False
    return True

def loeschian_by_search(n):
    """Independent, brute-force witness for is_loeschian on small n."""
    r = 0
    while r * r <= 4 * n:
        r += 1
    for A in range(-r, r + 1):
        for B in range(-r, r + 1):
            if A * A + A * B + B * B == n:
                return True
    return False

# ---------------------------------------------------------------------------
# The construction of Section 2.
# ---------------------------------------------------------------------------
def exponent_matrix(p, ex, ey):
    """H[i][j] = zeta_6^{E[i][j]} for H = [[X, Y], [Y*, -X*]], equation (expmat)."""
    n = 2 * p
    E = [[0] * n for _ in range(n)]
    for i in range(p):
        for j in range(p):
            E[i][j] = ex[(j - i) % p]
            E[i][p + j] = ey[(j - i) % p]
            E[p + i][j] = (-ey[(i - j) % p]) % 6
            E[p + i][p + j] = (3 - ex[(i - j) % p]) % 6
    return E

def block_matrix_over_Z6(p, ex, ey):
    """The same H, built instead as the literal block matrix [[X, Y], [Y*, -X*]] with
    entries as elements of Z[zeta_6].  Used only to confirm equation (expmat)."""
    n = 2 * p
    X = [[ZPOW[ex[(c - r) % p]] for c in range(p)] for r in range(p)]
    Y = [[ZPOW[ey[(c - r) % p]] for c in range(p)] for r in range(p)]
    Xs = [[conj(X[c][r]) for c in range(p)] for r in range(p)]      # X*
    Ys = [[conj(Y[c][r]) for c in range(p)] for r in range(p)]      # Y*
    H = [[None] * n for _ in range(n)]
    for r in range(p):
        for c in range(p):
            H[r][c] = X[r][c]
            H[r][p + c] = Y[r][c]
            H[p + r][c] = Ys[r][c]
            H[p + r][p + c] = neg(Xs[r][c])
    return H

def gamma(p, e, s):
    """gamma_s = sum_i x_i conj(x_{i+s}) = sum_i zeta^{e_i - e_{i+s}}."""
    cnt = [0] * 6
    for i in range(p):
        cnt[(e[i] - e[(i + s) % p]) % 6] += 1
    return reduce_counts(cnt)

def gamma_violations(p, ex, ey):
    return [s for s in range(1, p) if add(gamma(p, ex, s), gamma(p, ey, s)) != (0, 0)]

def gram_violations(E, n):
    """(bad, npairs) over all i <= j: <H_i, H_j> must be n on the diagonal, 0 off it."""
    bad, npairs = [], 0
    for i in range(n):
        Ei = E[i]
        for j in range(i, n):
            Ej = E[j]
            cnt = [0] * 6
            for t in range(n):
                cnt[(Ei[t] - Ej[t]) % 6] += 1
            v = reduce_counts(cnt)
            npairs += 1
            if v != ((n, 0) if i == j else (0, 0)):
                bad.append((i, j, v))
    return bad, npairs

def rowsum_norm(e):
    cnt = [0] * 6
    for m in e:
        cnt[m] += 1
    return norm(reduce_counts(cnt))

def subgroup(p, g, k):
    return sorted({pow(g, k * t, p) for t in range(p)})

def cosets(p, g, G0):
    out, seen = [], set()
    for j in range(p - 1):
        base = pow(g, j, p)
        if base in seen:
            continue
        cs = sorted({(base * h) % p for h in G0})
        seen.update(cs)
        out.append(cs)
    return out

def expand_coset_code(p, g, G0, code):
    e = [None] * p
    e[0] = 0
    for idx, cs in enumerate(cosets(p, g, G0)):
        for m in cs:
            e[m] = code[idx]
    return e

def constant_on_cosets(p, e, G0):
    """e is simple of index k iff e[0] = 0 and e is constant on each coset of G0."""
    if e[0] != 0:
        return False
    for i in range(1, p):
        for h in G0:
            if e[(i * h) % p] != e[i]:
                return False
    return True

# ---------------------------------------------------------------------------
# The check recorder.
# ---------------------------------------------------------------------------
_LINES = []
_FAILED = 0

def check(ok, name, detail=''):
    global _FAILED
    if ok:
        _LINES.append('PASS %s%s' % (name, ('  ' + detail) if detail else ''))
    else:
        _FAILED += 1
        _LINES.append('FAIL %s%s' % (name, ('  ' + detail) if detail else ''))

def emit(text=''):
    _LINES.append(text)

# ---------------------------------------------------------------------------
# 1. Sanity of the arithmetic layer itself (a checker nobody checked is an assertion).
# ---------------------------------------------------------------------------
emit('--- the exact Z[zeta_6] layer ---')
check(all(mul(ZPOW[i], ZPOW[j]) == ZPOW[(i + j) % 6] for i in range(6) for j in range(6)),
      'zeta6-multiplication-table', 'z^i z^j = z^(i+j mod 6) for all 36 pairs')
check(all(mul(ZPOW[i], conj(ZPOW[i])) == (1, 0) for i in range(6)),
      'zeta6-unit-modulus', 'z^i conj(z^i) = 1 for i = 0..5')
check(mul((0, 1), (0, 1)) == (-1, 1) and norm((1, 1)) == 3,
      'zeta6-minimal-polynomial', 'z^2 = z - 1 and |1+z|^2 = 3')
check(all(reduce_counts([1 if m == i else 0 for m in range(6)]) == ZPOW[i] for i in range(6)),
      'reduce-counts-agrees-with-zeta-powers', 'equation (red) on the six unit vectors')
check(all(is_loeschian(n) == loeschian_by_search(n) for n in range(0, 200)),
      'loeschian-criterion-vs-brute-force', 'agree on every n in 0..199')

# ---------------------------------------------------------------------------
# 2. The three witnesses.
# ---------------------------------------------------------------------------
for W in WITNESSES:
    p, ex, ey, nm = W['p'], W['e_x'], W['e_y'], W['name']
    n = 2 * p
    emit('--- %s: %s, p = %d, k = %d ---' % (nm, W['claim'], p, W['k']))

    ok = (len(ex) == len(ey) == p and ex[0] == 0 and ey[0] == 0
          and all(0 <= v <= 5 for v in ex + ey))
    check(ok, 'shape/%s' % nm, 'both vectors have length %d, first entry 0, entries in Z_6' % p)

    E = exponent_matrix(p, ex, ey)
    check(all(0 <= E[i][j] <= 5 for i in range(n) for j in range(n)),
          'entries-in-Omega6/%s' % nm, 'all %d entries of H are sixth roots of unity' % (n * n))

    B = block_matrix_over_Z6(p, ex, ey)
    check(all(B[i][j] == ZPOW[E[i][j]] for i in range(n) for j in range(n)),
          'block-form-matches-exponent-matrix/%s' % nm,
          'equation (expmat) reproduces [[X,Y],[Y*,-X*]] entrywise')

    gv = gamma_violations(p, ex, ey)
    check(gv == [], 'autocorrelation/%s' % nm,
          'gamma_s(x) + gamma_s(y) = 0 for all s = 1..%d, %d violations' % (p - 1, len(gv)))

    bad, npairs = gram_violations(E, n)
    check(bad == [], 'gram/%s' % nm,
          'H H* = %d I_%d exactly: %d inner products with i <= j, %d violations'
          % (n, n, npairs, len(bad)))

    nx, ny = rowsum_norm(ex), rowsum_norm(ey)
    check((nx, ny) == (W['nx'], W['ny']) and nx + ny == n
          and is_loeschian(nx) and is_loeschian(ny),
          'row-sum-split/%s' % nm,
          '|sum x|^2 = %d, |sum y|^2 = %d, sum = %d = 2p, both Loeschian' % (nx, ny, nx + ny))

    check(constant_on_cosets(p, ex, W['G0']) and constant_on_cosets(p, ey, W['G0']),
          'simple-of-index-k/%s' % nm,
          'both first rows are constant on the cosets of G_0 = %s, so simple of index %d'
          % (W['G0'], W['k']))

    if W['g'] is not None:
        G0 = subgroup(p, W['g'], W['k'])
        check(G0 == sorted(W['G0']),
              'subgroup-rederived/%s' % nm,
              'G_0 = {h^%d mod %d} = %s, as printed' % (W['k'], p, G0))
        cs = cosets(p, W['g'], G0)
        check(len(cs) == W['k'] and sorted(x for c in cs for x in c) == list(range(1, p)),
              'cosets-partition/%s' % nm,
              '%d cosets of G_0 partition Z_%d^*' % (len(cs), p))
        check(expand_coset_code(p, W['g'], G0, W['a']) == ex
              and expand_coset_code(p, W['g'], G0, W['b']) == ey,
              'coset-code-expands/%s' % nm,
              'a and b expand over the cosets of G_0 to exactly the printed e_x, e_y')
        # Lemma coset, first half: gamma_s depends only on the coset s G_0
        okc = True
        for c in cs:
            gs = {gamma(p, ex, s) for s in c}
            if len(gs) != 1:
                okc = False
        check(okc, 'gamma-constant-on-cosets/%s' % nm,
              'gamma_s(x) takes one value on each of the %d cosets (Lemma coset)' % len(cs))

# ---------------------------------------------------------------------------
# 3. The gamma values quoted in the paper.
# ---------------------------------------------------------------------------
emit('--- the gamma values printed in Section 3 ---')
p = 23
gx = [gamma(p, W2['e_x'], s) for s in range(1, p)]
gy = [gamma(p, W2['e_y'], s) for s in range(1, p)]
check(all(v[1] == 0 for v in gx) and all(v[1] == 0 for v in gy),
      'gamma-real-at-k11', 'every gamma_s is a rational integer at p = 23, k = 11 '
                           '(-1 in G_0, Lemma coset)')
check(all(gx[s - 1] == gx[p - s - 1] for s in range(1, p)),
      'gamma-palindromic-at-k11', 'gamma_s = gamma_{-s} for all s = 1..22')
check([2 * gx[s - 1][0] for s in range(1, 12)] == W2_TWO_GAMMA_X,
      'W2-gamma-vector', '(2 gamma_s(x))_{s=1..11} = %s, as printed' % W2_TWO_GAMMA_X)
check([2 * gy[s - 1][0] for s in range(1, 12)] == [-v for v in W2_TWO_GAMMA_X],
      'W2-gamma-vector-negated', '(2 gamma_s(y))_{s=1..11} is its negative')
check(sum(W2_TWO_GAMMA_X) == W2['nx'] - p,
      'W2-lemma-sum-identity',
      '2 sum_{s=1..11} gamma_s(x) = %d = |sum x|^2 - p = %d - %d'
      % (sum(W2_TWO_GAMMA_X), W2['nx'], p))
g1 = gamma(37, W3['e_x'], 1)
g14 = gamma(37, W3['e_x'], 14)
check(g1 == W3_GAMMA_1 and g14 == W3_GAMMA_14,
      'W3-gamma-values', 'gamma_1(x) = -4 + z and gamma_14(x) = 6 - 7z, as printed')
check(any(gamma(37, W3['e_x'], s)[1] != 0 for s in range(1, 37)) and 36 % 2 == 0
      and (-1) % 37 not in W3['G0'],
      'W3-gamma-genuinely-non-real',
      '-1 = 36 is not in G_0 = %s at p = 37, and some gamma_s(x) has nonzero z-part'
      % W3['G0'])

# ---------------------------------------------------------------------------
# 4. The branch enumeration of Section 2.
# ---------------------------------------------------------------------------
emit('--- the Loeschian branch enumeration ---')
check(not is_loeschian(23) and not loeschian_by_search(23),
      'no-circulant-BH-23-6', '23 is prime and 2 mod 3, hence not Loeschian, so '
                              '|sum x|^2 = 23 is impossible')
check(is_loeschian(37) and loeschian_by_search(37),
      '37-is-Loeschian', '37 = 3^2 + 3*4 + 4^2, so the same argument does NOT apply at p = 37')
check([n for n in range(1, 46) if n % 2 == 1 and is_loeschian(n)] == ODD_LOESCHIAN_BELOW_46,
      'odd-Loeschian-below-46', 'exactly %s' % ODD_LOESCHIAN_BELOW_46)
check(sorted((u, 46 - u) for u in ODD_LOESCHIAN_BELOW_46
             if u < 46 - u and is_loeschian(46 - u) and (46 - u) % 2 == 1) == SPLITS_46,
      'admissible-splits-of-46', 'exactly the five branches %s' % (SPLITS_46,))
check(tuple(sorted((W1['nx'], W1['ny']))) in SPLITS_46
      and tuple(sorted((W2['nx'], W2['ny']))) in SPLITS_46
      and sorted((W1['nx'], W1['ny'])) != sorted((W2['nx'], W2['ny'])),
      'W1-and-W2-on-distinct-branches',
      '{%d,%d} and {%d,%d} are two different members of the five'
      % (W1['nx'], W1['ny'], W2['nx'], W2['ny']))
check(W1['e_x'] != W2['e_x'] and W1['e_x'][2] == 5 and W2['e_x'][2] == 2,
      'W1-and-W2-differ', 'the two first rows already differ at coordinate 2 (5 vs 2)')
check(sorted(k for k in range(1, 23) if 22 % k == 0) == [1, 2, 11, 22],
      'admissible-indices-p23', 'k | 22 gives k in {1, 2, 11, 22}')
check(36 % 12 == 0 and sorted(k for k in range(1, 47) if 46 % k == 0) == [1, 2, 23, 46],
      'admissible-indices-p37-p47', 'k = 12 divides 36; at p = 47, k | 46 gives {1, 2, 23, 46}')

# ---------------------------------------------------------------------------
# 5. Proposition (index 1 and index 2 are empty at p = 23 and p = 47).
# ---------------------------------------------------------------------------
emit('--- the exhaustive index-<=2 sweeps (Proposition) ---')

def index1_norms(p):
    out = []
    for a in range(6):
        e = [0] + [a] * (p - 1)
        out.append(rowsum_norm(e))
    return out

check(index1_norms(23) == INDEX1_NORMS_23, 'index1-norms-p23',
      '|1 + 22 z^a|^2 = %s, all > 46' % INDEX1_NORMS_23)
check(index1_norms(47) == INDEX1_NORMS_47, 'index1-norms-p47',
      '|1 + 46 z^a|^2 = %s, all > 94' % INDEX1_NORMS_47)
check(sorted({norm(add(ZPOW[i], ZPOW[j])) for i in range(6) for j in range(6)}) == [0, 1, 3, 4],
      'two-sixth-roots-moduli', '|z^a + z^b|^2 in {0, 1, 3, 4}, so the sum is 0 or has '
                                'modulus at least 1')

def index2_vectors(p):
    """All first rows simple of index 2 (this family contains the index-1 rows)."""
    qr = {(t * t) % p for t in range(1, p)}
    out = []
    for a1, a2 in product(range(6), repeat=2):
        e = [0] * p
        for i in range(1, p):
            e[i] = a1 if i in qr else a2
        out.append(e)
    return out

for p in (23, 47):
    vs = index2_vectors(p)
    gvecs = [tuple(gamma(p, e, s) for s in range(1, p)) for e in vs]
    hits = 0
    for i in range(len(gvecs)):
        for j in range(len(gvecs)):
            if all(add(gvecs[i][t], gvecs[j][t]) == (0, 0) for t in range(p - 1)):
                hits += 1
    nord = len(vs) * len(vs)
    nun = len(vs) * (len(vs) + 1) // 2
    check(hits == 0, 'index2-sweep-p%d' % p,
          'no two-circulant BH(%d,6) with both first rows simple of index <= 2: '
          '0 hits over %d ordered (%d unordered) pairs' % (2 * p, nord, nun))
    ns = sorted({rowsum_norm(e) for e in vs})
    bound = ((p - 1) // 2 - 1) ** 2
    check(all(v == 1 or v >= bound for v in ns), 'index2-norm-gap-p%d' % p,
          'every index-<=2 row-sum norm is 1 or at least %d; observed set %s'
          % (bound, ns))

check(not is_loeschian(45), '45-not-Loeschian',
      'so the split {1, 45} of 46 is impossible for any first rows at all')
check(is_loeschian(93) and 93 not in {rowsum_norm(e) for e in index2_vectors(47)},
      '93-Loeschian-but-unreachable',
      '93 = 4^2 + 4*7 + 7^2 IS Loeschian, so the split {1, 93} of 94 is not excluded by '
      'the Loeschian test; it is excluded because no index-<=2 row attains 93')

# ---------------------------------------------------------------------------
# 6. Controls, both polarities.
# ---------------------------------------------------------------------------
emit('--- controls ---')
DFT6 = [[(r * c) % 6 for c in range(6)] for r in range(6)]
bad, npairs = gram_violations(DFT6, 6)
check(bad == [], 'control-positive-DFT6',
      'the published BH(6,6) H[r][c] = z^{rc} passes the same Gram checker: '
      '%d inner products, 0 violations' % npairs)

BADDFT = [row[:] for row in DFT6]
BADDFT[0][0] = 1
bad2, _ = gram_violations(BADDFT, 6)
check(len(bad2) > 0 and ((0, 1) in [(i, j) for i, j, _ in bad2]),
      'control-negative-DFT6-corrupted',
      'corrupting H[0][0] to z breaks it: %d violations, including the (0,1) inner '
      'product, which by hand is z - 1 != 0' % len(bad2))
check([v for i, j, v in bad2 if (i, j) == (0, 1)] == [add((0, 1), (-1, 0))],
      'control-negative-value-is-z-minus-1',
      'and that inner product equals exactly z - 1 = %s' % (add((0, 1), (-1, 0)),))

pert = W1['e_x'][:]
pert[7] = (pert[7] + 1) % 6
gvp = gamma_violations(23, pert, W1['e_y'])
Ep = exponent_matrix(23, pert, W1['e_y'])
badp, _ = gram_violations(Ep, 46)
check(len(gvp) > 0 and len(badp) > 0,
      'control-negative-perturbed-witness',
      'moving one exponent of W1 by +1 produces %d autocorrelation violations and %d '
      'Gram violations, so neither detector is vacuous' % (len(gvp), len(badp)))
check((len(gvp) == 0) == (len(badp) == 0),
      'lemma-reduce-consistency',
      'the autocorrelation test and the full Gram test agree on the perturbed vector, '
      'as Lemma reduce requires')

# ---------------------------------------------------------------------------
# Scope, and the verdict.
# ---------------------------------------------------------------------------
npass = sum(1 for L in _LINES if L.startswith('PASS '))
nfail = sum(1 for L in _LINES if L.startswith('FAIL '))
for L in _LINES:
    print(L)
print()
print('NOT RE-RUN: the two lane censuses. This program does NOT re-derive the search '
      'programs\' complete-count claims -- 10,120 ordered / 5,060 unordered complementary '
      'pairs over all 6^11 = 362,797,056 first rows per side at (p,k) = (23,11), and '
      '5,904 ordered / 2,952 unordered over 6^12 = 2,176,782,336 at (37,12) -- nor the '
      'shard-cover accounting of the 128-thread p = 37 run, nor the searchers\' own '
      'forced-positive controls at p = 3,5,7,11,13,17 and p = 31 (k=6), nor the checksum '
      'of the C engine. Those rest on the search programs\' counters; their standard '
      'output was not preserved, and reproducing them needs 456 s single-threaded at '
      'p = 23 and 903 s on 128 threads at p = 37. Nothing above depends on them: '
      'existence needs one object, and the objects are the printed exponent vectors.')
print('NOT RE-RUN: minimality, or any classification of BH(46,6)/BH(74,6). Nothing here '
      'says these are the smallest or the only such matrices. The only separation '
      'established between W1 and W2 is that they sit on different Loeschian branches, so '
      'no symmetry of the two-circulant ansatz carries one to the other; whether H(W1) '
      'and H(W2) are inequivalent as Butson Hadamard matrices, under the far larger group '
      'of row/column permutations and sixth-root scalings, is NOT tested. The k = 22 lane '
      'at p = 23 (6^22 vectors) and the k = 23, 46 lanes at p = 47 were not searched, and '
      'BH(86,6) at p = 43 was not searched at all.')
print('NOT RE-RUN: the literature search. This program checks mathematics, not novelty; '
      'the unread bibliographic channels are named in the paper.')
print()
if nfail:
    print('VERDICT: %d of %d CHECKS FAILED' % (nfail, npass + nfail))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % npass)
sys.exit(0)
