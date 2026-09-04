#!/usr/bin/env python3
"""verify.py -- checker for

    "Four Pairwise Orthogonal Semigroup Operations on a Four-Element Set,
     and the 160 Extremal Families"

Python 3.9+, STANDARD LIBRARY ONLY (no numpy/sympy/networkx, no external data file).
Exact integer arithmetic throughout; no floating-point value enters any decision.

INPUT.  The four Cayley tables printed in Section 2 of the paper, in BOTH of the forms the paper
prints them (the 4x4 arrays and the flat row-concatenated strings), and the ten extremal-family
representatives printed in the table of Section 4.  Nothing else: no file, no artifact, no data.

WHAT IT RE-DERIVES.
 (A) the exhibited family: associativity, value-balance, all six pairwise orthogonalities, the four
     self-pairs (which must NOT be orthogonal), the Rees-Sushkevich type of each member, and the
     identification of f4 with the table published in the source paper;
 (B) the F_2 picture that makes the lower bound a hand theorem, including the classification of all
     sixteen linear maps into associative and non-associative, and the tightness of the rank-1
     idempotent mechanism;
 (C) the full unpruned census for n = 1,2,3,4: every associative table on a fixed n-set is generated
     from scratch with no lemma and no symmetry reduction, checked against OEIS A023814 (labelled)
     and A027851 (up to isomorphism); the orthogonality graph is built and its cliques counted twice
     at n = 4 -- once level by level, once by an algorithm-free scan of all C(48,4) and C(48,5)
     subsets -- so the load-bearing negative "there is no 5-clique" does not depend on any clique
     algorithm being correct;
 (D) the classification of the 160 extremal families into 10 orbits under the diagonal action of
     Sym(4), and every structural statement the paper makes about them.

The balance criterion (an orthogonal partner forces value-balance) and the equivalence
(balanced <=> completely simple) are not assumed: the balance criterion is instantiated as an
explicit certificate for each unbalanced table, and the equivalence is checked as a set identity in
both directions over the whole census, for n <= 4 only.

The closing SCOPE block states what is not covered.  Exit status 0 iff every check passed.
"""
import itertools
import random
import sys
import time
from collections import Counter

_OK = []
_BAD = []


def ck(name, cond, detail=''):
    """One check.  `PASS <name> <detail>` or `FAIL <name> <detail>`; no other line begins with
    either word, so the PASS count of the transcript is the number of checks."""
    (_OK if cond else _BAD).append(name)
    print('%s %s %s' % ('PASS' if cond else 'FAIL', name, detail))


def note(text):
    print('NOTE ' + text)


# ---------------------------------------------------------------------------
# 1. THE OBJECTS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------
# Section 2, flat row-concatenated form, 1-based symbols on the fixed set {1,2,3,4}.
PAPER_FLAT = (('f1', '1111222233334444'), ('f2', '1234123434123412'),
              ('f3', '1234214321431234'), ('f4', '1234432112344321'))
# Section 2, the same four operations as printed 4x4 arrays (row a lists a*1,...,a*4).
PAPER_ARRAY = {
    'f1': ((1, 1, 1, 1), (2, 2, 2, 2), (3, 3, 3, 3), (4, 4, 4, 4)),
    'f2': ((1, 2, 3, 4), (1, 2, 3, 4), (3, 4, 1, 2), (3, 4, 1, 2)),
    'f3': ((1, 2, 3, 4), (2, 1, 4, 3), (2, 1, 4, 3), (1, 2, 3, 4)),
    'f4': ((1, 2, 3, 4), (4, 3, 2, 1), (1, 2, 3, 4), (4, 3, 2, 1)),
}
FLAT_L = '1111222233334444'          # L(x,y) = x, the left-zero band
FLAT_R = '1234123412341234'          # R(x,y) = y, the right-zero band
FLAT_KLEIN = '1234214334124321'      # the Klein table 1234|2143|3412|4321

# Section 4: the ten extremal-family representatives, with the orbit size the paper prints.
PAPER_ORBITS = (
    (24, ('1133224411332244', '1212212134344343', '1234123443214321', '1234341234121234')),
    (24, ('1133224411332244', '1212212134344343', '1234123443214321', '3412123412343412')),
    (24, ('1133224411332244', '1221122143344334', '1234341234121234', '1414232332324141')),
    (24, ('1133224411332244', '1221211234434334', '1234123434123412', '1414232332324141')),
    (24, ('1133224411332244', '1221211234434334', '1234123434123412', '4141323223231414')),
    (8, ('1111222233334444', '1234123434123412', '1234214321431234', '1234432112344321')),
    (8, ('1111222233334444', '1234123434123412', '1234214321431234', '4321123443211234')),
    (8, ('1133224411332244', '1221122143344334', '1234214334124321', '1414323232321414')),
    (8, ('1133224433114422', '1221211234434334', '1234123412341234', '1414232332324141')),
    (8, ('1133224433114422', '1221211234434334', '1234123412341234', '4141323223231414')),
)

A023814 = (1, 1, 8, 113, 3492)       # labelled associative operations on an n-set, n = 0..4
A027851 = (1, 1, 5, 24, 188)         # semigroups of order n up to isomorphism, n = 0..4
EXPECT_BALANCED = {1: 1, 2: 4, 3: 5, 4: 48}
EXPECT_EDGES = {1: 0, 2: 5, 3: 7, 4: 354}
EXPECT_CLIQUES = {1: {1: 1}, 2: {1: 4, 2: 5, 3: 2, 4: 0},
                  3: {1: 5, 2: 7, 3: 3, 4: 0}, 4: {1: 48, 2: 354, 3: 542, 4: 160, 5: 0}}
EXPECT_MAX = {1: 1, 2: 3, 3: 3, 4: 4}
EXPECT_EXTREMAL = {1: 1, 2: 2, 3: 3, 4: 160}
EXPECT_ORBITS = {2: 1, 3: 1, 4: 10}
EXPECT_ORBIT_SIZES_4 = [8, 8, 8, 8, 8, 24, 24, 24, 24, 24]
EXPECT_TYPES_4 = {(1, 1, 4): 1, (2, 1, 2): 6, (4, 1, 1): 1,
                  (1, 2, 2): 12, (2, 2, 1): 12, (1, 4, 1): 16}


# ---------------------------------------------------------------------------
# 2. TABLES
# ---------------------------------------------------------------------------
def decode(flat):
    """Flat 1-based string of n*n symbols -> 0-based row-major tuple with T[a*n+b] = a*b."""
    n = 1
    while n * n < len(flat):
        n += 1
    assert n * n == len(flat)
    t = tuple(int(c) - 1 for c in flat)
    assert all(0 <= v < n for v in t)
    return t


def encode(t):
    return ''.join(str(v + 1) for v in t)


def assoc(t, n):
    for x in range(n):
        xn = x * n
        for y in range(n):
            u = t[xn + y] * n
            yn = y * n
            for z in range(n):
                if t[u + z] != t[xn + t[yn + z]]:
                    return False
    return True


def balanced(t, n):
    c = Counter(t)
    return len(c) == n and all(v == n for v in c.values())


def orth(t, u, n):
    return len(set(zip(t, u))) == n * n


def idempotents(t, n):
    return frozenset(a for a in range(n) if t[a * n + a] == a)


def is_quasigroup(t, n):
    for a in range(n):
        if len(set(t[a * n:a * n + n])) != n:
            return False
        if len(set(t[b * n + a] for b in range(n))) != n:
            return False
    return True


def simple(t, n):
    """S a S = S for every a in S.  For a FINITE semigroup, simple <=> completely simple."""
    full = frozenset(range(n))
    for a in range(n):
        two = set()
        for s in set(t[x * n + a] for x in range(n)):
            sn = s * n
            for y in range(n):
                two.add(t[sn + y])
        if frozenset(two) != full:
            return False
    return True


def classes(t, n, right):
    """Partition of S by principal right ideal (right=True) or principal left ideal."""
    key = {}
    for a in range(n):
        if right:
            s = frozenset(set(t[a * n + b] for b in range(n)) | {a})
        else:
            s = frozenset(set(t[b * n + a] for b in range(n)) | {a})
        key.setdefault(s, []).append(a)
    return list(key.values())


def rees_type(t, n):
    """(|I|, |G|, |Lambda|): |I| R-classes, |Lambda| L-classes, structure group of order
    n / (|I| |Lambda|).  Meaningful for a completely simple operation."""
    i = len(classes(t, n, True))
    lam = len(classes(t, n, False))
    return (i, n // (i * lam), lam)


def group_orders(t, n):
    """Sorted element orders of e S e for the least idempotent e -- the structure group."""
    e = min(idempotents(t, n))
    h = sorted(set(t[e * n + t[x * n + e]] for x in range(n)))
    out = []
    for x in h:
        y, k = x, 1
        while y != e and k <= n:
            y = t[y * n + x]
            k += 1
        out.append(k)
    return tuple(sorted(out))


def apply_perm(t, s, n):
    out = [0] * (n * n)
    for a in range(n):
        sa = s[a] * n
        an = a * n
        for b in range(n):
            out[sa + s[b]] = s[t[an + b]]
    return tuple(out)


# ---------------------------------------------------------------------------
# 3. THE UNPRUNED CENSUS
# ---------------------------------------------------------------------------
def all_assoc(n):
    """EVERY associative table on the fixed set {0,...,n-1}: plain backtracking over the n*n cells
    in row-major order with a FULL n^3 associativity rescan at every node.  No lemma, no balance
    pruning, no symmetry reduction -- the vertex set is the whole of S_n."""
    N = n * n
    t = [-1] * N
    out = []
    nodes = [0]
    rng = range(n)

    def consistent():
        for x in rng:
            xn = x * n
            for y in rng:
                u = t[xn + y]
                if u < 0:
                    continue
                un = u * n
                yn = y * n
                for z in rng:
                    v = t[yn + z]
                    if v < 0:
                        continue
                    w1 = t[un + z]
                    if w1 < 0:
                        continue
                    w2 = t[xn + v]
                    if w2 >= 0 and w1 != w2:
                        return False
        return True

    def rec(p):
        if p == N:
            out.append(tuple(t))
            return
        for val in rng:
            t[p] = val
            nodes[0] += 1
            if consistent():
                rec(p + 1)
        t[p] = -1

    rec(0)
    return out, nodes[0]


def brute_assoc(n):
    """The same set by unconditional brute force over all n^(n*n) tables -- a control on the
    generator at n = 2 and n = 3; hopeless at n = 4 (4^16 tables)."""
    return [t for t in itertools.product(range(n), repeat=n * n) if assoc(t, n)]


def orbit_count(tables, n):
    perms = list(itertools.permutations(range(n)))
    seen = set()
    reps = 0
    for t in tables:
        if t in seen:
            continue
        reps += 1
        for s in perms:
            seen.add(apply_perm(t, s, n))
    return reps


def cliques_by_level(adj, m):
    """All cliques level by level, size k -> k+1 until a level is empty.  No heuristic and no
    branch-and-bound, so an empty level is an exhaustive non-existence statement."""
    levels = {1: [(v,) for v in range(m)]}
    k = 1
    while levels[k]:
        nxt = []
        for c in levels[k]:
            common = adj[c[0]]
            for v in c[1:]:
                common &= adj[v]
            for v in range(c[-1] + 1, m):
                if common >> v & 1:
                    nxt.append(c + (v,))
        k += 1
        levels[k] = nxt
    return levels


def subset_scan(adj, m, k, keep=0):
    """The algorithm-free count: EVERY k-subset of the m vertices, tested pair by pair.
    Returns (cliques found, sample of them, NUMBER OF SUBSETS ACTUALLY EXAMINED) -- the third value
    is counted by the loop rather than quoted from a binomial coefficient, so the transcript's
    "all C(m,k) subsets" is a measurement and not an assertion."""
    total = 0
    found = []
    examined = 0
    for c in itertools.combinations(range(m), k):
        examined += 1
        good = True
        for i in range(k - 1):
            ai = adj[c[i]]
            for j in range(i + 1, k):
                if not ai >> c[j] & 1:
                    good = False
                    break
            if not good:
                break
        if good:
            total += 1
            if len(found) < keep:
                found.append(c)
    return total, found, examined


def binom(m, k):
    r = 1
    for i in range(k):
        r = r * (m - i) // (i + 1)
    return r


# ---------------------------------------------------------------------------
# 4. 2x2 MATRICES OVER GF(2)
# ---------------------------------------------------------------------------
def mat_apply(m, a):
    """m = ((m00,m01),(m10,m11)); a in 0..3 read as the column (a0,a1) with a = a0 + 2 a1."""
    a0, a1 = a & 1, (a >> 1) & 1
    return ((m[0][0] * a0 + m[0][1] * a1) & 1) + 2 * ((m[1][0] * a0 + m[1][1] * a1) & 1)


def mat_mul(m, p):
    return tuple(tuple(sum(m[i][k] * p[k][j] for k in range(2)) & 1 for j in range(2))
                 for i in range(2))


def mat_add(m, p):
    return tuple(tuple((m[i][j] + p[i][j]) & 1 for j in range(2)) for i in range(2))


def det2(m):
    return (m[0][0] * m[1][1] + m[0][1] * m[1][0]) & 1


def mat_rank(m):
    if m[0] == (0, 0) and m[1] == (0, 0):
        return 0
    if m[0] == (0, 0) or m[1] == (0, 0) or m[0] == m[1]:
        return 1
    return 2


ALL_MATS = [((a, b), (c, d)) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)]
M2 = ((0, 0), (0, 1))
M3 = ((1, 1), (0, 0))
M4 = ((1, 0), (1, 0))
IDENT = ((1, 0), (0, 1))
ZERO = ((0, 0), (0, 0))


def table_from_matrix(m):
    """f_M(a,b) = b XOR M(a) on {0,1,2,3} = (Z/2)^2, as a 0-based row-major table."""
    return tuple(b ^ mat_apply(m, a) for a in range(4) for b in range(4))


def gf4_mul(x, y):
    """GF(4) = F_2[w]/(w^2+w+1); 0,1,w,w+1 encoded 0,1,2,3 and addition = XOR."""
    r = 0
    for i in range(2):
        if y >> i & 1:
            r ^= x << i
    if r >> 2 & 1:
        r ^= 0b111
    return r & 3


# ===========================================================================
def main():
    t_start = time.time()
    print('verification of: four pairwise orthogonal semigroup operations on a four-element set,')
    print('and the 160 extremal families.  Source: arXiv:2608.25092v1, the seventh problem of')
    print('the closing list of problems, the cell n = 4.')
    print('python %d.%d, standard library only, exact integer arithmetic.' % sys.version_info[:2])
    print('')
    print('--- A. the family exhibited in Section 2 ---')
    n = 4
    fam = {}
    for name, flat in PAPER_FLAT:
        t = decode(flat)
        fam[name] = t
        ck('flat-decodes-' + name, len(t) == 16 and encode(t) == flat, flat)
    for name in ('f1', 'f2', 'f3', 'f4'):
        rows = PAPER_ARRAY[name]
        ck('array-agrees-with-flat-' + name,
           tuple(v - 1 for row in rows for v in row) == fam[name],
           'printed array = printed flat string')
    F = [fam[k] for k in ('f1', 'f2', 'f3', 'f4')]
    for name in ('f1', 'f2', 'f3', 'f4'):
        ck('associative-' + name, assoc(fam[name], n), 'all 64 triples')
    for name in ('f1', 'f2', 'f3', 'f4'):
        ck('balanced-' + name, balanced(fam[name], n), 'each value fills 4 of the 16 cells')
    for i, j in itertools.combinations(range(4), 2):
        ck('orthogonal-f%d-f%d' % (i + 1, j + 1), len(set(zip(F[i], F[j]))) == 16,
           'the joint map is a bijection of X^2')
    for i in range(4):
        ck('self-pair-not-orthogonal-f%d' % (i + 1), len(set(zip(F[i], F[i]))) == 4,
           'image lies in the diagonal, 4 < 16 points')
    ck('lower-bound-max4-at-least-4',
       len(set(F)) == 4 and all(assoc(t, n) for t in F)
       and all(orth(F[i], F[j], n) for i, j in itertools.combinations(range(4), 2)),
       'four distinct associative operations, pairwise orthogonal')
    ck('f1-is-the-left-zero-band',
       all(fam['f1'][a * n + b] == a for a in range(4) for b in range(4)), 'f1(x,y) = x')
    ck('f1-rees-type-I4-triv-L1',
       simple(fam['f1'], n) and rees_type(fam['f1'], n) == (4, 1, 1), '(|I|,|G|,|L|) = (4,1,1)')
    for name in ('f2', 'f3', 'f4'):
        ck('rees-type-I1-Z2-L2-' + name,
           simple(fam[name], n) and rees_type(fam[name], n) == (1, 2, 2)
           and group_orders(fam[name], n) == (1, 2), 'right group Z_2 x R_2')
    ck('f4-is-the-source-papers-published-table', encode(fam['f4']) == '1234432112344321',
       'the Mace4 example at tex 706-720')
    ck('f4-idempotents-are-1-and-3',
       sorted(a + 1 for a in idempotents(fam['f4'], n)) == [1, 3], 'E(S) = {1,3}, as printed there')

    print('')
    print('--- B. the F_2 picture of Section 3 ---')
    ck('xor-form-f1', all(fam['f1'][a * 4 + b] == a for a in range(4) for b in range(4)),
       'f1(a,b) = a')
    for name, m, mn in (('f2', M2, 'M2'), ('f3', M3, 'M3'), ('f4', M4, 'M4')):
        ck('xor-form-' + name, fam[name] == table_from_matrix(m),
           'f(a,b) = b + %s(a); %s(0,1,2,3) = %s'
           % (mn, mn, tuple(mat_apply(m, a) for a in range(4))))
    for m, mn in ((M2, 'M2'), (M3, 'M3'), (M4, 'M4')):
        ck('idempotent-' + mn, mat_mul(m, m) == m, 'M^2 = M')
        ck('rank-one-' + mn, mat_rank(m) == 1, 'rank 1')
    for (m, mn), (p, pn) in itertools.combinations(((M2, 'M2'), (M3, 'M3'), (M4, 'M4')), 2):
        ck('sum-invertible-%s-%s' % (mn, pn), det2(mat_add(m, p)) == 1, 'det = 1 over GF(2)')
    nid = sum(1 for m in ALL_MATS if mat_mul(m, m) == m)
    ck('associative-iff-idempotent-over-the-16-linear-maps',
       all(assoc(table_from_matrix(m), 4) == (mat_mul(m, m) == m) for m in ALL_MATS),
       'f_M is associative for exactly the %d idempotent M of 16' % nid)
    ck('orthogonal-iff-sum-invertible-over-all-linear-pairs',
       all(orth(table_from_matrix(m), table_from_matrix(p), 4) == (det2(mat_add(m, p)) == 1)
           for m in ALL_MATS for p in ALL_MATS), 'all 256 ordered pairs')
    ck('left-zero-band-orthogonal-to-every-linear-f',
       all(orth(fam['f1'], table_from_matrix(m), 4) for m in ALL_MATS), 'invertible M or not')
    rank1 = [m for m in ALL_MATS if mat_mul(m, m) == m and mat_rank(m) == 1]
    ck('exactly-six-rank-one-idempotents', len(rank1) == 6, 'over GF(2), 2x2')
    r_adj = [0] * 6
    r_edges = 0
    for i, j in itertools.combinations(range(6), 2):
        if det2(mat_add(rank1[i], rank1[j])) == 1:
            r_adj[i] |= 1 << j
            r_adj[j] |= 1 << i
            r_edges += 1
    ck('rank-one-graph-3-regular-9-edges',
       r_edges == 9 and sorted(bin(x).count('1') for x in r_adj) == [3] * 6, '9 edges, 3-regular')
    r_lv = cliques_by_level(r_adj, 6)
    ck('rank-one-graph-has-two-triangles-and-no-4-clique',
       len(r_lv[3]) == 2 and not r_lv[4], '1 + 3 = 4 is tight for this mechanism')
    ck('the-M2-M3-M4-triangle-is-one-of-them',
       frozenset((rank1.index(M2), rank1.index(M3), rank1.index(M4)))
       in {frozenset(c) for c in r_lv[3]}, '')
    ck('higher-rank-idempotents-are-not-excluded',
       mat_mul(ZERO, ZERO) == ZERO and mat_mul(IDENT, IDENT) == IDENT
       and det2(mat_add(ZERO, IDENT)) == 1, 'M = 0, N = I: idempotent, M + N invertible')
    sq = {}
    for a in (1, 2, 3):
        sq[a] = tuple(gf4_mul(a, x) ^ y for x in range(4) for y in range(4))
        ck('latin-square-a%d' % a, is_quasigroup(sq[a], 4), 'L_a(x,y) = a x + y over GF(4)')
    ck('three-mols-of-order-4',
       sum(1 for a, b in itertools.combinations((1, 2, 3), 2) if orth(sq[a], sq[b], 4)) == 3,
       'N(4) >= 3, so the classical ceiling N(4) + 2 = 5 is not vacuous')
    ck('only-one-of-the-three-mols-is-associative',
       sum(1 for a in (1, 2, 3) if assoc(sq[a], 4)) == 1, 'a^2 = a forces a = 1')

    print('')
    print('--- C. the unpruned census for n = 1,2,3,4 ---')
    census = {}
    for m in (1, 2, 3, 4):
        t0 = time.time()
        tables, nodes = all_assoc(m)
        census[m] = tables
        note('n=%d: %d tables from %d backtracking nodes, %d ms'
             % (m, len(tables), nodes, int(1000 * (time.time() - t0))))
        ck('A023814-n%d' % m, len(tables) == A023814[m], 'A023814(%d) = %d' % (m, A023814[m]))
        ck('census-all-associative-n%d' % m, all(assoc(t, m) for t in tables), 're-tested')
        ck('census-distinct-n%d' % m, len(set(tables)) == len(tables), 'no repeats')
        ck('A027851-n%d' % m, orbit_count(tables, m) == A027851[m],
           'A027851(%d) = %d, as the diagonal S_%d orbit count' % (m, A027851[m], m))
    for m in (2, 3):
        ck('backtracking-agrees-with-brute-force-n%d' % m,
           set(census[m]) == set(brute_assoc(m)), 'same set as a sweep of all %d^%d tables'
           % (m, m * m))

    graph = {}
    for m in (1, 2, 3, 4):
        tables = census[m]
        bal = sorted(t for t in tables if balanced(t, m))
        unbal = [t for t in tables if not balanced(t, m)]
        ck('balanced-count-n%d' % m, len(bal) == EXPECT_BALANCED[m],
           '%d of %d tables are value-balanced' % (len(bal), len(tables)))
        ck('every-unbalanced-table-has-an-overfull-value-n%d' % m,
           all(max(Counter(t).values()) > m for t in unbal),
           'each of the %d unbalanced tables has a value in > %d cells'
           % (len(unbal), m))
        ck('balanced-iff-completely-simple-n%d' % m,
           all(simple(t, m) == balanced(t, m) for t in tables),
           'set identity over all %d tables, both directions' % len(tables))
        adj = [0] * len(bal)
        e = 0
        for i, j in itertools.combinations(range(len(bal)), 2):
            if orth(bal[i], bal[j], m):
                adj[i] |= 1 << j
                adj[j] |= 1 << i
                e += 1
        ck('edges-n%d' % m, e == EXPECT_EDGES[m], '%d of C(%d,2) = %d pairs orthogonal'
           % (e, len(bal), len(bal) * (len(bal) - 1) // 2))
        ck('no-self-loops-n%d' % m, m == 1 or not any(orth(t, t, m) for t in bal),
           'no operation is orthogonal to itself for n >= 2')
        graph[m] = (bal, adj)

    rnd = random.Random(20260903)
    tab4 = census[4]
    set4 = set(graph[4][0])
    unb4 = [t for t in tab4 if t not in set4]
    hits = sum(1 for _ in range(100000)
               if orth(rnd.choice(unb4), rnd.choice(tab4), 4))
    ck('control-unbalanced-tables-are-isolated', hits == 0,
       '100000 seeded random pairs with an unbalanced member: %d orthogonal' % hits)

    print('')
    print('--- C2. cliques, level by level and by algorithm-free subset scan ---')
    for m in (1, 2, 3, 4):
        bal, adj = graph[m]
        levels = cliques_by_level(adj, len(bal))
        got = {k: len(v) for k, v in sorted(levels.items()) if k <= max(EXPECT_CLIQUES[m])}
        ck('cliques-by-size-n%d' % m, got == EXPECT_CLIQUES[m], '%s' % got)
        top = max(k for k, v in levels.items() if v)
        ck('max-n%d' % m, top == EXPECT_MAX[m], 'max(%d) = %d' % (m, top))
        ck('extremal-count-n%d' % m, len(levels[top]) == EXPECT_EXTREMAL[m],
           '%d labelled extremal families' % len(levels[top]))
    bal4, adj4 = graph[4]
    n4, found4, seen4 = subset_scan(adj4, len(bal4), 4, keep=200)
    ck('subset-scan-C48-4', n4 == 160 and seen4 == binom(48, 4) == 194580,
       '%d 4-subsets examined pair by pair (= C(48,4)): %d are cliques' % (seen4, n4))
    n5, _, seen5 = subset_scan(adj4, len(bal4), 5)
    ck('subset-scan-C48-5-is-empty', n5 == 0 and seen5 == binom(48, 5) == 1712304,
       '%d 5-subsets examined pair by pair (= C(48,5)): %d are cliques, so max(4) <= 4'
       % (seen5, n5))
    lv4 = cliques_by_level(adj4, len(bal4))
    ck('two-methods-agree-on-the-160',
       {frozenset(c) for c in lv4[4]} == {frozenset(c) for c in found4}, 'same 160 families')
    ck('max4-equals-4', len(lv4[4]) == 160 and not lv4[5] and not n5,
       'max(4) = 4 with exactly 160 labelled extremal families')

    types = Counter(rees_type(t, 4) for t in bal4)
    ck('rees-type-breakdown-of-the-48', dict(types) == EXPECT_TYPES_4,
       '%s' % sorted(types.items()))
    z4 = sum(1 for t in bal4 if rees_type(t, 4) == (1, 4, 1) and group_orders(t, 4) == (1, 2, 4, 4))
    v4 = sum(1 for t in bal4 if rees_type(t, 4) == (1, 4, 1) and group_orders(t, 4) == (1, 2, 2, 2))
    ck('sixteen-group-tables-twelve-cyclic-four-klein', (z4, v4) == (12, 4),
       '%d labelled Z_4 and %d labelled Klein tables' % (z4, v4))
    grp = [i for i, t in enumerate(bal4) if is_quasigroup(t, 4)]
    ck('group-operations-span-no-edge',
       len(grp) == 16 and not any(adj4[a] >> b & 1 for a, b in itertools.combinations(grp, 2)),
       'the 16 group operations on [4] are pairwise non-orthogonal')

    print('')
    print('--- D. the classification of Section 4 ---')
    fams = [frozenset(bal4[v] for v in c) for c in lv4[4]]
    ck('extremal-families-are-4-sets', all(len(f) == 4 for f in fams), '160 families of size 4')
    perms4 = list(itertools.permutations(range(4)))
    seen = set()
    orbits = []
    for f in fams:
        if f in seen:
            continue
        orb = set(frozenset(apply_perm(t, s, 4) for t in f) for s in perms4)
        orbits.append(orb)
        seen |= orb
    ck('orbit-count-n4', len(orbits) == EXPECT_ORBITS[4], '%d orbits under diagonal Sym(4)'
       % len(orbits))
    sizes = sorted(len(o) for o in orbits)
    ck('orbit-sizes-n4', sizes == EXPECT_ORBIT_SIZES_4, '%s' % sizes)
    ck('orbits-partition-the-160', sum(sizes) == 160 and len(seen) == 160, 'sum 160, cover 160')
    for m in (2, 3):
        balm, adjm = graph[m]
        lvm = cliques_by_level(adjm, len(balm))
        top = max(k for k, v in lvm.items() if v)
        fm = [frozenset(balm[v] for v in c) for c in lvm[top]]
        sm, cnt = set(), 0
        for f in fm:
            if f in sm:
                continue
            cnt += 1
            for s in itertools.permutations(range(m)):
                sm.add(frozenset(apply_perm(t, s, m) for t in f))
        ck('orbit-count-n%d' % m, cnt == EXPECT_ORBITS[m],
           '%d extremal families at n = %d in %d orbit(s)' % (len(fm), m, cnt))

    orbit_of = {}
    for idx, orb in enumerate(orbits):
        for g in orb:
            orbit_of[g] = idx
    hit = set()
    for k, (size, flats) in enumerate(PAPER_ORBITS, start=1):
        f = frozenset(decode(x) for x in flats)
        ok = f in orbit_of and len(f) == 4
        ck('paper-orbit-%d-is-an-extremal-family' % k, ok, 'the four printed tables')
        ck('paper-orbit-%d-has-size-%d' % (k, size),
           ok and len(orbits[orbit_of[f]]) == size, 'as printed')
        if ok:
            hit.add(orbit_of[f])
    ck('paper-lists-all-ten-orbits-once', len(hit) == 10,
       'the ten representatives meet %d distinct orbits' % len(hit))
    ck('section-2-family-is-the-sixth-printed-representative',
       frozenset(F) == frozenset(decode(x) for x in PAPER_ORBITS[5][1]), '')

    tL, tR = decode(FLAT_L), decode(FLAT_R)
    ck('L-and-R-are-vertices', tL in set4 and tR in set4, 'both balanced')
    ck('L-orth-R', orth(tL, tR, 4), 'their joint map is the identity of X^2')
    only_l = sum(1 for f in fams if tL in f and tR not in f)
    only_r = sum(1 for f in fams if tR in f and tL not in f)
    both = sum(1 for f in fams if tL in f and tR in f)
    neither = sum(1 for f in fams if tL not in f and tR not in f)
    ck('L-R-split-16-16-128-0', (only_l, only_r, neither, both) == (16, 16, 128, 0),
       'L only %d, R only %d, neither %d, both %d' % (only_l, only_r, neither, both))
    ck('no-extremal-family-contains-both-L-and-R', both == 0,
       'both-count over the 160 families is 0')
    withgrp = [f for f in fams if any(is_quasigroup(t, 4) for t in f)]
    ck('only-eight-families-contain-a-group-operation', len(withgrp) == 8,
       '%d of the 160' % len(withgrp))
    gops = set(t for f in withgrp for t in f if is_quasigroup(t, 4))
    ck('those-group-operations-are-the-four-klein-tables',
       len(gops) == 4 and all(group_orders(t, 4) == (1, 2, 2, 2) for t in gops),
       'no Z_4 table occurs in any extremal family')
    ck('the-printed-klein-table-occurs', decode(FLAT_KLEIN) in gops, FLAT_KLEIN)
    ck('the-group-containing-families-are-one-orbit-of-8',
       len({orbit_of[f] for f in withgrp}) == 1 and len(orbits[orbit_of[withgrp[0]]]) == 8, '')
    naive = frozenset({tL, tR, decode(FLAT_KLEIN)})
    ck('the-naive-3-family-is-not-extremal',
       len(naive) == 3 and all(orth(a, b, 4) for a, b in itertools.combinations(naive, 2))
       and not any(naive <= f for f in fams),
       '{L, R, Klein} is a 3-clique contained in no maximum family')

    print('')
    print('SCOPE -- what this program does NOT establish:')
    print('NOTE (a) it settles the cells n = 1,2,3,4 of the source Problem and nothing beyond them;')
    print('NOTE     no table of order above 4 is generated anywhere, so no value of max(n)')
    print('NOTE     for n >= 5 is computed or supported here.')
    print('NOTE (b) the ceiling max(n) <= N(n) + 2 quoted in the paper is classical (orthogonal')
    print('NOTE     arrays / MOLS). This program only exhibits three MOLS of order 4 to show the')
    print('NOTE     ceiling is not vacuous; it does not prove N(4) = 3, nor N(6) = 1 (Tarry).')
    print('NOTE (c) the attributions of Section 5 are textual and bibliographic facts about')
    print('NOTE     arXiv:2608.25092v1 and about the Belousov school; no program checks them. The')
    print('NOTE     identification of f4 with the table at tex 706-720 is checked here only at the')
    print('NOTE     level of the table itself and of E(S) = {1,3}.')
    print('NOTE (d) nothing here is a claim of novelty or of minimality; those are judgements.')
    print('')
    note('elapsed %d ms' % int(1000 * (time.time() - t_start)))
    if _BAD:
        print('VERDICT: %d of %d checks did NOT pass: %s'
              % (len(_BAD), len(_OK) + len(_BAD), _BAD[:10]))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % len(_OK))
    return 0


if __name__ == '__main__':
    sys.exit(main())
