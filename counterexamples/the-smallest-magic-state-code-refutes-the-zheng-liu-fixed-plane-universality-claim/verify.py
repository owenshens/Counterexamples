#!/usr/bin/env python3
"""verify.py -- re-derives every computational claim of paper.tex from the objects PRINTED IN THAT PAPER.

Python 3.9+, standard library only (no numpy, no sympy, no networkx), no external data file, no
network, no randomness.  Every stabilizer group, logical frame, spoiling product and coset
representative below is a literal copied from paper.tex; nothing is read from disk.

All arithmetic is exact: Pauli strings are handled as pairs of n-bit integers over GF(2) and every
count is an ordinary Python integer.  There is no floating point anywhere in a decision.  (Two
float numbers are PRINTED as commentary in the lambda table of Section 4; both are also checked in
exact rational form with Fraction, and the decision is taken on the Fraction.)

WHAT IS CHECKED
  1. the four named witness codes of the paper (its own [[3,1,1]] and [[4,1,1]], the n=2 minimum,
     and the 2016 Howard-Dawkins code that is itself an unrecognised counterexample): the group is
     abelian and independent, |N(S)| = 2^(n+k), the number of logical classes is 4^k, the printed
     spoiling products close as Pauli arithmetic and are Z-free, the printed class representatives
     are pairwise in distinct cosets and exhaust the classes, and -- independently of the printed
     certificate -- an exhaustive sweep of N(S) confirms the number of Z-clean classes.
  2. the k=2 certificate with exactly one clean class and no frame (the symplectic obstruction).
  3. the infinite family S(n,k), on every member with 2n <= 16: all 4^k classes are dirty, and the
     printed tensor-product representatives are Z-free and hit every class.
  4. the exhaustive census of ALL stabilizer groups (mod phase) on n <= 4 qubits with 1 <= n-k <= n-1:
     17478 codes, cell sizes against the closed form, 10232 for which z=0's condition is
     unsatisfiable, 5927 for which no coordinate plane's condition is satisfiable, and the full
     clean-class histogram at n=4, k=2 with its 664/864 split.
  5. the counting bound 4^k <= 3^n of Proposition 7, and the empirical lambda = (3/2)^n / 2^k
     table of Section 4 of the paper.
  6. the positive controls: the paper's own [[5,1,3]] and [[7,1,3]] with their published logical Z,
     the Z-type codes, and a code with no invariant coordinate plane in any direction.

Exit 0 if and only if every check passed.  Prints one `PASS <name> <detail>` line per check, then
`VERDICT: ALL <n> CHECKS PASS`.  What is NOT covered is printed as `NOT RE-RUN:` lines at the end.
"""

from fractions import Fraction
from itertools import combinations
import sys

# ---------------------------------------------------------------------------------------------
# 0.  Pauli arithmetic modulo phase.  A string over {I,X,Y,Z} of length n is the pair (a,b) of
#     n-bit integers with bit i of a set iff factor i is X or Y, and bit i of b set iff factor i
#     is Z or Y.  The product is the componentwise XOR; phases are discarded, which is sound
#     because the predicate under test reads only which LETTERS occur.
# ---------------------------------------------------------------------------------------------
LET = {'I': (0, 0), 'X': (1, 0), 'Z': (0, 1), 'Y': (1, 1)}
INV = {(0, 0): 'I', (1, 0): 'X', (0, 1): 'Z', (1, 1): 'Y'}


def enc(s):
    a = b = 0
    for i, c in enumerate(s):
        ai, bi = LET[c]
        a |= ai << i
        b |= bi << i
    return (a, b)


def dec(p, n):
    a, b = p
    return ''.join(INV[((a >> i) & 1, (b >> i) & 1)] for i in range(n))


def mul(p, q):
    return (p[0] ^ q[0], p[1] ^ q[1])


def commutes(p, q):
    return (bin(p[0] & q[1]).count('1') + bin(p[1] & q[0]).count('1')) % 2 == 0


def has_letter(p, L, n):
    """does the string contain at least one tensor factor equal to the letter L?"""
    a, b = p
    mk = (1 << n) - 1
    if L == 'Z':
        return ((~a) & b & mk) != 0
    if L == 'X':
        return (a & (~b) & mk) != 0
    if L == 'Y':
        return (a & b) != 0
    raise ValueError(L)


def span(gens):
    out = [(0, 0)]
    for g in gens:
        out = out + [mul(x, g) for x in out]
    return out


# ---------------------------------------------------------------------------------------------
# 1.  The check harness.
# ---------------------------------------------------------------------------------------------
_PASSES = []
_FAILS = []


def ck(name, cond, detail=''):
    if cond:
        _PASSES.append(name)
        print('PASS %s %s' % (name, detail))
    else:
        _FAILS.append(name)
        print('FAILED-CHECK %s %s' % (name, detail))


# ---------------------------------------------------------------------------------------------
# 2.  The analysis of one stabilizer code.
# ---------------------------------------------------------------------------------------------
def normaliser(n, gens):
    out = []
    for a in range(1 << n):
        for b in range(1 << n):
            p = (a, b)
            if all(commutes(p, g) for g in gens):
                out.append(p)
    return out


def classes_of(n, gens):
    """-> (S, N(S), [coset]) with each coset a list of its 2^(n-k) members."""
    S = span(gens)
    N = normaliser(n, gens)
    seen = set()
    cos = []
    for p in N:
        if p in seen:
            continue
        c = [mul(p, s) for s in S]
        seen.update(c)
        cos.append(c)
    return S, N, cos


def clean_classes(cos, L, n):
    """the classes that are CLEAN for the letter L: no member of the coset is L-free."""
    return [c for c in cos if all(has_letter(q, L, n) for q in c)]


def _independent_mod_S(sub, Sset):
    for r in range(1, len(sub) + 1):
        for cc in combinations(sub, r):
            p = (0, 0)
            for q in cc:
                p = mul(p, q)
            if p in Sset:
                return False
    return True


def frame_exists(k, clean, Sset):
    """is there a legitimate logical frame all of whose k classes are clean?  The k nominated
    operators must lie in N(S), pairwise commute, and be independent modulo S."""
    if k == 0:
        return True, None
    reps = [c[0] for c in clean]
    if len(reps) < k:
        return False, None
    for sub in combinations(reps, k):
        if all(commutes(sub[i], sub[j]) for i in range(k) for j in range(i + 1, k)) \
                and _independent_mod_S(sub, Sset):
            return True, sub
    return False, None


def strict_frame_exists(k, clean, Sset, S, n, L='Z'):
    """the STRICT reading: every nontrivial product of the nominated operators must also have a
    clean class.  For k = 1 it coincides with frame_exists."""
    if k == 0:
        return True
    reps = [c[0] for c in clean]
    if len(reps) < k:
        return False
    for sub in combinations(reps, k):
        if not (all(commutes(sub[i], sub[j]) for i in range(k) for j in range(i + 1, k))
                and _independent_mod_S(sub, Sset)):
            continue
        good = True
        for r in range(2, k + 1):
            for cc in combinations(sub, r):
                p = (0, 0)
                for q in cc:
                    p = mul(p, q)
                if any(not has_letter(mul(p, s), L, n) for s in S):
                    good = False
                    break
            if not good:
                break
        if good:
            return True
    return False


# ---------------------------------------------------------------------------------------------
# 3.  THE OBJECTS, copied literally out of paper.tex.
# ---------------------------------------------------------------------------------------------
# The named witness codes: Theorems 3 and 4 (Section 2) and the prior-art code of Section 5.  `spoil` holds the printed lines
# (class representative, stabilizer element, product).
WITNESSES = [
    dict(tag='W1', n=3, k=1,
         label='the target paper own [[3,1,1]], tex line 283',
         gens=['XZI', 'ZXX'],
         group=['III', 'XZI', 'ZXX', 'YYX'],
         ZL='IIX', XL='IZZ',
         spoil=[('III', 'III', 'III'),
                ('IIX', 'III', 'IIX'),
                ('IZZ', 'YYX', 'YXY'),
                ('IZY', 'XZI', 'XIY')],
         cleanZ=0, cleanY=1, cleanX=0),
    dict(tag='W2', n=4, k=1,
         label='the target paper own [[4,1,1]], tex line 298',
         gens=['XZII', 'ZXZX', 'IZXI'],
         group=['IIII', 'XZII', 'ZXZX', 'IZXI', 'YYZX', 'XIXI', 'ZYYX', 'YXYX'],
         ZL='IIIX', XL='IZIZ',
         spoil=[('IIII', 'IIII', 'IIII'),
                ('IIIX', 'IIII', 'IIIX'),
                ('IZIZ', 'YXYX', 'YYYY'),
                ('IZIY', 'XZII', 'XIIY')],
         cleanZ=0, cleanY=1, cleanX=0),
    dict(tag='W3', n=2, k=1,
         label='the global minimum n=2',
         gens=['XZ'],
         group=['II', 'XZ'],
         ZL='IZ', XL='ZX',
         spoil=[('II', 'II', 'II'),
                ('IZ', 'XZ', 'XI'),
                ('ZX', 'XZ', 'YY'),
                ('ZY', 'XZ', 'YX')],
         cleanZ=0, cleanY=1, cleanX=0),
    dict(tag='W4', n=3, k=1,
         label='Howard-Dawkins 2016 eq. (3a-gen), prior art, not ours',
         gens=['ZIZ', 'XZX'],
         group=['III', 'ZIZ', 'XZX', 'YZY'],
         ZL='XXY', XL='IXZ',
         spoil=[('III', 'III', 'III'),
                ('XXY', 'III', 'XXY'),
                ('IXZ', 'XZX', 'XYY'),
                ('XIX', 'III', 'XIX')],
         cleanZ=0, cleanY=1, cleanX=1),
]

# The k=2 certificate of Section 4: exactly one clean class, and still no frame.
CERT_K2 = dict(gens=['ZXXI', 'IZZZ'], n=4, k=2, cleanZ=1,
               clean_rep_coset=['IZZI', 'ZYYI', 'IIIZ', 'ZXXZ'])

# The n=2 atom's normaliser, coset partition and complete set of Z-free coset representatives
# (paper, Lemma 5 and its proof).
ATOM_REPS = ['II', 'XI', 'YY', 'YX']
ATOM_N_PRINTED = ['II', 'XZ', 'IZ', 'XI', 'ZX', 'YY', 'ZY', 'YX']
ATOM_COSETS_PRINTED = [['II', 'XZ'], ['IZ', 'XI'], ['ZX', 'YY'], ['ZY', 'YX']]

# The 16 printed Z-free coset representatives of S(4,2) (paper, display in Section 3).
FAMILY_42_REPS = ['IIII', 'IIXI', 'IIYY', 'IIYX',
                  'XIII', 'XIXI', 'XIYY', 'XIYX',
                  'YYII', 'YYXI', 'YYYY', 'YYYX',
                  'YXII', 'YXXI', 'YXYY', 'YXYX']

# Table 1 of the paper (Section 4), and its n=4, k=2 clean-class histogram.
CENSUS = [(2, 1, 1, 15, 4, 0),
          (3, 1, 2, 63, 0, 0),
          (3, 2, 1, 315, 168, 69),
          (4, 1, 3, 255, 0, 0),
          (4, 2, 2, 5355, 1528, 404),
          (4, 3, 1, 11475, 8532, 5454)]
CENSUS_TOTALS = (17478, 10232, 5927)
HIST_42 = {0: 208, 1: 456, 2: 1608, 3: 546, 4: 792, 6: 1172, 7: 573}
SPLIT_42 = (664, 864)

# Positive controls printed in the paper.
CTRL_513 = (5, ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ'], 'ZZZZZ')
CTRL_713 = (7, ['IIIXXXX', 'IXXIIXX', 'XIXIXIX', 'IIIZZZZ', 'IZZIIZZ', 'ZIZIZIZ'], 'ZZZZZZZ')
CTRL_NOPLANE = (3, ['ZZI', 'YXZ'])


def family_gens(n, k):
    """S(n,k) = < X_{2i-1} Z_{2i} : i = 1..k >  together with  < Z_j : j = 2k+1..n >."""
    g = []
    for i in range(1, k + 1):
        s = ['I'] * n
        s[2 * i - 2] = 'X'
        s[2 * i - 1] = 'Z'
        g.append(''.join(s))
    for j in range(2 * k + 1, n + 1):
        s = ['I'] * n
        s[j - 1] = 'Z'
        g.append(''.join(s))
    return g


def cell_size(n, m):
    """the number of m-dimensional isotropic subspaces of the symplectic space F_2^{2n}."""
    v = 1
    for i in range(1, m + 1):
        v = v * (2 ** (2 * (n - i + 1)) - 1) // (2 ** i - 1)
    return v


def isotropic_subspaces(n, m):
    """every m-dimensional isotropic subspace of F_2^{2n}, exactly once, as a list of m
    generators in reduced row echelon form.  A subspace is enumerated by its pivot columns and
    the free entries to the right of each pivot, which is a bijection onto RREF matrices."""
    d = 2 * n
    for piv in combinations(range(d), m):
        pivset = set(piv)
        free = [[c for c in range(piv[i] + 1, d) if c not in pivset] for i in range(m)]
        cnt = [len(f) for f in free]
        tot = sum(cnt)
        for bits in range(1 << tot):
            rows = []
            off = 0
            for i in range(m):
                v = 1 << piv[i]
                for j, c in enumerate(free[i]):
                    if (bits >> (off + j)) & 1:
                        v |= 1 << c
                off += cnt[i]
                # coordinates 0..n-1 are the b-bits, n..2n-1 the a-bits; which half is which is
                # immaterial, the enumeration is over subspaces either way.
                rows.append(((v >> n) & ((1 << n) - 1), v & ((1 << n) - 1)))
            ok = True
            for i in range(m):
                for j in range(i + 1, m):
                    if not commutes(rows[i], rows[j]):
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                yield rows


# ---------------------------------------------------------------------------------------------
# 4.  The checks.
# ---------------------------------------------------------------------------------------------
def check_witnesses():
    print('--- 1. the named witness codes of Section 2 and Section 5 ------------------------------------------')
    for w in WITNESSES:
        n, k, tag = w['n'], w['k'], w['tag']
        gens = [enc(g) for g in w['gens']]
        ck('%s.abelian' % tag, all(commutes(g, h) for g in gens for h in gens),
           'S = <%s> is abelian (%s)' % (','.join(w['gens']), w['label']))
        S = span(gens)
        ck('%s.independent' % tag, len(set(S)) == 2 ** len(gens),
           '|S| = 2^%d = %d, so the %d generators are independent'
           % (len(gens), 2 ** len(gens), len(gens)))
        ck('%s.group-printed' % tag, set(dec(p, n) for p in S) == set(w['group']),
           'the printed element list %s is exactly the closure of the generators' % (w['group'],))
        ck('%s.dimension' % tag, n - len(gens) == k, 'n - dim S = %d - %d = k = %d' % (n, len(gens), k))
        S, N, cos = classes_of(n, gens)
        Sset = set(S)
        ck('%s.normaliser' % tag, len(N) == 2 ** (n + k),
           '|N(S)| = %d = 2^(n+k) = 2^%d' % (len(N), n + k))
        ck('%s.class-count' % tag, len(cos) == 4 ** k,
           '|N(S)/S| = %d = 4^k = 4^%d' % (len(cos), k))
        # the printed logical frame
        ZL, XL = enc(w['ZL']), enc(w['XL'])
        ck('%s.ZL-in-N' % tag, all(commutes(ZL, g) for g in gens) and ZL not in Sset,
           'the frame operator Z_L = %s of Section 2 lies in N(S) and outside S' % w['ZL'])
        ck('%s.XL-in-N' % tag, all(commutes(XL, g) for g in gens) and XL not in Sset,
           'the frame operator X_L = %s of Section 2 lies in N(S) and outside S' % w['XL'])
        ck('%s.anticommute' % tag, not commutes(ZL, XL),
           'Z_L = %s and X_L = %s anticommute, so they are a genuine conjugate logical pair'
           % (w['ZL'], w['XL']))
        YL = mul(ZL, XL)
        ck('%s.YL' % tag, dec(YL, n) == w['spoil'][3][0],
           'Z_L * X_L = %s, the third nontrivial class representative printed alongside it'
           % dec(YL, n))
        # the printed spoiling lines
        keys = []
        for (q, s, p) in w['spoil']:
            pq, ps, pp = enc(q), enc(s), enc(p)
            ck('%s.arith[%s*%s]' % (tag, q, s), mul(pq, ps) == pp,
               '%s * %s = %s closes as Pauli arithmetic' % (q, s, p))
            ck('%s.multiplier[%s]' % (tag, s), ps in Sset,
               'the multiplier %s is a genuine element of S' % s)
            ck('%s.inN[%s]' % (tag, q), all(commutes(pq, g) for g in gens),
               'the class representative %s lies in N(S)' % q)
            ck('%s.zfree[%s]' % (tag, p), not has_letter(pp, 'Z', n),
               '%s has no tensor factor Z, so the class of %s is dirty' % (p, q))
            keys.append(min(mul(pq, s2) for s2 in S))
        ck('%s.distinct-cosets' % tag, len(set(keys)) == len(keys),
           'the %d printed representatives lie in pairwise distinct cosets of S' % len(keys))
        ck('%s.exhausts' % tag, len(set(keys)) == 4 ** k,
           '%d distinct dirty classes = 4^k = %d: every logical class is accounted for'
           % (len(set(keys)), 4 ** k))
        # and now, independently of the printed certificate, an exhaustive sweep
        for L, want in (('Z', w['cleanZ']), ('Y', w['cleanY']), ('X', w['cleanX'])):
            cl = clean_classes(cos, L, n)
            ck('%s.sweep-%s' % (tag, L), len(cl) == want,
               'exhaustive sweep of all %d elements of N(S): %d classes are %s-clean (paper says %d)'
               % (len(N), len(cl), L, want))
        clZ = clean_classes(cos, 'Z', n)
        okf, _ = frame_exists(k, clZ, Sset)
        ck('%s.no-frame' % tag, not okf,
           'no legitimate logical Z frame has all its classes Z-clean: Proposition 3 hypothesis '
           'is unsatisfiable for this code')
        ck('%s.no-strict-frame' % tag, not strict_frame_exists(k, clZ, Sset, S, n),
           'and unsatisfiable under the stricter all-products reading too')
        clY = clean_classes(cos, 'Y', n)
        oky, _ = frame_exists(k, clY, Sset)
        ck('%s.y-plane' % tag, oky == (w['cleanY'] >= k),
           'the analogous condition in the y direction is %s for this code'
           % ('satisfiable' if oky else 'unsatisfiable'))


def check_cert_k2():
    print('--- 2. the k=2 certificate with one clean class and no frame ------------------------')
    n, k = CERT_K2['n'], CERT_K2['k']
    gens = [enc(g) for g in CERT_K2['gens']]
    S, N, cos = classes_of(n, gens)
    Sset = set(S)
    ck('cert2.abelian', all(commutes(g, h) for g in gens for h in gens),
       'S = <%s> is abelian and n=%d, k=%d' % (','.join(CERT_K2['gens']), n, k))
    ck('cert2.normaliser', len(N) == 2 ** (n + k) and len(cos) == 4 ** k,
       '|N(S)| = %d = 2^%d and |N(S)/S| = %d = 4^%d' % (len(N), n + k, len(cos), k))
    clZ = clean_classes(cos, 'Z', n)
    ck('cert2.one-clean', len(clZ) == CERT_K2['cleanZ'],
       'exactly %d class is Z-clean' % len(clZ))
    printed = set(CERT_K2['clean_rep_coset'])
    ck('cert2.clean-coset', len(clZ) == 1 and set(dec(p, n) for p in clZ[0]) == printed,
       'the clean class is the coset %s printed in the paper' % sorted(printed))
    ck('cert2.clean-has-Z', all(has_letter(enc(s), 'Z', n) for s in CERT_K2['clean_rep_coset']),
       'every one of its %d members carries a tensor factor Z' % len(printed))
    okf, _ = frame_exists(k, clZ, Sset)
    ck('cert2.no-frame', not okf,
       'and still no frame exists: one clean class cannot furnish k=2 commuting independent ones')


def check_family():
    print('--- 3. the infinite family S(n,k) --------------------------------------------------')
    ck('family.atom-reps', all(not has_letter(enc(s), 'Z', 2) for s in ATOM_REPS),
       'the four printed representatives %s of the atom <XZ> in Lemma 5 are all Z-free' % ATOM_REPS)
    a = [enc('XZ')]
    S2, N2, cos2 = classes_of(2, a)
    ck('family.atom-normaliser',
       sorted(dec(p, 2) for p in N2) == sorted(ATOM_N_PRINTED),
       'N(<XZ>) is exactly the eight strings %s printed in the proof of Lemma 5'
       % sorted(ATOM_N_PRINTED))
    ck('family.atom-coset-partition',
       sorted(sorted(dec(q, 2) for q in c) for c in cos2)
       == sorted(sorted(c) for c in ATOM_COSETS_PRINTED),
       'and it splits into exactly the four printed cosets %s' % (ATOM_COSETS_PRINTED,))
    ck('family.atom-pad', [dec(p, 1) for p in normaliser(1, [enc('Z')])] == ['I', 'Z'],
       'the normaliser of the one-qubit pad <Z> is {I,Z}, a single coset containing the Z-free I')
    keys = [min(mul(enc(s), t) for t in S2) for s in ATOM_REPS]
    ck('family.atom-cosets', len(set(keys)) == 4 == len(cos2),
       'and they lie in 4 = 4^1 pairwise distinct cosets, so all four classes of <XZ> are dirty')
    ck('family.42-reps-zfree', all(not has_letter(enc(s), 'Z', 4) for s in FAMILY_42_REPS),
       'all %d printed representatives of S(4,2) in the display of Section 3 are Z-free' % len(FAMILY_42_REPS))
    g42 = [enc(g) for g in family_gens(4, 2)]
    S42 = span(g42)
    keys = [min(mul(enc(s), t) for t in S42) for s in FAMILY_42_REPS]
    ck('family.42-reps-distinct', len(set(keys)) == 16,
       'they lie in 16 = 4^2 pairwise distinct cosets of S(4,2), hence every class is dirty')
    ck('family.42-reps-inN', all(all(commutes(enc(s), g) for g in g42) for s in FAMILY_42_REPS),
       'and every one of them lies in N(S(4,2))')
    for k in (1, 2, 3):
        for n in range(2 * k, 2 * k + 3):
            if 2 * n > 16:
                continue
            gens = family_gens(n, k)
            pg = [enc(g) for g in gens]
            S, N, cos = classes_of(n, pg)
            clZ = clean_classes(cos, 'Z', n)
            ck('family.n%dk%d' % (n, k),
               len(N) == 2 ** (n + k) and len(cos) == 4 ** k and len(clZ) == 0,
               'S(%d,%d) = <%s>: |N(S)| = 2^%d, %d = 4^%d classes, ALL dirty, 0 clean'
               % (n, k, ','.join(gens), n + k, len(cos), k))


def check_census():
    print('--- 4. the exhaustive census over n <= 4 -------------------------------------------')
    tot_codes = tot_z = tot_any = 0
    hist42 = {}
    split42 = [0, 0]
    for (n, m, k, want_codes, want_z, want_any) in CENSUS:
        codes = zfail = anyfail = 0
        for gens in isotropic_subspaces(n, m):
            codes += 1
            S, N, cos = classes_of(n, gens)
            Sset = set(S)
            res = {}
            for L in 'ZYX':
                cl = clean_classes(cos, L, n)
                res[L] = (len(cl), frame_exists(k, cl, Sset)[0])
            if not res['Z'][1]:
                zfail += 1
            if not (res['Z'][1] or res['Y'][1] or res['X'][1]):
                anyfail += 1
            if (n, m) == (4, 2):
                nc = res['Z'][0]
                hist42[nc] = hist42.get(nc, 0) + 1
                if not res['Z'][1]:
                    split42[0 if nc < k else 1] += 1
        ck('census.n%dm%d.count' % (n, m), codes == want_codes == cell_size(n, m),
           'n=%d, dim S=%d, k=%d: %d stabilizer groups enumerated = closed form %d = paper %d'
           % (n, m, k, codes, cell_size(n, m), want_codes))
        ck('census.n%dm%d.zfail' % (n, m), zfail == want_z,
           'of those, %d admit NO Z-clean logical frame (paper says %d)' % (zfail, want_z))
        ck('census.n%dm%d.anyfail' % (n, m), anyfail == want_any,
           'and %d admit no clean frame in ANY of the three coordinate directions (paper says %d)'
           % (anyfail, want_any))
        tot_codes += codes
        tot_z += zfail
        tot_any += anyfail
    ck('census.totals', (tot_codes, tot_z, tot_any) == CENSUS_TOTALS,
       'totals over the six cells: %d codes, %d with z=0 unsatisfiable, %d with no plane at all '
       '(paper says %s)' % (tot_codes, tot_z, tot_any, CENSUS_TOTALS))
    ck('census.hist42', hist42 == HIST_42,
       'the n=4, k=2 clean-class histogram is %s, exactly as printed'
       % dict(sorted(hist42.items())))
    ck('census.hist42-sums', sum(hist42.values()) == 5355,
       'and its buckets sum to 5355, the whole cell')
    ck('census.split42', tuple(split42) == SPLIT_42 and sum(split42) == 1528,
       '%d of the 1528 counterexamples have fewer than k=2 clean classes and %d have two or more '
       'clean classes yet no commuting independent pair' % (split42[0], split42[1]))
    ck('census.minimal-k2', CENSUS[1][4] == 0 and CENSUS[4][4] > 0,
       'n=3, k=2 is counterexample-free (0 of 63), so the least n for a k=2 counterexample is 4')


def check_distinctness():
    """Section 5 of the paper: the 2016 Howard-Dawkins code is not a relabelling of the target's."""
    print('--- 4b. the two three-qubit codes are genuinely different ---------------------------')
    A = [enc(g) for g in WITNESSES[0]['gens']]      # <XZI,ZXX>, the target's [[3,1,1]]
    B = [enc(g) for g in WITNESSES[3]['gens']]      # <ZIZ,XZX>, Howard-Dawkins eq. (3a-gen)
    def i_bearing(gens):
        out = []
        for p in span(gens):
            s = dec(p, 3)
            if 'I' in s and set(s) != {'I'}:
                out.append(s)
        return out
    ia, ib = i_bearing(A), i_bearing(B)
    ck('distinct.unique-I-element', len(ia) == 1 and len(ib) == 1,
       'each group has exactly one non-identity element containing a factor I: %s and %s'
       % (ia[0], ib[0]))
    ck('distinct.multisets', sorted(ia[0]) != sorted(ib[0]),
       'their letter multisets differ (%s vs %s), and a qubit permutation preserves the multiset, '
       'so no permutation of qubits carries one group to the other'
       % (sorted(ia[0]), sorted(ib[0])))
    # and the two groups really are distinct as sets, permutation or not
    ck('distinct.as-sets', set(dec(p, 3) for p in span(A)) != set(dec(p, 3) for p in span(B)),
       'the two stabilizer groups are distinct as sets of strings')


def check_bounds():
    print('--- 5. the counting bound and the lambda table --------------------------------------')
    for n in (2, 3, 4, 5):
        cnt = sum(1 for a in range(1 << n) for b in range(1 << n)
                  if not has_letter((a, b), 'Z', n))
        ck('bound.zfree-count-n%d' % n, cnt == 3 ** n,
           'the Z-free strings of length %d number %d = 3^%d, the bound used in Proposition 7'
           % (n, cnt, n))
    # 4^k <= 3^n  is exactly  k/n <= (log_2 3)/2, and the paper prints that constant as 0.79248.
    # The bracket below pins it to five decimals in exact integer arithmetic, with no logarithm.
    ck('bound.threshold', 4 ** 79248 <= 3 ** 100000 < 4 ** 79249,
       'the collapse inequality 4^k <= 3^n is exactly k <= c*n with c = (log_2 3)/2, and '
       '4^79248 <= 3^100000 < 4^79249 pins c to 0.79248 <= c < 0.79249 without any logarithm')
    for (n, m, k, codes, zf, af) in CENSUS:
        collapse_possible = 4 ** k <= 3 ** n
        ck('bound.n%dk%d' % (n, k), (zf == 0) or collapse_possible,
           '4^%d = %d vs 3^%d = %d: total collapse is %s, and %d of %d codes fail'
           % (k, 4 ** k, n, 3 ** n, 'not excluded' if collapse_possible else 'EXCLUDED', zf, codes))
    for (n, m, k, codes, zf, af) in CENSUS:
        lam = Fraction(3, 2) ** n / Fraction(2) ** k
        ck('lambda.n%dk%d' % (n, k), (lam > 1) == (zf > 0),
           'lambda = (3/2)^%d / 2^%d = %s (%.4f) and the failure rate is %d/%d = %.2f%%: the sign '
           'of lambda - 1 agrees with the sign of the count'
           % (n, k, lam, float(lam), zf, codes, 100.0 * zf / codes))
    ck('bound.family-consistent', all(4 ** k <= 3 ** n for k in (1, 2, 3) for n in (2 * k,)),
       'the family lives at k <= n/2, strictly inside the collapse-permitted region k <= 0.79248 n')


def check_controls():
    print('--- 6. positive controls ------------------------------------------------------------')
    for (n, gens, ZL) in (CTRL_513, CTRL_713):
        pg = [enc(g) for g in gens]
        S, N, cos = classes_of(n, pg)
        k = n - len(gens)
        ck('ctrl.n%d.abelian' % n, all(commutes(g, h) for g in pg for h in pg),
           'the [[%d,%d,3]] code of the target paper is abelian with |N(S)| = %d = 2^%d'
           % (n, k, len(N), n + k))
        pz = enc(ZL)
        key = min(mul(pz, s) for s in S)
        clZ = clean_classes(cos, 'Z', n)
        ck('ctrl.n%d.published-ZL-clean' % n,
           any(min(mul(c[0], s) for s in S) == key for c in clZ),
           'its published logical Z = %s has a CLEAN class: every element of Z_L*S carries a Z, '
           'so the paper condition holds here exactly as the paper says' % ZL)
        ck('ctrl.n%d.clean-count' % n, len(clZ) == 1,
           'and it is the only clean class of the %d' % len(cos))
    n, gens = CTRL_NOPLANE
    pg = [enc(g) for g in gens]
    S, N, cos = classes_of(n, pg)
    Sset = set(S)
    k = n - len(gens)
    res = {}
    for L in 'ZYX':
        cl = clean_classes(cos, L, n)
        res[L] = frame_exists(k, cl, Sset)[0]
    ck('ctrl.noplane', not any(res.values()),
       'S = <%s> at n=%d, k=%d admits no clean frame in ANY coordinate direction '
       '(z:%s y:%s x:%s)' % (','.join(gens), n, k, res['Z'], res['Y'], res['X']))
    # Z-type codes always satisfy the condition
    bad = 0
    zonly = 0
    for n in (2, 3, 4):
        for m in range(1, n):
            for gens in isotropic_subspaces(n, m):
                if any(a != 0 for (a, b) in gens):
                    continue                      # not a Z-type code
                zonly += 1
                S, N, cos = classes_of(n, gens)
                cl = clean_classes(cos, 'Z', n)
                if not frame_exists(n - m, cl, set(S))[0]:
                    bad += 1
    ck('ctrl.z-type', bad == 0 and zonly > 0,
       'all %d Z-type stabilizer codes with n <= 4 satisfy the condition, so the conjecture is not '
       'vacuous and an unstructured survey of protocols could miss the failure' % zonly)


def main():
    print('verify.py -- the Zheng-Liu fixed-plane universality claim, refuted')
    print('objects: read from the literals of paper.tex; no external data file, exact arithmetic')
    print('')
    check_witnesses()
    check_cert_k2()
    check_family()
    check_census()
    check_distinctness()
    check_bounds()
    check_controls()
    print('')
    print('NOT RE-RUN: the census here is exhaustive for n <= 4 only (all 17478 stabilizer groups '
          'with 1 <= dim S <= n-1). The n = 5 cells (1023, 86955, 782595, 11475-analogue 782595) '
          'and every n >= 6 row quoted in the discovery record are NOT recomputed here; a miss '
          'there would be inconclusive rather than negative.')
    print('NOT RE-RUN: the target paper [[15,1,3]] and [[14,2,2]] protocols, and its [[6,1,2]] '
          'code, are not decided here -- only its [[3,1,1]], [[4,1,1]], [[5,1,3]] and [[7,1,3]] '
          'are. The paper own numerics are therefore neither confirmed nor contradicted by this '
          'program beyond those four codes.')
    print('NOT RE-RUN: nothing here computes the dynamical map of the target paper. The claim '
          'checked is exactly that Proposition 3 hypothesis is unsatisfiable; whether z = 0 is '
          'nonetheless an invariant plane for these codes is a different question, and the '
          'y-direction checks above show these codes DO have an invariant coordinate plane.')
    print('NOT RE-RUN: the lambda = (3/2)^n / 2^k rule is checked only for sign agreement on the '
          'six cells above. It is an empirical heuristic fitted to census data and no proof of it '
          'is claimed or attempted.')
    print('')
    if _FAILS:
        print('%d CHECKS FAILED: %s' % (len(_FAILS), ', '.join(_FAILS)))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % len(_PASSES))
    return 0


if __name__ == '__main__':
    sys.exit(main())
