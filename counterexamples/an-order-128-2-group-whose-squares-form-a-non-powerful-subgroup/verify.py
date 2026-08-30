#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification program for

    "An Order-128 2-Group whose Squares Form a Non-Powerful Subgroup"

Python 3.9+, STANDARD LIBRARY ONLY (itertools, sys).  No third-party package, no external
data file, no floating point anywhere: every quantity below is an exact integer or an exact
permutation of {1,...,8}.

WHAT IT READS.  Its whole input is the object printed in the paper, transcribed verbatim into
the block PAPER_* below:

  * the seven permutations f1,...,f7 of {1,...,8}, in the cycle notation of Section 2;
  * the 28 relators of the presentation of Section 2 (seven squares, 21 commutators);
  * the sixteen elements of the square set with a square root for each, from Table 1;
  * every numerical invariant the paper asserts.

Nothing else is supplied.  In particular NO group-theory library, no SmallGroups library and
no GAP: the group is rebuilt from the printed permutations by closure under multiplication,
and the presentation is handled by a collection procedure that applies the printed relators
and nothing else.

CONVENTIONS, matching the paper.  Permutations act on the right and are composed left to
right: (pq)(x) = q(p(x)).  Commutators are [x,y] = x^{-1} y^{-1} x y.

Output: one `PASS <name>` line per check, the `NOT RE-RUN:` statements of what is out of
scope, then `VERDICT: ALL <n> CHECKS PASS`.  Exit status 0 iff every check passed.
"""

import itertools
import sys

# =============================================================================
# THE OBJECT, TRANSCRIBED FROM THE PAPER
# =============================================================================

DEGREE = 8

# Section 2, the seven generators as permutations of {1,...,8}.
PAPER_PERMS = {
    'f1': '(1 5)(2 6)(3 7)(4 8)',
    'f2': '(3 4)',
    'f3': '(5 7)(6 8)',
    'f4': '(3 4)(7 8)',
    'f5': '(1 3)(2 4)(5 7)(6 8)',
    'f6': '(5 6)(7 8)',
    'f7': '(1 2)(3 4)(5 6)(7 8)',
}

# Section 2, the six non-trivial commutator relators, [f_a, f_b] = w with a > b.
PAPER_COMMUTATORS = {
    (2, 1): 'f4',
    (3, 1): 'f5',
    (6, 1): 'f7',
    (5, 2): 'f6f7',
    (4, 3): 'f6',
    (5, 4): 'f7',
}
# Every other pair of generators commutes; every generator squares to the identity.
PAPER_RELATOR_COUNT = 28

# Section 2, the generators of the iterated wreath product C2 wr C2 wr C2 acting on the
# complete binary tree with leaves 1,...,8 (level generators, bottom to top).
PAPER_WREATH_GENS = ['(1 2)', '(3 4)', '(5 6)', '(7 8)',
                     '(1 3)(2 4)', '(5 7)(6 8)',
                     '(1 5)(2 6)(3 7)(4 8)']

# Table 1: the sixteen elements of S = {g^2 : g in G}, each with a square root in G and that
# root's order.  Words are read left to right in the generators f1,...,f7.
PAPER_ROOTS = [
    ('1',        '1',          1),
    ('f4',       'f1f2',       4),
    ('f5',       'f1f3',       4),
    ('f6',       'f3f4',       4),
    ('f7',       'f1f6',       4),
    ('f4f5',     'f1f2f3f5',   8),
    ('f4f6',     'f1f2f4f5',   4),
    ('f4f7',     'f1f2f6',     4),
    ('f5f6',     'f1f3f4',     4),
    ('f5f7',     'f1f3f6',     4),
    ('f6f7',     'f2f5',       4),
    ('f4f5f6',   'f1f2f3f6',   8),
    ('f4f5f7',   'f1f2f3f4',   8),
    ('f4f6f7',   'f1f2f5',     4),
    ('f5f6f7',   'f1f3f4f5',   4),
    ('f4f5f6f7', 'f1f2f3',     8),
]

# Every numeric assertion of the paper, in one place.
PAPER_CLAIMS = {
    'order_G': 128,
    'order_G_is_2_to_the': 7,
    'exponent_G': 8,
    'order_dist_G': {1: 1, 2: 43, 4: 68, 8: 16},
    'factorial_8': 40320,
    'odd_part_factorial_8': 315,
    'order_S': 16,
    'exponent_S': 4,
    'order_dist_S': {1: 1, 2: 11, 4: 4},
    'involutions_in_S': 11,
    'involutions_in_C2_x_Q8': 3,
    'order_S_derived': 2,
    'order_S_squared': 2,
    'order_S_fourth': 1,
    'order_centre_G': 2,
    'lower_central_orders': [16, 4, 2, 1],   # gamma_2, gamma_3, gamma_4, gamma_5
    'order_commutator_S_G': 4,
    'order_D8_factor': 8,
    'order_f4f5': 4,
}

# =============================================================================
# EXACT PERMUTATION ARITHMETIC
# =============================================================================

IDENT = tuple(range(1, DEGREE + 1))


def parse_perm(text):
    """Cycle notation -> image tuple.  Raises on anything that is not a permutation."""
    img = list(range(1, DEGREE + 1))
    body = text.strip()
    if body in ('()', '1'):
        return tuple(img)
    moved = set()
    i = 0
    while i < len(body):
        if body[i] != '(':
            raise ValueError('bad cycle notation %r' % text)
        j = body.index(')', i)
        pts = [int(tok) for tok in body[i + 1:j].replace(',', ' ').split()]
        if len(pts) < 2:
            raise ValueError('cycle of length < 2 in %r' % text)
        for p in pts:
            if not (1 <= p <= DEGREE):
                raise ValueError('point %d out of range in %r' % (p, text))
            if p in moved:
                raise ValueError('point %d repeated in %r' % (p, text))
            moved.add(p)
        for k, p in enumerate(pts):
            img[p - 1] = pts[(k + 1) % len(pts)]
        i = j + 1
    return tuple(img)


def mul(p, q):
    """First p, then q."""
    return tuple(q[p[i] - 1] for i in range(DEGREE))


def inv(p):
    r = [0] * DEGREE
    for i in range(DEGREE):
        r[p[i] - 1] = i + 1
    return tuple(r)


def comm(x, y):
    """[x,y] = x^{-1} y^{-1} x y."""
    return mul(mul(inv(x), inv(y)), mul(x, y))


def order(p):
    o, q = 1, p
    while q != IDENT:
        q = mul(q, p)
        o += 1
    return o


def closure(gens):
    """The subgroup generated by `gens`, as a set of image tuples."""
    H = {IDENT}
    frontier = [IDENT]
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = mul(a, g)
                if b not in H:
                    H.add(b)
                    nxt.append(b)
        frontier = nxt
    return H


def cycle_str(p):
    seen, out = set(), ''
    for i in range(1, DEGREE + 1):
        if i in seen or p[i - 1] == i:
            continue
        cyc, j = [i], p[i - 1]
        seen.add(i)
        while j != i:
            cyc.append(j)
            seen.add(j)
            j = p[j - 1]
        out += '(' + ' '.join(str(x) for x in cyc) + ')'
    return out or '()'


F = {name: parse_perm(txt) for name, txt in PAPER_PERMS.items()}
GEN = [F['f%d' % i] for i in range(1, 8)]


def word(w):
    """A word such as 'f1f2f3' in the printed generators -> permutation.  '1' is empty."""
    if w in ('1', ''):
        return IDENT
    q = IDENT
    i = 0
    while i < len(w):
        if w[i] != 'f':
            raise ValueError('bad word %r' % w)
        q = mul(q, F[w[i:i + 2]])
        i += 2
    return q


def letters(w):
    """'f1f2' -> [1, 2]."""
    return [] if w in ('1', '') else [int(w[i + 1]) for i in range(0, len(w), 2)]


# =============================================================================
# CHECK BOOKKEEPING
# =============================================================================

_passes, _fails = [], []


def check(name, ok, detail=''):
    if ok:
        _passes.append(name)
        print('PASS %s%s' % (name, (' ' + detail) if detail else ''))
    else:
        _fails.append(name)
        print('FAIL %s%s' % (name, (' ' + detail) if detail else ''))


# =============================================================================
# A. THE PRINTED PERMUTATIONS SATISFY THE PRINTED PRESENTATION
# =============================================================================

check('perms-parse-as-degree-8-permutations',
      all(sorted(F['f%d' % i]) == list(range(1, DEGREE + 1)) for i in range(1, 8))
      and all(cycle_str(F['f%d' % i]) == PAPER_PERMS['f%d' % i].replace(',', ' ')
              for i in range(1, 8)),
      'seven permutations of {1..8}, cycle notation round-trips')

check('seven-square-relators',
      all(mul(F['f%d' % i], F['f%d' % i]) == IDENT for i in range(1, 8)),
      'f_i^2 = 1 for i = 1..7')

for (a, b), rhs in sorted(PAPER_COMMUTATORS.items()):
    got = comm(F['f%d' % a], F['f%d' % b])
    check('commutator-relator-f%d-f%d' % (a, b), got == word(rhs),
          '[f%d,f%d] = %s = %s' % (a, b, rhs, cycle_str(got)))

_trivial = [(a, b) for a in range(1, 8) for b in range(1, a)
            if (a, b) not in PAPER_COMMUTATORS]
check('remaining-pairs-commute',
      len(_trivial) == 15
      and all(comm(F['f%d' % a], F['f%d' % b]) == IDENT for a, b in _trivial),
      'the other %d of the 21 generator pairs commute' % len(_trivial))

check('relator-count-is-28', 7 + 21 == PAPER_RELATOR_COUNT,
      '7 squares + C(7,2) = 21 commutators = 28 relators, 6 of them non-trivial')

# =============================================================================
# B. THE GROUP THE PERMUTATIONS GENERATE
# =============================================================================

G = closure(GEN)
Gs = sorted(G)

check('group-order-128', len(G) == PAPER_CLAIMS['order_G'],
      '|<f1,...,f7>| = %d' % len(G))

check('group-is-a-2-group',
      len(G) == 2 ** PAPER_CLAIMS['order_G_is_2_to_the'],
      '128 = 2^7, so G is a finite 2-group')

_wreath = closure([parse_perm(t) for t in PAPER_WREATH_GENS])
check('equals-iterated-wreath-product-C2-wr-C2-wr-C2',
      _wreath == G,
      'the closure of the seven binary-tree level generators is the same set of %d '
      'permutations' % len(_wreath))

check('sylow-2-subgroup-of-S8',
      PAPER_CLAIMS['factorial_8'] == 2 * 3 * 4 * 5 * 6 * 7 * 8
      and PAPER_CLAIMS['factorial_8'] == 128 * PAPER_CLAIMS['odd_part_factorial_8']
      and PAPER_CLAIMS['odd_part_factorial_8'] % 2 == 1,
      '8! = 40320 = 128 * 315 with 315 odd, so a subgroup of S_8 of order 128 is a full '
      'Sylow 2-subgroup')

check('transitive-on-8-points',
      len({g[0] for g in G}) == DEGREE,
      'the orbit of the point 1 is all of {1..8}')

_dist = {}
for g in G:
    _dist[order(g)] = _dist.get(order(g), 0) + 1
check('order-distribution-of-G', _dist == PAPER_CLAIMS['order_dist_G'],
      'element orders 1/2/4/8 occur %s times, summing to %d'
      % ('/'.join(str(_dist[k]) for k in sorted(_dist)), sum(_dist.values())))

check('exponent-of-G-is-8',
      max(_dist) == PAPER_CLAIMS['exponent_G'],
      'exp G = %d exactly: an element of order 8 exists and none has order 16' % max(_dist))

# =============================================================================
# C. THE PRESENTATION DETERMINES G (COLLECTION FROM THE PRINTED RELATORS ONLY)
# =============================================================================

NORMAL = {}
for bits in itertools.product((0, 1), repeat=7):
    q = IDENT
    for k, b in enumerate(bits):
        if b:
            q = mul(q, GEN[k])
    NORMAL[bits] = q

check('128-normal-words-are-distinct',
      len(set(NORMAL.values())) == 128 and set(NORMAL.values()) == G,
      'the 128 words f1^a1...f7^a7 with a_i in {0,1} take 128 distinct values, exactly G')


def collect(word_letters, cap=100000):
    """Rewrite a word to normal form using ONLY the printed relators.

    Two rewrites are used, each a direct instance of a relator:
      (i)  f_i f_i -> 1                            from f_i^2 = 1;
      (ii) f_a f_b -> f_b f_a w   for a > b        from f_a f_b = f_b f_a [f_a,f_b],
           where w is the printed right-hand side of [f_a,f_b] (empty when the pair
           commutes).  Every such w involves only generators of index > a.
    Returns the exponent-bit tuple of the resulting strictly ascending repetition-free word.
    """
    w = list(word_letters)
    steps = 0
    while True:
        steps += 1
        if steps > cap:
            raise RuntimeError('collection did not terminate within %d steps' % cap)
        for i in range(len(w) - 1):
            if w[i] == w[i + 1]:
                del w[i:i + 2]
                break
            if w[i] > w[i + 1]:
                a, b = w[i], w[i + 1]
                rhs = letters(PAPER_COMMUTATORS.get((a, b), '1'))
                w[i:i + 2] = [b, a] + rhs
                break
        else:
            bits = [0] * 7
            for x in w:
                bits[x - 1] = 1
            return tuple(bits)


_bad, _maxsteps = 0, 0
for u, v in itertools.product(sorted(NORMAL), repeat=2):
    lu = [k + 1 for k, b in enumerate(u) if b]
    lv = [k + 1 for k, b in enumerate(v) if b]
    got = collect(lu + lv)
    if NORMAL[got] != mul(NORMAL[u], NORMAL[v]):
        _bad += 1
check('collection-closes-all-16384-products', _bad == 0,
      'each of the 128*128 = 16384 products of normal words collects, by relator rewrites '
      'alone, to a normal word equal to the product in the permutation model')

check('presentation-presents-exactly-G', _bad == 0 and len(set(NORMAL.values())) == 128,
      'collection bounds the presented group by 128 while the permutation model realises '
      '128, so the 28 relators present G and no larger group')

# =============================================================================
# D. THE SET OF SQUARES IS A NON-ABELIAN SUBGROUP
# =============================================================================

S = {mul(g, g) for g in G}
Ss = sorted(S)

check('square-set-has-16-elements', len(S) == PAPER_CLAIMS['order_S'],
      '|{g^2 : g in G}| = %d' % len(S))

N4 = closure([F['f4'], F['f5'], F['f6'], F['f7']])
check('square-set-equals-N-generated-by-f4-f5-f6-f7', N4 == S,
      'S = N = <f4,f5,f6,f7>, order %d' % len(N4))

check('square-set-closed-under-multiplication',
      all(mul(a, b) in S for a in Ss for b in Ss),
      'all %d ordered pairs tested directly, with no subgroup machinery' % (len(Ss) ** 2))

check('square-set-closed-under-inverses',
      all(inv(a) in S for a in Ss),
      'all 16 inverses lie in S')

check('generated-subgroup-equals-square-set', closure(Ss) == S,
      '|<S>| = |S| = 16, so the SET of squares is a subgroup')

_witness = [(a, b) for a in Ss for b in Ss if mul(a, b) != mul(b, a)]
check('square-subgroup-is-non-abelian', len(_witness) > 0,
      '%d ordered non-commuting pairs in S' % len(_witness))

check('the-two-exhibited-squares',
      mul(word('f1f2'), word('f1f2')) == F['f4']
      and mul(word('f1f3'), word('f1f3')) == F['f5'],
      '(f1f2)^2 = f4 and (f1f3)^2 = f5, so f4 and f5 are squares')

check('the-exhibited-non-commuting-pair',
      comm(F['f4'], F['f5']) == F['f7'] and F['f7'] != IDENT
      and order(F['f7']) == 2,
      '[f4,f5] = f7 = %s, of order 2' % cycle_str(F['f7']))

_rootbad = []
for tgt, rt, ordr in PAPER_ROOTS:
    t, r = word(tgt), word(rt)
    if mul(r, r) != t or order(r) != ordr or t not in S:
        _rootbad.append(tgt)
check('all-16-square-roots-from-table-1', not _rootbad,
      'for each of the 16 elements s of S the printed root r satisfies r^2 = s')

check('table-1-covers-S-exactly',
      {word(t) for t, _r, _o in PAPER_ROOTS} == S and len(PAPER_ROOTS) == 16,
      'the 16 targets of Table 1 are pairwise distinct and exhaust S')

check('table-1-root-orders',
      all(order(word(rt)) == ordr for _t, rt, ordr in PAPER_ROOTS),
      'the printed order of every root is correct; four roots have order 8')

_sdist = {}
for g in Ss:
    _sdist[order(g)] = _sdist.get(order(g), 0) + 1
check('order-distribution-of-S', _sdist == PAPER_CLAIMS['order_dist_S'],
      '1 identity, %d involutions, %d elements of order 4'
      % (_sdist[2], _sdist[4]))

check('exponent-of-S-is-4', max(_sdist) == PAPER_CLAIMS['exponent_S'],
      'exp S = 4')

# =============================================================================
# E. S IS C2 x D8, NOT C2 x Q8
# =============================================================================

D = closure([F['f4'], F['f5']])
_r = mul(F['f4'], F['f5'])
check('f4-f5-generate-a-dihedral-group-of-order-8',
      len(D) == PAPER_CLAIMS['order_D8_factor']
      and order(_r) == PAPER_CLAIMS['order_f4f5']
      and mul(_r, _r) == F['f7']
      and mul(mul(F['f4'], _r), F['f4']) == inv(_r)
      and F['f4'] not in closure([_r]),
      '<f4,f5> has order 8, r = f4f5 has order 4 with r^2 = f7, and f4 is an involution '
      'outside <r> inverting r: the dihedral presentation of D8')

check('S-is-the-direct-product-of-that-D8-with-<f6>',
      all(comm(F['f6'], x) == IDENT for x in Ss)
      and closure([F['f6']]) & D == {IDENT}
      and len(D) * 2 == len(S)
      and closure(sorted(D) + [F['f6']]) == S,
      'f6 is central in S, <f6> meets <f4,f5> trivially, and the product is all of S: '
      'S = D8 x C2 = C2 x D8')

# Q8 built from scratch, as the eight unit quaternions with exact integer arithmetic.
_Q8 = [(s, i) for s in (1, -1) for i in range(4)]      # (sign, one of 1, i, j, k)


def _qmul(x, y):
    """Quaternion product on {+-1, +-i, +-j, +-k}, exact integers."""
    a = [0, 0, 0, 0]
    a[x[1]] = x[0]
    b = [0, 0, 0, 0]
    b[y[1]] = y[0]
    w = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
    i = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2]
    j = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1]
    k = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]
    v = [w, i, j, k]
    nz = [t for t in range(4) if v[t] != 0]
    assert len(nz) == 1 and abs(v[nz[0]]) == 1, v
    return (v[nz[0]], nz[0])


_C2Q8 = [(c, q) for c in (0, 1) for q in _Q8]


def _c2q8_mul(x, y):
    return ((x[0] + y[0]) % 2, _qmul(x[1], y[1]))


_e = (0, (1, 0))
_inv_count = 0
for x in _C2Q8:
    if x != _e and _c2q8_mul(x, x) == _e:
        _inv_count += 1
check('C2-x-Q8-has-only-3-involutions',
      len(_C2Q8) == 16 and _inv_count == PAPER_CLAIMS['involutions_in_C2_x_Q8'],
      'the other group of order 16 with a D8-like commutator structure has %d involutions, '
      'so S (with %d) is not isomorphic to it'
      % (_inv_count, PAPER_CLAIMS['involutions_in_S']))

# =============================================================================
# F. S IS NOT POWERFUL, ON EVERY READING IN PRINT
# =============================================================================

S_der = closure([comm(a, b) for a in Ss for b in Ss])
S_sq = closure([mul(a, a) for a in Ss])
S_4 = closure([mul(mul(a, a), mul(a, a)) for a in Ss])

check('order-of-S-derived-subgroup-is-2', len(S_der) == PAPER_CLAIMS['order_S_derived'],
      "|S'| = %d, generated by [f4,f5] = f7" % len(S_der))
check('order-of-S-squared-is-2', len(S_sq) == PAPER_CLAIMS['order_S_squared'],
      '|S^2| = |<s^2 : s in S>| = %d' % len(S_sq))
check('S-fourth-power-subgroup-is-trivial', len(S_4) == PAPER_CLAIMS['order_S_fourth'],
      '|S^4| = %d, because exp S = 4' % len(S_4))

check('S-IS-NOT-POWERFUL-published-p-2-convention',
      not (S_der <= S_4),
      "S' has order 2 and S^4 = 1, so S' is NOT contained in S^4 = P_2(S): S is not powerful")

G_der = closure([comm(a, b) for a in Gs for b in Gs])
check('derived-subgroup-of-G-equals-S', G_der == S,
      "G' = S, order %d" % len(G_der))

G_sq = closure([mul(g, g) for g in Gs])
check('subgroup-generated-by-squares-equals-S', G_sq == S,
      'G^2 = <g^2> = S, so here the SET of squares already equals the subgroup it generates')

Phi = closure([comm(a, b) for a in Gs for b in Gs] + [mul(g, g) for g in Gs])
check('frattini-subgroup-of-G-equals-S', Phi == S,
      "Phi(G) = G' G^2 = S, order %d; hence G/S = C2^3 and every square lies in S" % len(Phi))

check('S-is-normal-in-G',
      all(mul(mul(inv(g), s), g) in S for g in Gs for s in Ss),
      'g^{-1} S g = S for all 128 g and all 16 s')

_cosets = {frozenset(mul(g, s) for s in Ss) for g in Gs}
check('quotient-G-mod-S-is-elementary-abelian-of-order-8',
      len(_cosets) == 8
      and all(mul(g, g) in S for g in Gs)
      and all(comm(g, h) in S for g in Gs for h in Gs),
      'S has index 8; every g^2 lies in S so G/S has exponent 2, and every commutator lies '
      'in S so G/S is abelian: G/S = C2^3, and (gS)^2 = S for every g')

Z = [g for g in Gs if all(mul(g, h) == mul(h, g) for h in Gs)]
check('centre-of-G-has-order-2', len(Z) == PAPER_CLAIMS['order_centre_G'],
      '|Z(G)| = %d' % len(Z))

_lcs, _prev = [], set(G)
while True:
    _nxt = closure([comm(a, b) for a in sorted(_prev) for b in Gs])
    _lcs.append(len(_nxt))
    if _nxt == _prev or len(_nxt) == 1:
        break
    _prev = _nxt
check('lower-central-series-of-G', _lcs == PAPER_CLAIMS['lower_central_orders'],
      'gamma_2,...,gamma_5 have orders %s, so G has nilpotency class 4'
      % '/'.join(str(x) for x in _lcs))

SG = closure([comm(s, g) for s in Ss for g in Gs])
check('order-of-commutator-S-G-is-4', len(SG) == PAPER_CLAIMS['order_commutator_S_G'],
      '|[S,G]| = %d = |gamma_3(G)|' % len(SG))

check('S-IS-NOT-POWERFULLY-EMBEDDED-published-convention',
      not (SG <= S_4),
      '[S,G] has order 4 and S^4 = 1, so [S,G] is not contained in S^4')

check('S-IS-NOT-POWERFULLY-EMBEDDED-odd-p-convention',
      not (SG <= S_sq),
      '[S,G] has order 4 and S^2 has order 2, so [S,G] is not contained in S^2 either; '
      'the refutation does not depend on which of the two embedding conventions is used')

check('the-only-surviving-reading-is-vacuous-here',
      S_der <= S_sq,
      "S' <= S^2 does hold, which is why the NON-standard reading S' <= S^2 cannot be the "
      'intended one')

_seen, _vac_bad = set(), 0
for a, b in itertools.combinations(Gs, 2):
    H = frozenset(closure([a, b]))
    if H in _seen:
        continue
    _seen.add(H)
    Hsq = closure([mul(h, h) for h in H])
    if not all(comm(x, y) in Hsq for x in H for y in H):
        _vac_bad += 1
check('vacuity-of-that-reading-on-410-subgroups',
      _vac_bad == 0 and len(_seen) == 410,
      "H' <= H^2 for all %d distinct 2-generated subgroups H of G, none of them a "
      'counterexample: the reading is empty of content' % len(_seen))

# =============================================================================
# G. THE HYPOTHESES OF THE QUESTION, AND THE ANSWER
# =============================================================================

check('hypotheses-of-kourovka-21-137-are-met',
      len(G) == 2 ** 7 and max(_dist) == 8 and closure(Ss) == S
      and all(mul(a, b) in S for a in Ss for b in Ss),
      'G is a finite 2-group (order 2^7) of exponent 8, and the set of its squares is closed '
      'under multiplication, so the squares form a subgroup')

check('answer-to-the-third-sentence-is-NO',
      len(_witness) > 0,
      'that subgroup is not abelian')

check('answer-to-the-headline-sentence-is-NO-at-p-2',
      not (S_der <= S_4),
      'that subgroup is not powerful in the published p = 2 sense')

# =============================================================================
# SCOPE
# =============================================================================

print()
print('NOT RE-RUN: the CENSUS claims of the accompanying record are outside this program. It '
      'does not verify that 128 is the least order at which the phenomenon occurs, that '
      'exactly ten groups of order 128 are witnesses, or the counts 34 at order 256 and '
      '2094 at order 512; those need the SmallGroups library. Nothing here asserts that G '
      'is the unique witness of its order.')
print('NOT RE-RUN: the library label SmallGroup(128,928). This program uses no group-theory '
      'library, so it can neither confirm nor deny that catalogue number; the object is '
      'named here only by its permutations, by Syl_2(S_8) = C2 wr C2 wr C2, and by its 28 '
      'relators.')
print('NOT RE-RUN: the general lemma that S\' <= S^2 holds in EVERY finite group. That is '
      'proved in one line in the paper; the program checks it only on the 410 distinct '
      '2-generated subgroups of this G.')
print('NOT RE-RUN: the second sentence of Kourovka 21.137, the odd-prime branch. It is a '
      'different question, no object here bears on it, and it remains open.')
print('NOT RE-RUN: every bibliographic and priority claim. Whether this observation is new '
      'is not a computation and is not tested here.')
print()

if _fails:
    print('VERDICT: %d of %d CHECKS FAILED' % (len(_fails), len(_passes) + len(_fails)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(_passes))
sys.exit(0)
