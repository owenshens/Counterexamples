#!/usr/bin/env python3
"""verify.py -- checks every computational claim of

    "Four Factorizations of K_35 Realizing beta in {1,2,4,6} in the Hamilton-Waterloo
     Problem HWP(35;5,7)"

against the objects PRINTED IN THE PAPER.  The twenty-eight generating objects below are the
cycle lists of Section 4, transcribed character for character; nothing is read from a data file
and nothing is recomputed by a solver.

Python 3.9+, standard library only.  All arithmetic is exact integer arithmetic on the vertex
set {0,...,34}; there is no floating point and no tolerance anywhere in this program.

WHAT IT CHECKS
  A. the two prescribed automorphisms rho (p=7) and sigma (p=5): that each is a permutation of
     the 35 vertices of the stated order, and that its induced action on E(K_35) has exactly the
     orbit census the paper's bookkeeping uses (85 x 7 = 595 and 119 x 5 = 595), with every
     intra-class orbit a p-cycle and every inter-class orbit a perfect matching;
  B. each printed object: that its cycles are vertex-disjoint, cover all 35 vertices, have the
     declared length, and so give a 2-factor of K_35 with 35 distinct edges;
  C. each of the four cells: that the objects declared invariant really are fixed setwise, that
     each base object's orbit closes after exactly p steps into p distinct factors, that the
     expansion has exactly 17 factors, that every factor is 2-regular on all 35 vertices with
     the declared cycle type, that the 17 factors are pairwise edge-disjoint with union EQUAL to
     E(K_35) as a set, and that the counted (alpha,beta) is the claimed cell;
  D. a negative control: the SAME routine applied to a one-edge tamper of each cell must refuse;
  E. the ARITHMETIC of assembling the spectrum: that the four cells are exactly the beta in
     {1,2,4,6} left open at v=35 and that adjoining the other fourteen values of beta -- QUOTED
     from Wang-Lu-Cao, neither displayed nor verified here -- accounts for all eighteen pairs.
     The paper draws the conclusion that all eighteen pairs are realizable only as Corollary 2,
     under that explicit unverified hypothesis; this program decides the bookkeeping, not the
     hypothesis and so not the corollary.

Check A also carries one line (p7-lane-arithmetic-not-in-paper) that has NO counterpart in the
revised paper: see the comment at its call site and item (2) of the closing SCOPE NOTE.

Exit status 0 if and only if every check passes.
"""
import sys

N = 35
EK = set((a, b) for a in range(N) for b in range(a + 1, N))
assert len(EK) == 595

FAILED = []
NPASS = 0


def ck(name, cond, detail=''):
    """One check.  Prints exactly one PASS or FAIL line."""
    global NPASS
    if cond:
        NPASS += 1
        print('PASS %s%s' % (name, (' ' + detail) if detail else ''))
    else:
        FAILED.append(name)
        print('FAIL %s%s' % (name, (' ' + detail) if detail else ''))
    return bool(cond)


# ---------------------------------------------------------------------------
# the two lanes
# ---------------------------------------------------------------------------
def rho7(v):
    c, x = divmod(v, 7)
    return 7 * c + (x + 1) % 7


def sig5(v):
    j, y = divmod(v, 5)
    return 5 * j + (y + 1) % 5


LANE = {'p7': (rho7, 7, 5), 'p5': (sig5, 5, 7)}   # map, p, number of classes


def cls(v, p):
    return v // p


def e(a, b):
    return (a, b) if a < b else (b, a)


def image(F, f):
    return set(e(f(a), f(b)) for a, b in F)


def order_of(f):
    k, cur = 1, [f(v) for v in range(N)]
    while cur != list(range(N)):
        cur = [f(v) for v in cur]
        k += 1
        if k > 100:
            return None
    return k


def edge_orbits(f, p):
    seen, orbs = set(), []
    for x in sorted(EK):
        if x in seen:
            continue
        o, cur = [], x
        for _ in range(p):
            o.append(cur)
            seen.add(cur)
            cur = e(f(cur[0]), f(cur[1]))
        orbs.append(o)
    return orbs


def parse_cycles(s):
    """The printed cycle list -> list of cycles (lists of ints).  Format is checked, not assumed."""
    out, i = [], 0
    while True:
        i = s.find('(', i)
        if i < 0:
            break
        j = s.find(')', i)
        if j < 0:
            raise ValueError('unclosed cycle in %r' % s)
        c = [int(t) for t in s[i + 1:j].split(',')]
        for v in c:
            if not (0 <= v < N):
                raise ValueError('vertex %d out of range' % v)
        out.append(c)
        i = j + 1
    return out


def edges_of(cycles):
    E = set()
    for c in cycles:
        for k in range(len(c)):
            E.add(e(c[k], c[(k + 1) % len(c)]))
    return E


def cycle_type(F):
    """The cycle type of a 2-regular spanning subgraph, walked from the edge set itself.
    Returns None unless F is 2-regular on all 35 vertices."""
    adj = {}
    for a, b in F:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    if len(adj) != N or any(len(x) != 2 for x in adj.values()):
        return None
    seen, types = set(), []
    for s0 in range(N):
        if s0 in seen:
            continue
        n, prev, cur = 0, None, s0
        while True:
            seen.add(cur)
            n += 1
            nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
            prev, cur = cur, nxt
            if cur == s0:
                break
        types.append(n)
    return sorted(types)


WANT = {'C5': sorted([5] * 7), 'C7': sorted([7] * 5)}


def cycle_type_of_class(F, p):
    """Cycle lengths of F on its own support; None unless 2-regular there."""
    adj = {}
    for a, b in F:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    if any(len(x) != 2 for x in adj.values()):
        return None
    seen, types = set(), []
    for s0 in sorted(adj):
        if s0 in seen:
            continue
        n, prev, cur = 0, None, s0
        while True:
            seen.add(cur)
            n += 1
            nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
            prev, cur = cur, nxt
            if cur == s0:
                break
        types.append(n)
    return sorted(types)


def is_matching(F, p):
    """F is a perfect matching between two distinct classes of size p."""
    if len(F) != p:
        return False
    vs = [v for x in F for v in x]
    if len(set(vs)) != 2 * p:
        return False
    cs = sorted(set(cls(v, p) for v in vs))
    if len(cs) != 2:
        return False
    for c in cs:
        if sorted(v for v in vs if cls(v, p) == c) != sorted(range(c * p, c * p + p)):
            return False
    return all(cls(a, p) != cls(b, p) for a, b in F)


# ---------------------------------------------------------------------------
# THE OBJECTS, transcribed from Section 4 of the paper
# ---------------------------------------------------------------------------
CERTS = [
    ('16,1', 'p7', [
        ('C7', 'inv',
         '(0,2,4,6,1,3,5)(7,9,11,13,8,10,12)(14,16,18,20,15,17,19)(21,24,27,23,26,22,25)(28,31,34,30,33,29,32)'),
        ('C5', 'base',
         '(0,12,11,23,31)(1,22,18,28,33)(2,8,4,5,19)(3,10,13,26,16)(6,7,24,30,25)(9,15,27,21,32)(14,17,29,20,34)'),
        ('C5', 'base',
         '(0,22,2,5,34)(1,10,28,8,26)(3,13,21,30,29)(4,27,25,17,32)(6,14,15,24,18)(7,31,20,11,33)(9,19,12,16,23)'),
        ('C5', 'inv',
         '(0,16,11,27,30)(1,17,12,21,31)(2,18,13,22,32)(3,19,7,23,33)(4,20,8,24,34)(5,14,9,25,28)(6,15,10,26,29)'),
        ('C5', 'inv',
         '(0,18,10,31,24)(1,19,11,32,25)(2,20,12,33,26)(3,14,13,34,27)(4,15,7,28,21)(5,16,8,29,22)(6,17,9,30,23)'),
    ]),
    ('15,2', 'p7', [
        ('C5', 'base',
         '(0,23,31,21,28)(1,10,20,15,32)(2,6,17,33,12)(3,11,8,4,25)(5,18,27,16,30)(7,19,26,9,22)(13,24,14,34,29)'),
        ('C5', 'base',
         '(0,17,1,20,26)(2,8,24,10,16)(3,18,9,4,21)(5,12,32,7,31)(6,28,13,14,33)(11,15,27,19,30)(22,23,29,25,34)'),
        ('C5', 'inv',
         '(0,27,8,15,30)(1,21,9,16,31)(2,22,10,17,32)(3,23,11,18,33)(4,24,12,19,34)(5,25,13,20,28)(6,26,7,14,29)'),
        ('C7', 'inv',
         '(0,1,2,3,4,5,6)(7,9,11,13,8,10,12)(14,15,16,17,18,19,20)(21,24,27,23,26,22,25)(28,31,34,30,33,29,32)'),
        ('C7', 'inv',
         '(0,2,4,6,1,3,5)(7,8,9,10,11,12,13)(14,17,20,16,19,15,18)(21,23,25,27,22,24,26)(28,29,30,31,32,33,34)'),
    ]),
    ('13,4', 'p5', [
        ('C5', 'base',
         '(0,21,9,5,33)(1,15,16,27,26)(2,25,24,13,29)(3,18,11,10,20)(4,28,19,6,31)(7,14,23,22,17)(8,30,32,12,34)'),
        ('C5', 'base',
         '(0,20,29,34,26)(1,10,27,31,30)(2,18,7,25,32)(3,11,5,13,16)(4,9,24,6,22)(8,15,14,12,17)(19,23,28,21,33)'),
        ('C5', 'inv',
         '(0,1,2,3,4)(5,10,19,32,26)(6,11,15,33,27)(7,12,16,34,28)(8,13,17,30,29)(9,14,18,31,25)(20,22,24,21,23)'),
        ('C5', 'inv',
         '(0,2,4,1,3)(5,7,9,6,8)(10,29,21,19,31)(11,25,22,15,32)(12,26,23,16,33)(13,27,24,17,34)(14,28,20,18,30)'),
        ('C5', 'inv',
         '(0,6,10,33,24)(1,7,11,34,20)(2,8,12,30,21)(3,9,13,31,22)(4,5,14,32,23)(15,17,19,16,18)(25,27,29,26,28)'),
        ('C7', 'inv',
         '(0,7,26,19,22,30,11)(1,8,27,15,23,31,12)(2,9,28,16,24,32,13)(3,5,29,17,20,33,14)(4,6,25,18,21,34,10)'),
        ('C7', 'inv',
         '(0,8,32,21,15,25,12)(1,9,33,22,16,26,13)(2,5,34,23,17,27,14)(3,6,30,24,18,28,10)(4,7,31,20,19,29,11)'),
        ('C7', 'inv',
         '(0,9,26,11,23,33,17)(1,5,27,12,24,34,18)(2,6,28,13,20,30,19)(3,7,29,14,21,31,15)(4,8,25,10,22,32,16)'),
        ('C7', 'inv',
         '(0,10,23,9,29,16,31)(1,11,24,5,25,17,32)(2,12,20,6,26,18,33)(3,13,21,7,27,19,34)(4,14,22,8,28,15,30)'),
    ]),
    ('11,6', 'p5', [
        ('C5', 'base',
         '(0,1,4,16,15)(2,5,7,24,6)(3,20,18,33,23)(8,9,27,29,22)(10,11,13,34,21)(12,17,30,31,19)(14,25,26,32,28)'),
        ('C5', 'base',
         '(0,19,26,20,34)(1,32,2,9,33)(3,12,4,10,24)(5,11,8,18,14)(6,22,7,28,31)(13,21,23,30,25)(15,17,27,16,29)'),
        ('C5', 'inv',
         '(0,5,16,13,28)(1,6,17,14,29)(2,7,18,10,25)(3,8,19,11,26)(4,9,15,12,27)(20,21,22,23,24)(30,32,34,31,33)'),
        ('C7', 'inv',
         '(0,6,19,30,12,24,29)(1,7,15,31,13,20,25)(2,8,16,32,14,21,26)(3,9,17,33,10,22,27)(4,5,18,34,11,23,28)'),
        ('C7', 'inv',
         '(0,10,8,27,24,19,33)(1,11,9,28,20,15,34)(2,12,5,29,21,16,30)(3,13,6,25,22,17,31)(4,14,7,26,23,18,32)'),
        ('C7', 'inv',
         '(0,12,32,9,29,20,16)(1,13,33,5,25,21,17)(2,14,34,6,26,22,18)(3,10,30,7,27,23,19)(4,11,31,8,28,24,15)'),
        ('C7', 'inv',
         '(0,18,21,32,13,8,25)(1,19,22,33,14,9,26)(2,15,23,34,10,5,27)(3,16,24,30,11,6,28)(4,17,20,31,12,7,29)'),
        ('C7', 'inv',
         '(0,23,17,5,31,14,27)(1,24,18,6,32,10,28)(2,20,19,7,33,11,29)(3,21,15,8,34,12,25)(4,22,16,9,30,13,26)'),
        ('C7', 'inv',
         '(0,24,14,15,6,33,26)(1,20,10,16,7,34,27)(2,21,11,17,8,30,28)(3,22,12,18,9,31,29)(4,23,13,19,5,32,25)'),
    ]),
]


def probe_cell(cell, lane, objs):
    """The six conditions that make a printed certificate a factorization of the claimed cell.
    Returns [(condition, bool, detail)] and never prints, so the SAME code serves the positive
    checks and the negative control."""
    f, p, _ = LANE[lane]
    out = []
    objE = [(kind, role, edges_of(parse_cycles(txt))) for kind, role, txt in objs]

    inv = [(k, F) for k, r, F in objE if r == 'inv']
    bad = [k for k, F in inv if image(F, f) != F]
    out.append(('invariance', not bad,
                '%d object(s) declared invariant, all fixed setwise by the lane map' % len(inv)))

    factors, closes = [], True
    for kind, role, F in objE:
        if role == 'inv':
            factors.append((kind, F))
        else:
            G, orb = F, []
            for _ in range(p):
                orb.append(G)
                G = image(G, f)
            if G != F or len(set(frozenset(x) for x in orb)) != p:
                closes = False
            factors.extend((kind, x) for x in orb)
    nbase = sum(1 for _, r, _ in objE if r == 'base')
    out.append(('base-orbits', closes,
                '%d base object(s); each orbit closes after exactly %d steps into %d distinct factors'
                % (nbase, p, p)))

    out.append(('seventeen-factors', len(factors) == 17,
                '%d objects expand to %d 2-factors' % (len(objE), len(factors))))

    types = [cycle_type(F) for _, F in factors]
    ok_t = all(t == WANT[k] for (k, _), t in zip(factors, types))
    out.append(('cycle-types', ok_t,
                'all 17 factors 2-regular on 35 vertices, cycle type [5]*7 or [7]*5 as declared'))

    union = []
    for _, F in factors:
        union.extend(F)
    ok_u = (len(union) == 595 and set(union) == EK)
    out.append(('partition', ok_u,
                '595 edges, %d distinct, set equality with E(K_35): %s'
                % (len(set(union)), set(union) == EK)))

    alpha = sum(1 for k, _ in factors if k == 'C5')
    beta = sum(1 for k, _ in factors if k == 'C7')
    out.append(('cell', '%d,%d' % (alpha, beta) == cell and alpha + beta == 17,
                'counted (alpha,beta)=(%d,%d), claimed (%s), alpha+beta=%d' % (alpha, beta, cell, alpha + beta)))
    return out


def tamper(objs):
    """Replace ONE edge of the first object by an edge not in it, by rerouting a single cycle."""
    kind, role, txt = objs[0]
    C = parse_cycles(txt)
    c0 = list(C[0])
    c0[1], c0[2] = c0[2], c0[1]          # a transposition inside one cycle changes two edges
    C = [c0] + C[1:]
    txt2 = ''.join('(' + ','.join(str(x) for x in c) + ')' for c in C)
    return [(kind, role, txt2)] + list(objs[1:])


if __name__ == '__main__':
    print('verify.py -- HWP(35;5,7): four printed factorizations of K_35 realizing beta in {1,2,4,6},')
    print('re-derived from the paper')
    print('')
    print('=== A. the two prescribed automorphisms and the orbit census ===')
    for lane, nm, blk in (('p7', 'rho', 7), ('p5', 'sigma', 5)):
        f, p, ncl = LANE[lane]
        perm = sorted(f(v) for v in range(N))
        blocks = all(cls(f(v), p) == cls(v, p) for v in range(N))
        ck('%s-permutation' % nm,
           perm == list(range(N)) and order_of(f) == p and blocks,
           'a permutation of the 35 vertices, order %d, preserving the %d blocks of size %d'
           % (p, ncl, p))
    for lane, nm, n_intra, n_inter in (('p7', 'rho', 15, 70), ('p5', 'sigma', 14, 105)):
        f, p, ncl = LANE[lane]
        orbs = edge_orbits(f, p)
        intra = [o for o in orbs if cls(o[0][0], p) == cls(o[0][1], p)]
        inter = [o for o in orbs if cls(o[0][0], p) != cls(o[0][1], p)]
        ck('%s-edge-orbits' % nm,
           (len(orbs) == n_intra + n_inter and all(len(set(o)) == p for o in orbs)
            and len(intra) == n_intra and len(inter) == n_inter
            and (n_intra + n_inter) * p == 595),
           '%d orbits of size %d (%d intra-class + %d inter-class), %d x %d = 595 = |E(K_35)|'
           % (len(orbs), p, len(intra), len(inter), len(orbs), p))
        ok_i = all(cycle_type_of_class(set(o), p) == [p] for o in intra)
        ok_j = all(is_matching(set(o), p) for o in inter)
        ck('%s-orbit-shapes' % nm, ok_i and ok_j,
           'each of the %d intra-class orbits is a %d-cycle on its class; each of the %d '
           'inter-class orbits is a perfect matching between its two classes'
           % (len(intra), p, len(inter)))
    ck('necessity', (35 - 1) // 2 == 17 and len([1 for b in range(18)]) == 18,
       'v=35 is odd, so alpha+beta=(v-1)/2=17, and there are 18 pairs with alpha,beta>=0')
    # NO COUNTERPART IN THE PAPER.  The bound B7 <= 3 used below came from a p=7-lane lemma that
    # earlier drafts carried and called "Lemma 2" (a rho-invariant C_7-factor is a union of five
    # intra-class orbits).  The revised paper DELETES that lemma: it motivated a design choice and
    # was not needed for Theorem 1, and no claim of the paper depends on it.  The line is kept and
    # labelled -- not deleted -- so the transcript stays the same length as the recorded run; what
    # it establishes is arithmetic only, and it is not evidence for anything the paper states.
    reach7 = set(7 * a + b for a in range(3) for b in range(4) if 7 * a + b <= 17)
    ck('p7-lane-arithmetic-not-in-paper',
       4 not in reach7 and 6 not in reach7 and 1 in reach7 and 2 in reach7,
       'NO COUNTERPART IN THE PAPER: taking beta = 7*A7 + B7 with B7 <= 3, beta attains 1 and 2 '
       'but never 4, 5 or 6 -- which is why (13,4) and (11,6) were built in the p=5 lane.  The '
       'bound B7 <= 3 came from a p=7-lane lemma the revised paper deletes, so this line records '
       'the arithmetic alone and no claim of the paper rests on it')
    print('')
    print('=== B. the twenty-eight printed objects ===')
    for cell, lane, objs in CERTS:
        for i, (kind, role, txt) in enumerate(objs, start=1):
            C = parse_cycles(txt)
            L = 5 if kind == 'C5' else 7
            flat = [v for c in C for v in c]
            E = edges_of(C)
            ck('object(%s)/G%d' % (cell, i),
               (all(len(c) == L for c in C) and len(C) * L == N
                and sorted(flat) == list(range(N)) and len(E) == N),
               '%d disjoint %d-cycles spanning 0..34, %d distinct edges, %s'
               % (len(C), L, len(E), 'base' if role == 'base' else 'invariant'))
    print('')
    print('=== C. the four cells ===')
    for cell, lane, objs in CERTS:
        for cond, ok, detail in probe_cell(cell, lane, objs):
            ck('(%s)/%s' % (cell, cond), ok, '[lane %s] %s' % (lane, detail))
    print('')
    print('=== D. negative control: the same routine on a one-edge-per-cell tamper ===')
    for cell, lane, objs in CERTS:
        res = probe_cell(cell, lane, tamper(objs))
        broke = [c for c, ok, _ in res if not ok]
        ck('tamper-rejected(%s)' % cell, bool(broke),
           'the tampered certificate fails: %s' % ', '.join(broke))
    print('')
    print('=== E. assembling the spectrum of HWP(35;5,7): arithmetic on a QUOTED hypothesis ===')
    ours = sorted(int(cell.split(',')[1]) for cell, _, _ in CERTS)
    ck('our-four-cells', ours == [1, 2, 4, 6],
       'the four cells settled here have beta = %s, exactly the four values Wang-Lu-Cao except at v=35'
       % ours)
    theirs = sorted(set(range(18)) - set(ours))
    ck('spectrum-complete-under-quoted-hypothesis',
       sorted(set(ours) | set(theirs)) == list(range(18)) and len(theirs) == 14,
       'their 14 values %s together with our %s give every beta in 0..17.  This is the arithmetic '
       'of putting the two halves together, on a list of fourteen values QUOTED from Wang-Lu-Cao '
       'and not verified here; so what this line supports is the paper\'s Corollary 2 -- all 18 '
       'pairs with alpha+beta=17 are realizable IF the Wang-Lu-Cao paraphrase of Section 1 holds '
       '-- and that hypothesis is not decided by this program' % (theirs, ours))
    print('')
    print('SCOPE NOTE -- what this program does and does not decide.')
    print('  DECIDED HERE, in full: that the four printed certificates are factorizations of K_35')
    print('  into 17 2-factors realizing (alpha,beta) = (16,1), (15,2), (13,4) and (11,6), i.e.')
    print('  that those four cells lie in HWP(35;5,7).  Each is a finite, self-contained check on')
    print('  the objects printed in the paper.')
    print('  NOT DECIDED HERE (1): the fourteen remaining pairs with alpha+beta=17.  Those are')
    print('  Wang-Lu-Cao\'s theorem, quoted, not re-proved; check E only does the arithmetic of')
    print('  putting the two halves together.')
    print('  NOT IN THE PAPER AT ALL (2): the p=7-lane lemma earlier drafts called Lemma 2 (a')
    print('  rho-invariant C_7-factor is a union of five intra-class orbits, whence B7 <= 3).  The')
    print('  revised paper DELETES that lemma and no claim of the paper depends on it, so the line')
    print('  named p7-lane-arithmetic-not-in-paper corresponds to nothing in the paper: it records')
    print('  an arithmetic consequence of a bound the paper no longer states.  The finite orbit')
    print('  census it once consumed is checked on its own merits in check A and is used by the')
    print('  paper\'s bookkeeping in Section 2.')
    print('  NOT DECIDED HERE (3): anything about which cells the PRINTED journal version of')
    print('  Burgess-Danziger-Traetta excepts.  That is a bibliographic question; the CONDITIONAL')
    print('  completion of the paper\'s Corollary 2 rests on Wang-Lu-Cao as quoted plus these four')
    print('  objects, and does not depend on it.')
    print('  NOT CLAIMED ANYWHERE: uniqueness, enumeration, classification, minimality or')
    print('  nonexistence for these factorizations, and no priority for the four cells -- the')
    print('  paper claims none of these, and nothing above is evidence for any of them.')
    print('')
    if FAILED:
        print('VERDICT: %d CHECK(S) FAILED: %s' % (len(FAILED), ', '.join(FAILED)))
        sys.exit(1)
    print('VERDICT: ALL %d CHECKS PASS' % NPASS)
    sys.exit(0)
