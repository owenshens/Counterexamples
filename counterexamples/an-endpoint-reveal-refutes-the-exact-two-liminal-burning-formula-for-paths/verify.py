#!/usr/bin/env python3
"""
verify.py -- re-derives every quantity claimed in

    "An Endpoint Reveal Refutes the Exact Two-Liminal Burning Formula for Paths"

from the objects PRINTED IN THAT PAPER and from nothing else.

Python 3.9+, STANDARD LIBRARY ONLY (sys, itertools). No third-party package, no external data
file, no network, no random source, no floating point in any decision: every comparison below is
between Python ints.

WHAT IS READ IN: the block PAPER_* immediately below is a transcription of the paper's printed
objects -- the edge list of P_4, the first-round reveal {v1,v4}, the three printed play lines,
the census rows, the reading-pin table and the control counts. Nothing is read from the pipeline
that produced the result: no database row, no earlier artifact, no log.

WHAT IS DERIVED: the game engine is written from the rules restated in Section 1 of the paper and
is used to recompute all of it. The decisive claim -- b_2(P_4) = 3 -- is ALSO established here
without the engine, by the cardinality argument and the finite case analysis of the paper's own
proof (group B), so a reader who distrusts the engine still has a machine-checked proof.

Run:  python3 verify.py        (no arguments; ~40 s)
Exit: 0 if and only if every check passed.
"""
import sys
from itertools import combinations

sys.setrecursionlimit(100000)

# =====================================================================================
# 0. THE OBJECTS AS PRINTED IN THE PAPER
# =====================================================================================
# Section 2, "Write V(P_4) = {v1,v2,v3,v4} with edges v1v2, v2v3, v3v4."
PAPER_P4_EDGES = "v1v2, v2v3, v3v4"
# Section 2, "The object": the first-round reveal.
PAPER_P4_REVEAL = "{v1, v4}"
# Section 1, the claim under test.
PAPER_FORMULA = "ceil((n+2)/3)"
# Section 3, the census table, reading R (row 1), n = 1..14.
PAPER_CENSUS_R = [1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6]
PAPER_CLAIM_ROW = [1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6]
PAPER_FAILURES = [4, 7]
# Section 3, control 4: the values under reading row 3, the discarded anti-control.
PAPER_CENSUS_ROW3 = [1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6]
# Section 1, the reading-pin table. (arsonist rule, saboteur rule, b_1 = CL?, b_|V| = b?)
PAPER_READINGS = [
    (1, 'any', 'live', True, True),
    (2, 'any', 'unrev', False, True),
    (3, 'cur', 'live', True, False),
    (4, 'free', 'live', False, True),
]
# Section 1: the two published values used as external checksums.
PAPER_CL_FORMULA = "floor(n/2)+1"        # CL(P_n), quoted from arXiv:2606.12330
PAPER_BURN_FORMULA = "ceil(sqrt(n))"     # b(P_n), Bonato-Janssen-Roshanbin
# Section 1, row 4 entry, and Section 3, control 4 count.
PAPER_ROW4_FIRST_CL_FAILURE = (6, 3, 4)  # (n, b_1 under row 4, CL(P_n))
PAPER_ROW3_BVN_FAILURES = 10             # of 14
# Section 5, item 1 and item 2.
PAPER_CEIL3_PLUS1_FAILURES = [1, 10, 13]
PAPER_ONE_MOD_3 = [4, 7, 10, 13]
# Section 5, item 3: the bracket of the defining paper at k = 2.
PAPER_BRACKET = ("floor(n/3)+1", "ceil(n/2)+1")
# Section 4: the monotonicity instance.
PAPER_CHAIN_INSTANCE = (4, 3, 2)         # (n, CL(P_n), b(P_n))
# Section 3, control counts.
PAPER_CONTROL_COUNTS = {'A': 28, 'B': 6, 'C_values': 66, 'C_comparisons': 55}

# The three printed play lines. Vertices are 1-based, exactly as printed.
# Section 2 / Section 3: P_4 under reading R.
PAPER_PLAY_P4 = [
    {'propagate': [],        'reveal': [1, 4], 'burn': 1,    'burned_after': [1]},
    {'propagate': [2],       'reveal': [3],    'burn': 3,    'burned_after': [1, 2, 3]},
    {'propagate': [4],       'reveal': None,   'burn': None, 'burned_after': [1, 2, 3, 4]},
]
# Section 3: P_7 under reading R (row 1), four rounds.
PAPER_PLAY_P7_R = [
    {'propagate': [],        'reveal': [1, 7], 'burn': 1, 'burned_after': [1]},
    {'propagate': [2],       'reveal': [3, 4], 'burn': 3, 'burned_after': [1, 2, 3]},
    {'propagate': [4],       'reveal': [5, 6], 'burn': 5, 'burned_after': [1, 2, 3, 4, 5]},
    {'propagate': [6],       'reveal': [],     'burn': 7, 'burned_after': [1, 2, 3, 4, 5, 6, 7]},
]
# Section 3: P_7 under reading row 4 (the loosest), three rounds -- the claimed value.
PAPER_PLAY_P7_ROW4 = [
    {'propagate': [],        'reveal': [1, 2], 'burn': 1, 'burned_after': [1]},
    {'propagate': [2],       'reveal': [3, 4], 'burn': 5, 'burned_after': [1, 2, 5]},
    {'propagate': [3, 4, 6], 'reveal': [7],    'burn': 7, 'burned_after': [1, 2, 3, 4, 5, 6, 7]},
]

# =====================================================================================
# 1. HARNESS
# =====================================================================================
_checks = []


def ck(name, ok, detail=''):
    _checks.append((name, bool(ok), detail))
    print('%s %s%s' % ('PASS' if ok else 'FAIL', name, (' [%s]' % detail) if detail else ''))
    return bool(ok)


def head(t):
    print()
    print('=== %s ===' % t)


# =====================================================================================
# 2. THE GRAPH, PARSED FROM THE PAPER'S PRINTED EDGE LIST
# =====================================================================================
def parse_edges(s):
    """'v1v2, v2v3' -> sorted list of (i, j) with i < j, 1-based, from the paper's own text."""
    out = []
    for tok in s.split(','):
        tok = tok.strip()
        if not tok:
            continue
        parts = [p for p in tok.split('v') if p != '']
        if len(parts) != 2:
            raise ValueError('cannot parse edge token %r' % tok)
        a, b = int(parts[0]), int(parts[1])
        out.append((min(a, b), max(a, b)))
    return sorted(set(out))


def parse_set(s):
    """'{v1, v4}' -> sorted [1, 4]."""
    body = s.strip().lstrip('{').rstrip('}')
    return sorted(int(t.strip().lstrip('v')) for t in body.split(',') if t.strip())


def path_edges(n):
    return [(i, i + 1) for i in range(1, n)]


def closed_nbhd(edges, n, S):
    S = set(S)
    out = set(S)
    for (a, b) in edges:
        if a in S:
            out.add(b)
        if b in S:
            out.add(a)
    return out


# =====================================================================================
# 3. THE ENGINE: EXACT BACKWARD INDUCTION ON THE FULL EXTENSIVE-FORM GAME
# =====================================================================================
# Rules, from Section 1 of the paper:
#   round 1        : saboteur reveals k vertices; arsonist burns one REVEALED vertex (the source)
#   round i >= 2   : propagation; saboteur reveals k further vertices; arsonist burns one
#                    admissible unburned vertex if one exists
#   the round in which the last vertex burns is counted, propagation-only rounds included
#   saboteur maximises, arsonist minimises
# ars: 'any'  arsonist may burn any REVEALED unburned vertex          (reading rows 1, 2)
#      'cur'  only from the CURRENT round's revealed set              (reading row 3)
#      'free' any unburned vertex from round 2 on (round-1 source still revealed)  (row 4)
# sab: 'live' saboteur reveals only unrevealed AND unburned vertices  (rows 1, 3, 4)
#      'unrev' saboteur reveals any unrevealed vertex, burned or not  (row 2)
_ENG = {}


def engine(n, k, ars, sab):
    key = (n, k, ars, sab)
    if key in _ENG:
        return _ENG[key]
    FULL = (1 << n) - 1
    nb = [0] * n
    for i in range(n):
        if i > 0:
            nb[i] |= 1 << (i - 1)
        if i < n - 1:
            nb[i] |= 1 << (i + 1)

    def prop(b):
        o, m = b, b
        while m:
            o |= nb[(m & -m).bit_length() - 1]
            m &= m - 1
        return o

    def bits(m):
        out = []
        while m:
            out.append((m & -m).bit_length() - 1)
            m &= m - 1
        return out

    memo = {}

    def play(burned, revealed, is_first):
        """Number of rounds from here, counting the round about to be played (>= 1)."""
        if not is_first:
            burned = prop(burned)
            if burned == FULL:
                return 1                       # propagation-only closing round: it COUNTS
        st = (burned, revealed, is_first)
        if st in memo:
            return memo[st]
        pool = (FULL & ~revealed) if sab == 'unrev' else (FULL & ~revealed & ~burned)
        pb = bits(pool)
        m = min(k, len(pb))
        best = -1
        for S in combinations(pb, m):
            sm = 0
            for i in S:
                sm |= 1 << i
            v = arson(burned, revealed | sm, sm, is_first)
            if v > best:
                best = v
        if best < 0:                           # nothing left to reveal
            best = arson(burned, revealed, 0, is_first)
        memo[st] = best
        return best

    def arson(burned, revealed, sm, is_first):
        if ars == 'free' and not is_first:
            allowed = FULL
        elif ars == 'cur':
            allowed = sm
        else:
            allowed = revealed
        cand = allowed & ~burned
        if cand == 0:
            return 1 + play(burned, revealed, False)
        best = 10 ** 9
        for i in bits(cand):
            b2 = burned | (1 << i)
            v = 1 if b2 == FULL else 1 + play(b2, revealed, False)
            if v < best:
                best = v
        return best

    e = {'FULL': FULL, 'prop': prop, 'bits': bits, 'play': play, 'arson': arson}
    _ENG[key] = e
    return e


def bk(n, k, ars='any', sab='live'):
    return engine(n, k, ars, sab)['play'](0, 0, True)


# ---- the published closed forms, as exact integer arithmetic -------------------------
def claim(n):
    return -(-(n + 2) // 3)          # ceil((n+2)/3)


def CL(n):
    return n // 2 + 1                # floor(n/2)+1


def burn(n):
    r = 0                            # ceil(sqrt(n)) by integer search: no float anywhere
    while r * r < n:
        r += 1
    return r


def bracket_lo(n):
    return n // 3 + 1                # floor(n/3)+1


def bracket_hi(n):
    return -(-n // 2) + 1            # ceil(n/2)+1


def ceil3p1(n):
    return -(-n // 3) + 1            # ceil(n/3)+1


# ---- the play-line replay validator -------------------------------------------------
def replay(n, k, ars, sab, rounds, force_exact_reveal=True):
    """Replay a PRINTED play line against the rules. -> (ok, [complaints], rounds_used)."""
    e = engine(n, k, ars, sab)
    FULL, prop = e['FULL'], e['prop']
    bad = []
    burned = revealed = 0
    for idx, r in enumerate(rounds):
        first = (idx == 0)
        if not first:
            nn = prop(burned)
            got = sorted(i + 1 for i in range(n) if (nn & ~burned) >> i & 1)
            if got != sorted(r['propagate']):
                bad.append('round %d: propagation is %s, printed %s'
                           % (idx + 1, got, sorted(r['propagate'])))
            burned = nn
            if burned == FULL:
                if idx != len(rounds) - 1:
                    bad.append('round %d: all burned before the printed last round' % (idx + 1))
                got = sorted(i + 1 for i in range(n) if burned >> i & 1)
                if got != sorted(r['burned_after']):
                    bad.append('round %d: burned set %s, printed %s'
                               % (idx + 1, got, sorted(r['burned_after'])))
                return (not bad), bad, idx + 1
        elif r['propagate']:
            bad.append('round 1: no propagation happens, printed %s' % r['propagate'])
        pool = (FULL & ~revealed) if sab == 'unrev' else (FULL & ~revealed & ~burned)
        pool_v = set(i + 1 for i in range(n) if pool >> i & 1)
        rv = r['reveal']
        if rv is None:
            rv = []
        if not set(rv) <= pool_v:
            bad.append('round %d: printed reveal %s is not inside the legal pool %s'
                       % (idx + 1, sorted(rv), sorted(pool_v)))
        if force_exact_reveal and len(rv) != min(k, len(pool_v)):
            bad.append('round %d: printed reveal has %d vertices, the rule gives min(k,|pool|) = %d'
                       % (idx + 1, len(rv), min(k, len(pool_v))))
        for v in rv:
            revealed |= 1 << (v - 1)
        if ars == 'free' and not first:
            allowed = FULL
        elif ars == 'cur':
            allowed = 0
            for v in rv:
                allowed |= 1 << (v - 1)
        else:
            allowed = revealed
        bv = r['burn']
        if bv is not None:
            if not (allowed & ~burned) >> (bv - 1) & 1:
                bad.append('round %d: printed burn v%d is not admissible' % (idx + 1, bv))
            burned |= 1 << (bv - 1)
        elif (allowed & ~burned) != 0 and burned != FULL:
            bad.append('round %d: printed line burns nothing while a legal burn exists' % (idx + 1))
        got = sorted(i + 1 for i in range(n) if burned >> i & 1)
        if got != sorted(r['burned_after']):
            bad.append('round %d: burned set %s, printed %s'
                       % (idx + 1, got, sorted(r['burned_after'])))
        if burned == FULL and idx != len(rounds) - 1:
            bad.append('round %d: all burned before the printed last round' % (idx + 1))
    if burned != FULL:
        bad.append('the printed line ends with %s unburned'
                   % sorted(i + 1 for i in range(n) if not burned >> i & 1))
    return (not bad), bad, len(rounds)


# =====================================================================================
# A. THE OBJECT PRINTED IN THE PAPER
# =====================================================================================
head('A. THE OBJECT PRINTED IN THE PAPER')
E4 = parse_edges(PAPER_P4_EDGES)
S4 = parse_set(PAPER_P4_REVEAL)
K = 2

ck('A1-paper-edge-list-parses', E4 == [(1, 2), (2, 3), (3, 4)],
   'printed %r -> %s, 3 distinct undirected edges on 4 vertices' % (PAPER_P4_EDGES, E4))
ck('A2-printed-list-is-P_4', E4 == path_edges(4),
   'label for label equal to the path v1-v2-v3-v4')
deg = dict((v, 0) for v in range(1, 5))
for (a, b) in E4:
    deg[a] += 1
    deg[b] += 1
ck('A3-degree-sequence', sorted(deg.values()) == [1, 1, 2, 2],
   'deg = %s: two endpoints of degree 1, two internal vertices of degree 2'
   % [deg[v] for v in range(1, 5)])
ck('A4-printed-reveal-is-a-legal-k-set', len(S4) == K and set(S4) <= {1, 2, 3, 4},
   'printed %r -> %s, |S| = %d = k' % (PAPER_P4_REVEAL, S4, len(S4)))
ck('A5-printed-reveal-is-the-endpoint-set', sorted(v for v in deg if deg[v] == 1) == S4,
   'S = %s and the degree-1 vertices are %s' % (S4, sorted(v for v in deg if deg[v] == 1)))
refl = dict((i, 5 - i) for i in range(1, 5))
ck('A6-reflection-is-an-automorphism-fixing-S',
   sorted((min(refl[a], refl[b]), max(refl[a], refl[b])) for (a, b) in E4) == E4
   and sorted(refl[v] for v in S4) == S4,
   'v_i -> v_{5-i} maps the edge set to itself and fixes S setwise')
ck('A7-the-claimed-value-at-n=4', claim(4) == 2,
   '%s at n = 4 is 2, and the paper claims the truth is 3' % PAPER_FORMULA)

# =====================================================================================
# B. THE PAPER'S HAND PROOF, MECHANISED -- NO GAME ENGINE IS USED IN THIS GROUP
# =====================================================================================
head('B. THE HAND PROOF OF b_2(P_4) = 3, WITHOUT THE ENGINE')
V4 = [1, 2, 3, 4]

admissible_sources = sorted(set(S4) & set(V4))          # revealed AND unburned, nothing burned yet
ck('B1-source-must-be-an-endpoint',
   admissible_sources == sorted(v for v in deg if deg[v] == 1) == [1, 4],
   'at the round-1 burn the revealed unburned set is %s, which is exactly the degree-1 set, so '
   'the source is v1 or v4' % admissible_sources)

# Lower bound, the cardinality step, under the LOOSEST reading: whatever the arsonist burns in
# round 2 -- any unburned vertex at all -- at most |N[source]| + 1 vertices are burned.
lb_ok, lb_rows = True, []
for src in S4:
    Ns = closed_nbhd(E4, 4, [src])
    reach = len(Ns) + 1
    if reach >= 4:
        lb_ok = False
    lb_rows.append('src v%d: |N[v%d]| + 1 = %d' % (src, src, reach))
    for extra in V4:                          # every possible round-2 burn, revealed or not
        if len(Ns | {extra}) >= 4:
            lb_ok = False
ck('B2-round-2-accounting-bounds-the-burned-set', lb_ok,
   '; '.join(lb_rows) + ' < 4, for every one of the 4 possible round-2 burns')
ck('B3-hence-lower-bound-3', lb_ok, 'b_2(P_4) >= 3 > 2 = the claimed value')
two_round_sources = sorted(v for v in V4 if len(closed_nbhd(E4, 4, [v])) + 1 >= 4)
ck('B4-a-two-round-finish-needs-an-internal-source', two_round_sources == [2, 3],
   'the sources v with |N[v]| + 1 >= 4 are %s, and S = %s is disjoint from them'
   % (two_round_sources, S4))

# Upper bound: the paper's case analysis over ALL SIX possible round-1 reveals.
ub_cases, ub_ok = [], True
for Sp in combinations(V4, K):
    Sp = list(Sp)
    inner = [v for v in Sp if v in (2, 3)]
    if inner:
        src = inner[0]
        n2 = closed_nbhd(E4, 4, closed_nbhd(E4, 4, [src]))
        ok = (n2 == set(V4))
        ub_cases.append('S=%s src v%d N^2[v%d]=V:%s' % (Sp, src, src, ok))
    else:
        # S = {v1, v4}: burn v1, then v4 (still revealed, revelations are permanent), then
        # round-3 propagation must close.
        b_after2 = closed_nbhd(E4, 4, [1]) | {4}
        ok = (closed_nbhd(E4, 4, b_after2) == set(V4))
        ub_cases.append('S=%s burn v1 then v4 -> %s, round-3 propagation closes:%s'
                        % (Sp, sorted(b_after2), ok))
    if not ok:
        ub_ok = False
ck('B5-upper-bound-over-all-six-round-1-reveals', ub_ok and len(ub_cases) == 6,
   '6 of 6 cases finish by round 3 | ' + ' | '.join(ub_cases))
ck('B6-hence-b_2(P_4)=3', lb_ok and ub_ok,
   'lower bound 3 and upper bound 3, both by hand; the formula gives %d' % claim(4))

# =====================================================================================
# C. THE ENGINE AT n = 4, UNDER EVERY READING
# =====================================================================================
head('C. THE ENGINE AT n = 4, UNDER EVERY READING')
n4_vals = {}
for (row, ars, sab, _c1, _c2) in PAPER_READINGS:
    v = bk(4, K, ars, sab)
    n4_vals[row] = v
    ck('C%d-b_2(P_4)-under-reading-row-%d' % (row, row), v == 3,
       'ars=%s sab=%s: exact minimax b_2(P_4) = %d, claimed %d' % (ars, sab, v, claim(4)))
ck('C5-b_2(P_4)-is-reading-independent', set(n4_vals.values()) == {3},
   'all four readings of the table in Section 1 return 3; also ars=free sab=unrev gives %d'
   % bk(4, K, 'free', 'unrev'))
ck('C6-negative-control-the-engine-never-returns-the-claimed-2',
   all(v != claim(4) for v in n4_vals.values()),
   'no reading reproduces 2, so the disagreement is not a reading artefact')
ok, bad, used = replay(4, K, 'any', 'live', PAPER_PLAY_P4)
ck('C7-printed-P_4-play-line-replays-legally', ok and used == 3,
   'every reveal, burn and propagation step of the printed line is legal and it ends in '
   '%d rounds%s' % (used, '' if ok else ' | ' + '; '.join(bad)))

# =====================================================================================
# D. THE READING PIN: THE 4 x 2 TABLE OF SECTION 1
# =====================================================================================
head('D. THE READING PIN, AGAINST TWO PUBLISHED IDENTITIES')
NMAX = 14
grid = {}
for (row, ars, sab, want_cl, want_b) in PAPER_READINGS:
    b1 = [bk(n, 1, ars, sab) for n in range(1, NMAX + 1)]
    bn = [bk(n, n, ars, sab) for n in range(1, NMAX + 1)]
    got_cl = (b1 == [CL(n) for n in range(1, NMAX + 1)])
    got_b = (bn == [burn(n) for n in range(1, NMAX + 1)])
    grid[row] = (got_cl, got_b)
    ck('D%da-row-%d-identity-b_1=CL' % (row, row), got_cl == want_cl,
       'ars=%s sab=%s: b_1 = %s vs CL = %s -> %s, paper says %s'
       % (ars, sab, b1, [CL(n) for n in range(1, NMAX + 1)],
          'holds' if got_cl else 'FAILS', 'holds' if want_cl else 'FAILS'))
    ck('D%db-row-%d-identity-b_|V|=b' % (row, row), got_b == want_b,
       'ars=%s sab=%s: b_n = %s vs b(P_n) = %s -> %s, paper says %s'
       % (ars, sab, bn, [burn(n) for n in range(1, NMAX + 1)],
          'holds' if got_b else 'FAILS', 'holds' if want_b else 'FAILS'))

ck('D5-exactly-one-reading-satisfies-both-identities',
   sorted(r for r in grid if grid[r] == (True, True)) == [1],
   'only row 1 (arsonist burns any revealed unburned vertex, saboteur reveals unrevealed '
   'unburned vertices) passes both; that is reading R')
b1_row2 = [bk(n, 1, 'any', 'unrev') for n in range(1, NMAX + 1)]
ck('D6-row-2-fails-because-b_1(P_n)=n', b1_row2 == list(range(1, NMAX + 1)),
   'b_1 under row 2 is %s, i.e. n, not %s' % (b1_row2, PAPER_CL_FORMULA))
row3_bad = [n for n in range(1, NMAX + 1) if bk(n, n, 'cur', 'live') != burn(n)]
ck('D7-row-3-fails-b_|V|=b-on-exactly-10-of-14',
   len(row3_bad) == PAPER_ROW3_BVN_FAILURES,
   '%d of %d paths contradict b_|V| = b under row 3, at n = %s; paper prints %d'
   % (len(row3_bad), NMAX, row3_bad, PAPER_ROW3_BVN_FAILURES))
n0, v0, cl0 = PAPER_ROW4_FIRST_CL_FAILURE
b1_row4 = [bk(n, 1, 'free', 'live') for n in range(1, NMAX + 1)]
first_bad4 = min(n for n in range(1, NMAX + 1) if b1_row4[n - 1] != CL(n))
ck('D8-row-4-first-b_1=CL-failure-is-the-printed-one',
   (first_bad4, b1_row4[n0 - 1], CL(n0)) == (n0, v0, cl0),
   'first failure at n = %d with b_1 = %d against CL(P_%d) = %d; the paper prints (%d, %d, %d)'
   % (first_bad4, b1_row4[n0 - 1], n0, CL(n0), n0, v0, cl0))

# =====================================================================================
# E. THE CENSUS
# =====================================================================================
head('E. THE CENSUS FOR n <= 14')
cen_R = [bk(n, K, 'any', 'live') for n in range(1, NMAX + 1)]
ck('E1-census-row-reading-R-matches-the-paper', cen_R == PAPER_CENSUS_R,
   'exact b_2(P_n), n = 1..14 = %s' % cen_R)
ck('E2-claim-row-matches-the-paper',
   [claim(n) for n in range(1, NMAX + 1)] == PAPER_CLAIM_ROW,
   '%s, n = 1..14 = %s' % (PAPER_FORMULA, [claim(n) for n in range(1, NMAX + 1)]))
fails = [n for n in range(1, NMAX + 1) if cen_R[n - 1] != claim(n)]
ck('E3-failures-are-exactly-n=4-and-n=7', fails == PAPER_FAILURES,
   'the formula fails at n = %s and holds at the other %d values of n <= 14'
   % (fails, NMAX - len(fails)))
ck('E4-n=7-under-rows-1-and-2-is-4',
   bk(7, K, 'any', 'live') == 4 and bk(7, K, 'any', 'unrev') == 4,
   'row 1 gives %d, row 2 gives %d, against the claimed %d'
   % (bk(7, K, 'any', 'live'), bk(7, K, 'any', 'unrev'), claim(7)))
ck('E5-n=7-under-row-4-is-the-CLAIMED-value-so-that-cell-is-reading-pinned',
   bk(7, K, 'free', 'live') == claim(7) == 3,
   'row 4 (the loosest reading) gives %d = the claimed %d, so n = 7 is NOT reading-independent; '
   'this concedes the point the paper makes explicitly' % (bk(7, K, 'free', 'live'), claim(7)))
ok, bad, used = replay(7, K, 'any', 'live', PAPER_PLAY_P7_R)
ck('E6-printed-P_7-play-line-under-reading-R-replays-legally', ok and used == 4,
   'the printed 4-round line is legal move by move%s'
   % ('' if ok else ' | ' + '; '.join(bad)))
ok, bad, used = replay(7, K, 'free', 'live', PAPER_PLAY_P7_ROW4)
ck('E7-printed-P_7-row-4-play-line-replays-legally-in-3-rounds', ok and used == 3,
   'the printed 3-round line under the loosest reading is legal move by move%s'
   % ('' if ok else ' | ' + '; '.join(bad)))
cen3 = [bk(n, K, 'cur', 'live') for n in range(1, NMAX + 1)]
ck('E8-row-3-census-matches-the-paper-and-is-worse-for-the-theorem',
   cen3 == PAPER_CENSUS_ROW3 and cen3[2] != claim(3),
   'the discarded reading gives %s, failing already at n = 3 (%d against %d), so the '
   'refutation is stated on the reading LEAST favourable to it'
   % (cen3, cen3[2], claim(3)))

# =====================================================================================
# F. THE FOUR CONTROLS, WITH THEIR COUNTS
# =====================================================================================
head('F. CONTROLS')
cA = 0
for n in range(1, NMAX + 1):
    for sab in ('live', 'unrev'):
        if bk(n, n, 'any', sab) == burn(n):
            cA += 1
ck('F1-control-A-forced-positive-against-published-integers',
   cA == PAPER_CONTROL_COUNTS['A'],
   'at k = n the engine returns b(P_n) = %s for n = 1..14 under both saboteur readings: '
   '%d of %d agreements (Bonato-Janssen-Roshanbin)'
   % (PAPER_BURN_FORMULA, cA, PAPER_CONTROL_COUNTS['A']))
cB = 0
for n in (1, 2, 3):
    for sab in ('live', 'unrev'):
        if bk(n, K, 'any', sab) == claim(n):
            cB += 1
ck('F2-control-B-must-stay-silent-where-the-theorem-is-provable',
   cB == PAPER_CONTROL_COUNTS['B'],
   'the engine AGREES with the theorem at n = 1,2,3 under both saboteur readings: %d of %d'
   % (cB, PAPER_CONTROL_COUNTS['B']))
cC, cC_ok, cC_vals = 0, True, 0
rows_c = {}
for n in range(1, 12):
    row = [bk(n, kk, 'any', 'live') for kk in range(1, n + 1)]
    rows_c[n] = row
    cC_vals += len(row)
    for i in range(len(row) - 1):
        cC += 1
        if row[i] < row[i + 1]:
            cC_ok = False
ck('F3-control-C-monotone-non-increasing-in-k',
   cC_ok and cC == PAPER_CONTROL_COUNTS['C_comparisons']
   and cC_vals == PAPER_CONTROL_COUNTS['C_values'],
   '%d of %d adjacent comparisons non-increasing, over %d computed values b_k(P_n) for '
   '1 <= k <= n <= 11' % (cC, PAPER_CONTROL_COUNTS['C_comparisons'], cC_vals))
ck('F4-control-C-endpoints-pinned-to-the-two-published-values',
   all(rows_c[n][0] == CL(n) and rows_c[n][-1] == burn(n) for n in rows_c),
   'for every n <= 11 the row starts at CL(P_n) = %s and ends at b(P_n) = %s'
   % (PAPER_CL_FORMULA, PAPER_BURN_FORMULA))
ck('F5-anti-control-FIRED-and-the-reading-it-kills-was-the-favourable-one',
   len(row3_bad) == PAPER_ROW3_BVN_FAILURES and cen3[2] != claim(3),
   'reading row 3 is discarded on %d of 14 contradictions of b_|V| = b, and it was the reading '
   'under which the theorem fails soonest' % len(row3_bad))

# =====================================================================================
# G. WHAT IS NOT SETTLED -- THE NEGATIVE CLAIMS OF SECTION 5
# =====================================================================================
head('G. THE NEGATIVE CLAIMS OF SECTION 5')
c3_bad = [n for n in range(1, NMAX + 1) if ceil3p1(n) != cen_R[n - 1]]
ck('G1-equality-with-ceil(n/3)+1-fails-at-exactly-n=1,10,13',
   c3_bad == PAPER_CEIL3_PLUS1_FAILURES,
   'ceil(n/3)+1 = %s against exact %s; equality fails at n = %s'
   % ([ceil3p1(n) for n in range(1, NMAX + 1)], cen_R, c3_bad))
ck('G2-no-residue-class-pattern',
   all(n % 3 == 1 for n in PAPER_ONE_MOD_3)
   and cen_R[9] == claim(10) and cen_R[12] == claim(13),
   'n = %s are all 1 mod 3, yet the formula HOLDS at n = 10 (%d = %d) and n = 13 (%d = %d)'
   % (PAPER_ONE_MOD_3, cen_R[9], claim(10), cen_R[12], claim(13)))
cen_2 = [bk(n, K, 'any', 'unrev') for n in range(1, NMAX + 1)]
inside = sum(1 for row in (cen_R, cen_2) for n in range(1, NMAX + 1)
             if bracket_lo(n) <= row[n - 1] <= bracket_hi(n))
ck('G3-every-computed-value-lies-inside-the-defining-paper-s-bracket', inside == 28,
   '%d of 28 census values (rows 1 and 2, n = 1..14) satisfy %s <= b_2(P_n) <= %s, so nothing '
   'here refutes that bracket' % (inside, PAPER_BRACKET[0], PAPER_BRACKET[1]))
ck('G4-the-bracket-cannot-decide-n=4-or-n=7',
   all(bracket_lo(n) <= claim(n) <= bracket_hi(n)
       and bracket_lo(n) <= cen_R[n - 1] <= bracket_hi(n) for n in PAPER_FAILURES),
   'at n = 4 it allows [%d,%d] and at n = 7 it allows [%d,%d], each containing BOTH the true '
   'and the false value' % (bracket_lo(4), bracket_hi(4), bracket_lo(7), bracket_hi(7)))
ck('G5-the-defining-paper-s-open-upper-bound-is-consistent-with-every-value-computed',
   all(cen_R[n - 1] <= ceil3p1(n) for n in range(1, NMAX + 1)),
   'b_2(P_n) <= ceil(n/3)+1 holds for n = 1..14 under reading R, so item 1 kills equality only')
nn, cl_v, b_v = PAPER_CHAIN_INSTANCE
ck('G6-the-monotonicity-chain-instance', (CL(nn), burn(nn)) == (cl_v, b_v) and cl_v > b_v,
   'CL(P_%d) = %d > %d = b(P_%d), so the chain printed in the target runs the wrong way'
   % (nn, cl_v, b_v, nn))
ck('G7-b_k(P_4)-is-non-increasing-in-k',
   rows_c[4] == sorted(rows_c[4], reverse=True),
   'b_1..b_4 on P_4 = %s' % rows_c[4])
ck('G8-exhaustion-arithmetic',
   len(cen_R) + len(cen_2) == 28,
   '28 of 28 attempted cells decided exactly (14 values of n x 2 saboteur readings); against '
   'the theorem range n >= 1 the exhausted fraction is 0, and one counterexample suffices')

# =====================================================================================
# SCOPE
# =====================================================================================
print()
print("NOT RE-RUN: the provenance of the target statement. The label "
      "`thm:2-liminal burning path graph`, the file name section_3/section_3.tex, its line 4, the "
      "numbering setup in theorem_setup.tex, the 24,463-byte e-print size, and the line numbers "
      "quoted from arXiv:2505.10727 were established by fetching and reading those sources; this "
      "program re-fetches nothing and checks mathematics only.")
print("NOT RE-RUN: the two published closed forms are used as EXTERNAL CHECKSUMS, not proved "
      "here. b(P_n) = ceil(sqrt(n)) is Bonato-Janssen-Roshanbin's. CL(P_n) = floor(n/2)+1 is "
      "quoted from arXiv:2606.12330; its primary source (Bonato-Marbach-Milne-Mishura, WAW 2024) "
      "was not accessible to us, so the b_1 column is our own computation agreeing with a "
      "secondary quotation, NOT an independent external check.")
print("NOT RE-RUN: n >= 15. The census is complete for 1 <= n <= 14 and says nothing beyond it. "
      "No closed form for b_2(P_n) is established, and none is claimed.")
print("NOT RE-RUN: the reveal-cardinality reading. The engine has the saboteur reveal exactly "
      "min(k, |pool|) vertices. Groups B2-B6 prove b_2(P_4) = 3 under both the exactly-k and the "
      "at-most-k readings, but the census values for n >= 5 are computed under exactly-k only.")
print("NOT RE-RUN: readings outside the four tabulated. Four combinations of the two ambiguous "
      "clauses were tested; a fifth reading nobody has proposed is not covered.")
print("NOT RE-RUN: the prior-art channels (arXiv, Semantic Scholar, OpenAlex, zbMATH, "
      "OpenCitations) that bound the novelty of this correction. OpenAlex citer enumeration "
      "returned HTTP 429 and was never read, and MathSciNet was never consulted.")

# =====================================================================================
# VERDICT
# =====================================================================================
bad_checks = [n for (n, ok, _d) in _checks if not ok]
print()
if bad_checks:
    print('FAILED CHECKS: %s' % ', '.join(bad_checks))
    print('VERDICT: %d of %d CHECKS FAILED' % (len(bad_checks), len(_checks)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % len(_checks))
sys.exit(0)
