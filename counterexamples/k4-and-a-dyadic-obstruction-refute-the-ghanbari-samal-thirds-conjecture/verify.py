#!/usr/bin/env python3
"""Verification program for

    "A Dyadic Obstruction, and K_4, against the Thirds Conjecture for Cubic
     Graphs under Uniform Random Embedding"

Python 3.9+, STANDARD LIBRARY ONLY, no external data file, no network. All
decisions are made in exact integer / Fraction arithmetic; floats appear only
inside format strings.

WHAT THIS PROGRAM DOES

  1. It implements the face-tracing rule of the paper's equation (2) from the
     definitions of Ghanbari & Samal (main.tex lines 117-120) and nothing else:
     Phi = T o R on the 4m arrival states (f, y, t), with
        R: (f,y,t) -> depart y along g = pi_y^t(f), keeping t
        T: -> (g, z, t*lambda(g)),  z the far end of g
        rev(f,u,t) = (f,v,-t*lambda(f)),  f = (#Phi-orbits)/2
     and, for an edge f=uv with A=(f,u,+1), B=(f,u,-1),
        f good singular <=> A ~ B,  f bad singular <=> A ~ rev(B),  else regular.

  2. It CARRIES A 64-ROW TABLE OF K_4 SIGNATURE DATA AS DATA in the string
     K4_TABLE below, and compares it ROW BY ROW against its own trace --
     signature, face count, Euler genus, six-letter class string and
     (bad,good,reg) triple. Every quantity asserted about K_4 is re-derived
     from that table or from the trace, never assumed.

  3. It exhausts 100% of the signed-rotation space Omega of theta, K_4, K_3,3
     and the 3-prism (2^(n+m) systems each), and 100% of the 2^m signatures at
     each of two fixed rotations of Q_3, Wagner V_8 and Petersen. The second
     form is licensed by the switching lemma, which this program tests rather
     than assumes: see check "switching_invariance_lemma_2_of_the_paper".

  4. It carries both polarities of control, reproduces two published integers
     computed by other authors (Gross & Furst 1987; Bender & Richmond 1990),
     re-checks the source paper's own assertion at main.tex line 119 on every
     orientable embedding, and verifies the face-count-delta form of the
     paper's Theorem 4 exhaustively on theta and K_4.

Output contract: one `PASS <name> [detail]` line per check, closing with
`VERDICT: ALL <n> CHECKS PASS`; exit 0 iff every check passed.

Runtime: a few tens of seconds (the Petersen sweep dominates).
"""

import re
import sys
from fractions import Fraction

# ---------------------------------------------------------------------------
# 0.  THE 64-ROW TABLE OF K_4 SIGNATURE DATA, carried in this program as data
#     and checked row by row against this program's own trace below.
# ---------------------------------------------------------------------------
K4_TABLE = """
  lambda(e0..e5) faces Eulergenus  classes(e0..e5)  (bad,good,reg)
     ++++++        2     2          RBRRBR        (2,0,4)
     -+++++        1     3          GBGGBG        (2,4,0)
     +-++++        2     2          RBRRGR        (1,1,4)
     --++++        1     3          GBGBGB        (3,3,0)
     ++-+++        1     3          GBGGBG        (2,4,0)
     -+-+++        2     2          RBRGRG        (1,2,3)
     +--+++        1     3          GBGBGB        (3,3,0)
     ---+++        2     2          RBRBRB        (3,0,3)
     +++-++        1     3          GBGGBG        (2,4,0)
     -++-++        2     2          RRGRBG        (1,2,3)
     +-+-++        1     3          BBBGGG        (3,3,0)
     --+-++        1     3          BGGBBG        (3,3,0)
     ++--++        2     2          GRRRRG        (0,2,4)
     -+--++        3     1          RRRRRG        (0,1,5)
     +---++        1     3          GGBBGG        (2,4,0)
     ----++        2     2          RGRBRG        (1,2,3)
     ++++-+        2     2          RGRRBR        (1,1,4)
     -+++-+        1     3          GGBGBB        (3,3,0)
     +-++-+        3     1          RRRRRR        (0,0,6)
     --++-+        2     2          GRRRRB        (1,1,4)
     ++-+-+        1     3          BGGBBG        (3,3,0)
     -+-+-+        1     3          BBBGGG        (3,3,0)
     +--+-+        2     2          RRGBRR        (1,1,4)
     ---+-+        1     3          GBGBGB        (3,3,0)
     +++--+        1     3          GGBGBB        (3,3,0)
     -++--+        2     2          RRBRBB        (3,0,3)
     +-+--+        2     2          RRBGRR        (1,1,4)
     --+--+        1     3          GGBGBB        (3,3,0)
     ++---+        1     3          GGBBGG        (2,4,0)
     -+---+        2     2          RRBRGG        (1,2,3)
     +----+        2     2          RRBBRR        (2,0,4)
     -----+        1     3          GGBBGG        (2,4,0)
     +++++-        1     3          GBGGBG        (2,4,0)
     -++++-        2     2          RRGGRR        (0,2,4)
     +-+++-        1     3          BBBGGG        (3,3,0)
     --+++-        1     3          BGGGGB        (2,4,0)
     ++-++-        2     2          GRRGBR        (1,2,3)
     -+-++-        3     1          RRRGRR        (0,1,5)
     +--++-        1     3          GGBGBB        (3,3,0)
     ---++-        2     2          RGRGRB        (1,2,3)
     +++-+-        2     2          GBGRRR        (1,2,3)
     -++-+-        3     1          RRGRRR        (0,1,5)
     +-+-+-        2     2          BBBRRR        (3,0,3)
     --+-+-        2     2          BGGRRR        (1,2,3)
     ++--+-        3     1          GRRRRR        (0,1,5)
     -+--+-        4     0          RRRRRR        (0,0,6)
     +---+-        2     2          GGBRRR        (1,2,3)
     ----+-        3     1          RGRRRR        (0,1,5)
     ++++--        1     3          BGGBBG        (3,3,0)
     -+++--        1     3          BGGGGB        (2,4,0)
     +-++--        2     2          BRRRRG        (1,1,4)
     --++--        2     2          BRRRRB        (2,0,4)
     ++-+--        2     2          BRRBBR        (3,0,3)
     -+-+--        2     2          BRRGGR        (1,2,3)
     +--+--        1     3          BGGBBG        (3,3,0)
     ---+--        1     3          BGGGGB        (2,4,0)
     +++---        1     3          GBGBGB        (3,3,0)
     -++---        2     2          RRGRGB        (1,2,3)
     +-+---        1     3          BBBGGG        (3,3,0)
     --+---        1     3          BGGGGB        (2,4,0)
     ++----        2     2          GRRBGR        (1,2,3)
     -+----        3     1          RRRRGR        (0,1,5)
     +-----        1     3          GGBBGG        (2,4,0)
     ------        2     2          RGRRGR        (0,2,4)
  COLUMN SUMS over the 64 signatures at pi_0: bad=108 good=138 regular=138
"""

ROW_RE = re.compile(r'^\s*([+-]{6})\s+(\d+)\s+(\d+)\s+([BGR]{6})\s+\((\d+),(\d+),(\d+)\)\s*$')

# The eight orientable (coboundary) signatures of K_4 at pi_0, carried here as
# data and recomputed from the 16 vertex subsets below.
K4_COBOUNDARIES = ['++++++', '---+++', '-++--+', '+-+-+-',
                         '++-+--', '+----+', '-+--+-', '--++--']

# ---------------------------------------------------------------------------
# 1.  THE GRAPHS OF THE LADDER
# ---------------------------------------------------------------------------
K4_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
THETA_EDGES = [(0, 1), (0, 1), (0, 1)]
K33_EDGES = [(a, b) for a in (0, 1, 2) for b in (3, 4, 5)]
PRISM_EDGES = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (0, 3), (1, 4), (2, 5)]
Q3_EDGES = sorted((v, v ^ (1 << b)) for v in range(8) for b in range(3) if v < (v ^ (1 << b)))
V8_EDGES = [(i, (i + 1) % 8) for i in range(8)] + [(i, i + 4) for i in range(4)]
PETERSEN_EDGES = ([(i, (i + 1) % 5) for i in range(5)]
                  + [(i, i + 5) for i in range(5)]
                  + [(5 + i, 5 + ((i + 2) % 5)) for i in range(5)])


def _norm(edges):
    return [(min(u, v), max(u, v)) for (u, v) in edges]


# name, n, edges, mode, extra rotations sampled in 'fixed' mode
CELLS = [
    ('theta',    2, THETA_EDGES,    'full',  ()),
    ('K4',       4, K4_EDGES,       'full',  ()),
    ('K33',      6, K33_EDGES,      'full',  ()),
    ('prism',    6, PRISM_EDGES,    'full',  ()),
    ('Q3',       8, Q3_EDGES,       'fixed', (0, 0b10101010)),
    ('V8',       8, V8_EDGES,       'fixed', (0, 0b11111111)),
    ('petersen', 10, PETERSEN_EDGES, 'fixed', (0, 0b1111111111)),
]

# The expected values for the ladder, as exact rationals / integers.
EXPECT = {
    #            E[bad]                E[good]=E[reg]        per-edge P[bad] multiset with multiplicity
    'theta':    (Fraction(3, 4),      Fraction(9, 8),      {Fraction(1, 4): 3}),
    'K4':       (Fraction(27, 16),    Fraction(69, 32),    {Fraction(9, 32): 6}),
    'K33':      (Fraction(351, 128),  Fraction(801, 256),  {Fraction(39, 128): 9}),
    'prism':    (Fraction(693, 256),  Fraction(1611, 512), {Fraction(69, 256): 6, Fraction(93, 256): 3}),
    'Q3':       (Fraction(939, 256),  Fraction(2133, 512), {Fraction(313, 1024): 12}),
    'V8':       (Fraction(1887, 512), Fraction(4257, 1024),
                 {Fraction(653, 2048): 8, Fraction(581, 2048): 4}),
    'petersen': (Fraction(18945, 4096), Fraction(42495, 8192), {Fraction(1263, 4096): 15}),
}
EDGE_TRANSITIVE = {'theta', 'K4', 'K33', 'Q3', 'petersen'}
# relative deviation of part (1) (percent, 3 dp)
EXPECT_RELDEV = {'theta': '-25.000', 'K4': '-15.625', 'K33': '-8.594', 'prism': '-9.766',
                 'Q3': '-8.301', 'V8': '-7.861', 'petersen': '-7.495'}
# E[#singular] = E[bad] + E[good]
EXPECT_SINGULAR = {'K4': Fraction(123, 32), 'prism': Fraction(2997, 512),
                   'petersen': Fraction(80385, 8192)}
# aggregate E[bad]/m, to 6 dp
EXPECT_AGG = {'theta': '0.250000', 'K4': '0.281250', 'K33': '0.304688', 'prism': '0.300781',
              'Q3': '0.305664', 'V8': '0.307129', 'petersen': '0.308350'}


# ---------------------------------------------------------------------------
# 2.  THE TRACER
# ---------------------------------------------------------------------------
class Cell(object):
    """Pre-computed dart combinatorics for one graph. State index of (e, y, tb)
    is ((e*2 + y)*2 + tb) with tb = 0 for t = +1 and tb = 1 for t = -1."""

    def __init__(self, n, edges):
        self.n = n
        self.edges = _norm(edges)
        self.m = m = len(self.edges)
        self.NS = 4 * m
        inc = [[] for _ in range(n)]
        for i, (u, v) in enumerate(self.edges):
            inc[u].append((i, 0))
            inc[v].append((i, 1))
        for w in range(n):
            if len(inc[w]) != 3:
                raise SystemExit('vertex %d has degree %d; the ladder is cubic' % (w, len(inc[w])))
        pos = {}
        for w in range(n):
            for k, it in enumerate(inc[w]):
                pos[it] = k
        self.inc = inc
        # vertex of each state, and the two departure tables (rotation bit 0 / 1)
        self.vs = [0] * self.NS
        self.base = [[0] * self.NS, [0] * self.NS]   # arrival state with tb = 0
        self.gof = [[0] * self.NS, [0] * self.NS]    # edge traversed
        for e in range(m):
            for y in (0, 1):
                w = self.edges[e][y]
                k = pos[(e, y)]
                for tb in (0, 1):
                    t = 1 if tb == 0 else -1
                    s = ((e * 2 + y) * 2) + tb
                    self.vs[s] = w
                    for bit in (0, 1):
                        step = t if bit == 0 else -t
                        g, j = inc[w][(k + step) % 3]
                        self.base[bit][s] = ((g * 2 + (1 - j)) * 2)
                        self.gof[bit][s] = g

    def rot_tables(self, rotbits):
        """(BASE, G) specialised to one rotation system."""
        BASE = [0] * self.NS
        G = [0] * self.NS
        b0, b1 = self.base
        g0, g1 = self.gof
        for s in range(self.NS):
            if (rotbits >> self.vs[s]) & 1:
                BASE[s] = b1[s]
                G[s] = g1[s]
            else:
                BASE[s] = b0[s]
                G[s] = g0[s]
        return BASE, G

    def trace(self, BASE, G, lbit):
        """-> (faces, class string). lbit[e] = 0 for lambda(e)=+1, 1 for -1."""
        NS = self.NS
        nxt = [0] * NS
        for s in range(NS):
            nxt[s] = BASE[s] + ((s & 1) ^ lbit[G[s]])
        orb = [-1] * NS
        c = 0
        for s in range(NS):
            if orb[s] < 0:
                x = s
                while orb[x] < 0:
                    orb[x] = c
                    x = nxt[x]
                c += 1
        out = []
        for f in range(self.m):
            a = (f * 4)                      # (f, 0, +1)
            b = (f * 4) + 1                  # (f, 0, -1)
            rb = (f * 4) + 2 + lbit[f]       # rev(B) = (f, 1, lambda(f))
            oa = orb[a]
            if oa == orb[b]:
                out.append('G')
            elif oa == orb[rb]:
                out.append('B')
            else:
                out.append('R')
        return c // 2, ''.join(out)

    def coboundary_lbits(self):
        """The 2^n vertex cuts, as lbit tuples (the orientable signatures)."""
        out = set()
        for S in range(1 << self.n):
            out.add(tuple(1 if (((S >> u) & 1) != ((S >> v) & 1)) else 0 for (u, v) in self.edges))
        return out


# ---------------------------------------------------------------------------
# 3.  CHECK BOOKKEEPING
# ---------------------------------------------------------------------------
_N = [0]
_F = [0]
_LINES = []


def check(name, ok, detail=''):
    _N[0] += 1
    if not ok:
        _F[0] += 1
    _LINES.append('%s %s%s' % ('PASS' if ok else 'FAIL', name, (' [%s]' % detail) if detail else ''))
    print(_LINES[-1])
    sys.stdout.flush()


def pct(fr):
    return '%.3f' % (float(fr) * 100.0)


# ---------------------------------------------------------------------------
# 4.  THE CENSUS
# ---------------------------------------------------------------------------
def census(cell, rotbits_list):
    """Exhaust all 2^m signatures at each rotation in rotbits_list.

    -> dict with per-rotation edge counters and violation counters.
    """
    m, n = cell.m, cell.n
    cob = cell.coboundary_lbits()
    per_rot = []
    viol_sum = viol_euler = viol_l119 = 0
    genus0 = 0
    good_seen = 0
    for rotbits in rotbits_list:
        BASE, G = cell.rot_tables(rotbits)
        bad = [0] * m
        good = [0] * m
        reg = [0] * m
        for k in range(1 << m):
            lbit = [(k >> i) & 1 for i in range(m)]
            faces, cs = cell.trace(BASE, G, lbit)
            nb = cs.count('B')
            ng = cs.count('G')
            nr = cs.count('R')
            if nb + ng + nr != m:
                viol_sum += 1
            eg = 2 - (n - m + faces)
            if eg < 0 or faces != 2 - eg - n + m:
                viol_euler += 1
            if tuple(lbit) in cob:
                if ng != 0:
                    viol_l119 += 1
                if eg % 2 != 0:
                    viol_euler += 1
            if eg == 0:
                genus0 += 1
            if ng:
                good_seen += 1
            for i, ch in enumerate(cs):
                if ch == 'B':
                    bad[i] += 1
                elif ch == 'G':
                    good[i] += 1
                else:
                    reg[i] += 1
        per_rot.append((bad, good, reg))
    return {'per_rot': per_rot, 'viol_sum': viol_sum, 'viol_euler': viol_euler,
            'viol_l119': viol_l119, 'genus0': genus0, 'good_seen': good_seen,
            'n_traced': len(rotbits_list) << m}


def switching_violations(cell, rotbits_list, cap):
    """Switch at EVERY vertex of the first `cap` signatures of each listed
    rotation and count class changes. This is the paper's Lemma 2 tested rather
    than assumed."""
    m = cell.m
    at = [[] for _ in range(cell.n)]
    for i, (u, v) in enumerate(cell.edges):
        at[u].append(i)
        at[v].append(i)
    bad = 0
    tested = 0
    lim = min(cap, 1 << m)
    for rotbits in rotbits_list:
        BASE, G = cell.rot_tables(rotbits)
        for v in range(cell.n):
            B2, G2 = cell.rot_tables(rotbits ^ (1 << v))
            for k in range(lim):
                lbit = [(k >> i) & 1 for i in range(m)]
                _, cs = cell.trace(BASE, G, lbit)
                lb2 = list(lbit)
                for i in at[v]:
                    lb2[i] ^= 1
                _, cs2 = cell.trace(B2, G2, lb2)
                tested += 1
                if cs2 != cs:
                    bad += 1
    return bad, tested


def flip_violations(cell, rotbits_list):
    """The paper's Theorem 4 as a per-embedding statement: flipping lambda(e)
    exchanges regular with good singular, fixes bad singular, and moves the
    face count by -1 / +1 / 0 respectively."""
    m = cell.m
    bad = 0
    tested = 0
    want = {('R', 'G'): -1, ('G', 'R'): +1, ('B', 'B'): 0}
    for rotbits in rotbits_list:
        BASE, G = cell.rot_tables(rotbits)
        for k in range(1 << m):
            lbit = [(k >> i) & 1 for i in range(m)]
            f0, cs0 = cell.trace(BASE, G, lbit)
            for e in range(m):
                lb2 = list(lbit)
                lb2[e] ^= 1
                f1, cs1 = cell.trace(BASE, G, lb2)
                tested += 1
                key = (cs0[e], cs1[e])
                if key not in want or (f1 - f0) != want[key]:
                    bad += 1
    return bad, tested


# ---------------------------------------------------------------------------
# 5.  MAIN
# ---------------------------------------------------------------------------
def main():
    print('verification of "A Dyadic Obstruction, and K_4, against the Thirds '
          'Conjecture for Cubic Graphs under Uniform Random Embedding"')
    print('python %s; exact integer and Fraction arithmetic only'
          % sys.version.split()[0])
    print('')

    # ---- 5.1  the printed table, read as data ---------------------------
    tab = []
    for line in K4_TABLE.splitlines():
        mm = ROW_RE.match(line)
        if mm:
            sig, f, eg, cs, b, g, r = mm.groups()
            tab.append((sig, int(f), int(eg), cs, (int(b), int(g), int(r))))
    check('carried_table_has_64_rows', len(tab) == 64, '%d rows parsed' % len(tab))
    check('carried_table_has_64_distinct_signatures', len({t[0] for t in tab}) == 64,
          '%d distinct' % len({t[0] for t in tab}))
    bs = sum(1 for t in tab if (t[3].count('B'), t[3].count('G'), t[3].count('R')) != t[4])
    check('carried_table_class_string_agrees_with_its_triple', bs == 0, '%d mismatches' % bs)
    bs = sum(1 for t in tab if sum(t[4]) != 6)
    check('carried_table_bad_plus_good_plus_reg_equals_m_6', bs == 0, '%d violations' % bs)
    bs = sum(1 for t in tab if t[1] != 4 - t[2])
    check('carried_table_euler_f_equals_4_minus_eulergenus', bs == 0, '%d violations' % bs)
    TB = sum(t[4][0] for t in tab)
    TG = sum(t[4][1] for t in tab)
    TR = sum(t[4][2] for t in tab)
    check('carried_table_column_sums_108_138_138', (TB, TG, TR) == (108, 138, 138),
          '%d/%d/%d' % (TB, TG, TR))
    check('carried_table_total_equals_64_times_m', TB + TG + TR == 64 * 6, '%d' % (TB + TG + TR))

    # ---- 5.2  the ladder ------------------------------------------------
    cells = {}
    results = {}
    for name, n, edges, mode, extra in CELLS:
        c = Cell(n, edges)
        cells[name] = c
        rots = list(range(1 << n)) if mode == 'full' else list(extra)
        res = census(c, rots)
        results[name] = res
        m = c.m
        Ntot = 1 << m                                    # signatures per rotation
        bad0, good0, reg0 = res['per_rot'][-1]
        # gauge / reduction: every rotation must give the same per-edge counters
        same = all(pr == res['per_rot'][0] for pr in res['per_rot'])
        check('%s_per_rotation_counters_identical' % name, same,
              '%d rotation system(s) traced, %s' % (len(rots), 'all equal' if same else 'DIFFER'))
        check('%s_coverage' % name, res['n_traced'] == len(rots) * Ntot,
              '%d embeddings traced = %d of the 2^%d rotations x all 2^%d signatures (%s)'
              % (res['n_traced'], len(rots), n, m,
                 ('100%% of Omega = 2^%d' % (n + m)) if mode == 'full'
                 else 'exact by Lemma 2, tested below'))
        check('%s_bad_plus_good_plus_reg_equals_m_every_embedding' % name,
              res['viol_sum'] == 0, '%d violations' % res['viol_sum'])
        check('%s_euler_and_orientable_parity_every_embedding' % name,
              res['viol_euler'] == 0, '%d violations' % res['viol_euler'])
        check('%s_main_tex_L119_good_equals_0_on_every_orientable_embedding' % name,
              res['viol_l119'] == 0, '%d violations' % res['viol_l119'])

        Ebad = sum(Fraction(x, Ntot) for x in bad0)
        Egood = sum(Fraction(x, Ntot) for x in good0)
        Ereg = sum(Fraction(x, Ntot) for x in reg0)
        third = Fraction(m, 3)
        eb, eg_, orb = EXPECT[name]
        check('%s_E_bad_equals_%s' % (name, str(eb).replace('/', '_over_')), Ebad == eb,
              '%s = %.6f' % (Ebad, Ebad))
        check('%s_E_good_equals_E_reg_equals_%s' % (name, str(eg_).replace('/', '_over_')),
              Egood == eg_ and Ereg == eg_, '%s and %s' % (Egood, Ereg))
        check('%s_all_three_parts_differ_from_m_over_3' % name,
              Ebad != third and Egood != third and Ereg != third,
              'm/3 = %s vs (%s, %s, %s)' % (third, Ebad, Egood, Ereg))
        check('%s_m_over_3_not_a_member_of_the_three_expectations' % name,
              third not in {Ebad, Egood, Ereg}, 'no relabelling rescues any part')
        check('%s_E_bad_plus_E_good_plus_E_reg_equals_m' % name, Ebad + Egood + Ereg == m,
              '= %s' % (Ebad + Egood + Ereg))
        pe = {}
        for i in range(m):
            pe[Fraction(bad0[i], Ntot)] = pe.get(Fraction(bad0[i], Ntot), 0) + 1
        check('%s_per_edge_P_bad_orbit_values' % name, pe == orb,
              ', '.join('%s x%d' % (k, v) for k, v in sorted(pe.items())))
        pgr = sum(1 for i in range(m) if good0[i] != reg0[i])
        check('%s_P_good_equals_P_reg_on_every_edge' % name, pgr == 0,
              '%d edges disagree (Theorem 4)' % pgr)
        rd = (Ebad - third) / third
        check('%s_relative_deviation_of_part_1' % name, pct(rd) == EXPECT_RELDEV[name],
              '%s%%' % pct(rd))
        check('%s_aggregate_E_bad_over_m' % name,
              ('%.6f' % float(Ebad / m)) == EXPECT_AGG[name], '%.6f' % float(Ebad / m))
        if name in EXPECT_SINGULAR:
            check('%s_E_singular_equals_%s' % (name, str(EXPECT_SINGULAR[name]).replace('/', '_over_')),
                  Ebad + Egood == EXPECT_SINGULAR[name],
                  '%s = %.6f vs m/3 = %s' % (Ebad + Egood, Ebad + Egood, third))
        if name in EDGE_TRANSITIVE:
            ps = set(pe)
            p = next(iter(ps))
            den = p.denominator
            check('%s_dyadic_obstruction_arithmetic' % name,
                  len(ps) == 1 and (den & (den - 1)) == 0 and p != Fraction(1, 3),
                  'single per-edge value %s, denominator 2^%d, 3 divides no power of 2'
                  % (p, den.bit_length() - 1))
        else:
            check('%s_has_two_edge_orbits_so_the_dyadic_theorem_does_not_reach_it' % name,
                  len(pe) == 2, '%d distinct per-edge values' % len(pe))
        results[name]['E'] = (Ebad, Egood, Ereg)

    # ---- 5.3  the printed table against the engine ----------------------
    K4 = cells['K4']
    BASE, G = K4.rot_tables(0)
    mism = []
    for sig, f, eg, cs, tri in tab:
        lbit = [1 if ch == '-' else 0 for ch in sig]
        f2, cs2 = K4.trace(BASE, G, lbit)
        eg2 = 2 - (4 - 6 + f2)
        tri2 = (cs2.count('B'), cs2.count('G'), cs2.count('R'))
        if (f2, eg2, cs2, tri2) != (f, eg, cs, tri):
            mism.append(sig)
    check('carried_table_reproduced_row_by_row_from_the_definitions', not mism,
          '64 rows x (faces, Euler genus, class string, triple), %d mismatches' % len(mism))

    # ---- 5.4  the orientable block, Gross-Furst, Bender-Richmond --------
    cob = K4.coboundary_lbits()
    cob_sigs = sorted(''.join('-' if b else '+' for b in t) for t in cob)
    check('K4_eight_coboundary_signatures_recomputed_from_the_16_vertex_subsets',
          cob_sigs == sorted(K4_COBOUNDARIES),
          '%d signatures, 2^n/|ker| = 16/2' % len(cob_sigs))
    byrow = {t[0]: t for t in tab}
    orient = [byrow[s] for s in cob_sigs]
    check('K4_good_equals_0_on_all_8_orientable_rows', all(t[4][1] == 0 for t in orient),
          'good values %s' % sorted({t[4][1] for t in orient}))
    check('K4_euler_genus_even_on_all_8_orientable_rows', all(t[2] % 2 == 0 for t in orient),
          '%s' % sorted(t[2] for t in orient))
    e0 = sum(1 for t in orient if t[2] == 0)
    e2 = sum(1 for t in orient if t[2] == 2)
    check('gross_furst_1987_K4_genus_2_to_14_appears_as_1_to_7_at_pi_0', (e0, e2) == (1, 7),
          'sphere %d : torus %d' % (e0, e2))
    # the same published integer over the FULL space: 16 rotations x 8 coboundaries
    g0 = g2 = 0
    for rb in range(16):
        Bb, Gg = K4.rot_tables(rb)
        for t in cob:
            f2, _ = K4.trace(Bb, Gg, list(t))
            eg2 = 2 - (4 - 6 + f2)
            if eg2 == 0:
                g0 += 1
            elif eg2 == 2:
                g2 += 1
    check('gross_furst_1987_over_the_full_orientable_space_16_to_112', (g0, g2) == (16, 112),
          '%d sphere : %d torus = 8 x (2 : 14)' % (g0, g2))
    bound = {0: 0, 1: 1}
    att = {}
    okbr = True
    for _sig, f, eg, cs, tri in tab:
        bd = bound.get(eg, 3 * eg - 3)
        sing = tri[0] + tri[1]
        if sing > bd:
            okbr = False
        a = att.setdefault(eg, [bd, 0, False])
        a[1] += 1
        if sing == bd:
            a[2] = True
    check('bender_richmond_1990_bound_holds_on_all_64_rows', okbr,
          '; '.join('eulergenus %d: %d rows, bound %d, attained %s'
                    % (k, v[1], v[0], v[2]) for k, v in sorted(att.items())))
    check('bender_richmond_1990_attained_in_all_four_surface_classes',
          all(v[2] for v in att.values()), '4 of 4')
    agg_s = sum(t[4][0] + t[4][1] for t in tab)
    agg_b = sum(att[t[2]][0] for t in tab)
    check('bender_richmond_1990_aggregate_246_le_259', (agg_s, agg_b) == (246, 259) and agg_s <= agg_b,
          '%d <= %d' % (agg_s, agg_b))

    # ---- 5.5  the orientable reading of "random embedding" --------------
    ob = Fraction(sum(t[4][0] for t in orient), 8)
    og = Fraction(sum(t[4][1] for t in orient), 8)
    orr = Fraction(sum(t[4][2] for t in orient), 8)
    check('K4_orientable_model_all_three_expectations_differ_from_m_over_3',
          (ob, og, orr) == (Fraction(9, 4), Fraction(0), Fraction(15, 4))
          and ob != 2 and og != 2 and orr != 2 and ob + og + orr == 6,
          'E[bad]=%s > 2, E[good]=%s, E[reg]=%s; the sign of the deviation FLIPS' % (ob, og, orr))

    # ---- 5.6  controls, both polarities ---------------------------------
    r4 = results['K4']
    check('whitney_K4_has_exactly_2_to_the_n_equals_16_sphere_embeddings',
          r4['genus0'] == 16, '%d of 1024 systems have Euler genus 0' % r4['genus0'])
    check('whitney_K33_is_non_planar_zero_sphere_embeddings',
          results['K33']['genus0'] == 0, '%d of 32768' % results['K33']['genus0'])
    check('opposite_polarity_the_detector_can_output_good_singular',
          r4['good_seen'] == 880,
          '%d of the 1024 K_4 systems carry a good singular edge, so no negative is vacuous'
          % r4['good_seen'])
    # forced positive: the planar tetrahedron 0:(1,2,3) 1:(0,3,2) 2:(0,1,3) 3:(0,2,1), lambda = +1
    Bt, Gt = K4.rot_tables((1 << 1) | (1 << 3))
    ft, cst = K4.trace(Bt, Gt, [0] * 6)
    check('forced_positive_control_planar_tetrahedron_is_all_regular',
          (ft, 2 - (4 - 6 + ft), cst) == (4, 0, 'RRRRRR'),
          'f=%d, Euler genus %d, classes %s -- and the m/3 detector fires (0 != 2, 6 != 2)'
          % (ft, 2 - (4 - 6 + ft), cst))

    # ---- 5.7  the two lemmas of the paper, tested -----------------------
    sw_bad = sw_tested = 0
    for name, n, edges, mode, extra in CELLS:
        if name in ('theta', 'K4'):
            rots, cap = list(range(1 << n)), (1 << cells[name].m)   # exhaustive over all of Omega
        else:
            rots, cap = (list(range(1 << n)) if mode == 'full' else list(extra))[:2], 256
        b, t = switching_violations(cells[name], rots, cap)
        sw_bad += b
        sw_tested += t
    check('switching_invariance_lemma_2_of_the_paper', sw_bad == 0,
          '%d single-vertex switchings tested over all 7 cells (exhaustive on theta and K_4, '
          '2 rotations x 256 signatures x every vertex elsewhere), %d class changes; this is what '
          'licenses averaging at one fixed rotation' % (sw_tested, sw_bad))
    fl_bad = fl_tested = 0
    for name in ('theta', 'K4'):
        rots = list(range(1 << cells[name].n))
        b, t = flip_violations(cells[name], rots)
        fl_bad += b
        fl_tested += t
    check('theorem_4_involution_and_face_count_delta', fl_bad == 0,
          '%d (embedding, edge) flips on theta and K_4, %d violations of '
          'regular<->good / bad fixed and delta f = -1/+1/0' % (fl_tested, fl_bad))

    # ---- 5.8  the computed expectations against m/3 ---------------------
    Eb, Eg, Er = results['K4']['E']
    check('K4_all_three_expectations_differ_from_m_over_3',
          (Eb, Eg, Er) == (Fraction(27, 16), Fraction(69, 32), Fraction(69, 32))
          and Fraction(2) not in {Eb, Eg, Er},
          'E[bad]=27/16, E[good]=E[reg]=69/32, m/3=2; refutation of the conjecture assumes, '
          'not verified here, that its three parts assert E[bad]=E[good]=E[reg]=m/3 for simple '
          'cubic 3-connected graphs; simplicity and 3-connectivity of the hard-coded K_4 edge '
          'list are NOT CHECKED')
    check('all_seven_cells_have_all_three_expectations_different_from_m_over_3',
          all(results[k]['E'][j] != Fraction(cells[k].m, 3) for k in results for j in (0, 1, 2)),
          '7 of 7 cells, 3 of 3 computed expectations, 0 exceptions; theta/K_4/K_3,3/prism over '
          'all of Omega, Q_3/V_8/Petersen at 2 fixed rotations (exact by Lemma 2); the '
          'identification of these three expectations with the three parts of the conjecture is '
          'NOT VERIFIED here')
    print('')
    print('NOTE SCOPE. Exhausted 100% of Omega for theta, K_4, K_3,3 and the 3-prism, and 100% of '
          'the 2^m signatures at two fixed rotations for Q_3, Wagner V_8 and Petersen; the latter '
          'is exact by Lemma 2 of the paper, which is tested above rather than assumed, but the '
          'switching test is exhaustive only on theta and K_4 (all embeddings, all vertices) and '
          'covers 2 rotations x 256 signatures x every vertex on the other five cells. '
          'NOT RE-RUN: any bridgeless '
          'cubic graph outside these seven; the infinitude of the connected cubic arc-transitive '
          'family, which the paper cites rather than verifies; Bender & Richmond 1990, '
          'whose bound is transcribed from the source '
          'authors\' own restatement and not recomputed from the 1990 paper; and the '
          'Gross & Furst 1987 distribution 2:14, which is transcribed and then reproduced here, '
          'not proved.')
    print('')
    if _F[0]:
        print('VERDICT: %d CHECK(S) FAILED of %d' % (_F[0], _N[0]))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _N[0])
    return 0


if __name__ == '__main__':
    sys.exit(main())
