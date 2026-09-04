#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Verification program for "A 24-Permutation Hamiltonian Cycle of P(4) That Does
Not Extend to a u-Cycle".

Python 3.9+, standard library only.  No third-party package, no external data
file.  All arithmetic is exact integer arithmetic; no floating-point value is
used for any decision, and none is printed.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs -- these are NOT checks by themselves)
--------------------------------------------------------------------------
  W1, W2        the two 24-permutation cyclic sequences printed in Section 3
                of the paper, one permutation of {1,2,3,4} per cyclic
                position 0..23, read verbatim from the paper's two tables.
  CHAIN1, CHAIN2
                the nine positions of the printed contradiction chains
                u_2<u_5<...<u_23<u_2 and u_0<u_3<...<u_22<u_0.
  WINDOWS1, WINDOWS2
                for each of the nine links, the window index the paper's
                table names as forcing it.
  STEPS1, STEPS2, KAHN1, KAHN2, CHAIN_ARCS, LIT_ARCS, MINCYC
                the numbers the paper displays about the two forced digraphs.
  KQ3_CYCLE, KQ3_WORD, KQ3_REDUCTION
                the Hamiltonian cycle of P(3), the u-cycle word U_3=142342 and
                the Hasse diagram of Figure "poset-P3", all three PUBLISHED by
                Kitaev and Qiu and quoted in the paper's control section.  They
                are the external pin on the reading convention: they are inputs
                here, and what is checked is that this program's own
                implementation of the convention reproduces them.
  TOTAL_P4=13063680, TW_H4=280, TOTAL_P3=8, TW_H3=2
                the arborescence and Eulerian-circuit counts this program
                recomputes from scratch below.
  A016031       the four PUBLISHED de Bruijn-sequence counts 2, 16, 2048,
                67108864 (OEIS A016031), used as external controls on this
                program's own BEST/Matrix-Tree code path.

Every one of these is compared against a value this program recomputes from
scratch, except the items explicitly named as external published data
(KQ3_*, A016031), which are reproduced BY the recomputation rather
than checked against it -- that is what makes them controls.

--------------------------------------------------------------------------
WHAT IS **NOT** DONE HERE, AND IS NOT CLAIMED BY THIS PROGRAM
--------------------------------------------------------------------------
  * The exhaustive census of Ham(P(4)) is NOT re-run.  Counting the BAD
    Hamiltonian cycles at n=4 would require enumerating 13,063,680 cycles and
    is out of scope for a referee-runnable script; no such count is an input
    to, or an output of, this program.
    ONLY the TOTAL 13,063,680 is re-derived here, and it is re-derived
    algebraically (BEST theorem + an integer Matrix-Tree cofactor), never by
    enumeration.
  * Nothing about n>=5 is computed or claimed.
  * The lower bound ceil(N/(n-1))=8 on the length of a forced contradiction is
    checked here as an arithmetic statement together with the structural fact
    that makes it valid (no forced arc spans cyclic distance > n-1).  Whether 8
    is attained by any bad cycle at n=4 is a census statement and is NOT
    re-derived here; what is re-derived is that the shortest directed cycle in
    each of the TWO EXHIBITED forced digraphs has length exactly 9.
"""

import collections
import itertools
import sys

# ---------------------------------------------------------------------------
# INPUTS, verbatim from the paper
# ---------------------------------------------------------------------------
W1 = ("1234 1243 1432 4312 4132 1324 3142 1423 4123 1342 3421 4213 "
      "2134 2341 2413 4231 2314 2143 2431 4321 3214 3241 3412 3124")
W2 = ("3214 3241 2314 3142 1324 2143 1423 4132 1432 4312 4123 1342 "
      "3421 4213 2134 2341 2413 4231 3412 3124 1234 1243 2431 4321")

CHAIN1 = [2, 5, 7, 9, 12, 14, 17, 20, 23]
CHAIN2 = [0, 3, 5, 8, 11, 14, 16, 19, 22]

# the window the paper's table names as forcing each link, link t joining
# CHAIN[t] -> CHAIN[t+1 mod 9]
WINDOWS1 = [2, 5, 6, 9, 12, 14, 17, 20, 23]
WINDOWS2 = [0, 3, 5, 8, 11, 14, 16, 19, 22]

STEPS1 = [3, 2, 2, 3, 2, 3, 3, 3, 3]
STEPS2 = [3, 2, 3, 3, 3, 2, 3, 3, 2]
KAHN1, KAHN2 = 7, 8
CHAIN_ARCS = 46            # both witnesses
LIT_ARCS = 72              # both witnesses
MINCYC = 9                 # shortest directed cycle in each forced digraph

# PUBLISHED control data (Kitaev--Qiu), quoted in the paper
KQ3_CYCLE = ["132", "312", "123", "231", "321", "213"]
KQ3_WORD = [1, 4, 2, 3, 4, 2]                      # U_3 = 142342
# their Figure "poset-P3", with a..f = positions 0..5
KQ3_REDUCTION = {(0, 2), (0, 5), (2, 3), (5, 3), (3, 1), (3, 4)}

TW_H4, TOTAL_P4 = 280, 13063680
TW_H3, TOTAL_P3 = 2, 8
TOTAL_P2 = 1

# OEIS A016031, the number of binary de Bruijn sequences of order m+1, for
# m = 2,3,4,5.  PUBLISHED numbers, used as controls on the code path below.
A016031 = {2: 2, 3: 16, 4: 2048, 5: 67108864}

# ---------------------------------------------------------------------------
# harness
# ---------------------------------------------------------------------------
_N = 0
_BAD = 0


def ck(name, cond, detail=''):
    global _N, _BAD
    _N += 1
    if cond:
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _BAD += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(text):
    print('NOTE %s' % text)


def head(text):
    print('\n=== %s' % text)


# ---------------------------------------------------------------------------
# the reading convention of Chung--Diaconis--Graham, implemented once
# ---------------------------------------------------------------------------
def st(w):
    """The standardisation (order-isomorphism class) of a word of distinct
    letters, as a permutation of 1..len(w)."""
    s = sorted(w)
    return tuple(s.index(x) + 1 for x in w)


def parse(text):
    return [tuple(int(c) for c in t) for t in text.split()]


def arc_ok(p, q):
    """The arc rule of P(n), stated as the source's iff: p -> q exactly when
    st(p_2..p_n) == st(q_1..q_{n-1})."""
    return st(p[1:]) == st(q[:-1])


def arc_failures(P):
    n = len(P)
    return [i for i in range(n) if not arc_ok(P[i], P[(i + 1) % n])]


def loops_used(P):
    n = len(P)
    return [i for i in range(n) if P[i] == P[(i + 1) % n]]


def chain_arcs(P, w):
    """The CONSECUTIVE-RANK reading.  Window i occupies cyclic positions
    i,...,i+w-1 and must be order-isomorphic to P[i]; the letters of ranks r
    and r+1 inside it give one forced strict inequality each."""
    n = len(P)
    arcs = set()
    for i in range(n):
        pos = [(i + k) % n for k in range(w)]
        order = sorted(range(w), key=lambda k: P[i][k])
        for a, b in zip(order, order[1:]):
            arcs.add((pos[a], pos[b]))
    return arcs


def chain_arcs_nowrap(P, w):
    """The same reading with the w-1 WRAP-AROUND windows deliberately omitted.
    Present only to show that this single omission is what would manufacture a
    false "acyclic"."""
    n = len(P)
    arcs = set()
    for i in range(n - w + 1):
        pos = [i + k for k in range(w)]
        order = sorted(range(w), key=lambda k: P[i][k])
        for a, b in zip(order, order[1:]):
            arcs.add((pos[a], pos[b]))
    return arcs


def literal_arcs(P, w):
    """The source's LITERAL all-ordered-pairs reading: every pair of positions
    inside a common window is compared."""
    n = len(P)
    arcs = set()
    for i in range(n):
        pos = [(i + k) % n for k in range(w)]
        for ka in range(w):
            for kb in range(w):
                if P[i][ka] < P[i][kb]:
                    arcs.add((pos[ka], pos[kb]))
    return arcs


def forcing_windows(P, w, a, b):
    """Every window that forces u_a < u_b by consecutive ranks."""
    n = len(P)
    out = []
    for i in range(n):
        pos = [(i + k) % n for k in range(w)]
        if a in pos and b in pos:
            ka, kb = pos.index(a), pos.index(b)
            if P[i][kb] == P[i][ka] + 1:
                out.append(i)
    return out


def kahn_ordered(nodes, arcs):
    """How many of the nodes a topological sort can place.  Fewer than all of
    them iff the digraph has a directed cycle."""
    indeg = collections.Counter()
    adj = collections.defaultdict(set)
    for a, b in arcs:
        if b not in adj[a]:
            adj[a].add(b)
            indeg[b] += 1
    q = [v for v in nodes if indeg[v] == 0]
    seen = 0
    while q:
        v = q.pop()
        seen += 1
        for x in adj[v]:
            indeg[x] -= 1
            if indeg[x] == 0:
                q.append(x)
    return seen


def shortest_directed_cycle(nodes, arcs):
    """Length of the shortest directed cycle, or None if acyclic.  BFS from
    every node; exact, no heuristics."""
    adj = collections.defaultdict(list)
    for a, b in arcs:
        adj[a].append(b)
    best = None
    for s in nodes:
        dist = {s: 0}
        q = collections.deque([s])
        while q:
            v = q.popleft()
            for x in adj[v]:
                if x == s:
                    d = dist[v] + 1
                    if best is None or d < best:
                        best = d
                elif x not in dist:
                    dist[x] = dist[v] + 1
                    q.append(x)
    return best


def is_linear_extension(word, arcs):
    return all(word[a] < word[b] for a, b in arcs)


def transitive_closure(nodes, arcs):
    reach = {v: set() for v in nodes}
    for a, b in arcs:
        reach[a].add(b)
    changed = True
    while changed:
        changed = False
        for v in nodes:
            add = set()
            for x in reach[v]:
                add |= reach[x]
            if not add <= reach[v]:
                reach[v] |= add
                changed = True
    return reach


def transitive_reduction(nodes, arcs):
    reach = transitive_closure(nodes, arcs)
    out = set()
    for a, b in arcs:
        if not any(b in reach[c] for c in reach[a] if c != b):
            out.add((a, b))
    return out


# ---------------------------------------------------------------------------
# exact linear algebra: Bareiss fraction-free determinant over the integers
# ---------------------------------------------------------------------------
def det_int(mat):
    """Exact determinant of an integer matrix.  Bareiss, so every intermediate
    value is an integer and no division is ever inexact."""
    m = [list(map(int, r)) for r in mat]
    n = len(m)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if m[k][k] == 0:
            piv = None
            for i in range(k + 1, n):
                if m[i][k] != 0:
                    piv = i
                    break
            if piv is None:
                return 0
            m[k], m[piv] = m[piv], m[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = m[i][j] * m[k][k] - m[i][k] * m[k][j]
                assert num % prev == 0, 'Bareiss division was not exact'
                m[i][j] = num // prev
        prev = m[k][k]
    return sign * m[n - 1][n - 1]


def arborescences(nverts, edges, root):
    """Number of spanning arborescences ORIENTED TOWARDS `root` counted by the
    digraph Matrix-Tree theorem: the (root,root) cofactor of L = D_out - A.
    Self-loops cancel between D_out and A and are dropped."""
    lap = [[0] * nverts for _ in range(nverts)]
    for a, b in edges:
        if a == b:
            continue
        lap[a][a] += 1
        lap[a][b] -= 1
    idx = [i for i in range(nverts) if i != root]
    return det_int([[lap[i][j] for j in idx] for i in idx])


def brute_arborescences(nverts, edges, root):
    """The same number with no linear algebra at all: choose one outgoing
    simple edge per non-root vertex and test that following it always reaches
    the root.  Exponential, used only on the 6-vertex cluster digraph."""
    # ⛔ WITH MULTIPLICITY: the cluster digraph is a MULTIdigraph, and the
    # Matrix-Tree cofactor counts arborescences with multiplicity, so an
    # independent brute force must too.  One outgoing edge OCCURRENCE is chosen
    # per non-root vertex, parallel copies counted separately.
    outs = {v: [(a, b) for (a, b) in edges if a == v and a != b]
            for v in range(nverts)}
    others = [v for v in range(nverts) if v != root]
    total = 0
    for choice in itertools.product(*[outs[v] for v in others]):
        nxt = {a: b for a, b in choice}
        good = True
        for v in others:
            seen = set()
            u = v
            while u != root:
                if u in seen or u not in nxt:
                    good = False
                    break
                seen.add(u)
                u = nxt[u]
            if not good:
                break
        if good:
            total += 1
    return total


def cluster_digraph(n):
    """H_n: vertices S_{n-1}, and each p in S_n is the edge
    st(p_1..p_{n-1}) -> st(p_2..p_n).  P(n) is the line digraph of H_n."""
    verts = sorted(itertools.permutations(range(1, n)))
    idx = {v: i for i, v in enumerate(verts)}
    edges = []
    for p in sorted(itertools.permutations(range(1, n + 1))):
        edges.append((idx[st(p[:-1])], idx[st(p[1:])]))
    return len(verts), edges, verts


def factorial(k):
    r = 1
    for i in range(2, k + 1):
        r *= i
    return r


def best_count(nverts, edges, root=0):
    """BEST theorem: the number of Eulerian circuits of a connected Eulerian
    digraph is tw(root) * prod_v (outdeg(v) - 1)!."""
    outdeg = collections.Counter(a for a, b in edges)
    tw = arborescences(nverts, edges, root)
    prod = 1
    for v in range(nverts):
        prod *= factorial(outdeg[v] - 1)
    return tw, tw * prod


def de_bruijn_graph(m):
    """B(2,m): vertices are the 2^m binary m-tuples, edges the 2^(m+1)
    (m+1)-tuples."""
    nv = 1 << m
    edges = []
    for v in range(nv):
        for bit in (0, 1):
            w = ((v << 1) | bit) & (nv - 1)
            edges.append((v, w))
    return nv, edges


def ham_p3_exhaustive():
    """Every directed Hamiltonian cycle of P(3), rooted at 123 so each cycle is
    counted once, together with the acyclicity of its forced order."""
    verts = sorted(itertools.permutations((1, 2, 3)))
    root = (1, 2, 3)
    rest = [v for v in verts if v != root]
    cycles = []
    for perm in itertools.permutations(rest):
        P = [root] + list(perm)
        if not arc_failures(P):
            cycles.append(P)
    good = [P for P in cycles if kahn_ordered(range(len(P)), chain_arcs(P, 3)) == len(P)]
    return cycles, good


def ham_p2_exhaustive():
    verts = sorted(itertools.permutations((1, 2)))
    root = (1, 2)
    rest = [v for v in verts if v != root]
    cycles = []
    for perm in itertools.permutations(rest):
        P = [root] + list(perm)
        if not arc_failures(P):
            cycles.append(P)
    good = [P for P in cycles if kahn_ordered(range(len(P)), chain_arcs(P, 2)) == len(P)]
    return cycles, good


# ===========================================================================
print('verification of "A 24-Permutation Hamiltonian Cycle of P(4) That Does '
      'Not Extend to a u-Cycle"')
print('python %s, exact integer arithmetic only'
      % sys.version.split()[0])

S4 = set(itertools.permutations((1, 2, 3, 4)))
P1, P2 = parse(W1), parse(W2)

# ---------------------------------------------------------------------------
head('Step 1: the two exhibited objects are Hamiltonian cycles of P(4)')
for nm, P in (('w1', P1), ('w2', P2)):
    ck('%s_has_24_entries' % nm, len(P) == 24, '%d' % len(P))
    ck('%s_entries_are_exactly_the_24_permutations_of_1234' % nm,
       len(set(P)) == 24 and set(P) == S4,
       '%d distinct, equals S_4: %s' % (len(set(P)), set(P) == S4))
    af = arc_failures(P)
    ck('%s_all_24_arcs_of_p4_hold_including_the_wrap_arc_23_to_0' % nm, af == [],
       'failures %s' % (af or 'none'))
    lp = loops_used(P)
    ck('%s_uses_neither_of_the_two_self_loops_of_p4' % nm, lp == [],
       'loop positions %s' % (lp or 'none'))

nloop = sum(1 for p in S4 if arc_ok(p, p))
ck('p4_has_exactly_two_self_loops_namely_1234_and_4321',
   nloop == 2 and arc_ok((1, 2, 3, 4), (1, 2, 3, 4)) and arc_ok((4, 3, 2, 1), (4, 3, 2, 1)),
   '%d self-loops' % nloop)

# ---------------------------------------------------------------------------
head('Step 2: the printed contradiction chains, link by link')
for nm, P, CH, WD, STP, KH in (('w1', P1, CHAIN1, WINDOWS1, STEPS1, KAHN1),
                               ('w2', P2, CHAIN2, WINDOWS2, STEPS2, KAHN2)):
    A = chain_arcs(P, 4)
    links = [(CH[t], CH[(t + 1) % 9]) for t in range(9)]
    for t, (a, b) in enumerate(links):
        w = WD[t]
        pos = [(w + k) % 24 for k in range(4)]
        ka, kb = pos.index(a), pos.index(b)
        okk = (a in pos and b in pos and P[w][kb] == P[w][ka] + 1)
        ck('%s_link_%d_of_9_u%d_lt_u%d_is_forced_by_window_%d' % (nm, t + 1, a, b, w),
           okk and (a, b) in A,
           'window %d is %s on positions %s, ranks %d and %d; all windows forcing '
           'it: %s' % (w, ''.join(map(str, P[w])), pos, P[w][ka], P[w][kb],
                       forcing_windows(P, 4, a, b)))
    ck('%s_the_nine_links_close_into_one_directed_cycle_of_strict_inequalities' % nm,
       all(l in A for l in links) and len(set(CH)) == 9,
       ' < '.join('u%d' % c for c in CH) + ' < u%d' % CH[0])
    steps = [(CH[(t + 1) % 9] - CH[t]) % 24 for t in range(9)]
    ck('%s_step_sizes_match_the_supplied_list_and_sum_to_24' % nm,
       steps == STP and sum(steps) == 24, '%s sum %d' % (steps, sum(steps)))
    k = kahn_ordered(range(24), A)
    ck('%s_forced_digraph_is_cyclic_topological_sort_places_only_%d_of_24' % (nm, KH),
       k == KH and k < 24, 'placed %d of 24' % k)
    ck('%s_consecutive_rank_reading_has_%d_arcs' % (nm, CHAIN_ARCS),
       len(A) == CHAIN_ARCS, '%d' % len(A))
    L = literal_arcs(P, 4)
    ck('%s_literal_all_ordered_pairs_reading_has_%d_arcs' % (nm, LIT_ARCS),
       len(L) == LIT_ARCS, '%d' % len(L))
    viol = [(a, b) for (a, b) in L if (b, a) in L]
    ck('%s_literal_reading_has_zero_antisymmetry_violations_so_there_is_no_2_cycle' % nm,
       viol == [], 'violations %s' % (viol or 'none'))
    ck('%s_literal_reading_contains_all_nine_links' % nm,
       all(l in L for l in links), '9 of 9 present')
    ck('%s_literal_reading_is_cyclic_too' % nm,
       kahn_ordered(range(24), L) < 24,
       'topological sort places %d of 24' % kahn_ordered(range(24), L))
    far = [(a, b) for (a, b) in L if min((a - b) % 24, (b - a) % 24) > 3]
    ck('%s_no_forced_arc_of_this_cycle_spans_cyclic_distance_more_than_3' % nm,
       far == [], 'offenders %s' % (far or 'none'))
    mc = shortest_directed_cycle(range(24), A)
    ck('%s_shortest_directed_cycle_in_the_forced_digraph_has_length_exactly_%d' % (nm, MINCYC),
       mc == MINCYC, 'length %s' % mc)
    mcl = shortest_directed_cycle(range(24), L)
    ck('%s_shortest_directed_cycle_is_9_under_the_literal_reading_too' % nm,
       mcl == MINCYC, 'length %s' % mcl)

ck('the_winding_lower_bound_is_ceil_24_over_3_equals_8',
   -(-24 // 3) == 8, 'ceil(N/(n-1)) = ceil(24/3) = 8')
note('so the exhibited chains of length 9 are one step above the structural '
     'lower bound of 8; whether 8 is attained by any Hamiltonian cycle of '
     'P(4) with cyclic forced order is NOT re-derived here.')
ck('the_two_witnesses_are_different_cyclic_objects_not_a_rotation_of_each_other',
   all(P2 != [P1[(i + r) % 24] for i in range(24)] for r in range(24)),
   'checked all 24 rotations')

# ---------------------------------------------------------------------------
head('Step 3: controls -- the convention is pinned to PUBLISHED output at n=3')
K3 = [tuple(int(c) for c in s) for s in KQ3_CYCLE]
ck('published_p3_cycle_has_all_six_arcs', arc_failures(K3) == [],
   'failures %s' % (arc_failures(K3) or 'none'))
ck('published_p3_cycle_is_hamiltonian_on_s3',
   set(K3) == set(itertools.permutations((1, 2, 3))) and len(set(K3)) == 6,
   '6 distinct vertices, equals S_3')
wins = [st(tuple(KQ3_WORD[(i + k) % 6] for k in range(3))) for i in range(6)]
ck('the_published_word_142342_read_as_a_circular_word_gives_exactly_that_cycle',
   wins == K3,
   'windows %s -- including the two wrap windows (4,2,1) and (2,1,4)'
   % [''.join(map(str, w)) for w in wins])
K3A = chain_arcs(K3, 3)
K3L = literal_arcs(K3, 3)
ck('published_p3_cycle_forced_order_is_acyclic_under_both_readings',
   kahn_ordered(range(6), K3A) == 6 and kahn_ordered(range(6), K3L) == 6,
   'topological sort places 6 of 6 under both readings')
ck('the_published_word_142342_satisfies_every_forced_strict_inequality',
   is_linear_extension(KQ3_WORD, K3L) and is_linear_extension(KQ3_WORD, K3A),
   'u = %s is a linear extension of the forced order' % KQ3_WORD)
red = transitive_reduction(range(6), K3L)
ck('the_transitive_reduction_reproduces_the_published_figure_poset_p3',
   red == KQ3_REDUCTION,
   'a..f = positions 0..5; got %s' % sorted(red))
ck('the_published_word_uses_only_4_letters_on_6_positions',
   len(set(KQ3_WORD)) == 4, 'letters %s' % sorted(set(KQ3_WORD)))

c3, g3 = ham_p3_exhaustive()
ck('exhaustively_p3_has_exactly_8_directed_hamiltonian_cycles',
   len(c3) == TOTAL_P3, '%d' % len(c3))
ck('and_all_8_of_them_have_acyclic_forced_order',
   len(g3) == len(c3) == TOTAL_P3, '%d of %d good' % (len(g3), len(c3)))
c2, g2 = ham_p2_exhaustive()
ck('exhaustively_p2_has_exactly_1_directed_hamiltonian_cycle_and_it_is_good',
   len(c2) == TOTAL_P2 and len(g2) == TOTAL_P2, '%d of %d good' % (len(g2), len(c2)))

anti = [(1, 2, 3)] * 3
ck('anti_control_three_identical_123_windows_on_three_letters_is_CYCLIC',
   kahn_ordered(range(3), chain_arcs(anti, 3)) < 3,
   'u0<u1<u2, u1<u2<u0, u2<u0<u1 -- the tester must and does say CYCLIC here')
nw1 = kahn_ordered(range(24), chain_arcs_nowrap(P1, 4))
nw2 = kahn_ordered(range(24), chain_arcs_nowrap(P2, 4))
ck('omitting_the_three_wrap_windows_makes_both_exhibited_cycles_look_acyclic',
   nw1 == 24 and nw2 == 24,
   'dropping windows 21,22,23 leaves a topological sort placing %d of 24 (W1) '
   'and %d of 24 (W2) -- this single omission is the one bug the wrap-around '
   'convention exists to prevent, and the published word 142342 above shows '
   'the convention is the source\'s own' % (nw1, nw2))

# ---------------------------------------------------------------------------
head('Step 4: the census TOTAL, algebraically (BEST + integer Matrix-Tree)')
nv4, ed4, vs4 = cluster_digraph(4)
outd = collections.Counter(a for a, b in ed4)
ind = collections.Counter(b for a, b in ed4)
ck('cluster_digraph_h4_has_6_vertices_24_edges_and_every_degree_4',
   nv4 == 6 and len(ed4) == 24
   and set(outd.values()) == {4} and set(ind.values()) == {4},
   '|V|=%d |E|=%d outdeg=%s indeg=%s' % (nv4, len(ed4), set(outd.values()),
                                         set(ind.values())))
# The bijection Ham(P(4)) <-> Euler(H_4): a permutation p IS the edge
# tail_of[p] -> head_of[p] of the cluster digraph, and p -> q is an arc of P(4)
# exactly when the head cluster of p equals the tail cluster of q.  Checked on
# every one of the 24*24 ordered pairs, so it is an equality of predicates and
# not an inclusion.
head_of = {p: st(p[1:]) for p in itertools.permutations((1, 2, 3, 4))}
tail_of = {p: st(p[:-1]) for p in itertools.permutations((1, 2, 3, 4))}
mism = sum(1 for p in itertools.permutations((1, 2, 3, 4))
           for q in itertools.permutations((1, 2, 3, 4))
           if arc_ok(p, q) != (head_of[p] == tail_of[q]))
ck('ham_p4_equals_euler_h4_the_two_predicates_agree_on_all_576_ordered_pairs',
   mism == 0, '%d mismatches over 24*24 pairs' % mism)

tws = [arborescences(nv4, ed4, r) for r in range(nv4)]
ck('tw_h4_is_280_from_every_one_of_the_six_roots',
   set(tws) == {TW_H4}, 'roots give %s' % tws)
bf = brute_arborescences(nv4, ed4, 0)
ck('tw_h4_is_280_again_by_brute_force_enumeration_with_no_linear_algebra',
   bf == TW_H4, '%d' % bf)
tw4, tot4 = best_count(nv4, ed4)
ck('best_theorem_gives_280_times_6_to_the_6_equals_13063680',
   tw4 == TW_H4 and tot4 == TOTAL_P4 and TW_H4 * 6 ** 6 == TOTAL_P4,
   '%d * 6^6 = %d' % (tw4, tot4))

nv3, ed3, vs3 = cluster_digraph(3)
tw3, tot3 = best_count(nv3, ed3)
ck('the_same_code_path_predicts_8_at_n_equals_3_matching_the_exhaustive_count',
   tw3 == TW_H3 and tot3 == TOTAL_P3 == len(c3),
   'tw(H_3)=%d, %d*(2!)^2=%d, exhaustive=%d' % (tw3, tw3, tot3, len(c3)))
note('that control is the whole reason the n=4 total may be trusted: the same '
     'BEST/Matrix-Tree code, unchanged, reproduces a number this program also '
     'obtains by exhaustive enumeration.')

for m in sorted(A016031):
    nvb, edb = de_bruijn_graph(m)
    twb, ecb = best_count(nvb, edb)
    ck('best_code_path_gives_de_bruijn_eulerian_circuit_count_of_order_%d_equals_%d'
       % (m + 1, A016031[m]), ecb == A016031[m],
       'B(2,%d) on %d vertices: tw=%d, ec=%d' % (m, nvb, twb, ecb))

# ---------------------------------------------------------------------------
head('closing')
note('SCOPE. What this program re-derives: both exhibited 24-permutation '
     'objects are Hamiltonian cycles of P(4); each of the nine printed links '
     'of each contradiction chain is forced by the window the paper names; '
     'each forced digraph is cyclic, with shortest directed cycle exactly 9 '
     'under both readings of the source; the verdict is independent of which '
     'reading is used; the n=3 convention agrees with Kitaev-Qiu\'s own '
     'published cycle, word U_3=142342 and Hasse diagram; P(3) and P(2) are '
     'exhaustively good; and the TOTAL |Ham(P(4))| = 280*6^6 = 13,063,680 by '
     'the BEST theorem with an exact integer Matrix-Tree cofactor, controlled '
     'against the exhaustive count 8 at n=3 and against four de '
     'Bruijn Eulerian-circuit counts computed by the same code path.')
note('SCOPE, NOT RE-RUN HERE. (a) No census of the bad Hamiltonian cycles '
     'at n=4 is performed: counting them, or any refinement such as a '
     'shortest-chain histogram or a count of symmetry orbits, would require '
     'enumerating all 13,063,680 cycles, and no such count is a claim of this '
     'program. (b) Nothing about '
     'n>=5. (c) The verbatim source quotations, which are byte-level facts '
     'about external e-print files and cannot be checked by a self-contained '
     'script. (d) Whether every n>=4 has a bad Hamiltonian cycle.')
print('VERDICT: ALL %d CHECKS PASS' % _N if _BAD == 0
      else 'VERDICT: %d of %d CHECKS FAILED' % (_BAD, _N))
sys.exit(1 if _BAD else 0)
