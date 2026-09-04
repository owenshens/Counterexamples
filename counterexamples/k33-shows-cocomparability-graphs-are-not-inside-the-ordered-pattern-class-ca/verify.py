#!/usr/bin/env python3
"""
Verification program for the paper

    K_{3,3} is a cocomparability graph with no P_a-avoiding vertex order

The objects of the note itself -- K_{3,3}, its complement, the order 031425 and the
permutation 456123 -- are transcribed from the note.  Other graphs, orders and permutations
used below are constructed by this program and are not displayed in the note.  Nothing is
read from an external file and nothing is downloaded.

Python 3.9+, standard library only (sys, itertools).  All arithmetic is on integers; there
is no floating point anywhere, so no decision depends on a rounding mode.

The program prints one `PASS <name> [detail]` line per check and closes with

    VERDICT: ALL <n> CHECKS PASS

exiting 0 if and only if every check passed.  It exits 1 on the first-and-any failure count.
"""

import sys
import itertools

PASSES = 0
FAILS = 0


def check(name, ok, detail=''):
    global PASSES, FAILS
    if ok:
        PASSES += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        FAILS += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(s):
    print('NOTE %s' % s)


# ---------------------------------------------------------------------------
# graphs
# ---------------------------------------------------------------------------
def graph(n, edge_list):
    """(n, adj-bitmask-per-vertex, frozenset of normalised edges)."""
    adj = [0] * n
    E = set()
    for (u, v) in edge_list:
        if u == v:
            raise ValueError('loop')
        adj[u] |= 1 << v
        adj[v] |= 1 << u
        E.add((min(u, v), max(u, v)))
    return (n, adj, frozenset(E))


def complement(G):
    n, adj, E = G
    return graph(n, [(u, v) for u in range(n) for v in range(u + 1, n) if (u, v) not in E])


def degrees(G):
    n, adj, E = G
    return [bin(adj[v]).count('1') for v in range(n)]


def triangles(G):
    n, adj, E = G
    t = 0
    for a in range(n):
        for b in range(a + 1, n):
            if not (adj[a] >> b) & 1:
                continue
            for c in range(b + 1, n):
                if (adj[a] >> c) & 1 and (adj[b] >> c) & 1:
                    t += 1
    return t


def edge_str(G):
    return ' '.join('%d%d' % (u, v) for (u, v) in sorted(G[2]))


def components(G):
    n, adj, E = G
    seen = 0
    comps = []
    for s in range(n):
        if (seen >> s) & 1:
            continue
        stack = [s]
        cur = 0
        while stack:
            v = stack.pop()
            if (cur >> v) & 1:
                continue
            cur |= 1 << v
            rest = adj[v] & ~cur
            while rest:
                b = rest & -rest
                stack.append(b.bit_length() - 1)
                rest ^= b
        seen |= cur
        comps.append(cur)
    return comps


def is_bipartite(G):
    n, adj, E = G
    colour = [None] * n
    for s in range(n):
        if colour[s] is not None:
            continue
        colour[s] = 0
        stack = [s]
        while stack:
            v = stack.pop()
            rest = adj[v]
            while rest:
                b = rest & -rest
                w = b.bit_length() - 1
                rest ^= b
                if colour[w] is None:
                    colour[w] = 1 - colour[v]
                    stack.append(w)
                elif colour[w] == colour[v]:
                    return False
    return True


def induced(G, keep):
    """keep: sorted list of vertices. Returns the induced subgraph relabelled 0..k-1."""
    n, adj, E = G
    ix = {v: i for i, v in enumerate(keep)}
    el = [(ix[u], ix[v]) for (u, v) in E if u in ix and v in ix]
    return graph(len(keep), el)


def complete_multipartite(parts):
    """Vertices are laid out consecutively part by part."""
    n = sum(parts)
    lab = []
    for t, sz in enumerate(parts):
        lab += [t] * sz
    el = [(u, v) for u in range(n) for v in range(u + 1, n) if lab[u] != lab[v]]
    verts = []
    off = 0
    for sz in parts:
        verts.append(list(range(off, off + sz)))
        off += sz
    return graph(n, el), verts


def cycle(n):
    return graph(n, [(i, (i + 1) % n) for i in range(n)])


def complete(n):
    return graph(n, [(u, v) for u in range(n) for v in range(u + 1, n)])


def partitions(n, maxpart=None):
    """Partitions of n as non-increasing tuples."""
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


# ---------------------------------------------------------------------------
# the ordered patterns P_S of Feuilloley-Habib
#
#   P_S is the pattern on positions 1<2<3<4 with MANDATORY EDGES (1,3) and (2,4), the pairs
#   in S forced to be NON-EDGES, and every other pair FREE.  The four optional pairs are
#   named a=(1,2), b=(2,3), c=(3,4), d=(1,4).
# ---------------------------------------------------------------------------
PAIR_OF = {'a': (0, 1), 'b': (1, 2), 'c': (2, 3), 'd': (0, 3)}


def occurrence_naive(G, order, S):
    """Prune-free scan of all C(n,4) position quadruples. Returns a quadruple or None."""
    n, adj, E = G
    for quad in itertools.combinations(range(len(order)), 4):
        w = [order[t] for t in quad]
        if not (adj[w[0]] >> w[2]) & 1:
            continue
        if not (adj[w[1]] >> w[3]) & 1:
            continue
        bad = False
        for s in S:
            i, j = PAIR_OF[s]
            if (adj[w[i]] >> w[j]) & 1:
                bad = True
                break
        if not bad:
            return quad
    return None


def occurrence_pa_fast(G, order):
    """Incremental decider for P_a only, independent of occurrence_naive.

    An order contains P_a iff some NON-EDGE {v_p, v_q} with p<q has a neighbour of v_p at
    some position k>q and a neighbour of v_q at some position l>k.  Working in position
    space, that is  min(A) < max(B)  with A the later neighbours of v_p and B those of v_q.
    """
    n, adj, E = G
    pos = [0] * n
    for i, v in enumerate(order):
        pos[v] = i
    padj = [0] * n
    for i, v in enumerate(order):
        m = 0
        rest = adj[v]
        while rest:
            b = rest & -rest
            m |= 1 << pos[b.bit_length() - 1]
            rest ^= b
        padj[i] = m
    full = (1 << n) - 1
    for p in range(n):
        for q in range(p + 1, n):
            if (padj[p] >> q) & 1:
                continue
            after = full ^ ((1 << (q + 1)) - 1)
            A = padj[p] & after
            if not A:
                continue
            B = padj[q] & after
            if not B:
                continue
            if ((A & -A).bit_length() - 1) < (B.bit_length() - 1):
                return (p, q)
    return None


def count_avoiding(G, S=('a',), fast=False):
    """Exhaustive count of vertex orders avoiding P_S, plus one witnessing order."""
    n = G[0]
    hit = occurrence_pa_fast if fast else (lambda g, o: occurrence_naive(g, o, S))
    cnt = 0
    wit = None
    for order in itertools.permutations(range(n)):
        if hit(G, order) is None:
            cnt += 1
            if wit is None:
                wit = order
    return cnt, wit


def first_avoiding(G, S=('a',), fast=False):
    n = G[0]
    hit = occurrence_pa_fast if fast else (lambda g, o: occurrence_naive(g, o, S))
    for order in itertools.permutations(range(n)):
        if hit(G, order) is None:
            return order
    return None


# ---------------------------------------------------------------------------
# cocomparability, by two structurally unrelated deciders
# ---------------------------------------------------------------------------
def umbrella_violations(G, order):
    """The source's pattern on 1<2<3: edge (1,3), non-edges (1,2) and (2,3)."""
    n, adj, E = G
    bad = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                a, b, c = order[i], order[j], order[k]
                if (adj[a] >> c) & 1 and not (adj[a] >> b) & 1 and not (adj[b] >> c) & 1:
                    bad.append((i, j, k))
    return bad


def umbrella_free_order(G):
    for order in itertools.permutations(range(G[0])):
        if not umbrella_violations(G, order):
            return order
    return None


def orientation_is_transitive(n, arcs):
    """(ok, reason). arcs is a list of ordered pairs."""
    out = [0] * n
    for (u, v) in arcs:
        out[u] |= 1 << v
    for (u, v) in arcs:
        if (out[v] >> u) & 1:
            return False, 'both %d->%d and %d->%d' % (u, v, v, u)
    for (u, v) in arcs:
        if out[v] & ~out[u]:
            return False, 'transitivity fails at %d->%d' % (u, v)
    return True, ''


def transitive_orientation(G):
    """Exhaustive search over all 2^m orientations. Returns arcs or None."""
    n, adj, E = G
    el = sorted(E)
    m = len(el)
    if m > 16:
        return 'TOO-LARGE'
    for bits in range(1 << m):
        out = [0] * n
        for k in range(m):
            u, v = el[k]
            if (bits >> k) & 1:
                out[u] |= 1 << v
            else:
                out[v] |= 1 << u
        good = True
        for k in range(m):
            u, v = el[k]
            a, b = (u, v) if (bits >> k) & 1 else (v, u)
            if out[b] & ~out[a]:
                good = False
                break
        if good:
            return [((u, v) if (bits >> k) & 1 else (v, u)) for k, (u, v) in enumerate(el)]
    return None


def is_cocomparability_by_search(G):
    r = transitive_orientation(complement(G))
    if r == 'TOO-LARGE':
        return None
    return r is not None


# ---------------------------------------------------------------------------
# permutation graphs
# ---------------------------------------------------------------------------
def inversion_graph(pi):
    n = len(pi)
    return graph(n, [(i, j) for i in range(n) for j in range(i + 1, n) if pi[i] > pi[j]])


def block_reversal(parts):
    """Reverse the order of the blocks, keeping the order inside each block."""
    n = sum(parts)
    pi = [0] * n
    p = 0
    for t, m in enumerate(parts):
        start = sum(parts[t + 1:])
        for k in range(m):
            pi[p] = start + k
            p += 1
    return pi


# ---------------------------------------------------------------------------
# Theorem 2 (the classification) and its explicit orders
# ---------------------------------------------------------------------------
def criterion(parts):
    """n_2 <= 2 and n_3 <= 1, with parts sorted non-increasingly and n_3 = 0 when p < 3."""
    s = sorted(parts, reverse=True)
    n2 = s[1] if len(s) > 1 else 0
    n3 = s[2] if len(s) > 2 else 0
    return n2 <= 2 and n3 <= 1


def sufficiency_order(parts):
    """An order built by this program from the part sizes; not an order printed anywhere."""
    G, verts = complete_multipartite(parts)
    order_of = sorted(range(len(parts)), key=lambda t: -parts[t])
    t1 = order_of[0]
    s = sorted(parts, reverse=True)
    n2 = s[1] if len(s) > 1 else 0
    if n2 <= 1:
        out = []
        for t in range(len(parts)):
            if t != t1:
                out += verts[t]
        return tuple(out + verts[t1])
    t2 = order_of[1]
    W = verts[t2]
    rest = []
    for t in range(len(parts)):
        if t not in (t1, t2):
            rest += verts[t]
    return tuple([W[0]] + rest + verts[t1] + [W[1]])


# ---------------------------------------------------------------------------
# the exhaustive census on at most five vertices
# ---------------------------------------------------------------------------
def every_labeled_graph_is_in_Ca(n):
    """(covered, total, avoiding-in-the-identity-order).

    A graph is in C_a iff SOME relabelling of it avoids P_a in the identity order.  So it
    is enough to test the identity order on all 2^C(n,2) edge masks and then close the
    resulting set under the action of S_n on the pairs; `covered == total` says every
    labeled graph on n vertices is in C_a.
    """
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    npair = len(pairs)
    idx = {p: k for k, p in enumerate(pairs)}
    pairmaps = []
    for p in itertools.permutations(range(n)):
        pm = [0] * npair
        for k, (i, j) in enumerate(pairs):
            a, b = p[i], p[j]
            pm[k] = idx[(min(a, b), max(a, b))]
        pairmaps.append(pm)
    ident = tuple(range(n))
    avoid = []
    for mask in range(1 << npair):
        el = [pairs[k] for k in range(npair) if (mask >> k) & 1]
        if occurrence_pa_fast(graph(n, el), ident) is None:
            avoid.append(mask)
    covered = set()
    for mask in avoid:
        for pm in pairmaps:
            r = 0
            mm = mask
            while mm:
                b = mm & -mm
                r |= 1 << pm[b.bit_length() - 1]
                mm ^= b
            covered.add(r)
    return len(covered), 1 << npair, len(avoid)


# ===========================================================================
def main():
    print('verification of "K_{3,3} is a cocomparability graph with no P_a-avoiding vertex order"')
    print('python %s, exact integer arithmetic only' % sys.version.split()[0])
    print('')

    # ------------------------------------------------------------------ 1
    print('=== Step 1: K_{3,3}, the witness exhibited in the note')
    K33_EDGES = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]
    A, B = [0, 1, 2], [3, 4, 5]
    K33 = graph(6, K33_EDGES)
    note('E(K_{3,3}) as printed = %s' % edge_str(K33))
    note('parts A = %s, B = %s' % (A, B))
    check('printed_edge_list_is_exactly_the_cross_pairs_of_the_two_printed_parts',
          set(K33[2]) == set((min(a, b), max(a, b)) for a in A for b in B),
          '9 cross pairs, no within-part pair')
    check('order_is_6_and_size_is_9', (K33[0], len(K33[2])) == (6, 9),
          'n=%d, m=%d' % (K33[0], len(K33[2])))
    check('k33_is_3_regular', degrees(K33) == [3] * 6, 'degrees %s' % degrees(K33))
    check('k33_is_triangle_free', triangles(K33) == 0, '0 triangles')
    cK33 = complement(K33)
    comps = components(cK33)
    check('complement_of_k33_is_two_disjoint_triangles',
          sorted(edge_str(cK33).split()) == ['01', '02', '12', '34', '35', '45']
          and len(comps) == 2 and all(bin(c).count('1') == 3 for c in comps),
          'complement = %s, 2 components of order 3' % edge_str(cK33))
    print('')

    # ------------------------------------------------------------------ 2
    print('=== Step 2: no vertex order of K_{3,3} avoids P_a  (Theorem 1, second half)')
    cnt_fast, _ = count_avoiding(K33, fast=True)
    check('no_order_of_k33_avoids_pa_incremental_decider', cnt_fast == 0,
          '0 of 720 orders avoid P_a')
    cnt_naive, _ = count_avoiding(K33, S=('a',), fast=False)
    check('no_order_of_k33_avoids_pa_prunefree_quadruple_scan', cnt_naive == 0,
          '0 of 720 orders, all C(6,4)=15 position quadruples examined in each')
    disagree = 0
    for order in itertools.permutations(range(6)):
        if (occurrence_pa_fast(K33, order) is None) != (occurrence_naive(K33, order, ('a',)) is None):
            disagree += 1
    check('the_two_pa_deciders_agree_on_all_720_orders_of_k33', disagree == 0,
          '0 disagreements between two independent implementations')

    spot = (0, 3, 1, 4, 2, 5)          # a1 b1 a2 b2 a3 b3, the order spot-checked in the paper
    occ = occurrence_naive(K33, spot, ('a',))
    check('the_pa_occurrence_printed_for_the_alternating_order_is_an_occurrence',
          occ is not None and occ == (0, 2, 3, 5),
          'order %s, positions %s carry the vertices %s: non-edge 0-1, edges 0-4 and 1-5'
          % (''.join(str(v) for v in spot), occ,
             tuple(spot[t] for t in occ) if occ else None))

    # the pigeonhole step of the hand proof, on every order
    bad = 0
    for order in itertools.permutations(range(6)):
        posn = {v: i for i, v in enumerate(order)}
        p2 = sorted(posn[v] for v in A)[1]
        q2 = sorted(posn[v] for v in B)[1]
        afterA = sum(1 for v in B if posn[v] > p2)
        afterB = sum(1 for v in A if posn[v] > q2)
        if not (afterA >= 2 or afterB >= 2):
            bad += 1
    check('pigeonhole_step_of_the_hand_proof_holds_for_every_one_of_the_720_orders', bad == 0,
          'in every order, two vertices of one part follow the second vertex of the other')
    print('')

    # ------------------------------------------------------------------ 3
    print('=== Step 3: K_{3,3} IS a cocomparability graph  (Theorem 1, first half)')
    umb = umbrella_violations(K33, (0, 1, 2, 3, 4, 5))
    check('the_order_012345_of_k33_is_umbrella_free', umb == [],
          '0 violations among the C(6,3)=20 triples')
    ARCS = [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)]      # supplied by this program
    ok, why = orientation_is_transitive(6, ARCS)
    check('this_orientation_of_the_complement_of_k33_is_a_transitive_orientation',
          ok and set((min(u, v), max(u, v)) for (u, v) in ARCS) == set(cK33[2])
          and len(ARCS) == len(cK33[2]),
          '%d arcs = |E(complement)| = %d, 0 reversed pairs, 0 transitivity violations%s'
          % (len(ARCS), len(cK33[2]), '' if ok else ' -- ' + why))
    check('complement_of_k33_admits_a_transitive_orientation_by_independent_exhaustive_search',
          is_cocomparability_by_search(K33) is True,
          'search over all 2^6 = 64 orientations of the complement')
    PI = [3, 4, 5, 0, 1, 2]                                      # the paper's 4 5 6 1 2 3
    check('k33_is_the_inversion_graph_of_the_permutation_456123',
          inversion_graph(PI)[2] == K33[2] and sorted(PI) == list(range(6)),
          'pi = 4 5 6 1 2 3 (1-based); its 9 inversions are exactly E(K_{3,3})')
    print('')

    # ------------------------------------------------------------------ 4
    print('=== Step 4: Theorem 2 -- membership of complete multipartite graphs')
    for n in range(2, 8):
        mism = []
        for parts in partitions(n):
            if len(parts) < 2:
                continue
            G, _ = complete_multipartite(parts)
            brute = first_avoiding(G, fast=True) is not None
            if brute != criterion(parts):
                mism.append(parts)
        tot = sum(1 for p in partitions(n) if len(p) >= 2)
        check('theorem_2_criterion_agrees_with_exhaustive_search_on_n_%d_vertices' % n,
              not mism, '%d complete multipartite graphs, 0 mismatches' % tot)

    bad = []
    for n in range(2, 11):
        for parts in partitions(n):
            if len(parts) < 2 or not criterion(parts):
                continue
            G, _ = complete_multipartite(parts)
            if occurrence_pa_fast(G, sufficiency_order(parts)) is not None:
                bad.append(parts)
    check('the_orders_built_by_sufficiency_order_avoid_pa_for_every_admissible_partition',
          not bad, 'all admissible partitions with 2 <= n <= 10')

    check('theorem_2_places_k33_and_the_octahedron_outside_ca',
          (not criterion((3, 3))) and (not criterion((2, 2, 2))),
          'K_{3,3}: n_2=3>2;  K_{2,2,2}: n_3=2>1')
    check('theorem_2_places_k4_inside_ca_as_the_source_observation_states',
          criterion((1, 1, 1, 1)) and count_avoiding(complete(4), fast=True)[0] == 24,
          'all 24 orders of K_4 avoid P_a')
    print('')

    # ------------------------------------------------------------------ 5
    print('=== Step 5: Corollary 3 -- every complete multipartite graph is a permutation graph')
    bad = []
    for n in range(2, 9):
        for parts in partitions(n):
            if len(parts) < 2:
                continue
            G, _ = complete_multipartite(parts)
            pi = block_reversal(parts)
            if sorted(pi) != list(range(n)) or inversion_graph(pi)[2] != G[2]:
                bad.append(parts)
    check('the_block_reversal_permutation_realises_every_complete_multipartite_graph',
          not bad, 'inversion graph equals the graph for every partition with 2 <= n <= 8')
    bad = []
    for n in range(2, 9):
        for parts in partitions(n):
            if len(parts) < 2:
                continue
            G, verts = complete_multipartite(parts)
            arcs = [(a, b) for grp in verts for i, a in enumerate(grp) for b in grp[i + 1:]]
            ok, why = orientation_is_transitive(n, arcs)
            if not ok or set((min(u, v), max(u, v)) for (u, v) in arcs) != set(complement(G)[2]):
                bad.append((parts, why))
    check('orienting_each_part_by_index_transitively_orients_the_complement',
          not bad, 'a disjoint union of transitive tournaments, 2 <= n <= 8')
    print('')

    # ------------------------------------------------------------------ 6
    print('=== Step 6: the minimum order of a graph outside C_a is 6')
    for n in (4, 5):
        cov, tot, av = every_labeled_graph_is_in_Ca(n)
        check('every_labeled_graph_on_%d_vertices_is_in_ca' % n, cov == tot,
              '%d of %d edge masks covered (%d avoid P_a in the identity order)'
              % (cov, tot, av))
    note('P_a has four vertices, so every graph on at most 3 vertices is trivially in C_a')

    PRISM = graph(6, [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5)])
    OCTA, _ = complete_multipartite((2, 2, 2))
    note('E(prism), as constructed by this program = %s' % edge_str(PRISM))
    note('E(octahedron K_{2,2,2}), as constructed by this program = %s' % edge_str(OCTA))
    for nm, H in (('prism', PRISM), ('octahedron', OCTA)):
        c = count_avoiding(H, fast=True)[0]
        cc = is_cocomparability_by_search(H)
        check('the_%s_is_outside_ca_and_is_a_cocomparability_graph' % nm, c == 0 and cc is True,
              '0 of 720 orders avoid P_a; complement %s is bipartite=%s and transitively orientable'
              % (edge_str(complement(H)), is_bipartite(complement(H))))
    ms = (len(K33[2]), len(PRISM[2]), len(OCTA[2]))
    ts = (triangles(K33), triangles(PRISM), triangles(OCTA))
    check('the_three_order_6_graphs_used_here_are_pairwise_non_isomorphic',
          ms == (9, 9, 12) and ts == (0, 2, 8),
          'sizes %s, triangle counts %s' % (ms, ts))
    print('')

    # ------------------------------------------------------------------ 7
    print('=== Step 7: controls, both polarities, including published ones')
    check('complete_graphs_are_in_ca_because_pa_requires_a_non_edge',
          all(first_avoiding(complete(n), fast=True) is not None for n in range(4, 8)),
          'K_4..K_7 all in C_a')
    check('edgeless_graphs_are_in_ca',
          all(first_avoiding(graph(n, []), fast=True) is not None for n in range(4, 8)),
          'n=4..7')
    check('all_24_orders_of_k4_avoid_both_pa_and_pb_as_the_source_observation_states',
          count_avoiding(complete(4), S=('a',))[0] == 24
          and count_avoiding(complete(4), S=('b',))[0] == 24,
          '24 and 24')
    C6 = cycle(6)
    check('c6_is_in_ca_via_the_cycle_order_which_the_source_parenthetical_needs',
          occurrence_pa_fast(C6, (0, 1, 2, 3, 4, 5)) is None,
          'the cycle order has no crossing pair of edges at all')
    for n in (5, 6, 7):
        r = is_cocomparability_by_search(cycle(n))
        check('c%d_is_not_a_cocomparability_graph_reproducing_the_gallai_control' % n, r is False,
              'no transitive orientation among the 2^%d orientations of the complement'
              % len(complement(cycle(n))[2]))
    K23, v23 = complete_multipartite((2, 3))
    check('k23_is_in_ca_via_the_order_0_2_3_4_1',
          occurrence_pa_fast(K23, (0, 2, 3, 4, 1)) is None,
          'order 0 2 3 4 1')
    check('k2m_is_in_ca_for_m_up_to_5_via_the_order_built_here',
          all(occurrence_pa_fast(complete_multipartite((2, m))[0],
                                 tuple([0] + list(range(2, 2 + m)) + [1])) is None
              for m in range(2, 6)),
          'first vertex of the 2-part, then the m-part, then the second vertex of the 2-part')
    K34, _ = complete_multipartite((3, 4))
    check('k34_as_constructed_here_has_no_pa_avoiding_order',
          count_avoiding(K34, fast=True)[0] == 0, '0 of 5040 orders')
    check('k4_and_k23_are_not_in_c_empty_matching_the_forbidden_minors_of_outerplanar_graphs',
          count_avoiding(complete(4), S=())[0] == 0 and count_avoiding(K23, S=())[0] == 0,
          '0 of 24 and 0 of 120 orders avoid P_emptyset')

    # a certificate built here: it pins which optional pair is a and which is b
    SRC = graph(8, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (0, 7),
                    (0, 4), (2, 6)])
    SRC_ORDER = (0, 1, 3, 4, 5, 6, 2, 7)
    note('E(8-vertex graph built here) = %s   (C_8 with the chords 04 and 26, m=%d)'
         % (edge_str(SRC), len(SRC[2])))
    check('the_order_01345627_of_this_8_vertex_graph_avoids_pab_and_contains_pa',
          occurrence_naive(SRC, SRC_ORDER, ('a', 'b')) is None
          and occurrence_naive(SRC, SRC_ORDER, ('a',)) is not None,
          'order 0 1 3 4 5 6 2 7: avoids P_ab, contains P_a -- this fails if a and b are swapped')

    dis = 0
    for mask in range(1 << 6):
        pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
        H = graph(4, [pairs[k] for k in range(6) if (mask >> k) & 1])
        if (umbrella_free_order(H) is not None) != (is_cocomparability_by_search(H) is True):
            dis += 1
    check('the_two_cocomparability_deciders_agree_on_all_64_labeled_graphs_on_4_vertices',
          dis == 0, 'umbrella-avoidance versus transitive orientation of the complement')
    print('')

    print('NOTE SCOPE / NOT RE-RUN. Theorem 2 is proved in the paper for all orders; the '
          'exhaustive cross-check above covers 2 <= n <= 7 only, and the explicit sufficiency '
          'orders 2 <= n <= 10. The claim that the umbrella characterisation of '
          'cocomparability agrees with "the complement has a transitive orientation" is '
          'checked exhaustively at n=4 and on every named graph here, and is otherwise '
          'quoted from the literature, not proved. The minimum-order statement is exhaustive '
          'over ALL labeled graphs on at most 5 vertices; at order 6 this program exhibits '
          'three pairwise non-isomorphic witnesses and does NOT claim they are the only ones. '
          'No claim is made or checked here about orders 7 and above, about a finite '
          'forbidden-subgraph characterisation of C_a, or about the inclusion '
          '"cocomparability graphs are interval filament graphs", which is quoted from the '
          'source and not re-derived.')

    if FAILS:
        print('VERDICT: %d OF %d CHECKS FAILED' % (FAILS, PASSES + FAILS))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % PASSES)
    return 0


if __name__ == '__main__':
    sys.exit(main())
