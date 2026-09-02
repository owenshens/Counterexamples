#!/usr/bin/env python3
"""
verify.py -- verification program for

    Two exact wheel values for cover-free families on graphs:
    t(W_11) = t(W_12) = 8

Python 3.9+, STANDARD LIBRARY ONLY (sys, os, time, itertools, multiprocessing).
No third-party package, no external data file, no network, no randomness.
All arithmetic is on integers and finite sets; no floating point value is ever
compared or branched on (wall clocks are printed, never tested).

WHAT THIS PROGRAM READS AS INPUT
    Only the objects PRINTED IN THE PAPER: the two families R10 and R11 of
    Section 2, transcribed below, and the two integers 8 and 8 that the source's
    Table 4 prints for t(W_11) and t(W_12).  Everything else it derives.

WHAT IT DERIVES
    * the upper bound, by checking the printed families against the definition
      of a G-CFF and then applying the hub construction of Lemma 3;
    * the lower bound, by re-running the two exhaustive censuses in full:
      no W_11-CFF(7,11) and no W_12-CFF(7,12) exists.

    The accompanying note reports only the witness checks and the two complete
    t=7 censuses; the auxiliary checker controls and small-case recomputations
    kept below are not reported by the note.

One `PASS <name> [detail]` line per check, then
    VERDICT: ALL <n> CHECKS PASS
and exit status 0 iff every check passed.

RUNTIME.  The two decisive censuses visit 14,002,112 and 14,015,792 search-tree
nodes.  They are sharded over 1024 disjoint prefixes and run through a
multiprocessing Pool, so wall clock is about 2 x 130 s on 13 cores and about
2 x 1000 s on one core.  Everything else in this program takes under a minute.
"""

import itertools
import os
import sys
import time
from multiprocessing import Pool

# ---------------------------------------------------------------------------
# 0.  THE OBJECTS PRINTED IN THE PAPER  (Section 2)
# ---------------------------------------------------------------------------
# R10: the C_10-CFF(7,10) of the paper, ground set {1,...,7}, vertex order
# v_1, ..., v_10 around the cycle.
R10 = [{1}, {2, 3}, {2, 4}, {2, 5}, {2, 6}, {6, 7}, {3, 6}, {3, 5}, {5, 7}, {4, 7}]

# R11: the C_11-CFF(7,11) of the paper, ground set {1,...,7}.
R11 = [{1}, {2, 3, 4}, {2, 3, 5}, {2, 3, 6}, {2, 7}, {2, 4, 6}, {4, 5, 6},
       {3, 4, 6}, {3, 4, 7}, {3, 5, 7}, {5, 6, 7}]

# The hub set of the paper's two wheel families; 8 is the new ground-set element.
HUB = {8}

# TRANSCRIBED, NOT DERIVED -- the two non-bold entries of Table 4 of
# arXiv:2605.12634v1 whose exactness this note settles.
TABLE4_NONBOLD_WHEEL = {11: 8, 12: 8}

# TRANSCRIBED, NOT DERIVED -- the six BOLD entries of the t(W_n) column of the
# same table.  These are the source's own published exact values and are used
# here as controls: this program recomputes each of them from scratch.
TABLE4_BOLD_WHEEL = {5: 5, 6: 6, 7: 6, 8: 7, 9: 7, 10: 7}

# TRANSCRIBED, NOT DERIVED -- the node counts and shard counts the paper reports
# for the two decisive censuses (Section 3).
REPORTED = {
    11: {"nodes": 14002112, "shards": 1024},
    12: {"nodes": 14015792, "shards": 1024},
}

# ---------------------------------------------------------------------------
# 1.  THE DEFINITION, IMPLEMENTED LITERALLY ON PYTHON SETS
# ---------------------------------------------------------------------------
# Parida-Moura, "Cover-free families on graphs", arXiv:2605.12634v1,
# the two definitions of a G-cover-free family:
#
#   G = ([1,n], E).  A family (B_1, ..., B_n) of subsets of [1,t], B_i on
#   vertex i, is a G-CFF(t,n) when
#     (S)  for every edge {a,b} in E:  B_a \ B_b != {}  and  B_b \ B_a != {};
#     (CF) for every edge {a,b} in E and every i0 in [1,n] \ {a,b}:
#            B_{i0} \ (B_a u B_b) != {}.
#   t(G) = min{ t : a G-CFF(t, |V(G)|) exists }.


def violations(fam, n, edges, t):
    """Every way `fam` fails to be a G-CFF(t,n).  Empty list == it is one."""
    bad = []
    ground = set(range(1, t + 1))
    if len(fam) != n:
        bad.append(("wrong-length", len(fam), n))
        return bad
    for i, B in enumerate(fam):
        if not B <= ground:
            bad.append(("outside-ground-set", i + 1, sorted(B)))
    for (a, b) in edges:                          # vertices are 0-indexed here
        A, Bb = fam[a], fam[b]
        if not (A - Bb):
            bad.append(("S-fail", a + 1, b + 1))
        if not (Bb - A):
            bad.append(("S-fail", b + 1, a + 1))
        U = A | Bb
        for j in range(n):
            if j in (a, b):
                continue
            if not (fam[j] - U):
                bad.append(("CF-fail", j + 1, (a + 1, b + 1)))
    return bad


def is_cff(fam, n, edges, t):
    return not violations(fam, n, edges, t)


def n_cf_checks(n, edges):
    """How many (edge, outside vertex) pairs the CF clause imposes."""
    return len(edges) * (n - 2)


def is_antichain(fam):
    return all(fam[i] - fam[j] for i in range(len(fam)) for j in range(len(fam)) if i != j)


# ---------------------------------------------------------------------------
# 2.  THE GRAPHS
# ---------------------------------------------------------------------------
def path(n):
    return [(i, i + 1) for i in range(n - 1)]


def cycle(n):
    return [(i, (i + 1) % n) for i in range(n)]


def wheel(n):
    """W_n on n vertices: vertex 0 is the hub, vertices 1..n-1 the rim C_{n-1}."""
    r = n - 1
    e = [(0, i) for i in range(1, n)]
    e += [(i, i + 1) for i in range(1, r)] + [(r, 1)]
    return e


def degrees(n, edges):
    d = [0] * n
    for (a, b) in edges:
        d[a] += 1
        d[b] += 1
    return d


# ---------------------------------------------------------------------------
# 3.  THE SEARCH  (bitmask DFS over ALL 2^t subsets per vertex)
# ---------------------------------------------------------------------------
# The only pruning is the definitional constraints themselves, each tested as
# soon as every vertex it mentions has been assigned.  No derived lemma is used
# as a prune, so an empty census is a statement about the definition alone.
#
# SYMMETRY.  The constraint system is built from set differences only, so it is
# invariant under the full symmetric group S_t acting on the ground set.  Hence
# the first vertex in the search order may be restricted to one representative
# per size class -- the masks 0, 1, 3, 7, ..., 2^t - 1 -- and that reduction is
# exact.  `canon=False` enumerates the unreduced space; check 6 uses it.

def canon_reps(t):
    return [0] + [(1 << k) - 1 for k in range(1, t + 1)]


def mask2set(m, t):
    return set(i + 1 for i in range(t) if (m >> i) & 1)


def search(n, edges, t, order=None, canon=True, limit=1, prefix=None):
    """(nsol, witnesses_as_sets, nodes).  limit=None counts ALL solutions."""
    if order is None:
        order = list(range(n))
    pos = {v: k for k, v in enumerate(order)}
    full = (1 << t) - 1
    S_at = [[] for _ in range(n)]
    CF_at = [[] for _ in range(n)]
    for (a, b) in edges:
        da, db = pos[a], pos[b]
        S_at[max(da, db)].append((da, db))
        for j in range(n):
            if j in (a, b):
                continue
            dj = pos[j]
            CF_at[max(da, db, dj)].append((da, db, dj))
    cands0 = list(range(1 << t))
    reps = canon_reps(t)
    B = [0] * n
    sols = []
    stats = {"nsol": 0, "nodes": 0}

    def rec(d):
        stats["nodes"] += 1
        if d == n:
            stats["nsol"] += 1
            if len(sols) < 200:
                sols.append(list(B))
            return stats["nsol"] == limit
        if prefix is not None and d < len(prefix):
            choices = [prefix[d]]
        elif d == 0 and canon:
            choices = reps
        else:
            choices = cands0
        for m in choices:
            B[d] = m
            ok = True
            for (da, db) in S_at[d]:
                if not (B[da] & ~B[db] & full) or not (B[db] & ~B[da] & full):
                    ok = False
                    break
            if ok:
                for (da, db, dj) in CF_at[d]:
                    if not (B[dj] & ~(B[da] | B[db]) & full):
                        ok = False
                        break
            if ok and rec(d + 1):
                return True
        return False

    sys.setrecursionlimit(10000)
    rec(0)
    out = []
    for s in sols:
        out.append([mask2set(s[pos[v]], t) for v in range(n)])
    return stats["nsol"], out, stats["nodes"]


def shard_prefixes(t):
    """The 1024 shards of the paper: every (canonical depth-0 mask, any depth-1
    mask).  They are pairwise disjoint and their union is the whole canonical
    space, because a shard is exactly the set of leaves below a fixed pair of
    first two assignments and every canonical leaf has such a pair."""
    return [(m1, m2) for m1 in canon_reps(t) for m2 in range(1 << t)]


def _task(arg):
    n, edges, t, order, pre = arg
    ns, w, nodes = search(n, edges, t, order=order, canon=True, limit=None, prefix=pre)
    return (ns, nodes)


def census(n, edges, t, order=None, nproc=None):
    """Full sharded census of the canonical space.  (nsol, nodes, nshards)."""
    prefixes = shard_prefixes(t)
    args = [(n, edges, t, order, p) for p in prefixes]
    np_ = nproc or max(1, (os.cpu_count() or 2) - 1)
    if np_ == 1:
        res = [_task(a) for a in args]
    else:
        with Pool(np_) as pool:
            res = pool.map(_task, args, chunksize=1)
    return sum(r[0] for r in res), sum(r[1] for r in res), len(res)


# ---------------------------------------------------------------------------
# 4.  THE HUB CONSTRUCTION  (Lemma 3 of the paper)
# ---------------------------------------------------------------------------
def hub_extend(rim, t):
    """(B_0, B_1, ..., B_m) with B_0 = {t+1}: the family of Lemma 3."""
    return [{t + 1}] + [set(B) for B in rim]


# ---------------------------------------------------------------------------
# 5.  THE CHECK HARNESS
# ---------------------------------------------------------------------------
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, (" [%s]" % detail) if detail else ""),
          flush=True)


def main():
    t0 = time.time()
    print("verification of: two exact wheel values for cover-free families on graphs")
    print("                 t(W_11) = t(W_12) = 8   (Parida-Moura, arXiv:2605.12634v1, Table 4)")
    print("checks the two printed witnesses against (S)-(CF), their hub extensions, and")
    print("re-runs the two complete t=7 censuses; the note reports only those checks,")
    print("the other checks below being auxiliary controls the note does not report")
    print("python %s, standard library only, exact set and integer arithmetic" % sys.version.split()[0])
    print("cores available: %s" % (os.cpu_count(),))
    print("")

    # -- 5.1  the definition checker itself ---------------------------------
    print("--- the definition checker ---")
    # (S) must fail, and ONLY (S), when one end of an edge contains the other.
    bad = violations([{1, 4}, {1}, {2}, {3}], 4, path(4), 5)
    check("checker_flags_an_S_violation",
          [v[0] for v in bad] == ["S-fail"] and bad[0][1:] == (2, 1),
          "on P_4 with B_2 = {1} inside B_1 = {1,4}: exactly %d violation(s), kinds %s, "
          "the offending ordered pair being (B_2, B_1)"
          % (len(bad), sorted(set(v[0] for v in bad))))

    # (CF) must fail, and ONLY (CF), when a set is swallowed by an edge's union.
    bad = violations([{1}, {2}, {3}, {1, 2}], 4, path(4), 3)
    check("checker_flags_a_CF_violation",
          bad and set(v[0] for v in bad) == {"CF-fail"},
          "on P_4 with B_4 = {1,2} covered by B_1 u B_2: %d violation(s), kinds %s, "
          "no (S) failure" % (len(bad), sorted(set(v[0] for v in bad))))

    # A hand-verified positive: on P_3, three singletons.
    fam = [{1}, {2}, {3}]
    check("checker_accepts_a_hand_verified_P3_CFF",
          is_cff(fam, 3, path(3), 3),
          "({1},{2},{3}) is a P_3-CFF(3,3): 2 edges, %d CF conditions, 0 violations"
          % n_cf_checks(3, path(3)))

    # The bitmask DFS must agree with the literal set checker, cell by cell, on
    # every assignment -- not just on the ones the DFS reaches.
    agree = []
    for (lbl, n, edges, t) in [("P_4 t=3", 4, path(4), 3), ("C_4 t=3", 4, cycle(4), 3),
                              ("C_4 t=4", 4, cycle(4), 4), ("W_5 t=3", 5, wheel(5), 3)]:
        brute = 0
        subsets = [mask2set(m, t) for m in range(1 << t)]
        for tup in itertools.product(subsets, repeat=n):
            if is_cff(list(tup), n, edges, t):
                brute += 1
        ns, _, _ = search(n, edges, t, canon=False, limit=None)
        agree.append((lbl, brute, ns))
    check("brute_force_and_bitmask_search_agree_on_the_four_small_instances",
          all(b == s for (_, b, s) in agree),
          "; ".join("%s: literal %d = DFS %d" % x for x in agree))

    # The canonical reduction: one representative per S_t-orbit of a subset of
    # [t], and orbits of subsets are exactly the size classes.
    reps = canon_reps(7)
    sizes = sorted(bin(m).count("1") for m in reps)
    covered = all(any(bin(m).count("1") == bin(r).count("1") for r in reps) for m in range(1 << 7))
    check("canonical_representatives_meet_every_subset_size_class_at_t_7",
          len(reps) == 8 and sizes == list(range(8)) and covered,
          "t=7: %d reps of sizes %s; all 128 subsets have a rep of their size" % (len(reps), sizes))

    # Emptiness must not depend on the reduction, nor on the vertex order.
    unred = [("W_9 t=6", 9, wheel(9), 6), ("W_11 t=6", 11, wheel(11), 6)]
    urows = []
    for (lbl, n, edges, t) in unred:
        nc, _, nodes_c = search(n, edges, t, canon=True, limit=None)
        nu, _, nodes_u = search(n, edges, t, canon=False, limit=None)
        urows.append((lbl, nc, nodes_c, nu, nodes_u))
    orders = []
    for (lbl, n) in [("W_11", 11), ("W_12", 12)]:
        rim = list(range(1, n))
        half = (len(rim) + 1) // 2
        inter = [0] + [v for pair in itertools.zip_longest(rim[:half], rim[half:])
                       for v in pair if v is not None]
        na, _, nodes_a = search(n, wheel(n), 6, limit=None)
        nb, _, nodes_b = search(n, wheel(n), 6, order=inter, limit=None)
        orders.append((lbl, na, nb, nodes_a, nodes_b))
    check("canonical_reduction_and_vertex_order_agree_on_the_cases_tested",
          all(r[1] == r[3] == 0 for r in urows) and all(o[1] == o[2] == 0 for o in orders),
          "; ".join("%s canonical %d sols/%d nodes vs unreduced %d sols/%d nodes"
                    % r for r in urows)
          + "; " + "; ".join("%s t=6 natural %d sols/%d nodes vs interleaved %d sols/%d nodes"
                             % (o[0], o[1], o[3], o[2], o[4]) for o in orders))

    # S_t-invariance, exercised on the paper's own two families.
    W11 = hub_extend(R10, 7)
    W12 = hub_extend(R11, 7)
    perms = list(itertools.permutations(range(1, 9)))[::401]
    inv = True
    for p in perms:
        pi = {i + 1: p[i] for i in range(8)}
        for (fam, n, t) in [(W11, 11, 8), (W12, 12, 8)]:
            img = [set(pi[x] for x in B) for B in fam]
            if not is_cff(img, n, wheel(n), t):
                inv = False
    check("validity_is_preserved_under_the_sampled_ground_set_relabellings",
          inv,
          "%d of the 40320 permutations of [8], applied to both printed families, "
          "every image still a CFF" % len(perms))

    # -- 5.2  the two graphs ------------------------------------------------
    print("")
    print("--- the graphs ---")
    for (n, rimlen, nedge) in [(11, 10, 20), (12, 11, 22)]:
        e = wheel(n)
        d = degrees(n, e)
        rimedges = [x for x in e if 0 not in x]
        ok = (len(e) == nedge and d[0] == n - 1 and set(d[1:]) == {3}
              and len(rimedges) == rimlen
              and sorted(tuple(sorted(x)) for x in rimedges)
              == sorted(tuple(sorted((a + 1, b + 1))) for (a, b) in cycle(rimlen)))
        check("W_%d_is_the_hub_plus_rim_C_%d" % (n, rimlen), ok,
              "%d vertices, %d edges, hub degree %d, all %d rim degrees 3, rim edges "
              "form C_%d" % (n, len(e), d[0], n - 1, rimlen))

    # -- 5.3  the printed objects and the upper bound -----------------------
    print("")
    print("--- the printed families, and the upper bound ---")
    for (nm, R, m) in [("R10", R10, 10), ("R11", R11, 11)]:
        ground = set().union(*R)
        ok = (len(R) == m and all(R) and len({frozenset(B) for B in R}) == m
              and ground <= set(range(1, 8)))
        check("%s_is_%d_distinct_nonempty_subsets_of_a_7_set" % (nm, m), ok,
              "%d sets, %d distinct, largest element %d, sizes %s"
              % (len(R), len({frozenset(B) for B in R}), max(ground),
                 sorted(set(len(B) for B in R))))
        check("%s_is_an_antichain" % nm, is_antichain(R),
              "all %d ordered pairs have B_i \\ B_j nonempty" % (m * (m - 1)))
        v = violations(R, m, cycle(m), 7)
        check("%s_is_a_C_%d_CFF_7_%d" % (nm, m, m), not v,
              "%d edges incl. the wrap edge {%d,1}, %d Sperner conditions, %d CF "
              "conditions, %d violations"
              % (len(cycle(m)), m, 2 * len(cycle(m)), n_cf_checks(m, cycle(m)), len(v)))

    for (n, fam, R) in [(11, W11, R10), (12, W12, R11)]:
        v = violations(fam, n, wheel(n), 8)
        ok = (not v) and fam[0] == HUB and fam[1:] == [set(B) for B in R]
        check("hub_construction_gives_a_W_%d_CFF_8_%d" % (n, n), ok,
              "B_0 = %s over the printed rim: %d edges, %d CF conditions, %d violations "
              "-> t(W_%d) <= 8" % (sorted(fam[0]), len(wheel(n)), n_cf_checks(n, wheel(n)),
                                   len(v), n))

    # Lemma 3 exercised on the small cycles, independently of R10 and R11.
    lem = []
    for (m, t) in [(4, 4), (5, 5), (6, 5), (7, 6), (8, 6), (9, 6)]:
        ns, w, _ = search(m, cycle(m), t, limit=1)
        if not ns:
            lem.append((m, t, False, "no C_%d-CFF(%d,%d) found" % (m, t, m)))
            continue
        R = w[0]
        anti = is_antichain(R)
        ext = is_cff(hub_extend(R, t), m + 1, wheel(m + 1), t + 1)
        lem.append((m, t, anti and ext, "C_%d-CFF(%d,%d) -> W_%d-CFF(%d,%d)"
                    % (m, t, m, m + 1, t + 1, m + 1)))
    check("hub_construction_lemma_holds_on_the_small_cycles_m_4_to_9",
          all(r[2] for r in lem),
          "; ".join(r[3] for r in lem))

    # -- 5.4  controls: the source's own bold wheel values ------------------
    print("")
    print("--- controls: every bold t(W_n) entry of Table 4, recomputed ---")
    wits = {}
    for n in sorted(TABLE4_BOLD_WHEEL):
        tstar = TABLE4_BOLD_WHEEL[n]
        below, _, nb = search(n, wheel(n), tstar - 1, limit=None)
        at, w, na = search(n, wheel(n), tstar, limit=1)
        good = bool(w) and is_cff(w[0], n, wheel(n), tstar)
        if w:
            wits[n] = (tstar, w[0])
        check("published_bold_wheel_value_t_W_%d_is_%d" % (n, tstar),
              below == 0 and at >= 1 and good,
              "census at t=%d: %d solutions in %d nodes; witness at t=%d verified by the "
              "literal checker: %s" % (tstar - 1, below, nb, tstar, good))

    # Monotonicity in t, the step that turns an empty census into a lower bound.
    mono = True
    rows = []
    src = dict(wits)
    src[11] = (8, W11)
    src[12] = (8, W12)
    for n in sorted(src):
        tstar, fam = src[n]
        if not is_cff(fam, n, wheel(n), tstar + 1):
            mono = False
        rows.append("W_%d: %d->%d" % (n, tstar, tstar + 1))
    check("monotonicity_in_t_holds_for_the_witnesses_collected_here", mono,
          "each witness re-verified over a ground set one larger: " + ", ".join(rows))

    # The differential control: the SAME code at the SAME t = 7 answers YES on
    # W_10 and NO on W_11, so the negatives below are about n, not about the
    # decider's ability to search a wheel at t = 7.
    at10, w10, nodes10 = search(10, wheel(10), 7, limit=1)
    check("differential_control_same_code_t_7_solvable_for_W_10",
          at10 >= 1 and is_cff(w10[0], 10, wheel(10), 7),
          "W_10-CFF(7,10) found in %d nodes and verified: %s"
          % (nodes10, [sorted(B) for B in w10[0]]))

    # -- 5.5  the decisive censuses -----------------------------------------
    print("")
    print("--- the censuses ---")
    for n in (11, 12):
        ns, _, nodes = search(n, wheel(n), 6, limit=None)
        check("W_%d_has_no_CFF_at_t_6" % n, ns == 0,
              "%d solutions in %d nodes" % (ns, nodes))

    pre = shard_prefixes(7)
    ok = (len(pre) == 1024 and len(set(pre)) == 1024
          and len(canon_reps(7)) * (1 << 7) == 1024
          and set(pre) == {(a, b) for a in canon_reps(7) for b in range(1 << 7)})
    check("the_1024_shards_are_a_disjoint_cover_of_the_canonical_space", ok,
          "8 canonical depth-0 masks x 128 depth-1 masks = %d pairwise distinct prefixes, "
          "and every canonical leaf lies under exactly one" % len(pre))

    nproc = max(1, (os.cpu_count() or 2) - 1)
    print("    running the two decisive censuses on %d process(es); "
          "expect a few minutes" % nproc, flush=True)
    for n in (11, 12):
        c0 = time.time()
        ns, nodes, nsh = census(n, wheel(n), 7, nproc=nproc)
        want = REPORTED[n]
        check("W_%d_census_at_t_7_is_empty" % n,
              ns == 0 and nsh == want["shards"] and nodes == want["nodes"],
              "SOLUTIONS=%d shards=%d/%d nodes=%d (paper reports 0, %d, %d) in %.1f s "
              "-> t(W_%d) >= 8" % (ns, nsh, want["shards"], nodes, want["shards"],
                                   want["nodes"], time.time() - c0, n))

    # -- 5.6  the conclusion ------------------------------------------------
    print("")
    print("--- the conclusion ---")
    by_name = {c[0]: c[1] for c in CHECKS}
    for n in (11, 12):
        lower = by_name.get("W_%d_census_at_t_7_is_empty" % n, False)
        upper = by_name.get("hub_construction_gives_a_W_%d_CFF_8_%d" % (n, n), False)
        check("t_W_%d_equals_8" % n, lower and upper,
              "lower bound from the empty t=7 census plus monotonicity in t; upper bound "
              "from the printed W_%d-CFF(8,%d)" % (n, n))
    by_name = {c[0]: c[1] for c in CHECKS}          # re-read: the two lines above added rows
    check("both_values_confirm_the_non_bold_table_entries_as_exact",
          TABLE4_NONBOLD_WHEEL == {11: 8, 12: 8}
          and by_name.get("t_W_11_equals_8") and by_name.get("t_W_12_equals_8"),
          "Table 4 prints t(W_11) <= 8 and t(W_12) <= 8 non-bold; both are equalities, so "
          "the source's bounds are confirmed and nothing published is contradicted")

    # -- 5.7  verdict --------------------------------------------------------
    npass = sum(1 for c in CHECKS if c[1])
    nfail = len(CHECKS) - npass
    print("")
    print("NOTE SCOPE: this program re-derives both claims of the accompanying note and every "
          "number that note prints for them. NOT RE-RUN: (a) the five remaining non-bold cells "
          "of Table 4 -- t(C_10), t(P_11), t(C_11), t(P_12), t(C_12), each printed there as 7 "
          "-- which are neither used nor rechecked here, so the statement 'Table 4 is exact "
          "throughout for n <= 12' is NOT verified by this program; "
          "(b) any t(W_n) with n >= 13; (c) any minimality or least-order claim "
          "-- nothing here says W_11 is the smallest wheel whose table entry was inexact; "
          "(d) the transcription of Table 4 itself and of the two integers 8, 8, which are read "
          "off arXiv:2605.12634v1 by hand and cannot be recomputed.")
    print("total wall clock %.1f s" % (time.time() - t0))
    if nfail:
        print("VERDICT: %d of %d CHECKS FAILED" % (nfail, len(CHECKS)))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % npass)
    return 0


if __name__ == "__main__":
    sys.exit(main())
