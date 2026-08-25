#!/usr/bin/env python3
"""
Independent verification of "a Steiner quadruple system of order 38 whose
derived Steiner triple system at every point has block-chromatic index 19"
(an mcDSQS(38)).

Standard library only; exact integer arithmetic throughout; no floats are
used for any decision.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs, transcribed verbatim, never trusted)
--------------------------------------------------------------------------
  * the point set X = Z_37 u {oo};
  * the group generators g_{a,b}(x) = 10^a x + b  (a in {0,1,2}, b in Z_37),
    fixing oo;
  * the 27 block-orbit representatives, in three families of sizes 6, 6, 15
    (A_REPS, S_REPS, L_REPS);
  * the stabiliser generators claimed for the first two families
    (STAB_A, STAB_S);
  * the two exhibited block colourings: table TAB_INF for the design derived
    at oo, table TAB_ZERO for the design derived at 0, each given as one set
    C_0 of ten triples and six sets G_1..G_6 of twelve triples;
  * the recipe "the 19 classes are C_0 together with mu^j(G_i) \\ C_0",
    mu(x) = 10x, 1 <= i <= 6, 0 <= j <= 2;
  * the recipe "translation by p carries the colouring of D_0 to D_p".

--------------------------------------------------------------------------
DERIVED HERE (nothing below is copied from the paper; all of it is computed)
--------------------------------------------------------------------------
  * that the 111 maps form a group of permutations of X and that 10 has
    multiplicative order 3 mod 37;
  * the exact stabiliser of every representative, hence every orbit length,
    hence |B| (claimed 2109);
  * that the 27 orbits are pairwise disjoint;
  * the SQS(38) property: every one of the C(38,3) triples of X lies in
    exactly one block  <-- hypothesis of the statement acted on;
  * that B is invariant under the whole group, that x -> x+1 is a 37-cycle
    on the finite points (the system is rotational) and that mu is a
    multiplier of order 3;
  * that the design derived at each of the 38 points is an STS(37) with 222
    blocks;
  * that both exhibited tables expand to 19 pairwise-disjoint classes of
    genuine blocks of the relevant derived design, that the classes are
    partial parallel classes (triples within a class pairwise disjoint),
    that they cover all 222 blocks, and their size profile;
  * that every entry of both printed tables is a 3-subset of X and that the
    rows have the stated sizes 10 and 12;
  * that the block set the printed table exhibits at 0 (not merely the one
    recomputed from B, which would make the statement automatic) is carried
    by each translation exactly onto the derived design at that point;
  * that transporting the colouring of D_0 by every translation gives a
    valid 19-class colouring of all 37 finite-point derived designs, so all
    38 derived designs are 19-colourable  <-- load-bearing conclusion;
  * the lower bound 19 by exact pigeonhole from the computed block count and
    the computed maximum size of a partial parallel class;
  * hence chi'(D_x) = 19 for all x, chi'(37) = 19, minimum colourability of
    every derived design, and q'_0(3,4,38) = 1 + 19 = 20.
"""
import sys
from itertools import combinations

CHECKS = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (" [" + detail + "]" if detail else ""))
    return bool(ok)


def finish():
    n = len(CHECKS)
    bad = [c for c, o in CHECKS if not o]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        sys.exit(1)
    print("VERDICT: ALL %d CHECKS PASS" % n)
    sys.exit(0)


# ------------------------------------------------------------------ inputs
N = 37
INF = 37                      # symbol for the point at infinity
X = tuple(range(N)) + (INF,)  # the 38 points

A_REPS = [(INF, 0, 1, 11), (INF, 0, 2, 22), (INF, 0, 4, 7),
          (INF, 0, 8, 14), (INF, 0, 16, 28), (INF, 0, 19, 32)]
S_REPS = [(0, 1, 27, 34), (0, 2, 17, 31), (0, 3, 7, 28),
          (0, 2, 19, 24), (0, 5, 23, 29), (0, 1, 12, 28)]
L_REPS = [(0, 1, 13, 14), (0, 1, 8, 21), (0, 1, 9, 20),
          (0, 1, 10, 31), (0, 1, 16, 23), (0, 1, 2, 7),
          (0, 1, 5, 33), (0, 1, 4, 18), (0, 1, 17, 30),
          (0, 1, 3, 15), (0, 1, 19, 35), (0, 2, 9, 15),
          (0, 2, 8, 33), (0, 2, 18, 23), (0, 2, 25, 32)]
STAB_A = [1, 2, 4, 8, 16, 32]        # claimed b with g_{1,b} stabilising
STAB_S = [27, 17, 7, 19, 29, 28]     # the matching representative

TAB_INF = {
    "C0": [(1, 10, 26), (2, 15, 20), (3, 11, 17), (5, 13, 19), (6, 18, 27),
           (7, 33, 34), (8, 16, 22), (9, 35, 36), (14, 29, 31), (21, 25, 28)],
    "G1": [(0, 30, 34), (2, 8, 31), (3, 11, 17), (4, 5, 15), (6, 19, 24),
           (7, 9, 29), (10, 12, 32), (13, 21, 27), (14, 22, 28),
           (16, 20, 23), (18, 33, 35), (25, 26, 36)],
    "G2": [(0, 5, 24), (1, 17, 29), (3, 8, 27), (4, 16, 25), (6, 10, 13),
           (7, 12, 31), (9, 35, 36), (11, 23, 32), (14, 18, 21),
           (15, 19, 22), (20, 28, 34), (26, 30, 33)],
    "G3": [(0, 23, 31), (1, 6, 25), (2, 28, 29), (3, 26, 34), (4, 19, 21),
           (7, 13, 36), (8, 16, 22), (9, 10, 20), (11, 15, 18),
           (12, 24, 33), (14, 27, 32), (17, 30, 35)],
    "G4": [(0, 1, 11), (2, 18, 30), (3, 16, 21), (4, 6, 26), (5, 31, 32),
           (8, 20, 29), (9, 22, 27), (10, 14, 17), (12, 13, 23),
           (15, 28, 33), (19, 34, 36), (24, 25, 35)],
    "G5": [(0, 16, 28), (1, 24, 32), (2, 7, 26), (3, 29, 30), (4, 12, 18),
           (5, 20, 22), (6, 21, 23), (8, 17, 33), (9, 11, 31),
           (10, 19, 35), (13, 25, 34), (15, 27, 36)],
    "G6": [(0, 2, 22), (1, 9, 15), (3, 18, 20), (4, 8, 11), (5, 14, 30),
           (6, 7, 17), (10, 25, 27), (12, 16, 19), (13, 26, 31),
           (23, 24, 34), (28, 32, 35), (29, 33, 36)],
}

TAB_ZERO = {
    "C0": [(1, 5, 33), (3, 4, 30), (6, 8, 23), (7, 19, 26), (9, 12, 16),
           (10, 13, 34), (11, 27, 36), (14, 21, 29), (17, 22, 35),
           (18, 24, 32)],
    "G1": [(1, 16, 23), (2, 14, 36), (4, 6, 12), (5, 18, 20), (7, 19, 26),
           (8, 13, 31), (9, 22, 28), (10, 24, 35), (11, 21, 25),
           (15, 27, 33), (17, 29, 32), (30, 34, INF)],
    "G2": [(1, 2, 7), (3, 5, 14), (4, 20, 26), (8, 18, 27), (9, 10, 21),
           (11, 24, 33), (12, 34, 35), (13, 19, 28), (15, 17, INF),
           (16, 29, 36), (22, 23, 25), (30, 31, 32)],
    "G3": [(1, 10, 31), (2, 8, 33), (3, 6, 22), (4, 16, 27), (5, 35, 36),
           (7, 18, 25), (11, 13, 23), (12, 15, 20), (14, 21, 29),
           (17, 28, 34), (19, 32, INF), (24, 26, 30)],
    "G4": [(1, 24, 25), (2, 6, 13), (3, 19, 34), (4, 21, 23), (5, 10, 17),
           (7, 12, 14), (8, 20, 32), (9, 18, 31), (15, 29, 35),
           (16, 28, INF), (22, 27, 30), (26, 33, 36)],
    "G5": [(1, 9, 20), (2, 27, 35), (3, 13, 16), (4, 10, 22), (5, 12, 32),
           (6, 29, INF), (7, 8, 24), (11, 18, 30), (14, 33, 34),
           (15, 19, 25), (21, 31, 36), (23, 26, 28)],
    "G6": [(1, 11, INF), (2, 3, 21), (4, 25, 34), (5, 26, 31), (6, 7, 16),
           (8, 22, 29), (9, 32, 35), (10, 23, 27), (12, 13, 36),
           (14, 15, 30), (17, 19, 33), (20, 24, 28)],
}
GNAMES = ("G1", "G2", "G3", "G4", "G5", "G6")


# ------------------------------------------------- group of order 111 on X
def pow10(a):
    """10**a mod 37 by exact integer arithmetic."""
    r = 1
    for _ in range(a):
        r = (r * 10) % N
    return r


def perm(a, b):
    """g_{a,b} as an explicit tuple-indexed permutation of the 38 points."""
    m = pow10(a)
    p = [(m * x + b) % N for x in range(N)]
    p.append(INF)
    return tuple(p)


GROUP = [perm(a, b) for a in range(3) for b in range(N)]
MU = perm(1, 0)          # the multiplier x -> 10x, fixing oo and 0
TRANS = [perm(0, p) for p in range(N)]


def act(p, blk):
    """Apply a permutation to a block, returning a canonical sorted tuple."""
    return tuple(sorted(p[y] for y in blk))


def orbit(blk):
    return set(act(p, blk) for p in GROUP)


def check_group():
    ok = pow10(3) == 1 and pow10(1) != 1 and pow10(2) != 1
    d = "10^1=%d 10^2=%d 10^3=%d mod 37" % (pow10(1), pow10(2), pow10(3))
    ck("mult_order_of_10_is_3", ok, d)
    distinct = set(GROUP)
    bij = all(len(set(p)) == len(X) and set(p) == set(X) for p in GROUP)
    ck("group_has_111_distinct_permutations",
       len(GROUP) == 111 and len(distinct) == 111 and bij,
       "|G|=%d distinct=%d all_bijections=%s" % (len(GROUP), len(distinct), bij))
    closed = True
    for p in GROUP:
        for q in GROUP:
            if tuple(p[q[i]] for i in range(len(X))) not in distinct:
                closed = False
                break
        if not closed:
            break
    fixinf = all(p[INF] == INF for p in GROUP)
    ck("group_closed_and_fixes_infinity", closed and fixinf,
       "closed=%s fixes_oo=%s" % (closed, fixinf))


def check_representatives():
    reps = A_REPS + S_REPS + L_REPS
    shapes = all(len(set(r)) == 4 and all(y in X for y in r) for r in reps)
    ck("27_representatives_are_4_subsets_of_X",
       len(reps) == 27 and len(A_REPS) == 6 and len(S_REPS) == 6
       and len(L_REPS) == 15 and shapes,
       "|A|=%d |S|=%d |L|=%d total=%d well_formed=%s"
       % (len(A_REPS), len(S_REPS), len(L_REPS), len(reps), shapes))
    a_inf = all(INF in r for r in A_REPS)
    fin = all(INF not in r for r in S_REPS + L_REPS)
    zero = all(0 in r for r in reps)
    ck("representative_families_have_stated_shape", a_inf and fin and zero,
       "A_all_contain_oo=%s S_and_L_finite=%s all_contain_0=%s"
       % (a_inf, fin, zero))


def check_stabilisers():
    bad = []
    for fam, claimed in (("A", STAB_A), ("S", STAB_S)):
        reps = A_REPS if fam == "A" else S_REPS
        for r, b in zip(reps, claimed):
            stab = [(a, bb) for a in range(3) for bb in range(N)
                    if act(perm(a, bb), r) == tuple(sorted(r))]
            gen = perm(1, b)
            gen_ok = act(gen, r) == tuple(sorted(r))
            gen_pow = set()
            q = perm(0, 0)
            for _ in range(3):
                q = tuple(gen[q[i]] for i in range(len(X)))
                gen_pow.add(q)
            spans = gen_pow == set(perm(a, bb) for (a, bb) in stab)
            if not (len(stab) == 3 and gen_ok and spans):
                bad.append("%s%s|stab|=%d gen_fixes=%s spans=%s"
                           % (fam, r, len(stab), gen_ok, spans))
    ck("A_and_S_stabilisers_have_order_3_with_stated_generator", not bad,
       "12 representatives checked" if not bad else "; ".join(bad[:3]))
    trivial = []
    for r in L_REPS:
        stab = [(a, bb) for a in range(3) for bb in range(N)
                if act(perm(a, bb), r) == tuple(sorted(r))]
        if len(stab) != 1:
            trivial.append("%s|stab|=%d" % (str(r), len(stab)))
    ck("L_stabilisers_are_trivial", not trivial,
       "15 representatives checked" if not trivial
       else "; ".join(trivial[:3]))


def build_blocks():
    orbits = [orbit(r) for r in A_REPS + S_REPS + L_REPS]
    blocks = set()
    for o in orbits:
        blocks |= o
    return orbits, blocks


ORBITS, BLOCKS, DERIVED = [], set(), {}


def build_all():
    """Build B and the 38 derived designs, after the shape checks have run.

    Deferred rather than done at import time so that a malformed input datum
    produces a reported failure and a verdict line, not a traceback.
    """
    global ORBITS, BLOCKS, DERIVED
    ORBITS, BLOCKS = build_blocks()
    DERIVED = dict((x, derived(x)) for x in X)


def check_orbit_lengths_and_count():
    lens = [len(o) for o in ORBITS]
    short_ok = all(L == N for L in lens[:12])
    long_ok = all(L == 111 for L in lens[12:])
    ck("orbit_lengths_are_37_twelve_times_and_111_fifteen_times",
       short_ok and long_ok and len(lens) == 27,
       "lengths=%s" % (sorted(set(lens)),))
    total = sum(lens)
    disjoint = (len(BLOCKS) == total)
    ck("27_orbits_are_pairwise_disjoint", disjoint,
       "sum_of_orbit_sizes=%d distinct_blocks=%d" % (total, len(BLOCKS)))
    expect = 12 * N + 15 * 111
    ck("block_count_is_2109", len(BLOCKS) == 2109 == expect,
       "|B|=%d and 12*37+15*111=%d" % (len(BLOCKS), expect))


def check_sqs_property():
    """Every 3-subset of X in exactly one block: the defining SQS condition."""
    cover = {}
    for b in BLOCKS:
        for t in combinations(b, 3):
            cover[t] = cover.get(t, 0) + 1
    nchoose3 = 38 * 37 * 36 // 6
    multi = [t for t, c in cover.items() if c != 1]
    ck("every_triple_of_X_covered_exactly_once",
       len(cover) == nchoose3 == 8436 and not multi,
       "triples_seen=%d C(38,3)=%d triples_with_multiplicity!=1=%d"
       % (len(cover), nchoose3, len(multi)))
    sizes = set(len(b) for b in BLOCKS)
    ck("all_blocks_are_quadruples_inside_X",
       sizes == {4} and all(all(y in X for y in b) for b in BLOCKS),
       "block sizes present=%s" % (sorted(sizes),))


def check_group_invariance_and_rotation():
    bad = 0
    for p in GROUP:
        if set(act(p, b) for b in BLOCKS) != BLOCKS:
            bad += 1
    ck("block_set_invariant_under_all_111_group_elements", bad == 0,
       "group elements not preserving B: %d of %d" % (bad, len(GROUP)))
    sigma = perm(0, 1)                     # x -> x+1
    seen, y, cyc = set(), 0, 0
    while y not in seen:
        seen.add(y)
        y = sigma[y]
        cyc += 1
    rot = (cyc == 37 and sigma[INF] == INF
           and set(act(sigma, b) for b in BLOCKS) == BLOCKS)
    ck("system_is_rotational_37_cycle_fixing_infinity", rot,
       "cycle length on finite points=%d, oo fixed=%s" % (cyc, sigma[INF] == INF))
    q, order = perm(0, 0), 0
    while True:
        q = tuple(MU[q[i]] for i in range(len(X)))
        order += 1
        if q == perm(0, 0):
            break
    mult = (order == 3 and MU[INF] == INF and MU[0] == 0
            and set(act(MU, b) for b in BLOCKS) == BLOCKS)
    ck("mu_is_a_nontrivial_multiplier_of_order_3", mult,
       "order of mu=%d, fixes oo and 0=%s" % (order, MU[INF] == INF and MU[0] == 0))


def derived(x):
    """Blocks of the design derived at x, as canonical sorted triples."""
    return set(tuple(sorted(y for y in b if y != x))
               for b in BLOCKS if x in b)


def check_table_shapes():
    """Every printed table entry is a 3-subset of X, rows of size 10 and 12."""
    bad = []
    for label, table in (("inf", TAB_INF), ("zero", TAB_ZERO)):
        for name in ("C0",) + GNAMES:
            rows = table[name]
            for t in rows:
                if len(t) != 3 or len(set(t)) != 3 or not all(y in X for y in t):
                    bad.append("%s/%s bad entry %s" % (label, name, str(t)))
            want = 10 if name == "C0" else 12
            if len(rows) != want:
                bad.append("%s/%s has %d entries, expected %d"
                           % (label, name, len(rows), want))
    ck("both_printed_tables_hold_triples_of_X_in_rows_of_10_and_12", not bad,
       "2 tables x 7 rows = 164 triples checked" if not bad
       else "; ".join(bad[:3]))


def check_derived_are_sts37():
    bad = []
    for x in X:
        d = DERIVED[x]
        pts = set(X) - set([x])
        pair = {}
        for t in d:
            for p in combinations(t, 2):
                pair[p] = pair.get(p, 0) + 1
        want = 37 * 36 // 2
        if not (len(d) == 222 and len(pair) == want
                and all(c == 1 for c in pair.values())
                and set().union(*d) == pts and all(len(t) == 3 for t in d)):
            bad.append("x=%s blocks=%d pairs=%d" % (x, len(d), len(pair)))
    ck("all_38_derived_designs_are_STS_37_with_222_blocks", not bad,
       "38 derived designs, 222 blocks and 666 pairs each" if not bad
       else "; ".join(bad[:3]))


def expand_colouring(table):
    """Build the 19 classes from a table: C_0 plus mu^j(G_i) minus C_0."""
    c0 = set(tuple(sorted(t)) for t in table["C0"])
    classes = [c0]
    powers = [perm(0, 0), MU, tuple(MU[MU[i]] for i in range(len(X)))]
    for name in GNAMES:
        for j in range(3):
            img = set(act(powers[j], t) for t in table[name])
            classes.append(img - c0)
    return classes


def verify_table_colouring(label, table, design):
    """Return (classes, list of complaints) after checking every requirement."""
    classes = expand_colouring(table)
    bad = []
    if len(classes) != 19:
        bad.append("class count %d" % len(classes))
    for k, cls in enumerate(classes):
        for t in cls:
            if t not in design:
                bad.append("class %d holds non-block %s" % (k, str(t)))
        pts = [y for t in cls for y in t]
        if len(pts) != len(set(pts)):
            bad.append("class %d is not a partial parallel class" % k)
    flat = [t for cls in classes for t in cls]
    if len(flat) != len(set(flat)):
        bad.append("classes overlap")
    if set(flat) != design:
        bad.append("union misses %d / adds %d blocks"
                   % (len(design - set(flat)), len(set(flat) - design)))
    ck("colouring_%s_is_a_19_class_partition_into_parallel_classes" % label,
       not bad, "19 classes covering all %d blocks, no repeats" % len(design)
       if not bad else "; ".join(bad[:3]))
    prof = sorted((len(c) for c in classes), reverse=True)
    want = [12] * 14 + [11] * 4 + [10]
    ck("colouring_%s_has_profile_12pow14_11pow4_10" % label, prof == want,
       "profile=%s total=%d" % (prof, sum(prof)))
    return classes


def check_colourings():
    verify_table_colouring("at_infinity", TAB_INF, DERIVED[INF])
    verify_table_colouring("at_zero", TAB_ZERO, DERIVED[0])


def check_replication_numbers():
    pt = dict((y, 0) for y in X)
    pr = {}
    for b in BLOCKS:
        for y in b:
            pt[y] += 1
        for p in combinations(b, 2):
            pr[p] = pr.get(p, 0) + 1
    npt = set(pt.values())
    npr = set(pr.values())
    ck("point_and_pair_replication_numbers_are_222_and_18",
       npt == {222} and npr == {18} and len(pr) == 38 * 37 // 2,
       "point degrees=%s pair degrees=%s pairs=%d"
       % (sorted(npt), sorted(npr), len(pr)))


def transported_colouring(p):
    """Carry the colouring of D_0 to D_p by the translation x -> x+p."""
    tp = TRANS[p]
    return [set(act(tp, t) for t in cls)
            for cls in expand_colouring(TAB_ZERO)]


def check_all_38_derived_designs_19_colourable():
    """Load-bearing: an explicit 19-colouring is verified for every point."""
    upper = {}
    bad = []
    # The block set carried across is the one the printed table exhibits, not
    # the one recomputed from B; otherwise the transport statement would be an
    # automatic consequence of B being a union of orbits and could never fail.
    tab0 = set(t for cls in expand_colouring(TAB_ZERO) for t in cls)
    for x in X:
        design = DERIVED[x]
        if x == INF:
            classes = expand_colouring(TAB_INF)
        else:
            classes = transported_colouring(x)
            if (set(act(TRANS[x], t) for t in DERIVED[0]) != design
                    or set(act(TRANS[x], t) for t in tab0) != design):
                bad.append("translation by %d does not carry D_0 to D_%d"
                           % (x, x))
        flat = [t for cls in classes for t in cls]
        good = (len(classes) == 19 and set(flat) == design
                and len(flat) == len(set(flat))
                and all(len(set(y for t in cls for y in t)) == 3 * len(cls)
                        for cls in classes))
        if not good:
            bad.append("x=%s colouring invalid" % x)
        # an invalid colouring proves no upper bound at all
        upper[x] = len(classes) if good else 10 ** 6
    ck("translations_carry_D_0_to_every_finite_derived_design",
       not [m for m in bad if m.startswith("translation")],
       "37 translations checked")
    ck("every_one_of_the_38_derived_designs_has_a_verified_19_colouring",
       not bad and set(upper.values()) == {19},
       "38 designs, 19 classes each, %d classes verified disjoint"
       % (38 * 19) if not bad else "; ".join(bad[:3]))
    return upper


def check_lower_bound():
    """Exact pigeonhole from the computed block count; no float arithmetic."""
    bounds = {}
    for x in X:
        d = DERIVED[x]
        pts = len(set(y for t in d for y in t))
        cap_x = pts // 3
        bounds[x] = -((-len(d)) // cap_x)
    d = DERIVED[INF]
    npoints = len(set(y for t in d for y in t))
    nblocks = len(d)
    cap = npoints // 3                       # 37 // 3 = 12
    lower = -((-nblocks) // cap)             # ceil(222 / 12) = 19
    ck("lower_bound_19_holds_for_all_38_derived_designs",
       set(bounds.values()) == {19} and len(bounds) == 38,
       "per-design pigeonhole bounds=%s" % sorted(set(bounds.values())))
    biggest = max(len(c) for c in expand_colouring(TAB_INF))
    ck("largest_partial_parallel_class_has_exactly_12_blocks",
       3 * (cap + 1) > npoints and biggest == cap == 12,
       "points=%d so at most %d disjoint triples; largest exhibited class=%d"
       % (npoints, cap, biggest))
    ck("pigeonhole_forces_at_least_19_classes",
       lower == 19 and (lower - 1) * cap < nblocks <= lower * cap,
       "ceil(%d/%d)=%d and 18*%d=%d < %d"
       % (nblocks, cap, lower, cap, 18 * cap, nblocks))
    return lower


def check_conclusions(upper, lower):
    """chi' of every derived design, chi'(37), mcDSQS status, q'_0(3,4,38)."""
    chi = dict((x, upper[x]) for x in X if upper[x] == lower)
    ck("chromatic_index_of_every_derived_design_equals_19",
       len(chi) == 38 and set(chi.values()) == {19},
       "upper=lower=19 on %d of 38 points" % len(chi))
    v = 37
    formula = (v + 1) // 2 if v % 6 == 1 and v not in (7, 13) else None
    ck("chi_prime_37_equals_19_and_matches_the_v_plus_1_over_2_value",
       lower == 19 and max(upper.values()) == 19 and formula == 19,
       "constructed upper bound %d, pigeonhole lower bound %d, (v+1)/2=%s"
       % (max(upper.values()), lower, formula))
    ck("system_is_an_mcDSQS_38",
       len(BLOCKS) == 2109 and len(chi) == 38
       and all(chi[x] == lower for x in X),
       "all 38 derived STS(37) attain the minimum chi'(37)=%d" % lower)
    worst = max(upper.values())
    ck("q0_prime_3_4_38_equals_20",
       len(chi) == 38 and worst == lower == 19 and 1 + worst == 20,
       "1 + max_x chi'(D_x) = 1 + %d = %d" % (worst, 1 + worst))


def check_paper_arithmetic():
    """Numbers the paper prints, each recomputed from the exhibited object.

    The left side of every comparison is measured on the constructed system,
    so a corrupted object changes it.
    """
    ntrip = len(set(t for b in BLOCKS for t in combinations(b, 3)))
    prof = sorted((len(c) for c in expand_colouring(TAB_INF)), reverse=True)
    items = [
        ("blocks of a derived design = 37*36/6",
         len(DERIVED[INF]), 37 * 36 // 6),
        ("distinct triples enumerated = C(38,3)", ntrip, 38 * 37 * 36 // 6),
        ("|B| = 12*37 + 15*111", len(BLOCKS), 12 * 37 + 15 * 111),
        ("|B| = C(38,3)/C(4,3)", len(BLOCKS), (38 * 37 * 36 // 6) // 4),
        ("orbits used = 6+6+15", len(ORBITS), 6 + 6 + 15),
        ("|G| = 3*37", len(GROUP), 3 * 37),
        ("largest class = floor(37/3)", prof[0], 37 // 3),
        ("classes = ceil(222/12)", len(prof),
         -((-len(DERIVED[INF])) // (37 // 3))),
        ("profile sums to the block count", sum(prof), len(DERIVED[INF])),
        ("profile is 14*12 + 4*11 + 10", sum(prof), 14 * 12 + 4 * 11 + 10),
    ]
    wrong = ["%s -> %d != %d" % (n, got, want)
             for n, got, want in items if got != want]
    ck("printed_numbers_recomputed_from_the_object", not wrong,
       "%d printed quantities reproduced" % len(items)
       if not wrong else "; ".join(wrong))


if __name__ == "__main__":
    try:
        check_group()
        check_representatives()
        check_table_shapes()
        check_stabilisers()
        build_all()
        check_orbit_lengths_and_count()
        check_sqs_property()
        check_group_invariance_and_rotation()
        check_derived_are_sts37()
        check_colourings()
        check_replication_numbers()
        UPPER = check_all_38_derived_designs_19_colourable()
        LOWER = check_lower_bound()
        check_conclusions(UPPER, LOWER)
        check_paper_arithmetic()
    except Exception as exc:
        ck("every_stage_ran_to_completion", False,
           "aborted on malformed input: %r" % (exc,))
        finish()
    print("NOTE: the search that found this system and these two colourings "
          "is not re-run; the exhibited object and both colourings are "
          "verified exactly, and the lower bound 19 is an exact pigeonhole "
          "valid for every STS(37), so the theorem is fully established. "
          "The identity q'_0(3,4,38) = 1 + max_x chi'(D_x) is quoted from "
          "the cited literature and is applied, not proved, here.")
    finish()
