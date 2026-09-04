#!/usr/bin/env python3
# census.py -- driver for the nine-vertex census behind
# "The Maximum Euler-Sombor Index at Diameter Three, and Both Extrema at Nine Vertices".
#
# This is the computation the shipped verify.py explicitly does NOT do: it enumerates the
# cell G_{9,3} and settles EXHAUSTIVENESS and UNIQUENESS of both extrema by measurement.
#
# It runs two independent enumerations, neither of which trusts the other:
#
#   PASS 1 (isomorphism classes).  `geng -q -c 9` from nauty emits one graph6 string per
#   isomorphism class of connected 9-vertex graph.  census.c reads them, computes each exact
#   diameter with integer bitmask arithmetic, and buckets the classes by (diameter, EU edge
#   profile).  Since EU(G) depends only on the multiset of degree pairs over the edges, every
#   bucket has ONE exact EU value, evaluated here in 80-significant-digit Decimal arithmetic.
#   The number of isomorphism classes attaining an extremum is therefore read off directly,
#   which is what makes the uniqueness claim measurable and the number of ties printable.
#
#   PASS 2 (labelled, no isomorphism engine at all).  census.c sweeps ALL 2^36 labelled
#   graphs on {0,...,8} -- no canonical form, no nauty, no pruning lemma -- reports the
#   labelled count per diameter, and reports every labelled graph in the diameter-3 cell
#   whose EU comes within 1e-6 of the extrema found in pass 1, plus a count of any graph that
#   BEATS them (which must be zero).  Its connected total is checked against the exponential
#   formula for the number of connected labelled graphs on 9 nodes, computed here in exact
#   integer arithmetic.
#
# Controls are confirmed BEFORE any headline is printed.  If a control fails the program
# prints the failure and exits non-zero without printing an extremum.
#
# usage: python3 census.py <path-to-geng> <path-to-compiled-census-binary> [--skip-labelled]

import os
import subprocess
import sys
import time
from decimal import Decimal, getcontext, ROUND_DOWN, ROUND_HALF_EVEN
from itertools import permutations
from math import comb

getcontext().prec = 80

GENG = sys.argv[1] if len(sys.argv) > 1 else "/usr/local/bin/geng"
BIN = sys.argv[2] if len(sys.argv) > 2 else "./census"
SKIP_LABELLED = "--skip-labelled" in sys.argv

CHECKS = []


def ck(name, ok, extra=""):
    CHECKS.append((name, bool(ok)))
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, (" [%s]" % extra) if extra else ""))
    return bool(ok)


def note(s):
    print("NOTE " + s)


def fatal(msg):
    print("")
    print("STOP: " + msg)
    print("VERDICT: CONTROL FAILED -- no extremum is asserted by this run")
    sys.exit(1)


def sq(n):
    return Decimal(n).sqrt()


def dstr(x, k):
    return str(x.quantize(Decimal(1).scaleb(-k), rounding=ROUND_HALF_EVEN))


def dtrunc(x, k):
    return str(x.quantize(Decimal(1).scaleb(-k), rounding=ROUND_DOWN))


def eu_of_profile(prof):
    """prof: dict {(a,b): count} with a<=b.  EU = sum c*sqrt(a^2+ab+b^2), exactly."""
    tot = Decimal(0)
    for (a, b), c in prof.items():
        tot += Decimal(c) * sq(a * a + a * b + b * b)
    return tot


def parse_profile(s):
    prof = {}
    for tok in s.split():
        pair, c = tok.split(":")
        a, b = pair.split(",")
        prof[(int(a), int(b))] = int(c)
    return prof


def prof_key(prof):
    return tuple(sorted(prof.items()))


def prof_str(prof):
    return " ".join("(%d,%d)^%d" % (a, b, c) for (a, b), c in sorted(prof.items()))


# ---------------------------------------------------------------- graph helpers (n=9)
PAIRS = [(i, j) for j in range(1, 9) for i in range(0, j)]


def g6_to_adj(s):
    assert len(s) == 7 and s[0] == "H", s
    adj = [0] * 9
    for k, (i, j) in enumerate(PAIRS):
        if ((ord(s[1 + k // 6]) - 63) >> (5 - k % 6)) & 1:
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    return adj


def edges_to_adj(edges, n=9):
    adj = [0] * n
    for i, j in edges:
        adj[i] |= 1 << j
        adj[j] |= 1 << i
    return adj


def degrees(adj):
    return [bin(a).count("1") for a in adj]


def diameter(adj):
    n = len(adj)
    full = (1 << n) - 1
    nb = [adj[v] | (1 << v) for v in range(n)]
    cur = list(nb)
    for k in range(1, n + 1):
        if all(c == full for c in cur):
            return k
        nxt = []
        for v in range(n):
            acc = cur[v]
            s = cur[v]
            while s:
                u = (s & -s).bit_length() - 1
                s &= s - 1
                acc |= nb[u]
            nxt.append(acc)
        if nxt == cur:
            return -1
        cur = nxt
    return -1


def profile_of(adj):
    d = degrees(adj)
    prof = {}
    for i, j in PAIRS:
        if (adj[i] >> j) & 1:
            a, b = sorted((d[i], d[j]))
            prof[(a, b)] = prof.get((a, b), 0) + 1
    return prof


def canon_edges(adj, perm):
    """edge set of adj relabelled by perm (perm[v] is the image of v)"""
    out = set()
    for i, j in PAIRS:
        if (adj[i] >> j) & 1:
            a, b = perm[i], perm[j]
            out.add((a, b) if a < b else (b, a))
    return frozenset(out)


def edge_set(adj):
    return frozenset((i, j) for i, j in PAIRS if (adj[i] >> j) & 1)


def isomorphic(adjA, adjB):
    """brute force over all 9! relabellings, pruned by degree class"""
    dA, dB = degrees(adjA), degrees(adjB)
    if sorted(dA) != sorted(dB):
        return False
    EB = edge_set(adjB)
    targets = [[w for w in range(9) if dB[w] == dA[v]] for v in range(9)]
    for perm in permutations(range(9)):
        ok = True
        for v in range(9):
            if perm[v] not in targets[v]:
                ok = False
                break
        if ok and canon_edges(adjA, perm) == EB:
            return True
    return False


def aut_size(adj):
    d = degrees(adj)
    E = edge_set(adj)
    tot = 0
    targets = [[w for w in range(9) if d[w] == d[v]] for v in range(9)]
    for perm in permutations(range(9)):
        ok = True
        for v in range(9):
            if perm[v] not in targets[v]:
                ok = False
                break
        if ok and canon_edges(adj, perm) == E:
            tot += 1
    return tot


# ---------------------------------------------------------------- the two named objects
def build_H(a, b, n=9):
    """H(a,b): K_{n-2} on W, u joined to A (|A|=a), v joined to B (|B|=b), u !~ v"""
    W = list(range(n - 2))
    u, v = n - 2, n - 1
    edges = [(i, j) for i in W for j in W if i < j]
    for i in W[:a]:
        edges.append((i, u))
    for i in W[a:a + b]:
        edges.append((i, v))
    return edges_to_adj([(min(e), max(e)) for e in edges], n)


def build_theta(l1, l2, l3):
    """Theta(l1,l2,l3): two degree-3 hubs joined by internally disjoint paths of those lengths"""
    assert l1 + l2 + l3 == 10
    u, v = 0, 1
    nxt = 2
    edges = []
    for L in (l1, l2, l3):
        prev = u
        for _ in range(L - 1):
            edges.append((prev, nxt))
            prev = nxt
            nxt += 1
        edges.append((prev, v))
    assert nxt == 9, nxt
    return edges_to_adj([(min(e), max(e)) for e in edges], 9)


def build_double_star(a, b):
    """D(a,b): centres joined, a leaves on one, b on the other"""
    edges = [(0, 1)]
    nxt = 2
    for _ in range(a):
        edges.append((0, nxt)); nxt += 1
    for _ in range(b):
        edges.append((1, nxt)); nxt += 1
    return edges_to_adj([(min(e), max(e)) for e in edges], 9)


# ---------------------------------------------------------------- published control values
A001349_9 = 261080        # connected graphs on 9 nodes, up to isomorphism
A000055_9 = 47            # unlabelled trees on 9 nodes
PAPER_STRATA = [1, 91518, 148229, 19320, 1818, 180, 13, 1]      # paper Section 5, D=1..8
PAPER_D3_BY_SIZE = [3, 17, 76, 369, 1483, 4381, 9621, 16547, 23013, 26173, 24441,
                    18828, 12040, 6477, 2980, 1185, 412, 131, 38, 11, 3]   # m=8..28


def connected_labelled(n):
    """exact number of connected labelled graphs on n nodes (exponential formula)"""
    C = [0] * (n + 1)
    for k in range(1, n + 1):
        tot = 2 ** comb(k, 2)
        for j in range(1, k):
            tot -= comb(k - 1, j - 1) * C[j] * 2 ** comb(k - j, 2)
        C[k] = tot
    return C[n]


# ================================================================== run
print("=" * 100)
print("census of the nine-vertex diameter-three cell G_{9,3} for the Euler-Sombor index")
print("EU(G) = sum over edges xy of sqrt(d(x)^2 + d(x)d(y) + d(y)^2)")
print("exact arithmetic: every EU value below is an 80-significant-digit Decimal, never a float")
print("=" * 100)

nproc = int(subprocess.run(["nproc"], capture_output=True, text=True).stdout.strip())
try:
    gengver = subprocess.run([GENG, "-help"], capture_output=True, text=True,
                             stdin=subprocess.DEVNULL, timeout=20)
    gv = [l for l in (gengver.stdout + gengver.stderr).splitlines() if "nauty" in l.lower()]
except Exception as exc:                                   # a banner is decoration, not evidence
    gv = ["version banner unavailable: %r" % (exc,)]
note("host nproc = %d" % nproc)
note("generator  = %s  (%s)" % (GENG, gv[0].strip() if gv else "version banner not printed"))
note("kernel     = %s" % BIN)

# ---------------------------------------------------------------- PASS 1
print("")
print("=== Pass 1: every isomorphism class of connected 9-vertex graph, from geng")
print("command: %s -q -c 9 | %s classes" % (GENG, BIN))
t0 = time.time()
p1 = subprocess.Popen([GENG, "-q", "-c", "9"], stdout=subprocess.PIPE)
p2 = subprocess.Popen([BIN, "classes"], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
p1.stdout.close()
out1, _ = p2.communicate()
rc1 = p2.returncode
t_pass1 = time.time() - t0
note("pass 1 wall clock = %.1f s, exit status %d" % (t_pass1, rc1))
if rc1 != 0:
    fatal("the class enumeration kernel exited %d" % rc1)

nread = nbad = ndisc = None
buckets = []   # (diam, m, count, [reps], profile dict)
saw_end = False
for line in out1.splitlines():
    f = line.split()
    if not f:
        continue
    if f[0] == "READ":
        nread = int(f[1])
    elif f[0] == "MALFORMED":
        nbad = int(f[1])
    elif f[0] == "DISCONNECTED":
        ndisc = int(f[1])
    elif f[0] == "END":
        saw_end = True
    elif f[0] == "BUCKET":
        diam, m, cnt, nrep = int(f[1]), int(f[2]), int(f[3]), int(f[4])
        reps = f[5:5 + nrep]
        # NB the LAST bar: graph6 alphabet is chr(63..126) and so contains '|' itself
        bar = line.rindex("|")
        buckets.append((diam, m, cnt, reps, parse_profile(line[bar + 1:])))

if not saw_end:
    fatal("the class enumeration output was truncated (no END marker)")

total_classes = sum(b[2] for b in buckets)
note("geng emitted %s graph6 lines; kernel bucketed %s classes into %d (diameter, profile) buckets"
     % (nread, total_classes, len(buckets)))

ck("CONTROL_geng_emits_A001349_9_equal_261080_connected_classes", nread == A001349_9, str(nread))
ck("CONTROL_no_malformed_graph6_line_and_no_disconnected_graph_from_geng_c",
   nbad == 0 and ndisc == 0, "malformed=%s disconnected=%s" % (nbad, ndisc))
ck("CONTROL_every_class_landed_in_exactly_one_bucket", total_classes == A001349_9, str(total_classes))

strata = {}
for diam, m, cnt, reps, prof in buckets:
    strata[diam] = strata.get(diam, 0) + cnt
note("diameter strata measured: " + ", ".join("D=%d:%d" % (d, strata[d]) for d in sorted(strata)))
note("diameter strata in paper : " + ", ".join("D=%d:%d" % (i + 1, v) for i, v in enumerate(PAPER_STRATA)))
measured_strata = [strata.get(d, 0) for d in range(1, 9)]
ck("CONTROL_the_eight_measured_diameter_strata_reproduce_the_papers_Section_5_table",
   measured_strata == PAPER_STRATA, str(measured_strata))
ck("CONTROL_the_D_3_cell_holds_exactly_148229_isomorphism_classes",
   strata.get(3, 0) == 148229, str(strata.get(3, 0)))

trees = sum(cnt for diam, m, cnt, reps, prof in buckets if m == 8)
ck("CONTROL_the_m_equals_8_classes_reproduce_A000055_9_equal_47_unlabelled_trees",
   trees == A000055_9, str(trees))

d3 = [b for b in buckets if b[0] == 3]
by_size = {}
for diam, m, cnt, reps, prof in d3:
    by_size[m] = by_size.get(m, 0) + cnt
sizes = sorted(by_size)
note("D=3 sizes present: m = %d..%d (%d values)" % (sizes[0], sizes[-1], len(sizes)))
note("D=3 class count by size: " + ", ".join("%d:%d" % (m, by_size[m]) for m in sizes))
ck("CONTROL_the_D_3_stratum_has_exactly_21_sizes_m_8_to_28",
   sizes == list(range(8, 29)), "%d sizes" % len(sizes))
ck("CONTROL_the_21_per_size_class_counts_reproduce_the_papers_list",
   [by_size[m] for m in range(8, 29)] == PAPER_D3_BY_SIZE,
   str([by_size[m] for m in range(8, 29)]))
ck("CONTROL_the_per_size_counts_sum_to_the_cell_size_148229",
   sum(by_size.values()) == 148229, str(sum(by_size.values())))

if not all(ok for _, ok in CHECKS):
    fatal("a control on the enumeration itself failed; the extremes this run would print "
          "are worthless and are withheld")

# ---------------------------------------------------------------- exact evaluation
print("")
print("=== Pass 1, step 2: exact 80-digit evaluation of EU on every profile of the cell")
note("EU(G) is a function of the edge profile alone, so the %d classes of the cell reduce to "
     "%d distinct exact values to be compared" % (148229, len({prof_key(b[4]) for b in d3})))

prof_tab = {}      # profile key -> [EU, classes, reps, m]
for diam, m, cnt, reps, prof in d3:
    k = prof_key(prof)
    if k not in prof_tab:
        prof_tab[k] = [eu_of_profile(prof), 0, [], m, prof]
    prof_tab[k][1] += cnt
    prof_tab[k][2].extend(reps[:2])
ordered = sorted(prof_tab.items(), key=lambda kv: kv[1][0])
note("distinct profiles in the cell = %d, covering %d classes"
     % (len(ordered), sum(v[1] for _, v in ordered)))
ck("CONTROL_the_profile_reduction_is_lossless_classes_still_sum_to_148229",
   sum(v[1] for _, v in ordered) == 148229, str(sum(v[1] for _, v in ordered)))

gaps = [ordered[i + 1][1][0] - ordered[i][1][0] for i in range(len(ordered) - 1)]
mingap = min(gaps)
note("smallest gap between two DISTINCT profile values anywhere in the cell = %s" % dstr(mingap, 40))
ck("no_two_distinct_profiles_of_the_cell_are_within_1e_30_so_no_comparison_here_is_delicate",
   mingap > Decimal("1e-30"), "min gap = %s" % dstr(mingap, 34))

lo_prof, lo = ordered[0]
hi_prof, hi = ordered[-1]
lo2, hi2 = ordered[1][1], ordered[-2][1]

print("")
print("--- the five cheapest and the five dearest profiles of the cell")
for kv in ordered[:5] + [(None, None)] + ordered[-5:]:
    if kv[0] is None:
        print("     ...")
        continue
    v = kv[1]
    print("     EU = %s   m=%2d   classes=%6d   profile: %s" % (dstr(v[0], 12), v[3], v[1], prof_str(v[4])))

print("")
print("--- MINIMUM over G_{9,3}")
note("EU_min to 60 digits = %s" % dstr(lo[0], 57))
note("classes attaining it = %d   (representative graph6: %s)" % (lo[1], " ".join(lo[2][:4])))
note("runner-up value      = %s   (margin %s)" % (dstr(lo2[0], 12), dstr(lo2[0] - lo[0], 12)))
min_cf = 8 * sq(3) + 6 * sq(19)
ck("the_measured_minimum_equals_8sqrt3_plus_6sqrt19_to_1e_60",
   abs(lo[0] - min_cf) < Decimal("1e-60"), "difference = %s" % (lo[0] - min_cf))
ck("the_measured_minimum_truncates_to_the_papers_printed_40_009800121795059",
   dtrunc(lo[0], 15) == "40.009800121795059", dtrunc(lo[0], 15))
ck("THE_MINIMUM_IS_ATTAINED_BY_EXACTLY_ONE_ISOMORPHISM_CLASS_zero_ties",
   lo[1] == 1, "classes at the minimum = %d, ties = %d" % (lo[1], lo[1] - 1))

theta334 = build_theta(3, 3, 4)
argmin = g6_to_adj(lo[2][0])
note("argmin: degrees %s, m=%d, diameter %d"
     % (sorted(degrees(argmin)), sum(degrees(argmin)) // 2, diameter(argmin)))
ck("the_argmin_returned_by_the_census_is_isomorphic_to_the_theta_graph_Theta_3_3_4",
   isomorphic(argmin, theta334), "graph6 %s" % lo[2][0])
ck("the_argmin_has_diameter_exactly_3_and_is_bicyclic_m_equals_10",
   diameter(argmin) == 3 and sum(degrees(argmin)) // 2 == 10)
aut_min = aut_size(argmin)
note("|Aut(argmin)| = %d, so the argmin accounts for 9!/|Aut| = %d labelled graphs"
     % (aut_min, 362880 // aut_min))

print("")
print("--- MAXIMUM over G_{9,3}")
note("EU_max to 60 digits = %s" % dstr(hi[0], 57))
note("classes attaining it = %d   (representative graph6: %s)" % (hi[1], " ".join(hi[2][:4])))
note("runner-up value      = %s   (margin %s)" % (dstr(hi2[0], 12), dstr(hi[0] - hi2[0], 12)))
max_cf = 147 * sq(3) + sq(57) + 6 * sq(127)
ck("the_measured_maximum_equals_147sqrt3_plus_sqrt57_plus_6sqrt127_to_1e_60",
   abs(hi[0] - max_cf) < Decimal("1e-60"), "difference = %s" % (hi[0] - max_cf))
ck("the_measured_maximum_truncates_to_the_papers_printed_329_777869165403581",
   dtrunc(hi[0], 15) == "329.777869165403581", dtrunc(hi[0], 15))
ck("THE_MAXIMUM_IS_ATTAINED_BY_EXACTLY_ONE_ISOMORPHISM_CLASS_zero_ties",
   hi[1] == 1, "classes at the maximum = %d, ties = %d" % (hi[1], hi[1] - 1))

H16 = build_H(1, 6)
argmax = g6_to_adj(hi[2][0])
note("argmax: degrees %s, m=%d, diameter %d"
     % (sorted(degrees(argmax)), sum(degrees(argmax)) // 2, diameter(argmax)))
ck("the_argmax_returned_by_the_census_is_isomorphic_to_H_1_6",
   isomorphic(argmax, H16), "graph6 %s" % hi[2][0])
ck("the_argmax_has_diameter_exactly_3_and_m_equals_28",
   diameter(argmax) == 3 and sum(degrees(argmax)) // 2 == 28)
aut_max = aut_size(argmax)
note("|Aut(argmax)| = %d, so the argmax accounts for 9!/|Aut| = %d labelled graphs"
     % (aut_max, 362880 // aut_max))
ck("the_runner_up_to_the_maximum_is_H_2_5_at_147sqrt3_plus_2sqrt67_plus_5sqrt109",
   abs(hi2[0] - (147 * sq(3) + 2 * sq(67) + 5 * sq(109))) < Decimal("1e-60"),
   dstr(hi2[0], 10))
ck("the_margin_of_the_maximum_over_the_runner_up_is_the_papers_6_5941623645",
   dstr(hi[0] - hi2[0], 10) == "6.5941623645", dstr(hi[0] - hi2[0], 10))

# ---------------------------------------------------------------- published values inside the cell
print("")
print("=== Pass 1, step 3: four more values the census must reproduce, three of them published")
by_m = {}
for k, v in prof_tab.items():
    by_m.setdefault(v[3], []).append(v)
min_m8 = min(v[0] for v in by_m[8])
min_m9 = min(v[0] for v in by_m[9])
max_m9 = max(v[0] for v in by_m[9])
note("min over the trees of the cell (m=8)      = %s" % dstr(min_m8, 10))
note("min over the unicyclic members (m=9)      = %s" % dstr(min_m9, 10))
note("max over the unicyclic members (m=9)      = %s" % dstr(max_m9, 10))
ck("the_best_tree_in_the_cell_is_D_3_4_at_sqrt61_plus_3sqrt21_plus_4sqrt31_43_8290342121",
   abs(min_m8 - (sq(61) + 3 * sq(21) + 4 * sq(31))) < Decimal("1e-60"), dstr(min_m8, 10))
ck("the_unicyclic_runner_up_is_4sqrt21_plus_8sqrt3_plus_4sqrt7_equal_42_7697144846",
   abs(min_m9 - (4 * sq(21) + 8 * sq(3) + 4 * sq(7))) < Decimal("1e-60"), dstr(min_m9, 10))
ck("PUBLISHED_the_unicyclic_maximum_of_the_cell_is_sharp_for_Sekar_et_al_Theorem_3_2_at_n_9",
   abs(max_m9 - (5 * sq(57) + sq(79) + sq(67) + sq(13) + sq(19))) < Decimal("1e-60"),
   dstr(max_m9, 10))
ck("the_strict_order_of_Theorem_A_bicyclic_below_unicyclic_below_best_tree",
   lo[0] < min_m9 < min_m8,
   "gaps %s and %s" % (dstr(min_m9 - lo[0], 10), dstr(min_m8 - min_m9, 10)))

d4 = [b for b in buckets if b[0] == 4]
d4_vals = {}
for diam, m, cnt, reps, prof in d4:
    k = prof_key(prof)
    if k not in d4_vals:
        d4_vals[k] = [eu_of_profile(prof), 0, reps]
    d4_vals[k][1] += cnt
d4_min = min(d4_vals.values(), key=lambda v: v[0])
note("min over G_{9,4} = %s attained by %d class(es), representative %s"
     % (dstr(d4_min[0], 10), d4_min[1], d4_min[2][0]))
ck("PUBLISHED_the_minimum_of_the_D_4_cell_is_Kizilirmak_Lemma_2_EU_C_9_equal_18sqrt3",
   abs(d4_min[0] - 18 * sq(3)) < Decimal("1e-60"), dstr(d4_min[0], 10))
ck("that_D_4_minimiser_is_the_9_cycle_itself",
   d4_min[1] == 1 and sorted(degrees(g6_to_adj(d4_min[2][0]))) == [2] * 9)

print("")
print("=== Pass 1, step 4: the sweep can report non-uniqueness -- the paper's two genuine ties")
for D in (5, 6, 7, 8):
    vals = {}
    for diam, m, cnt, reps, prof in buckets:
        if diam != D:
            continue
        k = prof_key(prof)
        if k not in vals:
            vals[k] = [eu_of_profile(prof), 0]
        vals[k][1] += cnt
    best = min(vals.values(), key=lambda v: v[0])
    note("min over G_{9,%d} = %s attained by %d isomorphism class(es)"
         % (D, dstr(best[0], 10), best[1]))
    if D in (6, 7):
        ck("the_reported_minimum_at_n_9_D_%d_is_a_genuine_TIE_of_two_classes_as_the_paper_states" % D,
           best[1] == 2, "classes = %d" % best[1])

# ---------------------------------------------------------------- PASS 2
if not SKIP_LABELLED:
    print("")
    print("=== Pass 2: independent sweep of all 2^36 labelled graphs on {0..8} -- no geng, "
          "no canonical form, no isomorphism test")
    emin = float(lo[0])
    emax = float(hi[0])
    win = 1e-6
    print("command: %s labelled %.17g %.17g %g" % (BIN, emin, emax, win))
    t0 = time.time()
    r = subprocess.run([BIN, "labelled", repr(emin), repr(emax), repr(win)],
                       capture_output=True, text=True)
    t_pass2 = time.time() - t0
    note("pass 2 wall clock = %.1f s, exit status %d" % (t_pass2, r.returncode))
    if r.returncode != 0:
        fatal("the labelled sweep exited %d" % r.returncode)
    L = {}
    lab_diam = {}
    lo_profs, hi_profs = [], []
    for line in r.stdout.splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] == "LABELLED_DIAM":
            lab_diam[int(f[1])] = int(f[2])
        elif f[0] in ("SWEPT", "DISCONNECTED", "BELOW_THRESHOLD", "ABOVE_THRESHOLD",
                      "WINDOW_LO_COUNT", "WINDOW_HI_COUNT", "THREADS"):
            L[f[0]] = int(f[1])
        elif f[0] in ("DOUBLE_MIN", "DOUBLE_MAX"):
            L[f[0]] = float(f[1])
        elif f[0] == "WINDOW_LO_PROFILE":
            lo_profs.append((int(f[1]), parse_profile(line[line.rindex("|") + 1:])))
        elif f[0] == "WINDOW_HI_PROFILE":
            hi_profs.append((int(f[1]), parse_profile(line[line.rindex("|") + 1:])))

    conn_lab = sum(lab_diam.values())
    exact_conn = connected_labelled(9)
    note("labelled graphs swept   = %d (= 2^36)" % L["SWEPT"])
    note("labelled disconnected   = %d" % L["DISCONNECTED"])
    note("labelled connected      = %d" % conn_lab)
    note("exponential formula     = %d connected labelled graphs on 9 nodes" % exact_conn)
    note("labelled by diameter    = " + ", ".join("D=%d:%d" % (d, lab_diam[d]) for d in sorted(lab_diam)))
    ck("CONTROL_the_labelled_sweep_saw_all_2_to_the_36_graphs",
       L["SWEPT"] == 2 ** 36 and L["SWEPT"] == L["DISCONNECTED"] + conn_lab)
    ck("CONTROL_the_labelled_connected_count_equals_the_exponential_formula_value",
       conn_lab == exact_conn, str(conn_lab))
    # the two passes are independent; the only arithmetic that must link them is the orbit
    # sandwich  classes(D) <= labelled(D) <= 9! * classes(D)  and the two strata whose
    # automorphism group is forced: D=1 is K_9 alone, D=8 is the path P_9 alone (|Aut|=2).
    sandwich = all(measured_strata[D - 1] <= lab_diam[D] <= 362880 * measured_strata[D - 1]
                   for D in range(1, 9))
    ck("CROSS_PASS_every_stratum_obeys_classes_at_most_labelled_at_most_9_factorial_times_classes",
       sandwich, "D=3: %d classes, %d labelled" % (measured_strata[2], lab_diam[3]))
    ck("CROSS_PASS_the_D_1_stratum_is_the_single_labelled_graph_K_9",
       lab_diam[1] == 1, str(lab_diam[1]))
    ck("CROSS_PASS_the_D_8_stratum_is_9_factorial_over_2_equal_181440_labellings_of_the_path_P_9",
       lab_diam[8] == 181440 and measured_strata[7] == 1, str(lab_diam[8]))
    ck("NOTHING_IN_THE_CELL_BEATS_THE_MINIMUM_found_by_the_class_pass",
       L["BELOW_THRESHOLD"] == 0, "labelled graphs below EU_min-1e-6 = %d" % L["BELOW_THRESHOLD"])
    ck("NOTHING_IN_THE_CELL_BEATS_THE_MAXIMUM_found_by_the_class_pass",
       L["ABOVE_THRESHOLD"] == 0, "labelled graphs above EU_max+1e-6 = %d" % L["ABOVE_THRESHOLD"])
    ck("the_labelled_sweeps_own_double_extrema_agree_with_the_exact_class_pass_values",
       abs(L["DOUBLE_MIN"] - emin) < 1e-9 and abs(L["DOUBLE_MAX"] - emax) < 1e-9,
       "min %.12f max %.12f" % (L["DOUBLE_MIN"], L["DOUBLE_MAX"]))
    note("labelled graphs within 1e-6 of the minimum = %d, in %d distinct profile(s)"
         % (L["WINDOW_LO_COUNT"], len(lo_profs)))
    note("labelled graphs within 1e-6 of the maximum = %d, in %d distinct profile(s)"
         % (L["WINDOW_HI_COUNT"], len(hi_profs)))
    ck("the_minimum_window_holds_exactly_one_profile_and_it_is_the_class_passs_argmin_profile",
       len(lo_profs) == 1 and prof_key(lo_profs[0][1]) == lo_prof,
       prof_str(lo_profs[0][1]) if lo_profs else "empty")
    ck("the_maximum_window_holds_exactly_one_profile_and_it_is_the_class_passs_argmax_profile",
       len(hi_profs) == 1 and prof_key(hi_profs[0][1]) == hi_prof,
       prof_str(hi_profs[0][1]) if hi_profs else "empty")
    ck("UNIQUENESS_CONFIRMED_LABELLED_the_minimum_window_holds_exactly_9_factorial_over_Aut_graphs",
       L["WINDOW_LO_COUNT"] == 362880 // aut_min,
       "%d = 362880/%d, so exactly one isomorphism class attains the minimum"
       % (L["WINDOW_LO_COUNT"], aut_min))
    ck("UNIQUENESS_CONFIRMED_LABELLED_the_maximum_window_holds_exactly_9_factorial_over_Aut_graphs",
       L["WINDOW_HI_COUNT"] == 362880 // aut_max,
       "%d = 362880/%d, so exactly one isomorphism class attains the maximum"
       % (L["WINDOW_HI_COUNT"], aut_max))

# ---------------------------------------------------------------- headline
npass = sum(1 for _, ok in CHECKS if ok)
nfail = len(CHECKS) - npass
print("")
print("=" * 100)
print("HEADLINE (measured, not transcribed):")
print("  the cell G_{9,3} has exactly %d isomorphism classes" % strata[3])
print("  max EU over G_{9,3} = %s" % dstr(hi[0], 30))
print("                      = 147*sqrt3 + sqrt57 + 6*sqrt127, attained by %d class: H(1,6)" % hi[1])
print("  min EU over G_{9,3} = %s" % dstr(lo[0], 30))
print("                      = 8*sqrt3 + 6*sqrt19, attained by %d class: Theta(3,3,4)" % lo[1])
print("  ties at the maximum: %d      ties at the minimum: %d" % (hi[1] - 1, lo[1] - 1))
print("=" * 100)
print("VERDICT: %d checks, %d PASS, %d FAIL -- %s" %
      (len(CHECKS), npass, nfail,
       "BOTH EXTREMA OF G_{9,3} ARE EXHAUSTIVELY CONFIRMED AND UNIQUE"
       if nfail == 0 else "SOMETHING FAILED; SEE THE FAIL LINES ABOVE"))
sys.exit(0 if nfail == 0 else 1)
