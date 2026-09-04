#!/usr/bin/env python3
"""
verify.py -- checks the computational claims of

    "A 25-Vertex Weak Tensor Product of alpha-Free Trees With an alpha-Valuation"

Python 3.9+, STANDARD LIBRARY ONLY.  No third-party package, no external data
file, no network.  All arithmetic is exact integer arithmetic; no floats are
used anywhere, so no decision in this program depends on rounding.

Everything it consumes is printed in the paper: the two factor trees, the two
near alpha-valuations gamma (eq. 1) and delta (eq. 2), the 25 labels of the
witness b of Section 3, and the 32 labels of the witness of Theorem 8.  The
graphs themselves are NOT read in: they are rebuilt from the weak tensor
product definition quoted in Section 1 by an adjacency predicate over all pairs
of vertices, so the program never trusts a stored edge list.

The program prints one `PASS <name> [detail]` line per check and closes with

    VERDICT: ALL <n> CHECKS PASS

exiting 0 if and only if every check passed.

A checker that cannot say NO is worthless, so this one is asked questions whose
true answer is NO and must give it: `gamma-is-NOT-alpha`, `delta-is-NOT-alpha`,
`sigma-is-NOT-alpha`, `b-corruption-refused`, and
`control-cycles-alpha-iff-0-mod-4`, the last of which reproduces a published
if-and-only-if in both directions by exhaustion.
"""

import itertools
import sys
from collections import Counter

# --------------------------------------------------------------------------
# check harness
# --------------------------------------------------------------------------

_N = 0
_FAILED = []


def check(ok, name, detail=""):
    global _N
    _N += 1
    if ok:
        print("PASS %-38s %s" % (name, detail))
    else:
        print("FAIL %-38s %s" % (name, detail))
        _FAILED.append(name)


# --------------------------------------------------------------------------
# generic labelling predicates, straight from the definitions quoted in S1
# --------------------------------------------------------------------------

def is_beta(lab, E):
    """beta-valuation (graceful labelling): injection into {0..q} whose edge
    differences are exactly {1,...,q}, q = |E|."""
    q = len(E)
    vals = list(lab.values())
    if len(set(vals)) != len(vals):
        return False
    if min(vals) < 0 or max(vals) > q:
        return False
    return sorted(abs(lab[u] - lab[v]) for u, v in E) == list(range(1, q + 1))


def is_near_alpha(lab, E, small, large):
    """every edge joins `small` to `large` with the smaller label on the
    `small` side."""
    small, large = set(small), set(large)
    for u, v in E:
        if u in small and v in large:
            a, b = u, v
        elif v in small and u in large:
            a, b = v, u
        else:
            return False
        if not lab[a] < lab[b]:
            return False
    return True


def alpha_boundaries(lab, E):
    """every x in {0..q} that straddles every edge (the alpha condition of
    S1).  Returns the full list, so 'no alpha-valuation' is an empty list."""
    q = len(E)
    out = []
    for x in range(q + 1):
        if all((lab[u] <= x < lab[v]) or (lab[v] <= x < lab[u]) for u, v in E):
            out.append(x)
    return out


def degrees(V, E):
    d = {v: 0 for v in V}
    for u, v in E:
        d[u] += 1
        d[v] += 1
    return d


def connected(V, E):
    adj = {v: set() for v in V}
    for u, v in E:
        adj[u].add(v)
        adj[v].add(u)
    seen = {V[0]}
    stack = [V[0]]
    while stack:
        z = stack.pop()
        for w in adj[z]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == len(V)


def weak_tensor_product(EG, Gsmall, Glarge, EH, Hsmall, Hlarge):
    """The weak tensor product of S1, verbatim: vertices
    (Gsmall x Hsmall) u (Glarge x Hlarge), and (x,y) ~ (u,v) iff
    {x,u} in E(G) and {y,v} in E(H).  Built by testing the predicate on every
    ordered pair, never from a stored edge list."""
    esG = set(map(frozenset, EG))
    esH = set(map(frozenset, EH))
    P = [(x, y) for x in Gsmall for y in Hsmall]
    Q = [(u, v) for u in Glarge for v in Hlarge]
    V = P + Q
    E = []
    for i, a in enumerate(V):
        for j, b in enumerate(V):
            if i < j and frozenset((a[0], b[0])) in esG \
                     and frozenset((a[1], b[1])) in esH:
                E.append((a, b))
    return V, P, Q, E


# --------------------------------------------------------------------------
# the two factor trees, as printed in S2 and S5
# --------------------------------------------------------------------------

# S_{3,2}: centre c, subdivision vertices m1,m2,m3, leaves l1,l2,l3
V32 = ["c", "l1", "l2", "l3", "m1", "m2", "m3"]
E32 = [("c", "m1"), ("c", "m2"), ("c", "m3"),
       ("l1", "m1"), ("l2", "m2"), ("l3", "m3")]
S32_SMALL = ["c", "l1", "l2", "l3"]
S32_LARGE = ["m1", "m2", "m3"]

# gamma, equation (1): gamma(c)=0, (gamma(m_i),gamma(l_i)) = (5,4),(3,1),(6,2)
GAMMA = {"c": 0, "m1": 5, "l1": 4, "m2": 3, "l2": 1, "m3": 6, "l3": 2}

# S_{4,2}: centre c, subdivision vertices a0..a3, leaves b0..b3
V42 = ["c"] + ["b%d" % j for j in range(4)] + ["a%d" % j for j in range(4)]
E42 = [("c", "a%d" % j) for j in range(4)] + \
      [("b%d" % j, "a%d" % j) for j in range(4)]
S42_SMALL = ["c"] + ["b%d" % j for j in range(4)]
S42_LARGE = ["a%d" % j for j in range(4)]

# delta, equation (2): delta(c)=0, (delta(a_j),delta(b_j)) = (8,4),(7,2),(6,5),(3,1)
DELTA = {"c": 0,
         "a0": 8, "b0": 4,
         "a1": 7, "b1": 2,
         "a2": 6, "b2": 5,
         "a3": 3, "b3": 1}

# --------------------------------------------------------------------------
# PART 1 -- the factor S_{3,2}, gamma, and alpha-freeness (paper S2)
# --------------------------------------------------------------------------

print("=" * 78)
print("PART 1  the factor S_{3,2} of Section 2, its near alpha-valuation gamma,")
print("        and Proposition 3 (S_{3,2} has no alpha-valuation)")
print("=" * 78)

d32 = degrees(V32, E32)
check(len(V32) == 7 and len(E32) == 6 and connected(V32, E32)
      and len(V32) == len(E32) + 1
      and sorted(d32.values()) == [1, 1, 1, 2, 2, 2, 3]
      and all((u in S32_SMALL) != (v in S32_SMALL) for u, v in E32),
      "s32-structure",
      "|V|=7 |E|=6 tree connected degrees 3,2,2,2,1,1,1; "
      "bipartition {c,l1,l2,l3}|{m1,m2,m3}")

check(is_beta(GAMMA, E32), "gamma-is-beta-valuation",
      "differences %s == {1..6}"
      % sorted(abs(GAMMA[u] - GAMMA[v]) for u, v in E32))

check(is_near_alpha(GAMMA, E32, S32_SMALL, S32_LARGE),
      "gamma-is-near-alpha",
      "V_small labels %s below V_large labels %s pointwise on edges"
      % ([GAMMA[v] for v in S32_SMALL], [GAMMA[v] for v in S32_LARGE]))

# ANTI-CONTROL: the true answer here is NO.
check(alpha_boundaries(GAMMA, E32) == [],
      "gamma-is-NOT-alpha (anti-control)",
      "no x in {0..6} straddles every edge; max gamma(V_small)=4 > "
      "3=min gamma(V_large)")

# Proposition 3, Case A and Case B, as arithmetic over every admissible b(c).
caseA = {bc: 2 * (4 + 5 + 6) - (3 * bc + (0 + 1 + 2 + 3 - bc)) for bc in range(4)}
check(all(s % 2 == 0 for s in caseA.values()) and 21 not in caseA.values(),
      "s32-prop3-caseA",
      "sum = 24-2b(c) in %s, all even, never 21"
      % [caseA[bc] for bc in sorted(caseA)])

caseB = {bc: (3 * bc + (3 + 4 + 5 + 6 - bc)) - 2 * (0 + 1 + 2)
         for bc in range(3, 7)}
check(all(s % 2 == 0 for s in caseB.values()) and 21 not in caseB.values(),
      "s32-prop3-caseB",
      "sum = 2b(c)+12 in %s, all even, never 21"
      % [caseB[bc] for bc in sorted(caseB)])


def alpha_free_by_structure(V, E, small, large):
    """Complete search for alpha-valuations of a tree, using Fact 2(ii): the
    low class is one of the two bipartition classes and receives exactly the
    initial segment of labels.  Returns (n_candidates, list of alpha found)."""
    n = len(V)
    total = 0
    found = []
    for low, high in ((small, large), (large, small)):
        lo_labels = list(range(len(low)))
        hi_labels = list(range(len(low), n))
        for pl in itertools.permutations(lo_labels):
            base = dict(zip(low, pl))
            for ph in itertools.permutations(hi_labels):
                lab = dict(base)
                lab.update(zip(high, ph))
                total += 1
                if is_beta(lab, E) and alpha_boundaries(lab, E):
                    found.append(dict(lab))
    return total, found


tot, found = alpha_free_by_structure(V32, E32, S32_SMALL, S32_LARGE)
check(tot == 288 and found == [], "s32-alpha-free-complete",
      "%d candidates permitted by Fact 2(ii) examined, %d alpha-valuations"
      % (tot, len(found)))


def census_all_injections(V, E, small, large):
    """Brute force over ALL bijections V -> {0..q} (q=|E|, |V|=q+1 for a
    tree): counts beta, near-alpha (either orientation) and alpha."""
    q = len(E)
    nb = nn = na = 0
    for p in itertools.permutations(range(q + 1)):
        lab = dict(zip(V, p))
        if not is_beta(lab, E):
            continue
        nb += 1
        if is_near_alpha(lab, E, small, large) \
           or is_near_alpha(lab, E, large, small):
            nn += 1
        if alpha_boundaries(lab, E):
            na += 1
    return nb, nn, na


nb, nn, na = census_all_injections(V32, E32, S32_SMALL, S32_LARGE)
check((nb, nn, na) == (60, 24, 0), "s32-census-all-5040-injections",
      "7!=5040 injections -> %d beta, %d near-alpha, %d alpha" % (nb, nn, na))

# --------------------------------------------------------------------------
# PART 2 -- the product K and the witness b (paper S3)
# --------------------------------------------------------------------------

print()
print("=" * 78)
print("PART 2  K = S_{3,2} wtp S_{3,2} rebuilt from the definition of Section 1,")
print("        and the alpha-valuation b printed in Section 3")
print("=" * 78)

VK, PK, QK, EK = weak_tensor_product(E32, S32_SMALL, S32_LARGE,
                                     E32, S32_SMALL, S32_LARGE)

# The intrinsic description of S3 ("a hub h of degree 9 joined to all nine
# vertices of a 3x3 grid ...") built independently, then matched to K.
VI = [("h",)] + [("R", i) for i in (1, 2, 3)] + [("C", j) for j in (1, 2, 3)] \
     + [("P", i, j) for i in (1, 2, 3) for j in (1, 2, 3)] \
     + [("Q", i, j) for i in (1, 2, 3) for j in (1, 2, 3)]
EI = []
for i in (1, 2, 3):
    for j in (1, 2, 3):
        EI.append((("h",), ("Q", i, j)))
        EI.append((("R", i), ("Q", i, j)))
        EI.append((("C", j), ("Q", i, j)))
        EI.append((("P", i, j), ("Q", i, j)))
ISO = {("h",): ("c", "c")}
for i in (1, 2, 3):
    ISO[("R", i)] = ("l%d" % i, "c")
    ISO[("C", i)] = ("c", "l%d" % i)
    for j in (1, 2, 3):
        ISO[("P", i, j)] = ("l%d" % i, "l%d" % j)
        ISO[("Q", i, j)] = ("m%d" % i, "m%d" % j)
check(len(VI) == len(VK)
      and set(ISO[v] for v in VI) == set(VK)
      and set(frozenset((ISO[u], ISO[v])) for u, v in EI)
          == set(frozenset(e) for e in EK),
      "K-intrinsic-equals-product",
      "the hub/grid/row/column/pendant graph of S3 maps edge-for-edge onto "
      "S_{3,2} wtp S_{3,2} under h->cc, R_i->l_ic, C_j->cl_j, "
      "P_ij->l_il_j, Q_ij->m_im_j")

check((len(VK), len(PK), len(QK), len(EK)) == (25, 16, 9, 36),
      "K-order-and-size",
      "|V(K)|=%d |P|=%d |Q|=%d |E(K)|=%d" % (len(VK), len(PK), len(QK), len(EK)))

check(connected(VK, EK), "K-connected",
      "so K has a unique bipartition, which is therefore P|Q")

check(all((u in PK) != (v in PK) for u, v in EK),
      "K-every-edge-runs-P-to-Q", "all 36 edges join P to Q")

dK = degrees(VK, EK)
check(sorted(Counter(dK.values()).items()) == [(1, 9), (3, 6), (4, 9), (9, 1)],
      "K-degree-sequence",
      "9^1 4^9 3^6 1^9 as stated in S3; hub cc has degree %d" % dK[("c", "c")])

# the witness b, transcribed from the two displays of S3
B = {("c", "l2"): 0, ("l1", "l2"): 1, ("l1", "c"): 2, ("c", "l1"): 3,
     ("l2", "l3"): 4, ("l2", "l2"): 5, ("l3", "l1"): 7, ("l3", "c"): 8,
     ("l1", "l1"): 9, ("c", "l3"): 15, ("l2", "c"): 17, ("l3", "l3"): 18,
     ("l3", "l2"): 20, ("l2", "l1"): 21, ("l1", "l3"): 23, ("c", "c"): 27,
     ("m1", "m1"): 28, ("m3", "m3"): 29, ("m3", "m2"): 30, ("m3", "m1"): 31,
     ("m2", "m2"): 32, ("m2", "m3"): 33, ("m2", "m1"): 34, ("m1", "m3"): 35,
     ("m1", "m2"): 36}

check(set(B) == set(VK), "b-domain-is-V(K)",
      "%d labels, one per vertex of the rebuilt K" % len(B))

unused = sorted(set(range(37)) - set(B.values()))
check(len(set(B.values())) == 25 and min(B.values()) == 0
      and max(B.values()) == 36
      and unused == [6, 10, 11, 12, 13, 14, 16, 19, 22, 24, 25, 26],
      "b-injective-in-0-to-36",
      "25 distinct labels in [0,36]; the 12 unused are %s" % unused)

diffs = sorted(abs(B[u] - B[v]) for u, v in EK)
check(diffs == list(range(1, 37)), "b-is-beta-valuation",
      "the 36 edge differences are exactly {1,...,36}")

# the four printed difference blocks of S3, reproduced exactly
Qlab = {(i, j): B[("m%d" % i, "m%d" % j)] for i in (1, 2, 3) for j in (1, 2, 3)}
blocks = {
    "HUB": [abs(B[("c", "c")] - Qlab[(i, j)])
            for i in (1, 2, 3) for j in (1, 2, 3)],
    "ROWS": [abs(B[("l%d" % i, "c")] - Qlab[(i, j)])
             for i in (1, 2, 3) for j in (1, 2, 3)],
    "COLS": [abs(B[("c", "l%d" % j)] - Qlab[(i, j)])
             for j in (1, 2, 3) for i in (1, 2, 3)],
    "PENDS": [abs(B[("l%d" % i, "l%d" % j)] - Qlab[(i, j)])
              for i in (1, 2, 3) for j in (1, 2, 3)],
}
printed = {
    "HUB": [1, 9, 8, 7, 5, 6, 4, 3, 2],
    "ROWS": [26, 34, 33, 17, 15, 16, 23, 22, 21],
    "COLS": [25, 31, 28, 36, 32, 30, 20, 18, 14],
    "PENDS": [19, 35, 12, 13, 27, 29, 24, 10, 11],
}
check(blocks == printed, "b-printed-difference-blocks",
      "all four blocks of S3 reproduce value-for-value in the printed order")
check(sorted(sum(printed.values(), [])) == list(range(1, 37)),
      "b-printed-blocks-tile-1-to-36",
      "the 36 printed values, sorted, are 1,2,...,36 with no gap or repeat")

bd = alpha_boundaries(B, EK)
maxP, minQ = max(B[v] for v in PK), min(B[v] for v in QK)
check(bd == [27] and maxP == 27 and minQ == 28,
      "b-is-alpha-boundary-27",
      "the only straddling x in {0..36} is 27; max_P=%d < %d=min_Q"
      % (maxP, minQ))

wid = sum(dK[v] * B[v] for v in QK) - sum(dK[v] * B[v] for v in PK)
check(wid == 36 * 37 // 2 == 666, "b-weighted-degree-identity",
      "sum_Q deg.b - sum_P deg.b = %d = 1+2+...+36" % wid)

e36 = [(u, v) for u, v in EK if abs(B[u] - B[v]) == 36]
check(e36 == [(("c", "l2"), ("m1", "m2"))],
      "b-difference-36-spot-check",
      "difference 36 is carried by cl_2=0 ~ m_1m_2=36, and those two ARE "
      "adjacent under the definition")

check(B[("c", "l1")] - B[("c", "l2")] == 3
      and B[("l1", "l1")] - B[("l1", "l2")] == 8,
      "b-is-not-separable",
      "b(cl_1)-b(cl_2)=3 but b(l_1l_1)-b(l_1l_2)=8; a separable "
      "phi(x,y)=p(x)+q(y) would force these equal (S4)")

# ANTI-CONTROL: corrupt the witness and demand a refusal.
Bbad = dict(B)
Bbad[("c", "l2")], Bbad[("l1", "l2")] = B[("l1", "l2")], B[("c", "l2")]
Bbad2 = dict(B)
Bbad2[("c", "l2")] = B[("l1", "l2")]          # destroys injectivity
check(not is_beta(Bbad, EK) and not is_beta(Bbad2, EK),
      "b-corruption-refused (anti-control)",
      "swapping b(cl_2) with b(l_1l_2) is not graceful, and a "
      "non-injective labelling is refused outright")

# --------------------------------------------------------------------------
# PART 3 -- the obstruction of paper S4: sigma can never witness
# --------------------------------------------------------------------------

print()
print("=" * 78)
print("PART 3  Lemma 5 / Corollary 6: the product labelling sigma of the")
print("        literature is graceful on K but has NO boundary")
print("=" * 78)


def sigma_labelling(EG, Gs, Gl, gam, EH, Hs, Hl, delt):
    """sigma(a,c) = m.delta(c)+gamma(a) on Gs x Hs,
       sigma(b,d) = m.(delta(d)-1)+gamma(b) on Gl x Hl,  m = |E(G)|."""
    m = len(EG)
    V, P, Q, E = weak_tensor_product(EG, Gs, Gl, EH, Hs, Hl)
    lab = {}
    for (x, y) in P:
        lab[(x, y)] = m * delt[y] + gam[x]
    for (u, v) in Q:
        lab[(u, v)] = m * (delt[v] - 1) + gam[u]
    return V, P, Q, E, lab


_, sP, sQ, sE, SIG = sigma_labelling(E32, S32_SMALL, S32_LARGE, GAMMA,
                                     E32, S32_SMALL, S32_LARGE, GAMMA)
check(sE == EK and is_beta(SIG, sE), "sigma-is-beta-valuation",
      "sigma is graceful on the same 36 edges of K")

sm, sn = max(SIG[v] for v in sP), min(SIG[v] for v in sQ)
check(alpha_boundaries(SIG, sE) == [] and sm == 28 and sn == 15,
      "sigma-is-NOT-alpha",
      "max_P sigma = %d > %d = min_Q sigma, so no x straddles: the published "
      "construction cannot witness the counterexample" % (sm, sn))

m_, gA, gB = len(E32), max(GAMMA[v] for v in S32_SMALL), \
             min(GAMMA[v] for v in S32_LARGE)
hC, hD = gA, gB       # here delta = gamma
lhs, rhs = m_ * (hC - hD + 1), gB - gA
check((lhs < rhs) is False and lhs == 12 and rhs == -1,
      "lemma5-criterion-agrees",
      "m(h_C-h_D+1) = %d and g_B-g_A = %d, so the criterion %d < %d is FALSE, "
      "matching the previous check" % (lhs, rhs, lhs, rhs))

# POSITIVE CONTROL for Lemma 5 / Corollary 6: two factors that DO have
# alpha-valuations must give an alpha sigma (Snevily's theorem).
P3V = ["u0", "u1", "u2"]
P3E = [("u0", "u1"), ("u1", "u2")]
P3S, P3L = ["u0", "u2"], ["u1"]
P3A = {"u0": 0, "u2": 1, "u1": 2}
check(is_beta(P3A, P3E) and alpha_boundaries(P3A, P3E) == [1],
      "control-P3-has-alpha", "P_3 labelled 0,2,1 is alpha with boundary 1")
_, cP, cQ, cE, CSIG = sigma_labelling(P3E, P3S, P3L, P3A,
                                      P3E, P3S, P3L, P3A)
cA = alpha_boundaries(CSIG, cE)
check(is_beta(CSIG, cE) and cA,
      "lemma5-recovers-snevily",
      "both factors alpha => sigma is an alpha-valuation of the product "
      "(|V|=%d, |E|=%d, boundary %s), which is Snevily's theorem"
      % (len(cP) + len(cQ), len(cE), cA))

# --------------------------------------------------------------------------
# PART 4 -- two DISTINCT graphs (paper S5)
# --------------------------------------------------------------------------

print()
print("=" * 78)
print("PART 4  the non-isomorphic pair S_{3,2}, S_{4,2}: Proposition 7 and the")
print("        32-vertex witness of Theorem 8")
print("=" * 78)

d42 = degrees(V42, E42)
check(len(V42) == 9 and len(E42) == 8 and connected(V42, E42)
      and len(V42) == len(E42) + 1
      and sorted(d42.values()) == [1, 1, 1, 1, 2, 2, 2, 2, 4]
      and all((u in S42_SMALL) != (v in S42_SMALL) for u, v in E42),
      "s42-structure",
      "|V|=9 |E|=8 tree; degrees 4,2,2,2,2,1,1,1,1; S_{4,2} is not "
      "isomorphic to S_{3,2} (different order)")

check(is_beta(DELTA, E42), "delta-is-beta-valuation",
      "differences %s == {1..8}"
      % sorted(abs(DELTA[u] - DELTA[v]) for u, v in E42))
check(is_near_alpha(DELTA, E42, S42_SMALL, S42_LARGE),
      "delta-is-near-alpha",
      "V_small = {c,b_0,...,b_3} with labels %s"
      % [DELTA[v] for v in S42_SMALL])
check(alpha_boundaries(DELTA, E42) == [], "delta-is-NOT-alpha",
      "5=delta(b_2) exceeds 3=delta(a_3), so no x straddles every edge")

# Proposition 7, Case A: b(c) is FORCED, and then differences 1 and 2 clash.
sumsA = {bc: 2 * (5 + 6 + 7 + 8) - (4 * bc + (0 + 1 + 2 + 3 + 4 - bc))
         for bc in range(5)}
forcedA = [bc for bc, s in sumsA.items() if s == 36]
pairs1A = [(m, l) for m in (5, 6, 7, 8) for l in (0, 1, 3, 4) if m - l == 1]
pairs2A = [(m, l) for m in (5, 6, 7, 8) for l in (0, 1, 3, 4) if m - l == 2]
clashA = all(set(p) & set(q) for p in pairs1A for q in pairs2A)
check(forcedA == [2] and sorted(sumsA.values()) == [30, 33, 36, 39, 42]
      and pairs1A == [(5, 4)] and sorted(pairs2A) == [(5, 3), (6, 4)]
      and clashA,
      "s42-prop7-caseA",
      "sum = 42-3b(c) forces b(c)=2; centre gives {3,4,5,6} so leaf edges owe "
      "{1,2,7,8}; difference 1 needs %s and every difference-2 pair %s meets "
      "it, so 1 and 2 cannot coexist" % (pairs1A, sorted(pairs2A)))

sumsB = {bc: (4 * bc + (4 + 5 + 6 + 7 + 8 - bc)) - 2 * (0 + 1 + 2 + 3)
         for bc in range(4, 9)}
forcedB = [bc for bc, s in sumsB.items() if s == 36]
pairs1B = [(l, m) for l in (4, 5, 7, 8) for m in (0, 1, 2, 3) if l - m == 1]
pairs2B = [(l, m) for l in (4, 5, 7, 8) for m in (0, 1, 2, 3) if l - m == 2]
clashB = all(set(p) & set(q) for p in pairs1B for q in pairs2B)
check(forcedB == [6] and sorted(sumsB.values()) == [30, 33, 36, 39, 42]
      and pairs1B == [(4, 3)] and sorted(pairs2B) == [(4, 2), (5, 3)]
      and clashB,
      "s42-prop7-caseB",
      "sum = 3b(c)+18 forces b(c)=6; centre gives {3,4,5,6} so leaf edges owe "
      "{1,2,7,8}; difference 1 needs %s and every difference-2 pair %s meets "
      "it" % (pairs1B, sorted(pairs2B)))

tot42, found42 = alpha_free_by_structure(V42, E42, S42_SMALL, S42_LARGE)
check(tot42 == 5760 and found42 == [], "s42-alpha-free-complete",
      "%d candidates permitted by Fact 2(ii) examined, %d alpha-valuations"
      % (tot42, len(found42)))

# The product L = S_{3,2} wtp S_{4,2}, in the paper's naming
# (S_{3,2} written with m_i -> a_{i-1}, l_i -> b_{i-1}).
E32b = [("c", "a0"), ("c", "a1"), ("c", "a2"),
        ("b0", "a0"), ("b1", "a1"), ("b2", "a2")]
G3S, G3L = ["c", "b0", "b1", "b2"], ["a0", "a1", "a2"]
VL, PL, QL, EL = weak_tensor_product(E32b, G3S, G3L, E42, S42_SMALL, S42_LARGE)

LP = {"b0b1": 0, "b0c": 1, "cb3": 2, "cb0": 3, "cb2": 4, "b2b0": 5,
      "b2c": 9, "b2b2": 10, "b1c": 17, "b1b3": 18, "b1b1": 19, "b1b2": 20,
      "b0b0": 22, "b1b0": 23, "b2b1": 24, "cb1": 25, "b2b3": 26,
      "b0b2": 27, "b0b3": 28, "cc": 36}
LQ = {"a0a0": 37, "a1a1": 38, "a2a2": 39, "a2a3": 40, "a2a1": 41,
      "a2a0": 42, "a1a3": 43, "a1a2": 44, "a1a0": 45, "a0a3": 46,
      "a0a2": 47, "a0a1": 48}
NAMES = set(G3S) | set(G3L) | set(S42_SMALL) | set(S42_LARGE)
BL = {}
for token, val in list(LP.items()) + list(LQ.items()):
    parts = [(token[:k], token[k:]) for k in (1, 2, 3)
             if token[:k] in NAMES and token[k:] in NAMES]
    if len(parts) != 1:
        raise SystemExit("ambiguous vertex name %r" % token)
    BL[parts[0]] = val

check((len(VL), len(PL), len(QL), len(EL)) == (32, 20, 12, 48)
      and connected(VL, EL)
      and all((u in PL) != (v in PL) for u, v in EL),
      "L-order-and-size",
      "|V(L)|=%d |P|=%d |Q|=%d |E(L)|=%d, connected, every edge P->Q"
      % (len(VL), len(PL), len(QL), len(EL)))

dL = degrees(VL, EL)
check(dL[("c", "c")] == 12
      and all(dL[("c", "b%d" % j)] == 3 for j in range(4))
      and all(dL[("b%d" % i, "c")] == 4 for i in range(3))
      and all(dL[("b%d" % i, "b%d" % j)] == 1
              for i in range(3) for j in range(4))
      and 12 + 4 * 3 + 3 * 4 + 12 == 48,
      "L-degree-profile",
      "deg cc = 12, deg cb_j = 3, deg b_ic = 4, deg b_ib_j = 1, and "
      "12+12+12+12 = 48 edges as stated in Theorem 8")

check(set(BL) == set(VL) and len(set(BL.values())) == 32
      and min(BL.values()) == 0 and max(BL.values()) == 48,
      "L-labels-injective-in-0-to-48",
      "32 distinct labels in [0,48] = [0,|E(L)|]")

dl = sorted(abs(BL[u] - BL[v]) for u, v in EL)
check(dl == list(range(1, 49)), "L-is-beta-valuation",
      "the 48 edge differences are exactly {1,...,48}")

bdL = alpha_boundaries(BL, EL)
mP, mQ = max(BL[v] for v in PL), min(BL[v] for v in QL)
check(bdL == [36] and mP == 36 and mQ == 37,
      "L-is-alpha-boundary-36",
      "the only straddling x in {0..48} is 36; max_P=%d < %d=min_Q"
      % (mP, mQ))

widL = sum(dL[v] * BL[v] for v in QL) - sum(dL[v] * BL[v] for v in PL)
check(widL == 48 * 49 // 2 == 1176, "L-weighted-degree-identity",
      "sum_Q deg.b - sum_P deg.b = %d = 1+2+...+48" % widL)

# --------------------------------------------------------------------------
# PART 5 -- controls in BOTH polarities on the checker itself
# --------------------------------------------------------------------------

print()
print("=" * 78)
print("PART 5  controls in both polarities: the same predicates, applied to")
print("        objects whose status is already published")
print("=" * 78)

# Positive: the standard alpha-valuation of the path P_{m+1} with m edges.
ok = True
detail = []
for m in range(2, 9):
    Vp = ["v%d" % i for i in range(m + 1)]
    Ep = [("v%d" % i, "v%d" % (i + 1)) for i in range(m)]
    lp = {}
    for i in range(m + 1):
        lp["v%d" % i] = i // 2 if i % 2 == 0 else m - (i - 1) // 2
    bs = alpha_boundaries(lp, Ep)
    ok = ok and is_beta(lp, Ep) and bool(bs)
    detail.append("m=%d:x=%d" % (m, bs[0] if bs else -1))
check(ok, "control-paths-alpha", "zigzag labelling verifies: " + " ".join(detail))

# Positive: Rosa's alpha-valuation of K_{p,q}.
ok = True
detail = []
for p in range(1, 5):
    for q in range(1, 5):
        Va = ["x%d" % i for i in range(p)]
        Vb = ["y%d" % j for j in range(q)]
        Ec = [(a, b) for a in Va for b in Vb]
        lc = {"x%d" % i: i for i in range(p)}
        lc.update({"y%d" % j: p * (j + 1) for j in range(q)})
        bs = alpha_boundaries(lc, Ec)
        ok = ok and is_beta(lc, Ec) and bs == [p - 1]
        detail.append("K_%d,%d" % (p, q))
check(ok, "control-Kpq-alpha",
      "Rosa's labelling of %s all verify with boundary p-1" % ",".join(detail))

# BOTH polarities from one published iff: C_m has an alpha-valuation exactly
# when m = 0 mod 4.  Exhaustive over every injection, so the NO answers are
# proofs, not search failures.
res = {}
for m in range(4, 8):
    Vc = ["z%d" % i for i in range(m)]
    Ec = [("z%d" % i, "z%d" % ((i + 1) % m)) for i in range(m)]
    hit = False
    for perm in itertools.permutations(range(m + 1), m):
        lab = dict(zip(Vc, perm))
        if is_beta(lab, Ec) and alpha_boundaries(lab, Ec):
            hit = True
            break
    res[m] = hit
check(res == {4: True, 5: False, 6: False, 7: False},
      "control-cycles-alpha-iff-0-mod-4 (anti-control)",
      "exhaustive over all injections: alpha exists for C_4 and for NONE of "
      "C_5, C_6, C_7, reproducing the published 'C_m has an alpha-valuation "
      "precisely if m = 0 mod 4' on 4 <= m <= 7")

# --------------------------------------------------------------------------

print()
print("SCOPE.  This program checks the objects PRINTED in the paper and the")
print("        finite claims made about them.  It does not search for further")
print("        counterexamples, does not claim K is of least order, and")
print("        establishes nothing about factors on more than 9 vertices.")
print()
if _FAILED:
    print("VERDICT: %d of %d CHECKS FAILED: %s"
          % (len(_FAILED), _N, ", ".join(_FAILED)))
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % _N)
sys.exit(0)
