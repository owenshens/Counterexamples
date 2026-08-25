#!/usr/bin/env python3
"""Independent verification of a five-vertex counterexample to the reverse
implication of Conjecture 5.8 of Curto-Geneson-Morrison (permitted motifs in
composite CTLNs).

TAKEN FROM THE PAPER (inputs, transcribed and never used as evidence for
themselves):
  * the CTLN convention  W_ii = 0,  W_ij = -1 + eps  if j -> i,
    W_ij = -1 - delta otherwise;  theta > 0, delta > 0, 0 < eps < delta/(delta+1);
    full support permitted  <=>  the solution of (I-W)x = theta*1 is > 0;
  * the skeleton  Ghat  on {A,B,C,D} with arcs  A<->D, B<->C, B->D, C->A;
  * the components  G_A = K1, G_B = K2, G_C = K1, G_D = K1;
  * the exhibited composite G on {a,b1,b2,c,d} with arcs
    a<->d, b1<->b2, b_i<->c, b_i->d (i=1,2), c->a;
  * the two displayed matrices I-What and I-W, with alpha = 1-eps, beta = 1+delta;
  * hhat = 2delta^2 - 2 delta eps + 6 delta - eps^2 + 2 eps  and
    h    = 3delta^2 - 3 delta eps + 9 delta - 2 eps^2 + 4 eps;
  * the claimed solution vectors and the determinants -eps^2*hhat, -eps^3*h;
  * the appendix lists of 12 distinct principal and 15 distinct Cramer
    determinants;
  * Table 1's six permitted three-vertex skeletons with their eta and Q data.

DERIVED HERE (nothing below is assumed):
  * the two matrices are rebuilt from the graphs and compared entrywise;
  * G is rebuilt from the skeleton and the components by the composite rule;
  * all determinants, solution vectors and minors are recomputed by exact
    symbolic elimination over Q[eps,delta] and compared to the paper;
  * every sign claim is certified over the WHOLE legal region by the change of
    variables  delta = d,  eps = v*d/((1+v)(1+d))  with d,v > 0, which is a
    bijection onto the legal region: a polynomial whose numerator has all
    coefficients of one sign there is strictly of that sign;
  * the load-bearing conclusion (skeleton permitted, all components permitted,
    G forbidden) is decided from the Cramer determinants of the REBUILT
    matrices, never from the paper's transcribed solution vectors -- those are
    only required to agree -- and is recomputed a second time by exact rational
    linear solves at a grid of legal (eps,delta,theta);
  * every index-encoded copy of the skeleton used by the numeric sweeps is
    derived from the transcribed arc set, so no second transcription can
    silently diverge from the exhibited object;
  * the permitted 1-, 2- and 3-vertex skeleton censuses are re-enumerated
    exhaustively from the CTLN criterion, not read off the cited figure, and
    Table 1's corner data (multiaffinity, eta, the sets Q) is recomputed;
  * the corner argument is re-run for every permitted 2- and 3-vertex skeleton,
    so the exclusion of 2 and 3 components does not rest on the cited theorem;
  * the reduction lemma is tested exhaustively and exactly on every composite of
    at most 4 vertices with 2 or 3 arbitrary permitted components;
  * all 4096 skeletons on 4 vertices are swept with components (K1,K2,K1,K1) to
    confirm the exhibited graph really is one of the counterexamples that
    family contains.

NOT ESTABLISHED HERE: the reduction lemma and the total-activity bound
beta^{-1} < 1^T(I-W_i)^{-1}1 < alpha^{-1} (quoted by the paper from Lemma 3.1
of the cited work) are verified only on the finite family just described and on
the components the counterexample uses, not for components of unbounded order;
consequently the minimality theorem (Theorem 3), whose exclusion of two and
three components quantifies over components of arbitrary order, rests here on
the paper's hand proof of that lemma; and the forward implication of the
conjecture (the cited Theorem 5.2) is not re-proved, being irrelevant to the
refutation.  All of this is also printed by the program, in the closing
NOT RE-RUN paragraph, so a referee reading only stdout sees it.
"""
from itertools import permutations, combinations
from fractions import Fraction

# ---------- exact multivariate polynomials over Q in (eps, delta, d1, d2, d3)
NV = 5
VN = ("eps", "delta", "d1", "d2", "d3")


def mono(idx, e=1, c=1):
    k = [0] * NV
    k[idx] = e
    return {tuple(k): Fraction(c)}


def const(c):
    c = Fraction(c)
    return {} if c == 0 else {tuple([0] * NV): c}


def padd(*ps):
    r = {}
    for p in ps:
        for k, v in p.items():
            r[k] = r.get(k, Fraction(0)) + v
    return {k: v for k, v in r.items() if v != 0}


def pneg(p):
    return {k: -v for k, v in p.items()}


def psub(a, b):
    return padd(a, pneg(b))


def pmul(a, b):
    r = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            k = tuple(x + y for x, y in zip(k1, k2))
            r[k] = r.get(k, Fraction(0)) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}


def psmul(c, p):
    c = Fraction(c)
    return {} if c == 0 else {k: v * c for k, v in p.items()}


def pprod(*ps):
    r = const(1)
    for p in ps:
        r = pmul(r, p)
    return r


def pshow(p):
    if not p:
        return "0"
    out = []
    for k in sorted(p, key=lambda k: (-sum(k), k)):
        s = str(p[k])
        for i, e in enumerate(k):
            if e:
                s += "*" + VN[i] + ("^%d" % e if e > 1 else "")
        out.append(s)
    return " + ".join(out)


def pev(p, vals):
    s = Fraction(0)
    for k, c in p.items():
        t = c
        for i, e in enumerate(k):
            if e:
                t *= vals[i] ** e
        s += t
    return s


def psubst(p, table):
    out = {}
    for k, c in p.items():
        t = const(c)
        for s, vp in table.items():
            for _ in range(k[s]):
                t = pmul(t, vp)
        base = list(k)
        for s in table:
            base[s] = 0
        out = padd(out, pmul({tuple(base): Fraction(1)}, t))
    return out


EPS = mono(0)
DEL = mono(1)
ALPHA = psub(const(1), EPS)   # 1 - eps
BETA = padd(const(1), DEL)    # 1 + delta


def pdet(M):
    n = len(M)
    tot = {}
    for perm in permutations(range(n)):
        seen = [False] * n
        sgn = 1
        for i in range(n):
            if not seen[i]:
                j, L = i, 0
                while not seen[j]:
                    seen[j] = True
                    j = perm[j]
                    L += 1
                if L % 2 == 0:
                    sgn = -sgn
        t = const(sgn)
        for i in range(n):
            t = pmul(t, M[i][perm[i]])
            if not t:
                break
        tot = padd(tot, t)
    return tot


def pcramer(M, i):
    N = [row[:] for row in M]
    for r in range(len(M)):
        N[r][i] = const(1)
    return pdet(N)

CHECKS = []
# facts about the SCOPE of what was actually exercised, recorded while the
# checks run so the closing NOT RE-RUN paragraph can print derived numbers
# instead of hardcoded ones
SCOPE = {}


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    if detail:
        print("%s %s [%s]" % (tag, name, detail))
    else:
        print("%s %s" % (tag, name))
    return bool(ok)


# ---------- sign certificates on the legal region
DV = mono(0)   # reused slots: slot 0 = d, slot 1 = v
VV = mono(1)
ONE_D = padd(const(1), DV)
ONE_V = padd(const(1), VV)
# the single definition of the change of variables used by every certificate:
#   delta = d,   eps = EPS_NUM(d,v) / EPS_DEN(d,v),   d > 0, v > 0.
EPS_NUM = pmul(VV, DV)
EPS_DEN = pmul(ONE_V, ONE_D)


def cert(p):
    """Sign of p on {delta>0, 0<eps<delta/(delta+1)}: 1, -1, 0, or None."""
    if not p:
        return 0
    for k in p:
        if any(k[i] for i in (2, 3, 4)):
            raise ValueError("free component variable in certificate")
    m = max(k[0] for k in p)
    tot = {}
    for k, c in p.items():
        i, j = k[0], k[1]
        t = const(c)
        for _ in range(i):
            t = pmul(t, EPS_NUM)
        for _ in range(j):
            t = pmul(t, DV)
        for _ in range(m - i):
            t = pmul(t, EPS_DEN)
        tot = padd(tot, t)
    if not tot:
        return 0
    vs = list(tot.values())
    if all(v > 0 for v in vs):
        return 1
    if all(v < 0 for v in vs):
        return -1
    return None


def legal_points():
    """Exact rational (eps, delta) strictly inside the legal region."""
    pts = []
    for dl in (Fraction(1, 100), Fraction(1, 10), Fraction(1, 2), Fraction(1),
               Fraction(3), Fraction(20), Fraction(500)):
        hi = dl / (dl + 1)
        for t in (Fraction(1, 1000), Fraction(1, 20), Fraction(1, 3),
                  Fraction(1, 2), Fraction(3, 4), Fraction(99, 100)):
            pts.append((t * hi, dl))
    return pts


def rsolve(M, rhs):
    """Exact Gaussian solve of M x = rhs over Q; returns None if singular."""
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(n)] + [Fraction(rhs[i])]
         for i in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return None
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [v / pv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]


# ---------- graphs
def imw_named(nodes, arcs):
    n = len(nodes)
    return [[const(1) if i == j else
             (ALPHA if (nodes[j], nodes[i]) in arcs else BETA)
             for j in range(n)] for i in range(n)]


def imw_idx(n, arcs):
    return [[const(1) if i == j else
             (ALPHA if (j, i) in arcs else BETA) for j in range(n)]
            for i in range(n)]


def imw_num(n, arcs, eps, dl):
    a, b = 1 - eps, 1 + dl
    return [[Fraction(1) if i == j else (a if (j, i) in arcs else b)
             for j in range(n)] for i in range(n)]


def composite(sk_nodes, sk_arcs, comps):
    """comps[v] = (label list, internal arc set). Returns (nodes, arcs)."""
    nodes = []
    for v in sk_nodes:
        nodes.extend(comps[v][0])
    arcs = set()
    for v in sk_nodes:
        arcs |= set(comps[v][1])
    for u in sk_nodes:
        for w in sk_nodes:
            if u != w and (u, w) in sk_arcs:
                for x in comps[u][0]:
                    for y in comps[w][0]:
                        arcs.add((x, y))
    return nodes, arcs


def iso_classes(n):
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    seen, out = {}, []
    for mask in range(1 << len(pairs)):
        arcs = frozenset(p for k, p in enumerate(pairs) if mask >> k & 1)
        if arcs in seen:
            continue
        for pi in permutations(range(n)):
            seen[frozenset((pi[a], pi[b]) for a, b in arcs)] = len(out)
        out.append(arcs)
    return out


def uniform_status(n, arcs):
    """PERMITTED / FORBIDDEN over the whole legal region, or UNRESOLVED."""
    M = imw_idx(n, arcs)
    d = pdet(M)
    nums = [pcramer(M, i) for i in range(n)]
    prod_signs = [cert(pmul(N, d)) for N in nums]
    if all(s == 1 for s in prod_signs):
        return "PERMITTED"
    if any(s == -1 for s in prod_signs):
        return "FORBIDDEN"
    ns = [cert(N) for N in nums]
    if any(s == 0 for s in ns):
        return "FORBIDDEN"
    if any(s == 1 for s in ns) and any(s == -1 for s in ns):
        return "FORBIDDEN"   # positivity fails whichever sign det has
    return "UNRESOLVED"


# ---------- the exhibited object, transcribed from the paper
SK_NODES = ["A", "B", "C", "D"]
SK_ARCS = {("A", "D"), ("D", "A"), ("B", "C"), ("C", "B"),
           ("B", "D"), ("C", "A")}
COMPS = {"A": (["a"], set()),
         "B": (["b1", "b2"], {("b1", "b2"), ("b2", "b1")}),
         "C": (["c"], set()),
         "D": (["d"], set())}
G_NODES = ["a", "b1", "b2", "c", "d"]
G_ARCS = {("a", "d"), ("d", "a"), ("b1", "b2"), ("b2", "b1"),
          ("b1", "c"), ("c", "b1"), ("b2", "c"), ("c", "b2"),
          ("b1", "d"), ("b2", "d"), ("c", "a")}
HHAT = padd(psmul(2, pmul(DEL, DEL)), psmul(-2, pmul(DEL, EPS)),
            psmul(6, DEL), psmul(-1, pmul(EPS, EPS)), psmul(2, EPS))
HFUL = padd(psmul(3, pmul(DEL, DEL)), psmul(-3, pmul(DEL, EPS)),
            psmul(9, DEL), psmul(-2, pmul(EPS, EPS)), psmul(4, EPS))
A_ = "alpha"
B_ = "beta"
PAPER_SK_MATRIX = [["1", B_, A_, A_],
                   [B_, "1", A_, B_],
                   [B_, A_, "1", B_],
                   [A_, A_, B_, "1"]]
PAPER_G_MATRIX = [["1", B_, B_, A_, A_],
                  [B_, "1", A_, A_, B_],
                  [B_, A_, "1", A_, B_],
                  [B_, A_, A_, "1", B_],
                  [A_, A_, A_, B_, "1"]]


def check_object():
    nodes, arcs = G_NODES, G_ARCS
    ok = (len(nodes) == 5 and len(set(nodes)) == 5 and len(arcs) == 11)
    ok = ok and all(u in nodes and w in nodes and u != w for u, w in arcs)
    deg = {v: (sum(1 for x, y in arcs if x == v),
               sum(1 for x, y in arcs if y == v)) for v in nodes}
    ck("exhibited_graph_wellformed", ok,
       "V=%s |V|=%d |E|=%d outdeg,indeg=%s"
       % (",".join(nodes), len(nodes), len(arcs),
          ";".join("%s:%d/%d" % (v, deg[v][0], deg[v][1]) for v in nodes)))
    sk_ok = (len(SK_NODES) == 4 and len(SK_ARCS) == 6 and
             all(u in SK_NODES and w in SK_NODES and u != w
                 for u, w in SK_ARCS))
    bidir = {(u, w) for u, w in SK_ARCS if (w, u) in SK_ARCS}
    ck("skeleton_wellformed", sk_ok,
       "|V|=4 |E|=%d bidirected pairs=%d single arcs=%d"
       % (len(SK_ARCS), len(bidir) // 2, len(SK_ARCS) - len(bidir)))
    sizes = [len(COMPS[v][0]) for v in SK_NODES]
    cnodes, carcs = composite(SK_NODES, SK_ARCS, COMPS)
    ck("composite_construction_reproduces_G",
       sorted(cnodes) == sorted(G_NODES) and carcs == G_ARCS and
       sum(sizes) == 5,
       "component orders %s sum=%d; rebuilt arcs=%d, symmetric difference=%d"
       % (sizes, sum(sizes), len(carcs), len(carcs ^ G_ARCS)))


def check_matrices():
    def sym(M, paper):
        for i in range(len(M)):
            for j in range(len(M)):
                want = (const(1) if paper[i][j] == "1"
                        else (ALPHA if paper[i][j] == A_ else BETA))
                if M[i][j] != want:
                    return False, "entry (%d,%d)" % (i + 1, j + 1)
        return True, "all %d entries" % (len(M) ** 2)
    ok1, d1 = sym(imw_named(SK_NODES, SK_ARCS), PAPER_SK_MATRIX)
    ck("skeleton_matrix_matches_paper", ok1, d1)
    ok2, d2 = sym(imw_named(G_NODES, G_ARCS), PAPER_G_MATRIX)
    ck("full_matrix_matches_paper", ok2, d2)


def check_determinants():
    dS = pdet(imw_named(SK_NODES, SK_ARCS))
    claim = pneg(pprod(EPS, EPS, HHAT))
    ck("det_skeleton_equals_minus_eps2_hhat", dS == claim,
       "computed %s" % pshow(dS))
    dG = pdet(imw_named(G_NODES, G_ARCS))
    claimG = pneg(pprod(EPS, EPS, EPS, HFUL))
    ck("det_full_equals_minus_eps3_h", dG == claimG,
       "computed %s" % pshow(dG))
    ck("hhat_and_h_strictly_positive",
       cert(HHAT) == 1 and cert(HFUL) == 1,
       "certificate signs hhat: %s, h: %s on the whole legal region"
       % (cert(HHAT), cert(HFUL)))
    r1 = padd(psmul(2, pmul(DEL, psub(DEL, EPS))), psmul(6, DEL),
              pmul(EPS, psub(const(2), EPS)))
    r2 = padd(psmul(3, pmul(DEL, psub(DEL, EPS))), psmul(9, DEL),
              psmul(2, pmul(EPS, psub(const(2), EPS))))
    ck("paper_rewritings_of_hhat_and_h_are_identities",
       r1 == HHAT and r2 == HFUL,
       "2delta(delta-eps)+6delta+eps(2-eps) == hhat and "
       "3delta(delta-eps)+9delta+2eps(2-eps) == h")
    ck("legal_region_implies_eps_lt_min_delta_1",
       cert(psub(DEL, EPS)) == 1 and cert(psub(const(1), EPS)) == 1,
       "certificates delta-eps=+1, 1-eps=+1")


def solved_signs(M):
    """Signs of the coordinates of (I-W)^{-1}1, recomputed from M alone.

    x_i = C_i/det, so sign(x_i) = sign(C_i * det).  Nothing transcribed from the
    paper enters here, so a corrupted graph changes the answer.
    """
    d = pdet(M)
    return [cert(pmul(pcramer(M, i), d)) for i in range(len(M))]


def solution_ratio_equals(M, nums, den):
    """x_i == nums[i]/den for all i, as identities of rational functions."""
    d = pdet(M)
    bad = []
    for i in range(len(M)):
        if pmul(pcramer(M, i), den) != pmul(nums[i], d):
            bad.append(i + 1)
    return (not bad), bad


def check_hypotheses():
    MS = imw_named(SK_NODES, SK_ARCS)
    numsS = [DEL, padd(psmul(2, DEL), EPS), padd(psmul(2, DEL), EPS), DEL]
    ok, bad = solution_ratio_equals(MS, numsS, HHAT)
    ck("skeleton_solution_matches_paper", ok,
       "x = (delta, 2delta+eps, 2delta+eps, delta)/hhat" if ok
       else "coordinates %s differ" % bad)
    signs = solved_signs(MS)
    claimed = [cert(pmul(n, HHAT)) for n in numsS]
    ck("hypothesis_skeleton_is_permitted",
       all(s == 1 for s in signs) and claimed == signs,
       "all four coordinates of (I-W)^{-1}1 recomputed from the rebuilt "
       "skeleton matrix and certified > 0 on the legal region, signs=%s "
       "(paper's transcribed vector agrees: %s)" % (signs, claimed))
    k1 = rsolve([[Fraction(1)]], [Fraction(1)])
    MK2 = [[const(1), ALPHA], [ALPHA, const(1)]]
    n2 = [pcramer(MK2, 0), pcramer(MK2, 1)]
    d2 = pdet(MK2)
    claim_den = psub(const(2), EPS)
    k2_ok = (pmul(n2[0], claim_den) == pmul(const(1), d2) and
             pmul(n2[1], claim_den) == pmul(const(1), d2))
    n_k1 = sum(1 for v in SK_NODES if COMPS[v][0] and not COMPS[v][1]
               and len(COMPS[v][0]) == 1)
    ck("hypothesis_all_components_permitted",
       k1 == [Fraction(1)] and k1[0] > 0 and k2_ok and
       cert(claim_den) == 1,
       "K1 solution=1>0 (%d copies, counted from the transcribed components); "
       "K2 solution=(1,1)/(2-eps) with 2-eps certified > 0" % n_k1)


def check_conclusion():
    MG = imw_named(G_NODES, G_ARCS)
    t = padd(psmul(2, DEL), EPS)
    numsG = [pneg(pmul(DEL, DEL)), pmul(EPS, t), pmul(EPS, t), pmul(EPS, t),
             padd(pmul(DEL, DEL), psmul(3, pmul(DEL, EPS)), pmul(EPS, EPS))]
    den = pmul(EPS, HFUL)
    ok, bad = solution_ratio_equals(MG, numsG, den)
    ck("full_solution_matches_paper", ok,
       "x = (-delta^2, eps(2delta+eps) three times, "
       "delta^2+3delta*eps+eps^2) / (eps*h)" if ok
       else "coordinates %s differ" % bad)
    signs = solved_signs(MG)
    claimed = [cert(pmul(n, den)) for n in numsG]
    ck("conclusion_full_support_solution_has_a_negative_coordinate",
       signs[0] == -1 and all(s == 1 for s in signs[1:]) and
       claimed == signs,
       "x_a certified < 0 and x_b1,x_b2,x_c,x_d certified > 0, from the "
       "Cramer determinants of the rebuilt I-W(G) itself; signs=%s (the "
       "paper's transcribed vector -delta^2, eps(2delta+eps), "
       "delta^2+3delta*eps+eps^2 over eps*h agrees: %s)" % (signs, claimed))
    ck("conclusion_G_is_forbidden", signs[0] == -1 and ok and claimed == signs,
       "the unique full-support fixed-point candidate is not strictly "
       "positive for any legal (eps,delta,theta), so the full vertex set "
       "of G is forbidden")


def check_grid():
    """Second, formula-free route: exact rational solves at legal points."""
    MSn = None
    bad = []
    npts = 0
    for eps, dl in legal_points():
        for th in (Fraction(1), Fraction(1, 7), Fraction(13),
                   Fraction(1000000)):
            npts += 1
            si = {v: i for i, v in enumerate(SK_NODES)}
            sa = {(si[u], si[w]) for u, w in SK_ARCS}
            xs = rsolve(imw_num(len(SK_NODES), sa, eps, dl),
                        [th] * len(SK_NODES))
            gi = {v: i for i, v in enumerate(G_NODES)}
            ga = {(gi[u], gi[w]) for u, w in G_ARCS}
            xg = rsolve(imw_num(5, ga, eps, dl), [th] * 5)
            xk2 = rsolve(imw_num(2, {(0, 1), (1, 0)}, eps, dl), [th] * 2)
            xk1 = rsolve(imw_num(1, set(), eps, dl), [th])
            good = (xs is not None and all(v > 0 for v in xs) and
                    xk2 is not None and all(v > 0 for v in xk2) and
                    xk1 == [th] and
                    xg is not None and xg[0] < 0 and
                    all(v > 0 for v in xg[1:]))
            if not good:
                bad.append((eps, dl, th))
    ck("independent_grid_solve_confirms_refutation", not bad,
       "%d legal (eps,delta,theta) triples: skeleton > 0, K1 and K2 > 0, "
       "x_a < 0 in G at every one" % npts if not bad
       else "failures at %s" % bad[:3])


def check_projection_pattern():
    idx = {v: i for i, v in enumerate(SK_NODES)}
    Aadj = [[1 if (SK_NODES[j], SK_NODES[i]) in SK_ARCS else 0
             for j in range(4)] for i in range(4)]
    rows = [[Fraction(x) for x in r] for r in Aadj]
    rank = 0
    for c in range(4):
        p = next((r for r in range(rank, 4) if rows[r][c] != 0), None)
        if p is None:
            continue
        rows[rank], rows[p] = rows[p], rows[rank]
        pv = rows[rank][c]
        rows[rank] = [v / pv for v in rows[rank]]
        for r in range(4):
            if r != rank and rows[r][c] != 0:
                f = rows[r][c]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
    bA = ("B", "A") in SK_ARCS
    bC = ("B", "C") in SK_ARCS
    ck("skeleton_projection_pattern_is_not_uniform",
       rank >= 2 and bA != bC,
       "intercomponent pattern has rank %d (> 1) and B->A=%s differs from "
       "B->C=%s, so vertices of the K2 component project differently to "
       "G_A and G_C" % (rank, int(bA), int(bC)))


def paper_principal_list():
    d2 = pmul(DEL, DEL)
    return [const(1),
            pneg(pmul(DEL, padd(DEL, const(2)))),
            pmul(EPS, psub(const(2), EPS)),
            psub(pmul(EPS, padd(DEL, const(1))), DEL),
            pneg(pprod(DEL, EPS, padd(DEL, pneg(EPS), const(3)))),
            pmul(pmul(EPS, EPS), psub(const(3), psmul(2, EPS))),
            pmul(EPS, psub(pmul(padd(psmul(2, DEL), const(1)), EPS),
                           psmul(2, DEL))),
            pneg(pprod(DEL, EPS, EPS, padd(DEL, psmul(-2, EPS), const(4)))),
            pneg(pmul(EPS, padd(psmul(2, d2), psmul(4, DEL), EPS))),
            pneg(pmul(pmul(EPS, EPS), padd(psmul(2, d2),
                                           pneg(pmul(DEL, EPS)),
                                           psmul(5, DEL), EPS))),
            pneg(pprod(EPS, EPS, HHAT)),
            pneg(pprod(EPS, EPS, EPS, HFUL))]


def paper_cramer_list():
    t = padd(psmul(2, DEL), EPS)
    d2 = pmul(DEL, DEL)
    return [const(1), EPS, pneg(DEL), pmul(EPS, EPS), pneg(pmul(DEL, EPS)),
            pneg(pprod(DEL, EPS, EPS)), pprod(DEL, DEL, EPS, EPS),
            pneg(pmul(DEL, padd(DEL, psmul(2, EPS)))),
            pneg(pmul(EPS, t)), pneg(pprod(EPS, EPS, t)),
            pneg(pprod(EPS, EPS, EPS, t)),
            padd(d2, pmul(DEL, EPS), pmul(EPS, EPS)),
            pmul(EPS, padd(psmul(2, d2), psmul(2, pmul(DEL, EPS)),
                           pmul(EPS, EPS))),
            pneg(pmul(EPS, padd(psmul(2, d2), psmul(4, pmul(DEL, EPS)),
                                pmul(EPS, EPS)))),
            pneg(pmul(pmul(EPS, EPS), padd(d2, psmul(3, pmul(DEL, EPS)),
                                           pmul(EPS, EPS))))]


def _canon(p):
    return tuple(sorted(p.items()))


def all_minor_data(M):
    """(distinct principal dets, distinct Cramer dets, #minors, #Cramer).

    The last two are the sizes of the enumeration actually performed, so the
    counts printed by the callers are derived rather than asserted.
    """
    n = len(M)
    prin, cram = {}, {}
    n_minors = n_cramer = 0
    for k in range(1, n + 1):
        for S in combinations(range(n), k):
            sub = [[M[i][j] for j in S] for i in S]
            d = pdet(sub)
            prin[_canon(d)] = d
            n_minors += 1
            for i in range(k):
                c = pcramer(sub, i)
                cram[_canon(c)] = c
                n_cramer += 1
    return prin, cram, n_minors, n_cramer


def check_nondegeneracy():
    MG = imw_named(G_NODES, G_ARCS)
    prin, cram, n_minors, n_cramer = all_minor_data(MG)
    want_p = {_canon(p): p for p in paper_principal_list()}
    want_c = {_canon(p): p for p in paper_cramer_list()}
    ck("principal_determinant_list_matches_appendix",
       set(prin) == set(want_p),
       "%d principal minors of I-W collapse to %d distinct polynomials, "
       "equal as a set to the paper's %d" % (n_minors, len(prin), len(want_p))
       if set(prin) == set(want_p) else
       "computed %d vs paper %d; missing from paper %s; extra in paper %s"
       % (len(prin), len(want_p),
          [pshow(prin[k]) for k in set(prin) - set(want_p)][:3],
          [pshow(want_p[k]) for k in set(want_p) - set(prin)][:3]))
    ck("cramer_determinant_list_matches_appendix",
       set(cram) == set(want_c),
       "%d Cramer determinants collapse to %d distinct polynomials, equal "
       "as a set to the paper's %d" % (n_cramer, len(cram), len(want_c))
       if set(cram) == set(want_c) else
       "computed %d vs paper %d; missing from paper %s; extra in paper %s"
       % (len(cram), len(want_c),
          [pshow(cram[k]) for k in set(cram) - set(want_c)][:3],
          [pshow(want_c[k]) for k in set(want_c) - set(cram)][:3]))
    allp = list(prin.values()) + list(cram.values())
    signs = [cert(p) for p in allp]
    zero = [pshow(p) for p, s in zip(allp, signs) if s != 1 and s != -1]
    ck("uniform_nondegeneracy_no_determinant_vanishes", not zero,
       "all %d distinct principal and Cramer determinants certified of "
       "constant nonzero sign on the whole legal region (%d positive, "
       "%d negative)" % (len(allp), signs.count(1), signs.count(-1))
       if not zero else "uncertified: %s" % zero[:3])
    prinS, cramS, _, _ = all_minor_data(imw_named(SK_NODES, SK_ARCS))
    ck("skeleton_determinant_lists_are_sublists",
       set(prinS) <= set(prin) and set(cramS) <= set(cram),
       "skeleton contributes %d principal and %d Cramer polynomials, all "
       "already among the graph's" % (len(prinS), len(cramS)))


def component_d(M):
    """d = 1/(1^T (I-W)^{-1} 1) as a polynomial, or None if not polynomial."""
    n = len(M)
    den = {}
    for i in range(n):
        den = padd(den, pcramer(M, i))
    num = pdet(M)
    for cand in (psub(const(1), psmul(Fraction(1, 2), EPS)), const(1)):
        if pmul(cand, den) == num:
            return cand
    return None


def check_reduction_lemma():
    # every component matrix, its d, and the block partition are read off the
    # transcribed components, not hardwired
    n = len(SK_NODES)
    ds = [component_d(imw_named(COMPS[v][0], COMPS[v][1])) for v in SK_NODES]
    gpos = {x: i for i, x in enumerate(G_NODES)}
    blocks = [[gpos[x] for x in COMPS[v][0]] for v in SK_NODES]
    okb = all(d is not None and cert(psub(d, ALPHA)) == 1 and
              cert(psub(BETA, d)) == 1 for d in ds)
    ck("lemma_component_reciprocals_lie_strictly_between_alpha_and_beta",
       okb, "component reciprocals d = %s, each certified > alpha and < beta"
       % ", ".join("d(G_%s)=%s" % (v, pshow(d) if d else d)
                   for v, d in zip(SK_NODES, ds)))
    if not okb:
        ds = [d if d is not None else const(1) for d in ds]
    A = [[ds[i] if i == j else
          (ALPHA if (SK_NODES[j], SK_NODES[i]) in SK_ARCS else BETA)
          for j in range(n)] for i in range(n)]
    dA = pdet(A)
    MG = imw_named(G_NODES, G_ARCS)
    dG = pdet(MG)
    bad = []
    for i, blk in enumerate(blocks):
        tot = {}
        for j in blk:
            tot = padd(tot, pcramer(MG, j))
        if pmul(pcramer(A, i), dG) != pmul(tot, dA):
            bad.append(i + 1)
    ck("lemma_reduction_reproduces_block_totals", not bad,
       "A(d)^{-1}1 equals the four block activity sums of (I-W_G)^{-1}1 as "
       "identities of rational functions" if not bad
       else "components %s differ" % bad)
    sg = [cert(pmul(pcramer(A, i), dA)) for i in range(n)]
    ck("lemma_applied_to_G_also_reports_forbidden",
       sg[0] == -1 and all(s == 1 for s in sg[1:]),
       "the 4x4 reduction alone gives S_A < 0, S_B,S_C,S_D > 0; signs=%s"
       % sg)


# ---------- Table 1 of the paper, transcribed
Q1 = padd(DEL, const(1))                                  # delta + 1
Q2 = padd(psmul(2, DEL), pneg(EPS), const(3))              # 2delta - eps + 3
Q3 = padd(DEL, pneg(EPS), const(2))                        # delta - eps + 2
Q4 = psub(const(1), EPS)                                   # 1 - eps
Q5 = padd(DEL, psmul(-2, EPS), const(3))                   # delta - 2eps + 3
PAPER_TABLE = [
    ("empty", frozenset(), 1, [Q1, Q2]),
    ("1<->2", frozenset({(0, 1), (1, 0)}), -1, [Q1, Q3, Q2]),
    ("1<->2,1->3", frozenset({(0, 1), (1, 0), (0, 2)}), -1, [Q4, Q1, Q3]),
    ("1<->2,1<->3", frozenset({(0, 1), (1, 0), (0, 2), (2, 0)}), -1,
     [Q4, Q3, Q5]),
    ("1->2->3->1", frozenset({(0, 1), (1, 2), (2, 0)}), 1, [Q3, Q5, Q2]),
    ("K3", frozenset({(0, 1), (1, 0), (0, 2), (2, 0), (1, 2), (2, 1)}), 1,
     [Q4, Q5]),
]
CORNERS = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]


def a_of_d(arcs):
    """A(d) for a three-vertex skeleton, with d1,d2,d3 free."""
    return [[mono(2 + i) if i == j else
             (ALPHA if (j, i) in arcs else BETA) for j in range(3)]
            for i in range(3)]


def check_census():
    st1 = [uniform_status(1, a) for a in iso_classes(1)]
    ck("census_one_vertex_skeleton", st1 == ["PERMITTED"],
       "the single 1-vertex skeleton is permitted; one component is "
       "tautological")
    cl2 = iso_classes(2)
    st2 = [uniform_status(2, a) for a in cl2]
    perm2 = sorted(sorted(a) for a, s in zip(cl2, st2) if s == "PERMITTED")
    ck("census_two_vertex_permitted_skeletons",
       len(cl2) == 3 and "UNRESOLVED" not in st2 and
       perm2 == [[], [(0, 1), (1, 0)]],
       "3 isomorphism classes on 2 vertices; permitted are exactly the "
       "independent set and the bidirected 2-clique: %s" % perm2)
    cl3 = iso_classes(3)
    st3 = [uniform_status(3, a) for a in cl3]
    perm3 = [a for a, s in zip(cl3, st3) if s == "PERMITTED"]
    orbits = []
    for _, arcs, _, _ in PAPER_TABLE:
        orbits.append({frozenset((pi[u], pi[v]) for u, v in arcs)
                       for pi in permutations(range(3))})
    matched = [sum(1 for a in perm3 if a in orb) for orb in orbits]
    ck("census_three_vertex_permitted_skeletons",
       len(cl3) == 16 and "UNRESOLVED" not in st3 and len(perm3) == 6 and
       matched == [1] * 6,
       "16 isomorphism classes on 3 vertices; exactly %d permitted, and "
       "each of the paper's 6 rows matches exactly one of them (%s); the "
       "other %d are forbidden for every legal parameter"
       % (len(perm3), matched, len(cl3) - len(perm3)))


def check_certificates_against_evaluation():
    """Every certified sign must match the exact value at every legal point."""
    MG = imw_named(G_NODES, G_ARCS)
    prin, cram, _, _ = all_minor_data(MG)
    pool = list(prin.values()) + list(cram.values()) + [HHAT, HFUL,
                                                        Q1, Q2, Q3, Q4, Q5]
    pool += [psub(DEL, EPS), psub(const(1), EPS),
             psub(pmul(EPS, BETA), DEL),
             psub(pmul(padd(psmul(2, DEL), const(1)), EPS), psmul(2, DEL))]
    bad = []
    tested = 0
    for p in pool:
        s = cert(p)
        if s not in (1, -1):
            continue
        for eps, dl in legal_points():
            v = pev(p, [eps, dl, 0, 0, 0])
            tested += 1
            if (v > 0) != (s == 1):
                bad.append((pshow(p), str(eps), str(dl), str(v)))
    ck("certificates_agree_with_exact_evaluation", not bad,
       "%d certified sign claims x %d legal points = %d exact evaluations, "
       "every one agreeing with its certificate"
       % (len(pool), len(legal_points()), tested) if not bad
       else "disagreements: %s" % bad[:2])


def check_table_corners():
    bad_ma, bad_num, bad_q, bad_nz = [], [], [], []
    dpe2 = pmul(padd(DEL, EPS), padd(DEL, EPS))
    for name, arcs, eta, qs in PAPER_TABLE:
        A = a_of_d(arcs)
        polys = [pdet(A)] + [pcramer(A, i) for i in range(3)]
        if any(k[s] > 1 for p in polys for k in p for s in (2, 3, 4)):
            bad_ma.append(name)
        seen_q, nz = set(), [False] * 4
        for corner in CORNERS:
            tab = {2 + i: (ALPHA if corner[i] == 0 else BETA)
                   for i in range(3)}
            vals = [psubst(p, tab) for p in polys]
            for idx, v in enumerate(vals):
                if v:
                    nz[idx] = True
            for v in vals[1:]:
                if v and v != psmul(eta, dpe2):
                    bad_num.append((name, pshow(v)))
            if vals[0]:
                hit = None
                for q in qs:
                    if vals[0] == psmul(eta, pmul(dpe2, q)):
                        hit = pshow(q)
                        break
                if hit is None:
                    bad_q.append((name, pshow(vals[0])))
                else:
                    seen_q.add(hit)
        if len(seen_q) != len(qs):
            bad_q.append((name, "Q realised %d of %d" % (len(seen_q),
                                                         len(qs))))
        if not all(nz):
            bad_nz.append(name)
    ck("table_polynomials_are_multiaffine_in_d", not bad_ma,
       "det A(d) and its 3 Cramer numerators have degree <= 1 in each d_i "
       "for all %d skeletons" % len(PAPER_TABLE)
       if not bad_ma else "violations: %s" % bad_ma)
    ck("table_corner_numerators_match_paper", not bad_num,
       "at all %d corners d_i in {alpha,beta}, every Cramer numerator is 0 "
       "or eta*(delta+eps)^2 with the paper's eta, for all %d skeletons"
       % (len(CORNERS), len(PAPER_TABLE))
       if not bad_num else "mismatches: %s" % bad_num[:3])
    ck("table_corner_determinants_match_Q_sets", not bad_q,
       "every corner determinant is 0 or eta*(delta+eps)^2*q with q in the "
       "paper's Q, and every listed q is attained" if not bad_q
       else "mismatches: %s" % bad_q[:3])
    ck("table_each_polynomial_nonzero_at_some_corner", not bad_nz,
       "for all %d skeletons the determinant and each of the 3 numerators is "
       "nonzero at at least one corner" % len(PAPER_TABLE) if not bad_nz
       else "all-zero polynomials for %s" % bad_nz)
    allq = [Q1, Q2, Q3, Q4, Q5]
    ck("table_Q_entries_are_positive_on_legal_region",
       all(cert(q) == 1 for q in allq),
       "all %d distinct q certified > 0: %s" % (len(allq),
                                                "; ".join(pshow(q)
                                                          for q in allq)))


def corner_argument(n, arcs):
    """Multiaffine corner argument for A(d) on the open box (alpha,beta)^n.

    Returns (multiaffine, all_nonzero_corner_values_share_one_sign, eta).
    If it returns (True, True, eta) then, because a multiaffine function is the
    convex combination of its corner values with strictly positive weights,
    det A(d) and every Cramer numerator have the strict sign eta at every
    interior point, so A(d)^{-1}1 > 0 there.
    """
    A = [[mono(2 + i) if i == j else
          (ALPHA if (j, i) in arcs else BETA) for j in range(n)]
         for i in range(n)]
    polys = [pdet(A)] + [pcramer(A, i) for i in range(n)]
    ma = not any(k[s] > 1 for p in polys for k in p
                 for s in range(2, 2 + n))
    signs, nz = set(), [False] * len(polys)
    corners = [tuple((m >> i) & 1 for i in range(n)) for m in range(1 << n)]
    for corner in corners:
        tab = {2 + i: (ALPHA if corner[i] == 0 else BETA) for i in range(n)}
        for idx, p in enumerate(polys):
            v = psubst(p, tab)
            if v:
                nz[idx] = True
                signs.add(cert(v))
    ok = (len(signs) == 1 and signs <= {1, -1} and all(nz))
    return ma, ok, (signs.pop() if len(signs) == 1 else None)


def check_corner_argument():
    rows = []
    for n in (2, 3):
        for arcs in iso_classes(n):
            if uniform_status(n, arcs) != "PERMITTED":
                continue
            ma, ok, eta = corner_argument(n, arcs)
            rows.append((n, sorted(arcs), ma, ok, eta))
    bad = [r for r in rows if not (r[2] and r[3] and r[4] in (1, -1))]
    ck("corner_argument_certifies_every_permitted_small_skeleton", not bad,
       "for each of the %d permitted skeletons on 2 or 3 vertices, det A(d) "
       "and all Cramer numerators are multiaffine and every nonzero corner "
       "value carries one certified sign eta (eta by row: %s), so "
       "A(d)^{-1}1 > 0 on the whole open box (alpha,beta)^n and no such "
       "skeleton can yield a counterexample"
       % (len(rows), [r[4] for r in rows]) if not bad
       else "inconclusive rows: %s" % bad[:2])


def check_box_positivity():
    """Direct: A(d)^{-1}1 > 0 for d in the open box, for the 6 skeletons."""
    ts = (Fraction(1, 5), Fraction(1, 2), Fraction(4, 5))
    bad = []
    total = 0
    for name, arcs, eta, _ in PAPER_TABLE:
        for eps, dl in legal_points():
            a, b = 1 - eps, 1 + dl
            for t1 in ts:
                for t2 in ts:
                    for t3 in ts:
                        d = [a + t * (b - a) for t in (t1, t2, t3)]
                        M = [[d[i] if i == j else
                              (a if (j, i) in arcs else b)
                              for j in range(3)] for i in range(3)]
                        S = rsolve(M, [Fraction(1)] * 3)
                        total += 1
                        if S is None or not all(v > 0 for v in S):
                            bad.append((name, eps, dl, t1, t2, t3))
    ck("three_component_reduction_always_positive", not bad,
       "%d exact solves (%d permitted skeletons x %d legal (eps,delta) x %d "
       "interior d): A(d)^{-1}1 > 0 every time, so no permitted 3-vertex "
       "skeleton yields a counterexample"
       % (total, len(PAPER_TABLE), len(legal_points()), len(ts) ** 3)
       if not bad
       else "failures at %s" % bad[:3])


def num_solve(n, arcs, eps, dl):
    """Exact solution of (I-W)x = 1 and whether it is strictly positive."""
    x = rsolve(imw_num(n, arcs, eps, dl), [Fraction(1)] * n)
    return x, (x is not None and all(v > 0 for v in x))


def digraphs(n):
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    for mask in range(1 << len(pairs)):
        yield {p for k, p in enumerate(pairs) if mask >> k & 1}


def size_vectors(total):
    out = []
    for nn in (2, 3):
        for sv in [v for v in _compositions(nn, total)]:
            out.append(sv)
    return out


def _compositions(parts, total):
    if parts == 1:
        for k in range(1, total + 1):
            yield (k,)
        return
    for k in range(1, total + 1):
        for rest in _compositions(parts - 1, total - k):
            if k + sum(rest) <= total:
                yield (k,) + rest


def check_reduction_lemma_finite_family():
    """Exhaustive exact test of the reduction lemma on small composites."""
    pts = [(Fraction(1, 5), Fraction(1)), (Fraction(1, 50), Fraction(1, 10)),
           (Fraction(9, 10), Fraction(19))]
    bad_bound, bad_crit, bad_sing = [], [], []
    count = 0
    orders = set()
    for sv in size_vectors(4):
        nn = len(sv)
        for sk in digraphs(nn):
            for combo in _component_combos(sv):
                for eps, dl in pts:
                    a, b = 1 - eps, 1 + dl
                    ds, okc = [], True
                    for m, carcs in zip(sv, combo):
                        x, pos = num_solve(m, carcs, eps, dl)
                        if x is None or not pos or sum(x) == 0:
                            okc = False
                            break
                        ds.append(1 / sum(x))
                    if not okc:
                        continue
                    count += 1
                    orders.add(max(sv))
                    if not all(a < d < b for d in ds):
                        bad_bound.append((sv, sorted(sk), eps, dl))
                    A = [[ds[i] if i == j else
                          (a if (j, i) in sk else b) for j in range(nn)]
                         for i in range(nn)]
                    S = rsolve(A, [Fraction(1)] * nn)
                    if S is None:
                        continue
                    nodes, garcs = _labelled_composite(sv, sk, combo)
                    xg, gpos = num_solve(len(nodes), garcs, eps, dl)
                    if xg is None:
                        bad_sing.append((sv, sorted(sk), eps, dl))
                        continue
                    if gpos != all(v > 0 for v in S):
                        bad_crit.append((sv, sorted(sk), eps, dl))
    SCOPE["lemma_family_max_component_order"] = max(orders) if orders else 0
    ck("reduction_lemma_holds_on_all_small_composites",
       not bad_bound and not bad_crit and not bad_sing,
       "%d exact instances (every composite of <= 4 vertices with 2 or 3 "
       "arbitrary permitted components, at %d legal parameter points; largest "
       "component order actually exercised = %d): "
       "alpha < d_i < beta always, I-W_G nonsingular whenever A(d) is, and "
       "G permitted exactly when A(d)^{-1}1 > 0"
       % (count, len(pts), SCOPE["lemma_family_max_component_order"])
       if not (bad_bound or bad_crit or bad_sing)
       else "bound %s / criterion %s / singular %s"
       % (bad_bound[:1], bad_crit[:1], bad_sing[:1]))


def _component_combos(sv):
    pools = [list(digraphs(m)) for m in sv]
    out = [[]]
    for pool in pools:
        out = [o + [g] for o in out for g in pool]
    return out


def _labelled_composite(sv, sk, combo):
    comps = {}
    nxt = 0
    for i, m in enumerate(sv):
        labels = list(range(nxt, nxt + m))
        nxt += m
        internal = {(labels[u], labels[v]) for u, v in combo[i]}
        comps[i] = (labels, internal)
    nodes, arcs = composite(list(range(len(sv))), sk, comps)
    return nodes, arcs


def check_five_vertex_family():
    """Sweep every 4-vertex skeleton with components (K1,K2,K1,K1)."""
    si = {v: i for i, v in enumerate(SK_NODES)}
    gi = {v: i for i, v in enumerate(G_NODES)}
    comps = {si[v]: ([gi[x] for x in COMPS[v][0]],
                     {(gi[x], gi[y]) for x, y in COMPS[v][1]})
             for v in SK_NODES}
    blocks = [comps[i][0] for i in range(len(SK_NODES))]
    target = frozenset((si[u], si[w]) for u, w in SK_ARCS)
    nsk = sum(1 for _ in digraphs(4))   # sweep size: counted, not assumed
    totals, hits, disagree = [], [], 0
    for eps, dl in [(Fraction(1, 5), Fraction(1)),
                    (Fraction(1, 50), Fraction(1, 10))]:
        found, npermit = [], 0
        a, b = 1 - eps, 1 + dl
        for sk in digraphs(4):
            xs, skok = num_solve(4, sk, eps, dl)
            if not skok:
                continue
            npermit += 1
            nodes, garcs = composite(list(range(len(SK_NODES))), sk, comps)
            pos = {v: k for k, v in enumerate(nodes)}
            xg, gok = num_solve(len(nodes),
                                {(pos[x], pos[y]) for x, y in garcs}, eps, dl)
            if xg is not None and not gok:
                found.append(frozenset(sk))
            ds = []
            for i in range(len(SK_NODES)):
                lab, ia = comps[i]
                r = {(lab.index(x), lab.index(y)) for x, y in ia}
                xc, _ = num_solve(len(lab), r, eps, dl)
                ds.append(1 / sum(xc))
            A = [[ds[i] if i == j else (a if (j, i) in sk else b)
                  for j in range(len(SK_NODES))]
                 for i in range(len(SK_NODES))]
            S = rsolve(A, [Fraction(1)] * len(SK_NODES))
            if xg is not None and S is not None:
                if [sum(xg[pos[j]] for j in blk) for blk in blocks] != S:
                    disagree += 1
        totals.append((npermit, len(found)))
        hits.append(target in found)
    ck("five_vertex_family_sweep_contains_the_exhibited_graph",
       all(hits) and disagree == 0 and all(t[1] > 0 for t in totals),
       "over all %d skeletons on 4 vertices with components "
       "(K1,K2,K1,K1): (permitted skeletons, forbidden composites) = %s at "
       "two legal points; the paper's skeleton is one of them at both, and "
       "the 4x4 reduction agrees with the direct 5x5 solve everywhere "
       "(%d disagreements)" % (nsk, totals, disagree))


def check_singleton_composites():
    bad = []
    count = 0
    for n in (2, 3, 4):
        nodes = list(range(n))
        comps = {v: ([v], set()) for v in nodes}
        pairs = [(i, j) for i in nodes for j in nodes if i != j]
        for mask in range(1 << len(pairs)):
            arcs = {p for k, p in enumerate(pairs) if mask >> k & 1}
            cn, ca = composite(nodes, arcs, comps)
            count += 1
            if cn != nodes or ca != arcs:
                bad.append((n, sorted(arcs)))
    orders = [len(COMPS[v][0]) for v in SK_NODES]
    nonsingleton = max(orders) >= 2
    ck("all_singleton_composite_equals_its_skeleton",
       not bad and nonsingleton,
       "%d skeletons on 2,3,4 vertices: replacing every vertex by K1 "
       "returns the skeleton unchanged, so a composite on <= 4 vertices "
       "with 4 nonempty components is its own permitted skeleton -- and the "
       "exhibited object escapes that case, its component orders being %s "
       "with a non-singleton" % (count, orders)
       if not bad and nonsingleton
       else "counterexamples: %s; exhibited component orders %s"
       % (bad[:3], orders))
    done = dict(CHECKS)
    need = ["census_one_vertex_skeleton", "census_two_vertex_permitted_skeletons",
            "census_three_vertex_permitted_skeletons",
            "corner_argument_certifies_every_permitted_small_skeleton",
            "three_component_reduction_always_positive",
            "reduction_lemma_holds_on_all_small_composites",
            "all_singleton_composite_equals_its_skeleton",
            "conclusion_G_is_forbidden"]
    missing = [k for k in need if not done.get(k)]
    ck("minimality_lower_bound_is_five",
       not missing and len(G_NODES) == 5 and sum(orders) == 5 and nonsingleton,
       "every prerequisite established above (%s): 1, 2 and 3 components are "
       "excluded, 4 components on <= 4 vertices force all-singleton hence a "
       "permitted composite, so >= 5 vertices; and the exhibited forbidden "
       "composite has exactly 5 -- but the exclusion of 2 and 3 components "
       "quantifies over components of ARBITRARY order and therefore leans on "
       "the paper's hand-proved reduction lemma, which is re-verified here "
       "only for components of order <= %s (see the NOT RE-RUN paragraph)"
       % (", ".join(need),
          SCOPE.get("lemma_family_max_component_order", "?"))
       if not missing else "prerequisites not established: %s" % missing)


def check_certificate_machinery():
    """The substitution must cover the legal region and only the legal region,
    and the certificate must be able to say 'negative' and 'unknown'."""
    into = []
    for dd in (Fraction(1, 7), Fraction(2), Fraction(50), Fraction(1000)):
        for vv in (Fraction(1, 900), Fraction(1, 9), Fraction(1),
                   Fraction(11), Fraction(4000)):
            e = pev(EPS_NUM, [dd, vv, 0, 0, 0]) / pev(EPS_DEN,
                                                      [dd, vv, 0, 0, 0])
            into.append(dd > 0 and 0 < e < dd / (dd + 1))
    onto = []
    for eps, dl in legal_points():
        t = eps * (dl + 1) / dl
        if not 0 < t < 1:
            onto.append(False)
            continue
        vv = t / (1 - t)
        e = pev(EPS_NUM, [dl, vv, 0, 0, 0]) / pev(EPS_DEN, [dl, vv, 0, 0, 0])
        onto.append(vv > 0 and e == eps)
    ck("certificate_substitution_covers_the_legal_region",
       all(into) and all(onto) and len(onto) == len(legal_points()),
       "(d,v) > 0 lands strictly inside the legal region at all %d tested "
       "pairs, and every one of the %d legal test points is hit exactly by "
       "d = delta, v = t/(1-t) with t = eps(delta+1)/delta"
       % (len(into), len(onto)))
    ck("certificate_reports_negative_and_unknown",
       cert(psub(EPS, DEL)) == -1 and cert(psub(DEL, const(1))) is None and
       cert(psub(pmul(EPS, BETA), DEL)) == -1 and cert({}) == 0,
       "cert(eps-delta)=-1, cert(delta-1)=None (indefinite, as it must be), "
       "cert(eps(delta+1)-delta)=-1, cert(0)=0")


def main():
    groups = [("object", check_object), ("matrices", check_matrices),
              ("determinants", check_determinants),
              ("hypotheses", check_hypotheses),
              ("conclusion", check_conclusion), ("grid", check_grid),
              ("projection", check_projection_pattern),
              ("nondegeneracy", check_nondegeneracy),
              ("reduction_lemma", check_reduction_lemma),
              ("census", check_census),
              ("certificate_cross_validation",
               check_certificates_against_evaluation),
              ("table_corners", check_table_corners),
              ("corner_argument", check_corner_argument),
              ("box_positivity", check_box_positivity),
              ("reduction_lemma_family", check_reduction_lemma_finite_family),
              ("five_vertex_family", check_five_vertex_family),
              ("singleton_composites", check_singleton_composites),
              ("certificate_machinery", check_certificate_machinery)]
    for name, fn in groups:
        try:
            fn()
        except Exception as exc:                      # a corrupted input
            ck("%s_section_completed" % name, False,
               "raised %s: %s" % (type(exc).__name__, exc))
    exhibited_order = max(len(COMPS[v][0]) for v in SK_NODES)
    fam_order = SCOPE.get("lemma_family_max_component_order")
    ord_txt = ("%d" % max(fam_order, exhibited_order)) if fam_order else \
        "unavailable, because the finite-family section did not complete"
    n_pl = len({_canon(p) for p in paper_principal_list()})
    n_cl = len({_canon(p) for p in paper_cramer_list()})
    print("NOT RE-RUN: no external catalogue or search is needed for this "
          "paper, and every census it reports (1-, 2- and 3-vertex skeletons, "
          "all corners of Table 1, the %d and %d determinant lists) is "
          "re-enumerated above. The reduction lemma's general form (the "
          "paper's Lemma 2), for components of unbounded order, is checked "
          "only on the finite family reported above: the largest component "
          "order exercised anywhere in this program is %s, and the exhibited "
          "counterexample's own components have order at most %d. Two "
          "consequences of that, spelled out because they are not visible in "
          "the check lines. (i) The minimality claim (Theorem 3: every "
          "counterexample has at least five vertices) is NOT fully "
          "machine-checked here. Its exclusion of two and three components "
          "quantifies over components of arbitrary order; what this program "
          "establishes is the censuses of permitted 1-, 2- and 3-vertex "
          "skeletons, the multiaffine corner argument giving A(d)^{-1}1 > 0 on "
          "the whole open box (alpha,beta)^n for every permitted 2- and "
          "3-vertex skeleton (the one-component case being tautological), "
          "and the reduction lemma itself only for components of order <= %s. "
          "So minimality rests, for components of larger order, on the paper's "
          "hand proof of that lemma. (ii) The total-activity bound "
          "beta^{-1} < 1^T(I-W_i)^{-1}1 < alpha^{-1}, which the paper quotes "
          "from Lemma 3.1 of the cited work in order to get "
          "alpha < d_i < beta, is likewise only re-verified here for "
          "components of order <= %s; it is not proved in general and this "
          "program does not check it beyond that finite family. The forward "
          "implication of the conjecture (the cited Theorem 5.2) is also not "
          "re-proved, being irrelevant to the refutation."
          % (n_pl, n_cl, ord_txt, exhibited_order, ord_txt, ord_txt))


def report():
    n = len(CHECKS)
    bad = [c for c, o in CHECKS if not o]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


if __name__ == "__main__":
    main()
    raise SystemExit(report())
