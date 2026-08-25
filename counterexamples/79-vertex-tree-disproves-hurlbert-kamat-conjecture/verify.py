#!/usr/bin/env python3
"""Verification of a 79-vertex tree refuting an 18-HK conjecture.

TAKEN FROM THE PAPER (inputs, transcribed):
  * the edge list of the exhibited tree T on 79 vertices;
  * the distinguished degree-three vertex and the claimed rank r = 18.

DERIVED HERE (computed by this program, nothing copied):
  * T is a tree (connected, 78 edges, acyclic), order 79;
  * no vertex of T has degree two;
  * independence polynomials of T and of vertex-deleted/closed-neighbourhood
    subtrees, by exact integer dynamic programming;
  * i_r(T, v) = number of independent r-sets containing v, for every vertex v;
  * the comparison at r = 18 between the degree-three vertex and every leaf;
  * each of the four asserted counts a SECOND time, by the paper's algebraic
    route: the closed coefficient formula for [z^m] S^q together with ordinary
    coefficient convolution, using binomial arithmetic only -- no graph, no
    traversal and none of the polynomial primitives the rooted DP uses;
  * mu(T), the minimum size of a maximal independent set of T, which fixes the
    rank window r <= mu(T)/2 that recurs in this literature and that this
    program cannot check the source conjecture against (see NOT RE-RUN).
"""

from fractions import Fraction

CHECKS = []

# ---------------------------------------------------------------- polynomials
# A polynomial is a list of exact integers, index = degree.


def p_trim(a):
    b = list(a)
    while len(b) > 1 and b[-1] == 0:
        b.pop()
    return b


def p_add(a, b):
    r = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        r[i] += x
    for i, x in enumerate(b):
        r[i] += x
    return p_trim(r)


def p_mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return p_trim(r)


def p_pow(a, k):
    r = [1]
    for _ in range(k):
        r = p_mul(r, a)
    return r


def cf(a, m):
    return a[m] if 0 <= m < len(a) else 0


def binom(n, k):
    if k < 0 or k > n or n < 0:
        return 0
    num = 1
    den = 1
    for i in range(k):
        num *= n - i
        den *= i + 1
    return num // den


# ------------------------------------------------- the object, as in the paper
# TAKEN FROM THE PAPER: branch sizes, the rank r, and the four asserted counts.
BRANCH_SIZES = (9, 8, 8)          # t_1, t_2, t_3
R = 18                            # the rank at which the conjecture fails
CLAIM_ORDER = 79                  # |V(T)|
CLAIM_CENTER = 1526033336050350   # |I_c^18(T)|
CLAIM_LEAF8 = 1524522945307798    # |I_{u_8}^18(T)|, asserted to be the leaf max
CLAIM_LEAF9 = 1521164183552880    # |I_{u_9}^18(T)|
CLAIM_GAP = 1510390742552         # asserted difference


def build_tree(ts):
    """Edge list of T: a center c, hubs h_i, supports s_ij, leaves l_ijk."""
    c = ("c",)
    verts = [c]
    edges = []
    for i, t in enumerate(ts, start=1):
        h = ("h", i)
        verts.append(h)
        edges.append((c, h))
        for j in range(1, t + 1):
            s = ("s", i, j)
            verts.append(s)
            edges.append((h, s))
            for k in (1, 2):
                lf = ("l", i, j, k)
                verts.append(lf)
                edges.append((s, lf))
    return verts, edges


def adjacency(verts, edges):
    adj = dict((v, set()) for v in verts)
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    return adj


def ind_poly_forest(verts, adj):
    """Independence polynomial of a forest, by exact rooted DP.

    Raises ValueError if the input has a cycle, so the DP can never be
    applied silently to a graph it does not handle.
    """
    total = [1]
    parent = {}
    for root in verts:
        if root in parent:
            continue
        parent[root] = None
        order = [root]
        stack = [root]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in parent:
                    parent[w] = v
                    order.append(w)
                    stack.append(w)
                elif parent[v] != w:
                    raise ValueError("cycle through %r-%r" % (v, w))
        out = {}
        inn = {}
        for v in reversed(order):
            po = [1]
            pi = [0, 1]
            for w in adj[v]:
                if parent[w] == v:
                    po = p_mul(po, p_add(out[w], inn[w]))
                    pi = p_mul(pi, out[w])
            out[v] = po
            inn[v] = pi
        total = p_mul(total, p_add(out[root], inn[root]))
    return total


def induced(verts, adj, removed):
    keep = [v for v in verts if v not in removed]
    ks = set(keep)
    sub = dict((v, set(w for w in adj[v] if w in ks)) for v in keep)
    return keep, sub


def brute_ind_profile(verts, adj):
    """Exhaustive enumeration: (independence polynomial, per-vertex counts).

    per_vertex[v][k] = number of independent k-sets containing v.
    Used only on a small graph, to validate the DP and the star identity.
    """
    n = len(verts)
    idx = dict((v, i) for i, v in enumerate(verts))
    block = [0] * n
    for v in verts:
        for w in adj[v]:
            block[idx[v]] |= 1 << idx[w]
    poly = [0] * (n + 1)
    per = [[0] * (n + 1) for _ in range(n)]
    chosen = []

    def rec(i, blocked):
        if i == n:
            poly[len(chosen)] += 1
            for a in chosen:
                per[a][len(chosen)] += 1
            return
        rec(i + 1, blocked)
        if not (blocked >> i) & 1:
            chosen.append(i)
            rec(i + 1, blocked | block[i])
            chosen.pop()

    rec(0, 0)
    return p_trim(poly), dict((v, per[idx[v]]) for v in verts)


def min_maximal_independent(verts, adj):
    """Minimum size of a maximal independent set of a forest, exactly.

    A maximal independent set is exactly an independent dominating set, so this
    is the independent domination number.  Rooted DP with three states per
    vertex v: (0) v is in the set; (1) v is not, and is dominated by one of its
    children; (2) v is not, and is not dominated inside its own subtree, so its
    parent must be in the set.  Arithmetic is integer throughout; INF = n + 1
    is a value no feasible solution can reach, so infeasible states can never
    win a minimum.  Raises ValueError on a cycle.
    """
    n = len(verts)
    INF = n + 1
    parent = {}
    total = 0
    for root in verts:
        if root in parent:
            continue
        parent[root] = None
        order = [root]
        stack = [root]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in parent:
                    parent[w] = v
                    order.append(w)
                    stack.append(w)
                elif parent[v] != w:
                    raise ValueError("cycle through %r-%r" % (v, w))
        f = {}
        for v in reversed(order):
            kids = [w for w in adj[v] if parent[w] == v]
            a = 1 + sum(min(f[w][1], f[w][2]) for w in kids)
            if kids:
                b = sum(min(f[w][0], f[w][1]) for w in kids)
                if not any(f[w][0] <= f[w][1] for w in kids):
                    # every child prefers state 1; force the cheapest into the
                    # set so that v is dominated.  Each delta is positive here.
                    b += min(f[w][0] - f[w][1] for w in kids)
            else:
                b = INF
            c = sum(f[w][1] for w in kids)
            f[v] = (min(a, INF), min(b, INF), min(c, INF))
        total += min(f[root][0], f[root][1])
    return total


def brute_min_maximal(verts, adj):
    """Minimum size of a maximal independent set by exhaustive enumeration.

    Used only on small graphs, to validate min_maximal_independent().
    """
    n = len(verts)
    idx = dict((v, i) for i, v in enumerate(verts))
    block = [0] * n
    for v in verts:
        for w in adj[v]:
            block[idx[v]] |= 1 << idx[w]
    closed = [block[i] | (1 << i) for i in range(n)]
    full = (1 << n) - 1
    best = [n + 1]
    chosen = []

    def rec(i, blocked, covered):
        if i == n:
            if covered == full and len(chosen) < best[0]:
                best[0] = len(chosen)
            return
        rec(i + 1, blocked, covered)
        if not (blocked >> i) & 1:
            chosen.append(i)
            rec(i + 1, blocked | block[i], covered | closed[i])
            chosen.pop()

    rec(0, 0, 0)
    return best[0]


def star_counts(verts, adj):
    """|I_v^r(G)| for every v and every r, via r-1 sets off N[v]."""
    out = {}
    for v in verts:
        closed = set(adj[v])
        closed.add(v)
        kv, ka = induced(verts, adj, closed)
        p = ind_poly_forest(kv, ka)
        out[v] = [0] + list(p)
    return out


STATE = {}


def T():
    if "T" not in STATE:
        verts, edges = build_tree(BRANCH_SIZES)
        STATE["T"] = (verts, edges, adjacency(verts, edges))
    return STATE["T"]


def stars():
    if "stars" not in STATE:
        verts, edges, adj = T()
        STATE["stars"] = star_counts(verts, adj)
    return STATE["stars"]


def leaves():
    verts, edges, adj = T()
    return [v for v in verts if len(adj[v]) == 1]


def mu_T():
    if "mu" not in STATE:
        verts, edges, adj = T()
        STATE["mu"] = min_maximal_independent(verts, adj)
    return STATE["mu"]


# ------------------------------------------------------------------ check 1-3
def chk_object():
    verts, edges, adj = T()
    n, m = len(verts), len(edges)
    if len(set(verts)) != n:
        return False, "repeated vertex label"
    if len(set(frozenset(e) for e in edges)) != m:
        return False, "repeated edge"
    try:
        ind_poly_forest(verts, adj)   # raises on any cycle
    except ValueError as exc:
        return False, "not acyclic: %s" % exc
    seen = set([verts[0]])
    stack = [verts[0]]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    ok = (n == CLAIM_ORDER and m == n - 1 and len(seen) == n)
    kinds = {}
    for v in verts:
        kinds[v[0]] = kinds.get(v[0], 0) + 1
    return ok, ("order=%d size=%d connected=%s centers=%d hubs=%d "
                "supports=%d leaves=%d") % (
        n, m, len(seen) == n, kinds.get("c", 0), kinds.get("h", 0),
        kinds.get("s", 0), kinds.get("l", 0))


def chk_no_degree_two():
    verts, edges, adj = T()
    spec = {}
    for v in verts:
        d = len(adj[v])
        spec[d] = spec.get(d, 0) + 1
    deg3 = 1 + sum(BRANCH_SIZES) + sum(1 for t in BRANCH_SIZES if t == 2)
    ok = 2 not in spec
    ok = ok and spec.get(1, 0) == 2 * sum(BRANCH_SIZES)
    ok = ok and spec.get(3, 0) == deg3
    hub_degs = sorted(len(adj[("h", i)]) for i in (1, 2, 3))
    ok = ok and hub_degs == sorted(t + 1 for t in BRANCH_SIZES)
    return ok, "degree spectrum %s, hub degrees %s" % (
        sorted(spec.items()), hub_degs)


def chk_center_degree():
    verts, edges, adj = T()
    d = len(adj[("c",)])
    return d == 3, "deg(c)=%d" % d


def chk_alpha():
    verts, edges, adj = T()
    p = ind_poly_forest(verts, adj)
    alpha = len(p) - 1
    return alpha >= R, "alpha(T)=%d, r=%d admissible (1<=r<=alpha)" % (alpha, R)


# -------------------------------------------- checks 5-6: validate the machine
def chk_dp_vs_brute():
    """The DP independence polynomial equals exhaustive enumeration."""
    for ts in ((2, 2, 1), (3, 0, 2), (1, 1, 1)):
        verts, edges = build_tree(ts)
        adj = adjacency(verts, edges)
        dp = ind_poly_forest(verts, adj)
        bf, _ = brute_ind_profile(verts, adj)
        if dp != bf:
            return False, "mismatch at t=%s" % (ts,)
    verts, edges = build_tree((2, 2, 1))
    adj = adjacency(verts, edges)
    return True, "n=%d member agrees, total independent sets %d" % (
        len(verts), sum(ind_poly_forest(verts, adj)))


def chk_star_identity():
    """|I_v^r(G)| = [z^{r-1}] P_{G-N[v]} for every vertex and every rank.

    Checked against exhaustive enumeration on a small member of the family.
    """
    verts, edges = build_tree((2, 2, 1))
    adj = adjacency(verts, edges)
    _, per = brute_ind_profile(verts, adj)
    pred = star_counts(verts, adj)
    tested = 0
    for v in verts:
        for r in range(1, len(verts) + 1):
            a = per[v][r] if r < len(per[v]) else 0
            b = pred[v][r] if r < len(pred[v]) else 0
            if a != b:
                return False, "v=%s r=%d brute=%d identity=%d" % (v, r, a, b)
            tested += 1
    return True, "%d (vertex, rank) pairs agree with enumeration" % tested


# ------------------------------- checks 7-10: the paper's algebraic ingredients
def path3_poly():
    vs = [0, 1, 2]
    aj = {0: set([1]), 1: set([0, 2]), 2: set([1])}
    return ind_poly_forest(vs, aj)


def branch_graph(t):
    """One hub, t supports, two leaves per support."""
    verts = [("h",)]
    edges = []
    for j in range(t):
        s = ("s", j)
        verts.append(s)
        edges.append((("h",), s))
        for k in (0, 1):
            lf = ("l", j, k)
            verts.append(lf)
            edges.append((s, lf))
    return verts, edges


def H_poly(t):
    return p_add(p_pow(path3_poly(), t),
                 p_mul([0, 1], p_pow([1, 1], 2 * t)))


def chk_S_polynomial():
    p = path3_poly()
    return p == [1, 3, 1], "P_{P_3}(z) coefficients %s" % (p,)


def chk_H_is_branch_poly():
    """H_t(z) = S(z)^t + z(1+z)^{2t} really is the branch's polynomial."""
    for t in range(0, 11):
        verts, edges = branch_graph(t)
        dp = ind_poly_forest(verts, adjacency(verts, edges))
        if dp != H_poly(t):
            return False, "t=%d: DP %s vs claimed %s" % (t, dp, H_poly(t))
    return True, "t=0..10 agree; H_9 has %d coefficients" % len(H_poly(9))


def chk_S_coeff_formula():
    """Equation (5) of the paper, against honest polynomial powers."""
    tested = 0
    for q in range(0, 26):
        pq = p_pow([1, 3, 1], q)
        for m in range(0, 2 * q + 3):
            s = 0
            j = max(0, m - q)
            while j <= m // 2:
                s += binom(q, j) * binom(q - j, m - 2 * j) * 3 ** (m - 2 * j)
                j += 1
            if s != cf(pq, m):
                return False, "q=%d m=%d formula=%d true=%d" % (
                    q, m, s, cf(pq, m))
            tested += 1
    return True, "%d coefficients of S^q, q<=25, match" % tested


def chk_H_coeff_formula():
    """[z^m] H_t = [z^m] S^t + C(2t, m-1)."""
    tested = 0
    for t in range(0, 11):
        ht = H_poly(t)
        st = p_pow([1, 3, 1], t)
        for m in range(0, 2 * t + 3):
            if cf(ht, m) != cf(st, m) + binom(2 * t, m - 1):
                return False, "t=%d m=%d" % (t, m)
            tested += 1
    return True, "%d coefficients match" % tested


def deleted_poly(v):
    verts, edges, adj = T()
    closed = set(adj[v])
    closed.add(v)
    kv, ka = induced(verts, adj, closed)
    return kv, ka, ind_poly_forest(kv, ka)


def components(verts, adj):
    seen = set()
    comps = []
    for r in verts:
        if r in seen:
            continue
        seen.add(r)
        stack = [r]
        comp = [r]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    comp.append(w)
                    stack.append(w)
        comps.append(comp)
    return comps


def chk_center_deletion():
    """Equation (3): T - N[c] is 25 disjoint 3-paths, P = S^25."""
    kv, ka, dp = deleted_poly(("c",))
    comps = components(kv, ka)
    sizes = sorted(len(c) for c in comps)
    ok = (len(kv) == CLAIM_ORDER - 4 and sizes == [3] * sum(BRANCH_SIZES)
          and dp == p_pow([1, 3, 1], sum(BRANCH_SIZES)))
    return ok, "%d vertices, %d components of size %s, P == S^%d: %s" % (
        len(kv), len(comps), set(sizes), sum(BRANCH_SIZES),
        dp == p_pow([1, 3, 1], sum(BRANCH_SIZES)))


def leaf_in_branch(size):
    for i, t in enumerate(BRANCH_SIZES, start=1):
        if t == size:
            return ("l", i, 1, 1)
    raise ValueError("no branch of size %d" % size)


def chk_leaf8_deletion():
    """Equation (4): P_{T-N[u_8]} = (1+z)(z S^24 + H_9 H_8 H_7)."""
    kv, ka, dp = deleted_poly(leaf_in_branch(8))
    S = [1, 3, 1]
    claimed = p_mul([1, 1], p_add(p_mul([0, 1], p_pow(S, 24)),
                                 p_mul(p_mul(H_poly(9), H_poly(8)),
                                       H_poly(7))))
    return dp == claimed, "graph on %d vertices; identity holds: %s" % (
        len(kv), dp == claimed)


def chk_leaf9_deletion():
    """Equation (4b): P_{T-N[u_9]} = (1+z)(z S^24 + H_8^3)."""
    kv, ka, dp = deleted_poly(leaf_in_branch(9))
    S = [1, 3, 1]
    claimed = p_mul([1, 1], p_add(p_mul([0, 1], p_pow(S, 24)),
                                 p_pow(H_poly(8), 3)))
    return dp == claimed, "graph on %d vertices; identity holds: %s" % (
        len(kv), dp == claimed)


# ----------------------------------------- checks 14-17: the numeric claims
def chk_center_count():
    val = stars()[("c",)][R]
    return val == CLAIM_CENTER, "|I_c^%d(T)| computed %d" % (R, val)


def chk_leaf_counts():
    st = stars()
    vals = {}
    for u in leaves():
        vals.setdefault(st[u][R], []).append(u)
    best = max(vals)
    ok = (len(vals) == 2 and best == CLAIM_LEAF8
          and min(vals) == CLAIM_LEAF9
          and len(vals[CLAIM_LEAF8]) == 2 * 8 * 2
          and len(vals[CLAIM_LEAF9]) == 2 * 9)
    return ok, "leaf orbits: %s" % sorted(
        (v, len(u)) for v, u in vals.items())


# ------------------------- the paper's algebraic route, computed independently
# These four functions evaluate the coefficients the paper evaluates, by the
# route the paper describes: the closed formula (5) for [z^m] S(z)^q, the
# closed formula for [z^m] H_t(z), and ordinary coefficient convolution.  They
# call binom() and nothing else -- no graph, no traversal, no p_mul/p_pow, no
# dynamic programming -- so agreement with star_counts() is a genuine
# two-route agreement rather than a restatement.


def s_coeff(q, m):
    """[z^m] S(z)^q from equation (5), binomial arithmetic only."""
    if m < 0:
        return 0
    s = 0
    j = max(0, m - q)
    while j <= m // 2:
        s += binom(q, j) * binom(q - j, m - 2 * j) * 3 ** (m - 2 * j)
        j += 1
    return s


def h_coeff(t, m):
    """[z^m] H_t(z) = [z^m] S(z)^t + C(2t, m-1), binomial arithmetic only."""
    if m < 0:
        return 0
    return s_coeff(t, m) + binom(2 * t, m - 1)


def alg_center_coeff():
    """[z^{R-1}] S(z)^{sum t_i}, i.e. equation (3) evaluated by formula."""
    return s_coeff(sum(BRANCH_SIZES), R - 1)


def alg_leaf_coeff(size):
    """[z^{R-1}] of (1+z)(z S^{q-1} + H_{a}H_{b}H_{c}), by formula.

    Deleting N[u] for a leaf u in a branch with `size` supports removes that
    leaf's support, so the three branch parameters are BRANCH_SIZES with one
    copy of `size` replaced by size-1: (9,8,7) for size 8 and (8,8,8) for
    size 9, matching equations (4) and (4b).
    """
    ts = list(BRANCH_SIZES)
    ts.remove(size)
    ts.append(size - 1)
    q = sum(BRANCH_SIZES) - 1

    def inner(m):
        """[z^m] (z S^q + H_{ts[0]} H_{ts[1]} H_{ts[2]})."""
        if m < 0:
            return 0
        v = s_coeff(q, m - 1)
        for i in range(m + 1):
            for j in range(m - i + 1):
                v += (h_coeff(ts[0], i) * h_coeff(ts[1], j)
                      * h_coeff(ts[2], m - i - j))
        return v

    # multiplication by (1+z)
    return inner(R - 1) + inner(R - 2)


def chk_alg_center():
    """Second route to |I_c^18|: formula (5) alone, against the rooted DP."""
    val = alg_center_coeff()
    dp = stars()[("c",)][R]
    ok = (val == CLAIM_CENTER and val == dp)
    return ok, ("[z^%d]S^%d by binomial formula = %d; rooted DP from the edge "
                "list = %d; paper = %d") % (
        R - 1, sum(BRANCH_SIZES), val, dp, CLAIM_CENTER)


def chk_alg_leaf8():
    """Second route to |I_{u_8}^18|: equation (4) by formula and convolution."""
    val = alg_leaf_coeff(8)
    dp = stars()[leaf_in_branch(8)][R]
    ok = (val == CLAIM_LEAF8 and val == dp)
    return ok, ("[z^%d](1+z)(zS^%d+H_9H_8H_7) by formula and convolution = %d; "
                "rooted DP = %d; paper = %d") % (
        R - 1, sum(BRANCH_SIZES) - 1, val, dp, CLAIM_LEAF8)


def chk_alg_leaf9():
    """Second route to |I_{u_9}^18|: equation (4b) by formula and convolution."""
    val = alg_leaf_coeff(9)
    dp = stars()[leaf_in_branch(9)][R]
    ok = (val == CLAIM_LEAF9 and val == dp)
    return ok, ("[z^%d](1+z)(zS^%d+H_8^3) by formula and convolution = %d; "
                "rooted DP = %d; paper = %d") % (
        R - 1, sum(BRANCH_SIZES) - 1, val, dp, CLAIM_LEAF9)


def chk_alg_gap():
    """The asserted gap, from the algebraic route only."""
    c = alg_center_coeff()
    l8 = alg_leaf_coeff(8)
    gap = c - l8
    ok = (gap == CLAIM_GAP and gap > 0)
    return ok, "gap by the algebraic route alone = %d - %d = %d (paper %d)" % (
        c, l8, gap, CLAIM_GAP)


def chk_mu_and_rank_window():
    """mu(T), the minimum size of a maximal independent set of T.

    The DP is validated against exhaustive enumeration on small members first.
    mu(T) is reported because the admissible-rank window of the conjecture is
    the one thing this program cannot check; see the NOT RE-RUN line.
    """
    for ts in ((2, 2, 1), (3, 0, 2), (1, 1, 1)):
        verts, edges = build_tree(ts)
        adj = adjacency(verts, edges)
        dp = min_maximal_independent(verts, adj)
        bf = brute_min_maximal(verts, adj)
        if dp != bf:
            return False, "t=%s: DP %d, enumeration %d" % (ts, dp, bf)
    mu = mu_T()
    ok = (mu == 1 + sum(BRANCH_SIZES))
    return ok, ("DP validated against enumeration on 3 small members; "
                "mu(T)=%d = |{c}| + #supports, mu(T)/2=%s, r=%d, so "
                "r > mu(T)/2: %s") % (
        mu, Fraction(mu, 2), R, Fraction(R) > Fraction(mu, 2))


def chk_refutation():
    """LOAD-BEARING: c beats every leaf at r = 18, computed from the graph."""
    st = stars()
    c = st[("c",)][R]
    leaf_vals = [st[u][R] for u in leaves()]
    best = max(leaf_vals)
    gap = c - best
    ok = (c > best and gap == CLAIM_GAP
          and all(c > x for x in leaf_vals) and len(leaf_vals) == 50)
    return ok, "c=%d  max leaf=%d  gap=%d over %d leaves" % (
        c, best, gap, len(leaf_vals))


def chk_no_leaf_is_max_center():
    """No leaf is a maximum 18-center: the argmax is a non-leaf."""
    verts, edges, adj = T()
    st = stars()
    top = max(st[v][R] for v in verts)
    winners = [v for v in verts if st[v][R] == top]
    leafset = set(leaves())
    ok = all(w not in leafset for w in winners)
    kinds = {}
    for w in winners:
        kinds[w[0]] = kinds.get(w[0], 0) + 1
    return ok, "max_v |I_v^%d| = %d attained at %s, degrees %s" % (
        R, top, sorted(kinds.items()),
        sorted(set(len(adj[w]) for w in winners)))


def chk_hubs_plus_leaves_independent():
    """The paper's witness for 18 <= alpha: the hubs and all 50 leaves form
    an independent set of size 53."""
    verts, edges, adj = T()
    wit = [v for v in verts if v[0] in ("h", "l")]
    bad = [(u, v) for u, v in edges if u in set(wit) and v in set(wit)]
    ok = (not bad) and len(wit) == 3 + 2 * sum(BRANCH_SIZES) and len(wit) >= R
    return ok, "independent set of size %d (3 hubs + %d leaves), edges inside: %d" % (
        len(wit), len(wit) - 3, len(bad))


def chk_hypotheses_of_conjecture():
    """Every hypothesis of the statement acted on: tree, no degree-2 vertex,
    r admissible, and the witness vertex has degree three."""
    verts, edges, adj = T()
    p = ind_poly_forest(verts, adj)
    conds = [
        ("tree", len(edges) == len(verts) - 1),
        ("no_degree_2", all(len(adj[v]) != 2 for v in verts)),
        ("r_admissible", 1 <= R <= len(p) - 1),
        ("witness_deg3", len(adj[("c",)]) == 3),
        ("witness_not_leaf", len(adj[("c",)]) > 1),
    ]
    ok = all(b for _, b in conds)
    return ok, ", ".join("%s=%s" % (a, b) for a, b in conds)


# ------------------------------------------- check 18: census over the family
def family_report(ts):
    """For one member of the family, the ranks r at which no leaf is a
    maximum r-center.  Only orbit representatives are evaluated."""
    verts, edges = build_tree(ts)
    adj = adjacency(verts, edges)
    reps = [("c",)]
    for i, t in enumerate(ts, start=1):
        reps.append(("h", i))
        if t:
            reps.append(("s", i, 1))
            reps.append(("l", i, 1, 1))
    poly = {}
    for v in reps:
        closed = set(adj[v])
        closed.add(v)
        kv, ka = induced(verts, adj, closed)
        poly[v] = [0] + list(ind_poly_forest(kv, ka))
    alpha = len(ind_poly_forest(verts, adj)) - 1
    bad = []
    for r in range(1, alpha + 1):
        vals = dict((v, poly[v][r] if r < len(poly[v]) else 0) for v in reps)
        top = max(vals.values())
        leafbest = max([vals[v] for v in reps if len(adj[v]) == 1] or [0])
        if leafbest < top:
            bad.append(r)
    return bad


def paper_member_bad_ranks():
    """The ranks at which the exhibited tree fails to be r-HK."""
    if "bad_ranks" not in STATE:
        STATE["bad_ranks"] = family_report(BRANCH_SIZES)
    return STATE["bad_ranks"]


def chk_family_census():
    """Search the two-level family for members that fail to be HK."""
    vals = [0] + list(range(3, 12))
    found = []
    total = 0
    for a in vals:
        for b in vals:
            if b > a:
                continue
            for c in vals:
                if c > b:
                    continue
                total += 1
                bad = family_report((a, b, c))
                if bad:
                    found.append((1 + 3 + 3 * (a + b + c), (a, b, c), bad))
    found.sort()
    hit = [f for f in found if sorted(f[1], reverse=True)
           == sorted(BRANCH_SIZES, reverse=True)]
    ok = (bool(hit) and R in hit[0][2] and bool(found))
    return ok, ("%d members with t_i in {0,3..11} searched, %d not HK; "
                "smallest order %d at t=%s ranks %s; paper's member found "
                "with failing ranks %s") % (
        total, len(found), found[0][0], found[0][1], found[0][2],
        hit[0][2] if hit else None)


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + str(detail) + "]"
    print(line)


def finish():
    n = len(CHECKS)
    bad = [c for c, o in CHECKS if not o]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


CHECK_LIST = [
    ("object_is_a_tree_on_79_vertices", chk_object),
    ("hypothesis_no_vertex_of_degree_two", chk_no_degree_two),
    ("hypothesis_witness_has_degree_three", chk_center_degree),
    ("hypothesis_rank_18_is_admissible", chk_alpha),
    ("hubs_and_leaves_form_an_independent_set", chk_hubs_plus_leaves_independent),
    ("all_hypotheses_of_the_conjecture", chk_hypotheses_of_conjecture),
    ("dp_agrees_with_exhaustive_enumeration", chk_dp_vs_brute),
    ("star_identity_agrees_with_enumeration", chk_star_identity),
    ("path3_independence_polynomial", chk_S_polynomial),
    ("branch_polynomial_H_t", chk_H_is_branch_poly),
    ("coefficient_formula_for_S_powers", chk_S_coeff_formula),
    ("coefficient_formula_for_H_t", chk_H_coeff_formula),
    ("center_deletion_identity", chk_center_deletion),
    ("leaf8_deletion_identity", chk_leaf8_deletion),
    ("leaf9_deletion_identity", chk_leaf9_deletion),
    ("center_18_star_count", chk_center_count),
    ("leaf_18_star_counts_both_orbits", chk_leaf_counts),
    ("second_route_center_18_by_formula_only", chk_alg_center),
    ("second_route_leaf8_18_by_formula_only", chk_alg_leaf8),
    ("second_route_leaf9_18_by_formula_only", chk_alg_leaf9),
    ("second_route_gap_by_formula_only", chk_alg_gap),
    ("refutation_center_beats_every_leaf", chk_refutation),
    ("no_leaf_is_a_maximum_18_center", chk_no_leaf_is_max_center),
    ("minimum_maximal_independent_set_of_T", chk_mu_and_rank_window),
    ("family_census_two_level_trees", chk_family_census),
]


def main():
    verts, edges, adj = T()
    print("object: two-level tree, branch sizes %s, %d vertices, %d edges"
          % (str(BRANCH_SIZES), len(verts), len(edges)))
    for name, fn in CHECK_LIST:
        try:
            ok, detail = fn()
        except Exception as exc:                       # noqa: BLE001
            ok, detail = False, "raised %s: %s" % (type(exc).__name__, exc)
        ck(name, ok, detail)
    print("scope: the exhibited tree, all its hypotheses, all four asserted "
          "counts and all three polynomial identities are recomputed from the "
          "edge list alone by rooted dynamic programming; each of the four "
          "asserted counts is then recomputed a second time by the paper's "
          "algebraic route -- formula (5) for [z^m]S^q, the closed form for "
          "[z^m]H_t, and coefficient convolution, binomial arithmetic only, "
          "touching neither the graph nor the polynomial primitives the DP "
          "uses -- and the two routes agree digit for digit.")
    try:
        derived = ("mu(T)=%d, so the exhibited failure rank r=%d exceeds "
                   "mu(T)/2=%s, and the ranks at which T fails to be r-HK are "
                   "exactly %s" % (mu_T(), R, Fraction(mu_T(), 2),
                                   paper_member_bad_ranks()))
    except Exception as exc:                               # noqa: BLE001
        derived = ("mu(T) and the failing ranks could not be derived (%s: %s)"
                   % (type(exc).__name__, exc))
    print("NOT RE-RUN: the statement being refuted. This program reads no "
          "literature, so it cannot confirm that the conjecture quoted in the "
          "paper -- every tree with no vertex of degree two is r-HK for every "
          "r with 1 <= r <= alpha -- is the source statement, and in "
          "particular it cannot rule out a restriction on r in the source. "
          "That gap is load-bearing here: " + derived + ". Under any reading "
          "that admits only r <= mu(T)/2, the threshold recurring in the "
          "Holroyd-Talbot literature, every count above remains correct and "
          "none of them contradicts the conjecture. Also not tested: "
          "minimality of the order 79 over all trees, which the paper does "
          "not claim, and the census sweep, which covers only branch sizes in "
          "{0,3,...,11}.")
    return finish()


if __name__ == "__main__":
    raise SystemExit(main())
