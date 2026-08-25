#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- referee verification program for

    "The Gaetz-Gao Closure Conjecture for Finite Coxeter Systems",
    Theorem 1 / Proposition 5, on Conjecture 1.3 of Gaetz-Gao.

Python 3 standard library only.  No floating point is used in any decision:
all arithmetic is exact in Z[phi] = Z[(1+sqrt5)/2], represented by integer
pairs (a,b) <-> a + b*phi with the single rule phi^2 = phi + 1.

--------------------------------------------------------------------------
VALUES TAKEN FROM THE PAPER (inputs; never compared to themselves)
--------------------------------------------------------------------------
  (i)  the Coxeter matrix of type H_4:  m_12 = 5, m_23 = m_34 = 3,
       m_ij = 2 otherwise  (Section 2);
  (ii) the claimed census of Proposition 5:
         |W| = 14400, max Coxeter length = 60,
         BP membership tests = 230400, positive memberships = 78498,
         unordered distinct BP pairs = 201991, closure failures = 0;
  (iii)the claimed |BP_W(w)| profile
         3:2492 4:62 5:8872 6:150 7:65 8:24 9:2627 10:68 12:24 16:16;
  (iv) the claimed Poincare polynomial factors [2]_q[12]_q[20]_q[30]_q.
Everything in (ii)-(iv) is a claim that this program RECOMPUTES from (i)
alone and then compares.  Only (i) -- three integers of the Coxeter
diagram -- is trusted.

--------------------------------------------------------------------------
DERIVED HERE (the checks)
--------------------------------------------------------------------------
  * the four 4x4 reflection matrices over Z[phi], and the exact braid
    relations (s_i s_j) has order exactly m_ij, verified power by power;
  * the whole group by breadth-first search: its order, its Cayley graph,
    left and right multiplication tables, Coxeter lengths, left/right
    descent sets, supports (by two independent routes);
  * the length distribution, compared with [2]_q[12]_q[20]_q[30]_q;
  * all 16 standard parabolic subgroups W_J, their orders by two
    independent routes, and |W^J| * |W_J| = |W|;
  * for all 14400 x 16 = 230400 pairs (w,J): the descent-stripping
    decomposition of Lemma 4, its independence of the choice of descent
    (three different choice rules), w^J w_J = w exactly, w^J in W^J,
    w_J in W_J, and l(w) = l(w^J) + l(w_J);
  * the BP test Supp(w^J) cap J subset D_L(w_J) for all 230400 pairs;
  * the SAME 230400 decompositions and the SAME BP families recomputed by a
    second, opposite-direction algorithm (each coset u W_J grown forwards
    from its minimal representative instead of stripping descents off w),
    with the positive count and the pair count recounted by different loops;
  * closure of every BP_W(w) under union and intersection over all
    unordered pairs of distinct members;
  * the paper's own arithmetic consistency: the three identities the proof
    asserts about its own table are evaluated on the paper's numbers and
    compared with the independently derived counts;
  * Lemma 2 (direct product) on an explicit instance I_2(5) x A_2;
  * Theorem 1 itself, verified outright for 25 further finite Coxeter
    systems: A_1..A_7, D_4, D_5, D_6, E_6, H_2, H_3, and the reducible
    systems A_1xH_4, A_2xH_4 and I_2(5)xH_4 (144000 elements, 27131854
    pairs), plus every dihedral I_2(m) with m <= 200 built in a second,
    purely combinatorial model that is cross-checked against the matrix
    model at m = 5;
  * three self-tests showing the load-bearing checks can fail: the derived
    order changes when the input diagram changes; the BP predicate is not
    constant AND is sensitive to the exact criterion (replacing D_L by D_R
    moves the |BP_W(w)| profile, which is what makes the profile check a
    real discriminator -- the positive count and the closure verdict alone
    are unmoved by that substitution); and the closure detector fires on a
    deliberately damaged family.

NOT verifiable here: Theorem 1 as a whole is NOT re-derived.  The
Gaetz-Gao theorem for finite Weyl groups, which the paper cites rather
than reproves, is spot-checked above but not re-derived in general, so the
general type-A/D/E cases and the components B_n, F_4, E_7, E_8 of Theorem 1
continue to rest on that citation; the types needing 2cos(pi/4) (B, F) fall
outside Z[phi], and E_7, E_8 and I_2(m) for m > 200 are out of budget.
What IS re-derived in full is the H_4 census of Proposition 5.

Exit status is 0 if and only if every check passes.
"""

import sys
from collections import deque

_RESULTS = []


def check(name, ok, detail=""):
    """Record and print one check.  Every call site must be falsifiable."""
    ok = bool(ok)
    _RESULTS.append((name, ok))
    line = ("PASS " if ok else "FAIL ") + name
    if detail:
        line += " [" + detail + "]"
    print(line)
    return ok


def note(msg):
    print("NOTE " + msg)


def verdict():
    n = len(_RESULTS)
    bad = [nm for nm, ok in _RESULTS if not ok]
    if bad:
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(bad), n))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % n)
    return 0


# =====================================================================
# INPUTS TAKEN FROM THE PAPER
# =====================================================================

# (i) Coxeter diagram of H_4, Section 2.  Indices 0..3 are s_1..s_4.
H4_M = [[1, 5, 2, 2],
        [5, 1, 3, 2],
        [2, 3, 1, 3],
        [2, 2, 3, 1]]

# (ii)/(iii) the claimed census of Proposition 5.
PAPER_ORDER = 14400
PAPER_MAXLEN = 60
PAPER_TESTS = 230400
PAPER_POSITIVE = 78498
PAPER_PAIRS = 201991
PAPER_FAILURES = 0
PAPER_PROFILE = {3: 2492, 4: 62, 5: 8872, 6: 150, 7: 65,
                 8: 24, 9: 2627, 10: 68, 12: 24, 16: 16}
# (iv) claimed Poincare polynomial [2]_q [12]_q [20]_q [30]_q.
PAPER_DEGREES = (2, 12, 20, 30)

# Scope knob for the extra (beyond-the-paper) dihedral sweep.
DIHEDRAL_MAX = 200


# =====================================================================
# EXACT ARITHMETIC IN Z[phi],  phi^2 = phi + 1
# =====================================================================

ZERO = (0, 0)
ONE = (1, 0)
PHI = (0, 1)


def zmul(x, y):
    """(a+b phi)(c+d phi) = (ac+bd) + (ad+bc+bd) phi, using phi^2 = phi+1."""
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


def zadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def two_cos_pi_over_m(m):
    """2*cos(pi/m) as an element of Z[phi]; defined for m in {2,3,5}."""
    if m == 2:
        return ZERO
    if m == 3:
        return ONE
    if m == 5:
        return PHI             # 2*cos(pi/5) = phi exactly
    raise ValueError("2cos(pi/%d) is not in Z[phi]" % m)


# =====================================================================
# n x n MATRICES OVER Z[phi], stored flat: entry (r,c) at 2*(r*n+c)
# =====================================================================

def mat_identity(n):
    out = [0] * (2 * n * n)
    for r in range(n):
        out[2 * (r * n + r)] = 1
    return tuple(out)


def mat_get(M, n, r, c):
    k = 2 * (r * n + c)
    return (M[k], M[k + 1])


def mat_mul(A, B, n):
    """Honest n x n product over Z[phi] via zmul/zadd.  Used for the braid
    relations and as the reference for the sparse fast paths."""
    out = [0] * (2 * n * n)
    for r in range(n):
        for c in range(n):
            acc = ZERO
            for k in range(n):
                acc = zadd(acc, zmul(mat_get(A, n, r, k),
                                     mat_get(B, n, k, c)))
            j = 2 * (r * n + c)
            out[j] = acc[0]
            out[j + 1] = acc[1]
    return tuple(out)


def gen_matrix(M_cox, n, i):
    """Matrix of the simple reflection s_i in the basis alpha_1..alpha_n:
       column c is the image of alpha_c."""
    out = [0] * (2 * n * n)
    for c in range(n):
        if c == i:
            out[2 * (i * n + i)] = -1          # s_i(alpha_i) = -alpha_i
        else:
            out[2 * (c * n + c)] = 1           # alpha_c
            a, b = two_cos_pi_over_m(M_cox[i][c])
            k = 2 * (i * n + c)
            out[k] = a                          # + 2cos(pi/m_ic) alpha_i
            out[k + 1] = b
    return tuple(out)


def right_mul_gen(M, n, i, coefs):
    """M * rho(s_i), exploiting the sparsity of rho(s_i)."""
    out = list(M)
    for r in range(n):
        base = 2 * r * n
        ki = base + 2 * i
        ai = M[ki]
        bi = M[ki + 1]
        out[ki] = -ai
        out[ki + 1] = -bi
        for c in range(n):
            if c == i:
                continue
            cc, dd = coefs[c]
            if cc == 0 and dd == 0:
                continue
            kc = base + 2 * c
            out[kc] = M[kc] + cc * ai + dd * bi
            out[kc + 1] = M[kc + 1] + cc * bi + dd * ai + dd * bi
    return tuple(out)


def left_mul_gen(M, n, i, coefs):
    """rho(s_i) * M, exploiting the sparsity of rho(s_i)."""
    out = list(M)
    base = 2 * i * n
    for c in range(n):
        kc = base + 2 * c
        a = -M[kc]
        b = -M[kc + 1]
        for k in range(n):
            if k == i:
                continue
            cc, dd = coefs[k]
            if cc == 0 and dd == 0:
                continue
            j = 2 * (k * n + c)
            x = M[j]
            y = M[j + 1]
            a += cc * x + dd * y
            b += cc * y + dd * x + dd * y
        out[kc] = a
        out[kc + 1] = b
    return tuple(out)


# =====================================================================
# BUILD A FINITE COXETER GROUP FROM ITS COXETER MATRIX
# =====================================================================

def build_coxeter(M_cox, cap=200000):
    """Breadth-first enumeration of the geometric representation.
    Returns a dict with multiplication tables, lengths, descents, supports.
    Raises if the group exceeds `cap` (guards against a non-finite input)."""
    n = len(M_cox)
    coefs = [[two_cos_pi_over_m(M_cox[i][c]) if c != i else ZERO
              for c in range(n)] for i in range(n)]
    gens = [gen_matrix(M_cox, n, i) for i in range(n)]
    ident = mat_identity(n)

    index = {ident: 0}
    elems = [ident]
    length = [0]
    supp = [0]
    right = [[-1] * n]
    q = deque([0])
    while q:
        w = q.popleft()
        Mw = elems[w]
        lw = length[w]
        for i in range(n):
            Mn = right_mul_gen(Mw, n, i, coefs[i])
            j = index.get(Mn)
            if j is None:
                j = len(elems)
                if j > cap:
                    raise RuntimeError("group larger than cap")
                index[Mn] = j
                elems.append(Mn)
                length.append(lw + 1)
                supp.append(supp[w] | (1 << i))
                right.append([-1] * n)
                q.append(j)
            right[w][i] = j

    N = len(elems)
    left = [[-1] * n for _ in range(N)]
    for w in range(N):
        Mw = elems[w]
        for i in range(n):
            left[w][i] = index[left_mul_gen(Mw, n, i, coefs[i])]

    dl = [0] * N
    dr = [0] * N
    for w in range(N):
        lw = length[w]
        a = 0
        b = 0
        for i in range(n):
            if length[right[w][i]] < lw:
                a |= 1 << i
            if length[left[w][i]] < lw:
                b |= 1 << i
        dr[w] = a
        dl[w] = b

    return {"n": n, "M": M_cox, "coefs": coefs, "gens": gens,
            "elems": elems, "index": index, "N": N,
            "length": length, "supp": supp,
            "right": right, "left": left, "dl": dl, "dr": dr}


# =====================================================================
# CHECKS ON THE REPRESENTATION AND THE CAYLEY GRAPH
# =====================================================================

def check_coxeter_relations(G, tag):
    n, M, gens = G["n"], G["M"], G["gens"]
    ident = mat_identity(n)
    ok_sq = all(mat_mul(gens[i], gens[i], n) == ident and gens[i] != ident
                for i in range(n))
    check("%s.gens_are_involutions" % tag, ok_sq)
    ok_ord = True
    detail = ""
    for i in range(n):
        for j in range(i + 1, n):
            m = M[i][j]
            P = mat_mul(gens[i], gens[j], n)
            X = P
            for k in range(1, m + 1):
                is_id = (X == ident)
                if k < m and is_id:
                    ok_ord = False
                    detail = "s%ds%d has order %d < %d" % (i + 1, j + 1, k, m)
                if k == m and not is_id:
                    ok_ord = False
                    detail = "s%ds%d order != %d" % (i + 1, j + 1, m)
                X = mat_mul(X, P, n)
    check("%s.braid_orders_exact" % tag, ok_ord, detail)


def check_sparse_matches_dense(G, tag, step=1):
    """The fast sparse generator multiplications must equal honest 4x4
    products with rho(s_i).  `step` = 1 means EVERY element is checked
    against the dense reference (the sparse path is the only engine that
    builds the group, so partial coverage would leave it unaudited)."""
    n, gens, elems = G["n"], G["gens"], G["elems"]
    coefs, index, right, left = G["coefs"], G["index"], G["right"], G["left"]
    bad = 0
    tested = 0
    for w in range(0, G["N"], step):
        Mw = elems[w]
        for i in range(n):
            tested += 1
            if right_mul_gen(Mw, n, i, coefs[i]) != mat_mul(Mw, gens[i], n):
                bad += 1
            if left_mul_gen(Mw, n, i, coefs[i]) != mat_mul(gens[i], Mw, n):
                bad += 1
            if elems[right[w][i]] != mat_mul(Mw, gens[i], n):
                bad += 1
            if elems[left[w][i]] != mat_mul(gens[i], Mw, n):
                bad += 1
    check("%s.sparse_products_match_dense" % tag, bad == 0 and tested > 0,
          "%d sampled elements, %d mismatches" % (tested, bad))


def check_cayley_graph(G, tag):
    n, N = G["n"], G["N"]
    length, right, left = G["length"], G["right"], G["left"]
    bad_inv = bad_step = bad_assoc = 0
    for w in range(N):
        lw = length[w]
        rw = right[w]
        lw_row = left[w]
        for i in range(n):
            u = rw[i]
            if right[u][i] != w:
                bad_inv += 1
            if abs(length[u] - lw) != 1:
                bad_step += 1
            v = lw_row[i]
            if left[v][i] != w:
                bad_inv += 1
            if abs(length[v] - lw) != 1:
                bad_step += 1
            for j in range(n):
                if left[right[w][j]][i] != right[left[w][i]][j]:
                    bad_assoc += 1
    check("%s.generators_act_as_involutions" % tag, bad_inv == 0,
          "%d violations" % bad_inv)
    # NOTE.  `bad_step == 0` alone cannot fail here: every generator matrix is
    # the identity with one row replaced, so det rho(s_i) = -1, the parity of
    # the BFS distance is forced, and a BFS edge always changes it by exactly
    # one.  The falsifiable content is bundled in below: only the identity may
    # have an empty descent set.  That is what makes descent stripping stop
    # exactly at w^J (in particular w^S = e), and a corrupted multiplication
    # table does break it.
    no_dr = [w for w in range(N) if G["dr"][w] == 0]
    no_dl = [w for w in range(N) if G["dl"][w] == 0]
    check("%s.cayley_graph_graded_and_descents_nonempty" % tag,
          bad_step == 0 and no_dr == [0] and no_dl == [0],
          "%d edges not changing length by 1; %d elements with empty D_R, "
          "%d with empty D_L (only the identity is allowed)"
          % (bad_step, len(no_dr), len(no_dl)))
    check("%s.left_and_right_actions_commute" % tag, bad_assoc == 0,
          "%d violations" % bad_assoc)
    # the unique longest element is the unique element with all descents
    mx = max(length)
    fullmask = (1 << n) - 1
    tops = [w for w in range(N) if length[w] == mx]
    allr = [w for w in range(N) if G["dr"][w] == fullmask]
    alll = [w for w in range(N) if G["dl"][w] == fullmask]
    check("%s.unique_longest_element_has_all_descents" % tag,
          len(tops) == 1 and allr == tops and alll == tops,
          "l=%d count %d, |D_R|=n count %d, |D_L|=n count %d"
          % (mx, len(tops), len(allr), len(alll)))


def q_int(k):
    """[k]_q = 1 + q + ... + q^{k-1} as a coefficient list."""
    return [1] * k


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


# =====================================================================
# PARABOLIC SUBGROUPS AND SUPPORTS
# =====================================================================

def parabolic_data(G, tag):
    """For every J subset S build W_J by BFS inside the Cayley graph, then
    cross-check its order against {w : Supp(w) subset J}, cross-check
    Supp itself against the smallest J containing w, and check
    |W^J| * |W_J| = |W|."""
    n, N = G["n"], G["N"]
    right, left, supp, dr = G["right"], G["left"], G["supp"], G["dr"]
    full = (1 << n) - 1
    members = {}
    for J in range(1 << n):
        seen = {0}
        q = deque([0])
        while q:
            w = q.popleft()
            for i in range(n):
                if not (J >> i) & 1:
                    continue
                u = right[w][i]
                if u not in seen:
                    seen.add(u)
                    q.append(u)
        members[J] = seen

    bad_supp_route = []
    bad_closed = 0
    bad_index = []
    for J in range(1 << n):
        seen = members[J]
        by_supp = sum(1 for w in range(N) if supp[w] & ~J == 0)
        if len(seen) != by_supp:
            bad_supp_route.append((J, len(seen), by_supp))
        for w in seen:
            for i in range(n):
                if (J >> i) & 1 and (right[w][i] not in seen
                                     or left[w][i] not in seen):
                    bad_closed += 1
        cosets = sum(1 for w in range(N) if dr[w] & J == 0)
        if cosets * len(seen) != N:
            bad_index.append((J, cosets, len(seen)))

    check("%s.parabolic_order_two_routes_agree" % tag, not bad_supp_route,
          "%d of %d subsets disagree" % (len(bad_supp_route), 1 << n))
    check("%s.parabolic_subgroups_closed" % tag, bad_closed == 0,
          "%d violations" % bad_closed)
    check("%s.minimal_coset_reps_times_parabolic_equals_order" % tag,
          not bad_index, "%d of %d subsets fail" % (len(bad_index), 1 << n))

    smallest = [full] * N
    for J in range(1 << n):
        for w in members[J]:
            smallest[w] &= J
    bad = sum(1 for w in range(N) if smallest[w] != supp[w])
    check("%s.support_two_routes_agree" % tag, bad == 0,
          "%d of %d elements disagree" % (bad, N))
    G["members"] = members
    return members


# =====================================================================
# LEMMA 4: DESCENT STRIPPING  w = w^J w_J
# =====================================================================

def strip(G, w, J, rule):
    """Repeatedly remove a right descent lying in J.  `rule` selects which
    descent: 0 = lowest index, 1 = highest index, 2 = deterministic
    pseudo-random.  Returns (w^J, w_J, list of chosen generators)."""
    right, left, dr = G["right"], G["left"], G["dr"]
    x = w
    p = 0
    word = []
    while True:
        d = dr[x] & J
        if d == 0:
            return x, p, word
        if rule == 0:
            i = (d & -d).bit_length() - 1
        elif rule == 1:
            i = d.bit_length() - 1
        else:
            bits = [b for b in range(G["n"]) if (d >> b) & 1]
            h = (w * 2654435761 + J * 40503 + len(word) * 97 + 12345) % 1000003
            i = bits[h % len(bits)]
        x = right[x][i]
        p = left[p][i]
        word.append(i)


def bp_scan(G, tag, check_alt_rules=True, emit=True):
    """For every (w,J): verify Lemma 4 and evaluate the BP condition
    Supp(w^J) cap J subset D_L(w_J).  Returns (bp, tests, positives)."""
    n, N = G["n"], G["N"]
    length, right, supp, dl, dr = (G["length"], G["right"], G["supp"],
                                   G["dl"], G["dr"])
    nsub = 1 << n
    bp = [0] * N
    tests = 0
    positives = 0
    bad_add = bad_lenp = bad_min = bad_par = bad_prod = bad_alt = 0
    for J in range(nsub):
        pairs = set()
        for w in range(N):
            tests += 1
            x, p, word = strip(G, w, J, 0)
            if length[x] + length[p] != length[w]:
                bad_add += 1
            if length[p] != len(word):
                bad_lenp += 1
            if dr[x] & J:
                bad_min += 1
            if supp[p] & ~J:
                bad_par += 1
            # reconstruct w^J * w_J from the reduced word of w_J
            cur = x
            for i in reversed(word):
                nxt = right[cur][i]
                if length[nxt] != length[cur] + 1:
                    bad_prod += 1
                cur = nxt
            if cur != w:
                bad_prod += 1
            if check_alt_rules:
                for rule in (1, 2):
                    x2, p2, _ = strip(G, w, J, rule)
                    if (x2, p2) != (x, p):
                        bad_alt += 1
            pairs.add((x, p))
            if supp[x] & J & ~dl[p] == 0:
                bp[w] |= 1 << J
                positives += 1
        if len(pairs) != N:
            bad_min += 1000000
    viol = bad_add + bad_lenp + bad_min + bad_par + bad_prod + bad_alt
    if emit:
        check("%s.length_additivity" % tag, bad_add == 0,
              "%d violations" % bad_add)
        check("%s.parabolic_factor_length_equals_steps" % tag, bad_lenp == 0,
              "%d violations" % bad_lenp)
        check("%s.wJ_is_minimal_coset_rep_and_bijective" % tag, bad_min == 0,
              "%d violations" % bad_min)
        check("%s.w_J_lies_in_W_J" % tag, bad_par == 0,
              "%d violations" % bad_par)
        check("%s.exact_product_wJ_times_w_J_equals_w" % tag, bad_prod == 0,
              "%d violations" % bad_prod)
        if check_alt_rules:
            check("%s.stripping_independent_of_descent_choice" % tag,
                  bad_alt == 0,
                  "%d violations over 2 alternative rules" % bad_alt)
    return bp, tests, positives, viol


def forward_coset_decomposition(G):
    """SECOND, INDEPENDENT route to the parabolic decomposition w = w^J w_J.

    `strip` works backwards: it peels right descents lying in J off w.  This
    routine works forwards: for each J it starts from every minimal coset
    representative u (an element with D_R(u) cap J = empty) and grows the coset
    u W_J by right multiplication by generators of J along strictly
    length-increasing edges, carrying the W_J-part along.  Every element v of
    u W_J is reached exactly once, with v = u t and l(v) = l(u) + l(t) true by
    construction.

    Returns a list indexed by J of pairs (rep, part) with rep[w] = w^J and
    part[w] = w_J, or -1 where the forward growth failed to reach w (which is
    itself a detectable defect).  No value produced by `strip`, `bp_scan` or
    `closure_scan` is read here."""
    n, N = G["n"], G["N"]
    right, length, dr = G["right"], G["length"], G["dr"]
    out = []
    for J in range(1 << n):
        Jb = [i for i in range(n) if (J >> i) & 1]
        rep = [-1] * N
        part = [-1] * N
        for u in range(N):
            if dr[u] & J:
                continue
            rep[u] = u
            part[u] = 0
            q = deque([(u, 0)])
            while q:
                v, t = q.popleft()
                lv = length[v]
                for i in Jb:
                    v2 = right[v][i]
                    if length[v2] == lv + 1 and rep[v2] < 0:
                        t2 = right[t][i]
                        rep[v2] = u
                        part[v2] = t2
                        q.append((v2, t2))
        out.append((rep, part))
    return out


# =====================================================================
# CLOSURE OF BP_W(w) UNDER UNION AND INTERSECTION
# =====================================================================

def closure_scan(G, bp, tag, emit=True):
    """Test every unordered pair of distinct members of every BP_W(w)."""
    n, N = G["n"], G["N"]
    full = (1 << n) - 1
    pairs = 0
    fail_union = 0
    fail_inter = 0
    wit_u = ""
    wit_i = ""
    no_empty = 0
    no_full = 0
    no_supp = 0
    supp = G["supp"]
    for w in range(N):
        fam = bp[w]
        if not (fam >> 0) & 1:
            no_empty += 1
        if not (fam >> full) & 1:
            no_full += 1
        # J = Supp(w) is the falsifiable companion of J = empty set: it forces
        # the stripping to reduce w completely inside W_{Supp(w)}, i.e.
        # w^{Supp(w)} = e.  (J = empty set on its own cannot fail for any
        # input, since the condition then reads  0 subset D_L(e).)
        if not (fam >> supp[w]) & 1:
            no_supp += 1
        mem = [J for J in range(1 << n) if (fam >> J) & 1]
        for a in range(len(mem)):
            Ja = mem[a]
            for b in range(a + 1, len(mem)):
                Jb = mem[b]
                pairs += 1
                if not (fam >> (Ja | Jb)) & 1:
                    fail_union += 1
                    if not wit_u:
                        wit_u = "first: w=%d J=%d K=%d" % (w, Ja, Jb)
                if not (fam >> (Ja & Jb)) & 1:
                    fail_inter += 1
                    if not wit_i:
                        wit_i = "first: w=%d J=%d K=%d" % (w, Ja, Jb)
    if emit:
        check("%s.empty_set_and_support_always_in_BP" % tag,
              no_empty == 0 and no_supp == 0,
              "%d elements lack the empty set, %d lack Supp(w)"
              % (no_empty, no_supp))
        check("%s.full_set_always_in_BP" % tag, no_full == 0,
              "%d elements lack S" % no_full)
        check("%s.BP_closed_under_union" % tag, fail_union == 0,
              ("%d failures, %s" % (fail_union, wit_u)) if fail_union
              else "%d pairs tested" % pairs)
        check("%s.BP_closed_under_intersection" % tag, fail_inter == 0,
              ("%d failures, %s" % (fail_inter, wit_i)) if fail_inter
              else "%d pairs tested" % pairs)
    return pairs, fail_union, fail_inter, no_empty + no_full + no_supp


# =====================================================================
# THE DIHEDRAL FAMILY I_2(m), BUILT COMBINATORIALLY (any m >= 2)
# =====================================================================

def build_dihedral(m):
    """I_2(m) with elements the alternating words; independent of the
    Z[phi] matrix model, so it also serves as a cross-check at m = 5."""
    def canon(x, k):
        if k == 0:
            return (0, 0)
        if k == m:
            return (0, m)
        return (x, k)

    labels = [(0, 0)]
    for k in range(1, m):
        labels.append((0, k))
        labels.append((1, k))
    labels.append((0, m))
    idx = dict((L, i) for i, L in enumerate(labels))
    N = len(labels)
    length = [L[1] for L in labels]
    supp = [0 if L[1] == 0 else ((1 << L[0]) if L[1] == 1 else 3)
            for L in labels]
    right = [[-1, -1] for _ in range(N)]
    left = [[-1, -1] for _ in range(N)]
    for L in labels:
        x, k = L
        w = idx[L]
        for i in (0, 1):
            if k == 0:
                right[w][i] = idx[canon(i, 1)]
                left[w][i] = idx[canon(i, 1)]
            elif k == m:
                xr = i if m % 2 == 1 else 1 - i
                right[w][i] = idx[canon(xr, m - 1)]
                left[w][i] = idx[canon(1 - i, m - 1)]
            else:
                last = x if k % 2 == 1 else 1 - x
                right[w][i] = idx[canon(x, k - 1) if i == last
                                  else canon(x, k + 1)]
                left[w][i] = idx[canon(1 - x, k - 1) if i == x
                                 else canon(i, k + 1)]
    dl = [0] * N
    dr = [0] * N
    for w in range(N):
        for i in (0, 1):
            if length[left[w][i]] < length[w]:
                dl[w] |= 1 << i
            if length[right[w][i]] < length[w]:
                dr[w] |= 1 << i
    return {"n": 2, "M": [[1, m], [m, 1]], "N": N, "labels": labels,
            "length": length, "supp": supp, "right": right, "left": left,
            "dl": dl, "dr": dr}


# =====================================================================
# LEMMA 2 (DIRECT PRODUCT) ON AN EXPLICIT INSTANCE
# =====================================================================

def embed_parabolic(G, sub, Gsub):
    """Simultaneous BFS identifying the standard parabolic of G on the
    generators `sub` (a list of indices of S) with the standalone group
    Gsub built on the induced Coxeter matrix.  Returns dict G-index -> Gsub-index."""
    m = {0: 0}
    q = deque([0])
    while q:
        w = q.popleft()
        for a, i in enumerate(sub):
            u = G["right"][w][i]
            v = Gsub["right"][m[w]][a]
            if u not in m:
                m[u] = v
                q.append(u)
            elif m[u] != v:
                return None
    return m


def check_product_lemma():
    """W = I_2(5) x A_2 : BP_W(w) must equal BP_{W_1}(w_1) x BP_{W_2}(w_2)."""
    Mprod = [[1, 5, 2, 2], [5, 1, 2, 2], [2, 2, 1, 3], [2, 2, 3, 1]]
    G = build_coxeter(Mprod)
    G1 = build_coxeter([[1, 5], [5, 1]])
    G2 = build_coxeter([[1, 3], [3, 1]])
    ok_order = (G["N"] == G1["N"] * G2["N"] == 60)
    check("product.order_factorises", ok_order,
          "|W|=%d |W1|=%d |W2|=%d" % (G["N"], G1["N"], G2["N"]))
    e1 = embed_parabolic(G, [0, 1], G1)
    e2 = embed_parabolic(G, [2, 3], G2)
    check("product.parabolic_embeddings_well_defined",
          e1 is not None and e2 is not None and len(e1) == G1["N"]
          and len(e2) == G2["N"])
    if e1 is None or e2 is None:
        return
    bp, _, _, v0 = bp_scan(G, "product", check_alt_rules=False, emit=False)
    bp1, _, _, v1 = bp_scan(G1, "product.f1", check_alt_rules=False, emit=False)
    bp2, _, _, v2 = bp_scan(G2, "product.f2", check_alt_rules=False, emit=False)
    check("product.lemma_4_holds_in_all_three_groups", v0 + v1 + v2 == 0,
          "%d violations" % (v0 + v1 + v2))
    bad = 0
    for w in range(G["N"]):
        x1, w1, _ = strip(G, w, 0b0011, 0)   # w1 = w_{S_1}, x1 = w^{S_1}
        if G["supp"][x1] & 0b0011:
            bad += 1
            continue
        a = e1[w1]
        b = e2[x1]
        for J in range(16):
            lhs = (bp[w] >> J) & 1
            rhs = ((bp1[a] >> (J & 3)) & 1) and ((bp2[b] >> (J >> 2)) & 1)
            if lhs != rhs:
                bad += 1
    check("product.BP_factorises_componentwise", bad == 0,
          "%d of %d (w,J) disagree" % (bad, 16 * G["N"]))


# =====================================================================
# POINCARE POLYNOMIAL
# =====================================================================

def check_poincare(G, tag, degrees):
    mx = max(G["length"])
    dist = [0] * (mx + 1)
    for lw in G["length"]:
        dist[lw] += 1
    poly = [1]
    for d in degrees:
        poly = poly_mul(poly, q_int(d))
    check("%s.length_distribution_equals_Poincare_polynomial" % tag,
          dist == poly,
          "derived degree %d vs %d, first mismatch %s"
          % (len(dist) - 1, len(poly) - 1,
             next((str(i) for i in range(min(len(dist), len(poly)))
                   if dist[i] != poly[i]), "none")))
    prod = 1
    for d in degrees:
        prod *= d
    check("%s.order_equals_product_of_degrees" % tag, G["N"] == prod,
          "derived %d vs %d" % (G["N"], prod))
    check("%s.max_length_equals_number_of_reflections" % tag,
          mx == sum(d - 1 for d in degrees),
          "derived %d vs %d" % (mx, sum(d - 1 for d in degrees)))
    return dist


# =====================================================================
# AGGREGATED SWEEP OVER THE DIHEDRAL FAMILY I_2(m)
# =====================================================================

def check_dihedral_family(mmax):
    bad_order = bad_struct = bad_lemma = bad_close = 0
    tot_pairs = 0
    for m in range(2, mmax + 1):
        D = build_dihedral(m)
        if D["N"] != 2 * m or max(D["length"]) != m:
            bad_order += 1
        for w in range(D["N"]):
            for i in (0, 1):
                if D["right"][D["right"][w][i]][i] != w:
                    bad_struct += 1
                if D["left"][D["left"][w][i]][i] != w:
                    bad_struct += 1
                if abs(D["length"][D["right"][w][i]] - D["length"][w]) != 1:
                    bad_struct += 1
                for j in (0, 1):
                    if (D["left"][D["right"][w][j]][i]
                            != D["right"][D["left"][w][i]][j]):
                        bad_struct += 1
        bp, _, _, v = bp_scan(D, "I2", check_alt_rules=True, emit=False)
        bad_lemma += v
        pr, fu, fi, fe = closure_scan(D, bp, "I2", emit=False)
        tot_pairs += pr
        bad_close += fu + fi + fe
    check("I2m.groups_have_order_2m_and_diameter_m", bad_order == 0,
          "%d of %d values of m fail" % (bad_order, mmax - 1))
    check("I2m.multiplication_tables_consistent", bad_struct == 0,
          "%d violations" % bad_struct)
    check("I2m.lemma_4_holds", bad_lemma == 0, "%d violations" % bad_lemma)
    check("I2m.BP_closed_for_all_m_up_to_%d" % mmax, bad_close == 0,
          "%d closure failures over %d pairs" % (bad_close, tot_pairs))


# =====================================================================
# OTHER FINITE COXETER TYPES REALISABLE OVER Z[phi]  (m in {2,3,5})
# =====================================================================
# Orders below come from the classification (Humphreys, Section 2.11),
# not from the paper under review; they are compared with derived values.
OTHER_TYPES = [
    ("A1", [[1]], 2),
    ("A1xA1", [[1, 2], [2, 1]], 4),
    ("A2", [[1, 3], [3, 1]], 6),
    ("H2=I2(5)", [[1, 5], [5, 1]], 10),
    ("A1xA1xA1", [[1, 2, 2], [2, 1, 2], [2, 2, 1]], 8),
    ("A2xA1", [[1, 3, 2], [3, 1, 2], [2, 2, 1]], 12),
    ("I2(5)xA1", [[1, 5, 2], [5, 1, 2], [2, 2, 1]], 20),
    ("A3", [[1, 3, 2], [3, 1, 3], [2, 3, 1]], 24),
    ("H3", [[1, 5, 2], [5, 1, 3], [2, 3, 1]], 120),
    ("A2xA2", [[1, 3, 2, 2], [3, 1, 2, 2], [2, 2, 1, 3], [2, 2, 3, 1]], 36),
    ("I2(5)xI2(5)", [[1, 5, 2, 2], [5, 1, 2, 2], [2, 2, 1, 5],
                     [2, 2, 5, 1]], 100),
    ("A4", [[1, 3, 2, 2], [3, 1, 3, 2], [2, 3, 1, 3], [2, 2, 3, 1]], 120),
    ("D4", [[1, 3, 3, 3], [3, 1, 2, 2], [3, 2, 1, 2], [3, 2, 2, 1]], 192),
    ("H3xA1", [[1, 5, 2, 2], [5, 1, 3, 2], [2, 3, 1, 2], [2, 2, 2, 1]], 240),
]


def check_other_finite_types():
    bad_order = []
    bad_lemma = 0
    bad_close = 0
    tot_pairs = 0
    tot_tests = 0
    for name, M, order in OTHER_TYPES:
        G = build_coxeter(M)
        if G["N"] != order:
            bad_order.append("%s:%d!=%d" % (name, G["N"], order))
        bp, tests, _, v = bp_scan(G, name, check_alt_rules=True, emit=False)
        bad_lemma += v
        tot_tests += tests
        pr, fu, fi, fe = closure_scan(G, bp, name, emit=False)
        tot_pairs += pr
        bad_close += fu + fi + fe
    check("other_types.orders_match_classification", not bad_order,
          ",".join(bad_order) if bad_order
          else "%d types" % len(OTHER_TYPES))
    check("other_types.lemma_4_holds", bad_lemma == 0,
          "%d violations over %d (w,J) pairs" % (bad_lemma, tot_tests))
    check("other_types.BP_closed_under_union_and_intersection",
          bad_close == 0,
          "%d failures over %d unordered pairs" % (bad_close, tot_pairs))


def diagram(n, edges):
    """Coxeter matrix on n generators; `edges` is a list of (i,j,m)."""
    M = [[1 if i == j else 2 for j in range(n)] for i in range(n)]
    for i, j, m in edges:
        M[i][j] = M[j][i] = m
    return M


# Larger finite Coxeter systems realisable over Z[phi] (all m in {2,3,5}).
# Orders are classification values (Humphreys), not paper values.
LARGE_TYPES = [
    ("A5", diagram(5, [(0, 1, 3), (1, 2, 3), (2, 3, 3), (3, 4, 3)]), 720),
    ("D5", diagram(5, [(0, 1, 3), (0, 2, 3), (0, 3, 3), (3, 4, 3)]), 1920),
    ("I2(5)xH3", diagram(5, [(0, 1, 5), (2, 3, 5), (3, 4, 3)]), 1200),
    ("A6", diagram(6, [(0, 1, 3), (1, 2, 3), (2, 3, 3), (3, 4, 3),
                       (4, 5, 3)]), 5040),
    ("H3xH3", diagram(6, [(0, 1, 5), (1, 2, 3), (3, 4, 5), (4, 5, 3)]), 14400),
    ("D6", diagram(6, [(0, 1, 3), (0, 2, 3), (0, 3, 3), (3, 4, 3),
                       (4, 5, 3)]), 23040),
    ("E6", diagram(6, [(0, 1, 3), (1, 2, 3), (2, 3, 3), (3, 4, 3),
                       (2, 5, 3)]), 51840),
    ("A7", diagram(7, [(0, 1, 3), (1, 2, 3), (2, 3, 3), (3, 4, 3),
                       (4, 5, 3), (5, 6, 3)]), 40320),
    ("A1xH4", diagram(5, [(0, 1, 5), (1, 2, 3), (2, 3, 3)]), 28800),
    ("A2xH4", diagram(6, [(0, 1, 3), (2, 3, 5), (3, 4, 3), (4, 5, 3)]), 86400),
    ("I2(5)xH4", diagram(6, [(0, 1, 5), (2, 3, 5), (3, 4, 3), (4, 5, 3)]),
     144000),
]


def check_large_types():
    """Direct verification of Theorem 1 for further finite types, including
    two reducible systems that actually contain the H_4 factor."""
    bad_order = []
    bad_lemma = 0
    bad_close = 0
    tot_pairs = 0
    tot_tests = 0
    for name, M, order in LARGE_TYPES:
        G = build_coxeter(M)
        if G["N"] != order:
            bad_order.append("%s:%d!=%d" % (name, G["N"], order))
        bp, tests, _, v = bp_scan(G, name, check_alt_rules=True, emit=False)
        pr, fu, fi, fe = closure_scan(G, bp, name, emit=False)
        bad_lemma += v
        bad_close += fu + fi + fe
        tot_tests += tests
        tot_pairs += pr
        print("OBJECT type=%s |W|=%d tests=%d pairs=%d failures=%d"
              % (name, G["N"], tests, pr, fu + fi + fe))
        del G, bp
    check("large_types.orders_match_classification", not bad_order,
          ",".join(bad_order) if bad_order else "%d types" % len(LARGE_TYPES))
    check("large_types.lemma_4_holds", bad_lemma == 0,
          "%d violations over %d (w,J) pairs" % (bad_lemma, tot_tests))
    check("large_types.BP_closed_under_union_and_intersection",
          bad_close == 0,
          "%d failures over %d unordered pairs" % (bad_close, tot_pairs))


def check_two_models_agree():
    """The combinatorial dihedral model and the Z[phi] matrix model must
    produce identical BP families for I_2(5)."""
    D = build_dihedral(5)
    G = build_coxeter([[1, 5], [5, 1]])
    bpD, _, _, _ = bp_scan(D, "x", check_alt_rules=False, emit=False)
    bpG, _, _, _ = bp_scan(G, "x", check_alt_rules=False, emit=False)
    kD = dict(((D["length"][w], D["dl"][w], D["dr"][w], D["supp"][w]), w)
              for w in range(D["N"]))
    kG = dict(((G["length"][w], G["dl"][w], G["dr"][w], G["supp"][w]), w)
              for w in range(G["N"]))
    ok_keys = (len(kD) == D["N"] == len(kG) == G["N"]
               and set(kD) == set(kG))
    bad = 0
    if ok_keys:
        for k in kD:
            if bpD[kD[k]] != bpG[kG[k]]:
                bad += 1
    check("crossmodel.I2_5_two_independent_models_agree",
          ok_keys and bad == 0,
          "keys ok=%s, %d elements disagree" % (ok_keys, bad))


# =====================================================================
# SELF-TESTS: the load-bearing checks must be able to FAIL
# =====================================================================

def selftests(G, bp, positives, tests):
    # (a) the pipeline is sensitive to the input diagram: turning m_12 = 5
    #     into 3 makes the diagram the A_4 chain, of order 120.
    Mbad = [row[:] for row in H4_M]
    Mbad[0][1] = Mbad[1][0] = 3
    Gbad = build_coxeter(Mbad)
    check("selftest.order_depends_on_the_input_diagram",
          Gbad["N"] == 120 and Gbad["N"] != G["N"],
          "m_12=3 gives |W|=%d (A_4), m_12=5 gives %d"
          % (Gbad["N"], G["N"]))
    # (b) the BP predicate is not constant, AND the program is sensitive to
    #     the exact criterion.  Replacing D_L(w_J) by D_R(w_J) -- the classic
    #     way to get a Coxeter-theoretic condition subtly wrong -- leaves both
    #     the positive count (78498) and the closure verdict (0 failures)
    #     untouched in H_4, so those two numbers alone do NOT pin the
    #     definition; only the |BP_W(w)| profile does.  Require here that the
    #     wrong criterion really does move the profile, so that
    #     census.bp_size_profile_matches_paper is a genuine discriminator.
    def _profile(fam):
        d = {}
        for w in range(G["N"]):
            k = bin(fam[w]).count("1")
            d[k] = d.get(k, 0) + 1
        return d

    Gswap = dict(G)
    Gswap["dl"] = list(G["dr"])
    bp_rd, _, pos_rd, _ = bp_scan(Gswap, "selftest.DR",
                                  check_alt_rules=False, emit=False)
    moved = _profile(bp_rd) != _profile(bp)
    check("selftest.BP_predicate_nonconstant_and_left_right_sensitive",
          0 < positives < tests and moved,
          "%d positives of %d tests; D_L -> D_R gives %d positives and a %s "
          "profile" % (positives, tests, pos_rd,
                       "different" if moved else "IDENTICAL"))
    # (c) the closure detector fires on a deliberately damaged family.
    full = (1 << G["n"]) - 1
    victim = next((w for w in range(G["N"])
                   if bp[w] == (1 << (full + 1)) - 1), None)
    fired = False
    if victim is not None:
        dmg = list(bp)
        dmg[victim] = bp[victim] & ~(1 << 0b0011)   # drop J = {s_1,s_2}
        _, fu, fi, _ = closure_scan(G, dmg, "damaged", emit=False)
        fired = fu > 0 and fi > 0
    check("selftest.closure_detector_fires_on_damaged_family", fired,
          "victim=%s" % victim)


# =====================================================================
# THE PAPER'S CENSUS
# =====================================================================

def check_census(G, bp, tests, positives, pairs, failures):
    n, N = G["n"], G["N"]
    mx = max(G["length"])
    profile = {}
    for w in range(N):
        k = bin(bp[w]).count("1")
        profile[k] = profile.get(k, 0) + 1
    print("OBJECT elements=%d max_length=%d tests=%d positives=%d "
          "pairs=%d failures=%d" % (N, mx, tests, positives, pairs, failures))
    print("OBJECT profile " + " ".join("%d:%d" % (k, profile[k])
                                       for k in sorted(profile)))
    check("census.elements_enumerated_is_14400", N == PAPER_ORDER,
          "derived %d, paper %d" % (N, PAPER_ORDER))
    check("census.maximum_coxeter_length_is_60", mx == PAPER_MAXLEN,
          "derived %d, paper %d" % (mx, PAPER_MAXLEN))
    check("census.bp_membership_tests_is_230400",
          tests == PAPER_TESTS and tests == N * (1 << n),
          "derived %d = %d x %d, paper %d"
          % (tests, N, 1 << n, PAPER_TESTS))
    check("census.positive_bp_memberships_is_78498",
          positives == PAPER_POSITIVE,
          "derived %d, paper %d" % (positives, PAPER_POSITIVE))
    check("census.unordered_distinct_bp_pairs_is_201991",
          pairs == PAPER_PAIRS,
          "derived %d, paper %d" % (pairs, PAPER_PAIRS))
    check("census.union_intersection_failures_is_zero",
          failures == PAPER_FAILURES,
          "derived %d, paper %d" % (failures, PAPER_FAILURES))
    check("census.bp_size_profile_matches_paper", profile == PAPER_PROFILE,
          "derived %s" % (" ".join("%d:%d" % (k, profile[k])
                                   for k in sorted(profile))))
    # The paper asserts three arithmetic identities about ITS table; they are
    # evaluated on the paper's numbers and compared with derived quantities.
    p0 = sum(PAPER_PROFILE.values())
    p1 = sum(k * c for k, c in PAPER_PROFILE.items())
    p2 = sum(k * (k - 1) // 2 * c for k, c in PAPER_PROFILE.items())
    check("census.paper_profile_sums_to_group_order", p0 == N,
          "paper table sums to %d, derived |W| = %d" % (p0, N))
    check("census.paper_profile_weighted_by_size_gives_positives",
          p1 == positives, "paper %d vs derived %d" % (p1, positives))
    check("census.paper_profile_weighted_by_binomial_gives_pair_count",
          p2 == pairs, "paper %d vs derived %d" % (p2, pairs))
    # A GENUINELY INDEPENDENT recomputation of the entire BP scan.
    #
    # (Comparing popcounts of `bp` with the running counter `positives` would
    # be a tautology: bp_scan sets bit J of bp[w] on exactly the iterations on
    # which it increments `positives`, and closure_scan increments `pairs` once
    # per unordered pair of set bits, so the two "accumulators" are the same
    # arithmetic and could not disagree for any input.  So instead:)
    #
    # Route 2 recomputes (w^J, w_J) by growing every coset u W_J forwards from
    # its minimal representative -- the opposite direction to the descent
    # stripping of Lemma 4 -- re-evaluates the BP predicate on those pairs,
    # and recounts the pairs with an ordered double loop and a halving instead
    # of a triangular loop.  It reads none of bp / positives / pairs.
    fwd = forward_coset_decomposition(G)
    supp, dl = G["supp"], G["dl"]
    bp2 = [0] * N
    pos2 = 0
    unassigned = 0
    for J in range(1 << n):
        rep, part = fwd[J]
        for w in range(N):
            if rep[w] < 0 or part[w] < 0:
                unassigned += 1
                continue
            if supp[rep[w]] & J & ~dl[part[w]] == 0:
                bp2[w] |= 1 << J
                pos2 += 1
    ordered = 0
    for w in range(N):
        mem = [J for J in range(1 << n) if (bp2[w] >> J) & 1]
        for Ja in mem:
            for Jb in mem:
                if Ja != Jb:
                    ordered += 1
    pairs2 = ordered // 2
    differ = sum(1 for w in range(N) if bp2[w] != bp[w])
    check("census.independent_accumulators_agree",
          unassigned == 0 and differ == 0 and pos2 == positives
          and pairs2 == pairs,
          "forward-coset route: %d positives, %d pairs, %d families differ, "
          "%d (w,J) unreached; stripping route: %d positives, %d pairs"
          % (pos2, pairs2, differ, unassigned, positives, pairs))
    return profile


# =====================================================================
# MAIN
# =====================================================================

def main():
    """Wrapper: an exception anywhere must still produce a FAIL and a
    VERDICT line, never a silent crash."""
    try:
        return _run()
    except Exception:
        import traceback
        traceback.print_exc()
        check("internal.verification_completed_without_error", False,
              "unhandled exception, see traceback on stderr")
        return verdict()


def _run():
    print("verify.py -- The Gaetz-Gao Closure Conjecture for Finite "
          "Coxeter Systems")
    note("input from the paper: H_4 Coxeter matrix m_12=5, m_23=m_34=3, "
         "m_ij=2 otherwise")
    note("exact arithmetic in Z[phi], phi^2 = phi+1; no floating point")

    G = build_coxeter(H4_M)
    check_coxeter_relations(G, "H4")
    check_sparse_matches_dense(G, "H4")
    check_cayley_graph(G, "H4")
    check_poincare(G, "H4", PAPER_DEGREES)
    parabolic_data(G, "H4")

    bp, tests, positives, _ = bp_scan(G, "H4", check_alt_rules=True)
    pairs, fu, fi, fe = closure_scan(G, bp, "H4")
    check_census(G, bp, tests, positives, pairs, fu + fi + fe)

    check_product_lemma()
    check_other_finite_types()
    check_two_models_agree()
    check_dihedral_family(DIHEDRAL_MAX)
    selftests(G, bp, positives, tests)
    check_large_types()

    extra_orders = dict((nm, od) for nm, _M, od in LARGE_TYPES)
    n_extra = len(OTHER_TYPES) + len(LARGE_TYPES)
    note("SCOPE: the paper's H_4 census (Proposition 5) is re-derived in "
         "full (%d elements, %d BP tests, %d pairs, %d failures)."
         % (G["N"], tests, pairs, fu + fi + fe))
    note("Theorem 1 is additionally verified outright for %d further "
         "finite Coxeter systems, among them I_2(5)xH_4 (%d elements) "
         "and A_2xH_4 (%d), plus I_2(m) for every m <= %d (which includes "
         "G_2 = I_2(6))."
         % (n_extra, extra_orders["I2(5)xH4"], extra_orders["A2xH4"],
            DIHEDRAL_MAX))
    note("NOT re-derived here: Theorem 1 as a whole.  Its general "
         "type-A/D/E cases and its components B_n, F_4, E_7, E_8 rest on "
         "the Gaetz-Gao theorem for finite Weyl groups, which the paper "
         "cites and this program does not reprove; only A_1..A_7, D_4..D_6, "
         "E_6 and products of these are re-verified.")
    note("NOT re-derived here: I_2(m) for m > %d, the types B/F needing "
         "2cos(pi/4) (outside Z[phi]), E_7 and E_8."
         % DIHEDRAL_MAX)
    return verdict()


if __name__ == "__main__":
    sys.exit(main())
