#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- checks every computational claim of

    "For odd primes p > q with q dividing p^2 - 1, every nonassociative right Bol
     loop of order pq has right multiplication group of order p^2 q"
     (an affirmative answer to Problem 7.1 of Kinyon,
     Nagy and Vojtechovsky, J. Algebra 473 (2017) 481-512)

against the objects printed in paper.tex.  Python 3.9+, STANDARD LIBRARY ONLY.
All arithmetic is exact integer arithmetic in F_p and F_{p^2}; permutation
groups are enumerated as byte strings.  There is no floating point anywhere and
no randomness anywhere, so the output is byte-stable.

The paper's decisive argument is a hand proof (Steps 1-5 of Section 3) valid for
all admissible (p,q) at once.  This program does NOT prove it; it checks it:

  * Section A re-derives, from the normal form printed in the paper, the
    admissible parameter vectors u for the small cells the paper exhibits, and
    checks they are exactly the vectors printed there.
  * Section B builds the multiplication tables and checks they are right Bol
    loops by the identity, not by fiat.
  * Section C runs the census over the ASSUMED theta-family (plus Z_{pq}) for
    pq <= 255, with |Mlt_r| obtained by brute-force closure of the right
    translations, independently of the proof.  That the family is exhaustive is
    imported from KNV Theorem 1.1(ii)-(iii) and is not checked here.
  * Section D checks each individual step of the proof as an identity on finitely
    many cells -- not symbolically for all p and q.
  * Section E runs the controls, including a FORCED POSITIVE: an object on which
    the decider must return the upper branch p^3 q, so that the 490 negatives
    mean something.
  * Section F sweeps the dimension formula far past the brute-force range. Its
    first two checks, F1 and F2, verify cardinality formulas for Gamma and for
    the nonassociative members per cell that the paper does NOT state; they are
    kept, and labelled as such, rather than deleted.

Run:   python3 verify.py
Exits 0 iff every check passed.
"""

import sys
from fractions import Fraction  # exactness guard only; no float decisions anywhere

# ---------------------------------------------------------------------------
# check harness
# ---------------------------------------------------------------------------
_PASSES = []
_FAILS = []


def check(name, ok, detail=''):
    if ok:
        _PASSES.append(name)
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _FAILS.append(name)
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))


def note(text):
    print('NOTE %s' % text)


# ---------------------------------------------------------------------------
# arithmetic
# ---------------------------------------------------------------------------
def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def inv_mod(a, p):
    a %= p
    if a == 0:
        raise ZeroDivisionError('inverse of 0 mod %d' % p)
    return pow(a, p - 2, p)


def nonresidue(p):
    """Least quadratic nonresidue mod p (p an odd prime). Deterministic."""
    for t in range(2, p):
        if pow(t, (p - 1) // 2, p) == p - 1:
            return t
    raise RuntimeError('no nonresidue mod %d' % p)


class F2(object):
    """F_{p^2} = F_p[z]/(z^2 - t), t the least nonresidue. Elements are pairs
    (a, b) meaning a + b z. Frobenius is (a, b) -> (a, -b) because
    z^p = t^((p-1)/2) z = -z."""

    def __init__(self, p):
        self.p = p
        self.t = nonresidue(p)

    def add(self, x, y):
        p = self.p
        return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)

    def sub(self, x, y):
        p = self.p
        return ((x[0] - y[0]) % p, (x[1] - y[1]) % p)

    def mul(self, x, y):
        p, t = self.p, self.t
        a, b = x
        c, d = y
        return ((a * c + t * b * d) % p, (a * d + b * c) % p)

    def one(self):
        return (1, 0)

    def zero(self):
        return (0, 0)

    def emb(self, a):
        return (a % self.p, 0)

    def pw(self, x, n):
        r = self.one()
        n %= (self.p * self.p - 1) if x != (0, 0) else 1
        while n:
            if n & 1:
                r = self.mul(r, x)
            x = self.mul(x, x)
            n >>= 1
        return r

    def inv(self, x):
        # (a+bz)^-1 = (a-bz)/(a^2 - t b^2)
        p, t = self.p, self.t
        a, b = x
        n = (a * a - t * b * b) % p
        ni = inv_mod(n, p)
        return ((a * ni) % p, (-b * ni) % p)

    def frob(self, x):
        return (x[0] % self.p, (-x[1]) % self.p)

    def order(self, x):
        n = self.p * self.p - 1
        o = n
        d = 2
        m = n
        fac = []
        while d * d <= m:
            if m % d == 0:
                fac.append(d)
                while m % d == 0:
                    m //= d
            d += 1
        if m > 1:
            fac.append(m)
        for f in fac:
            while o % f == 0 and self.pw(x, o // f) == self.one():
                o //= f
        return o


def root_of_unity(F, q):
    """The FIRST element of F_{p^2} of exact order q, in lexicographic order of
    (b, a) with x = a + b z. Deterministic; a different primitive q-th root only
    relabels the index set Z_q by i -> c i."""
    p = F.p
    n = p * p - 1
    assert n % q == 0
    for b in range(p):
        for a in range(p):
            x = (a, b)
            if x == (0, 0):
                continue
            y = F.pw(x, n // q)
            if y != F.one():
                # y^q = 1 and y != 1 with q prime, so y has exact order q
                return y
    raise RuntimeError('no q-th root of unity')


def root_in_Fp(p, q):
    """Least w in F_p^* of exact order q (used when q | p-1)."""
    for a in range(2, p):
        if pow(a, q, p) == 1:
            return a
    raise RuntimeError('no order-%d element mod %d' % (q, p))


# ---------------------------------------------------------------------------
# the normal form printed in paper.tex, Section 2
# ---------------------------------------------------------------------------
def branch(p, q):
    if (p - 1) % q == 0:
        return 'p-1'
    if (p + 1) % q == 0:
        return 'p+1'
    return None


def setup(p, q):
    """(F, w, winv, gammas) with gammas the raw admissible set (before the
    gamma ~ 1-gamma identification), each as an element of F_{p^2}."""
    F = F2(p)
    br = branch(p, q)
    if br == 'p-1':
        w0 = root_in_Fp(p, q)
        w = F.emb(w0)
        sub = set(pow(w0, i, p) for i in range(q))
        raw = []
        for g in range(1, p):                      # gamma != 0 (KNV)
            if (1 - inv_mod(g, p)) % p in sub:     # u_i = 0 for some i
                continue
            raw.append(F.emb(g))
    else:
        w = root_of_unity(F, q)
        a = inv_mod(2, p)                          # Tr(gamma) = 2a = 1
        sub = set(F.pw(w, i) for i in range(q))
        raw = []
        for b in range(p):
            g = (a, b)
            if F.sub(F.one(), F.inv(g)) in sub:
                continue
            raw.append(g)
    return F, w, F.inv(w), raw


def uvec(F, w, winv, g, q):
    """u_i = gamma w^i + (1-gamma) w^{-i}. Returns the tuple in F_p, or None if
    some u_i falls outside F_p."""
    p = F.p
    gg = F.sub(F.one(), g)
    out = []
    for i in range(q):
        x = F.add(F.mul(g, F.pw(w, i)), F.mul(gg, F.pw(winv, i)))
        if x[1] % p:
            return None
        out.append(x[0] % p)
    return tuple(out)


def admissible(u, p):
    """KNV's loop conditions on u, exactly as printed in paper.tex."""
    q = len(u)
    if u[0] % p != 1 % p:
        return False, 'u_0 != 1'
    for x in u:
        if x % p == 0:
            return False, 'some u_i = 0'
    for i in range(q):
        for j in range(q):
            if (u[i] + u[j]) % p == 0:
                return False, 'some u_i + u_j = 0'
    return True, ''


def members(p, q):
    """The Gamma-classes of the cell (p,q): a list of (gamma, u), one per
    isomorphism class in KNV's theta-family, gamma ~ 1-gamma identified.
    Deterministic order."""
    F, w, winv, raw = setup(p, q)
    seen = set()
    out = []
    for g in raw:
        u = uvec(F, w, winv, g, q)
        if u is None:
            continue
        rev = tuple(u[(-i) % q] for i in range(q))
        if u in seen or rev in seen:
            continue
        seen.add(u)
        seen.add(rev)
        out.append((g, u))
    return F, w, winv, out


# ---------------------------------------------------------------------------
# the multiplication table
# ---------------------------------------------------------------------------
def table_affine(p, q, u):
    """T[n][m] on Q = Z_q x F_p indexed n = i*p + j, built from the AFFINE form
    of Step 1: (i,j)(k,l) = (i+k, alpha_i j + m(1+alpha_i)),
    alpha_i = u_i/u_{i+k}, m = l u_k/(u_k + 1)."""
    n = p * q
    T = [[0] * n for _ in range(n)]
    for k in range(q):
        uk = u[k] % p
        mfac = (uk * inv_mod(uk + 1, p)) % p
        alpha = [(u[i] * inv_mod(u[(i + k) % q], p)) % p for i in range(q)]
        for l in range(p):
            m = (l * mfac) % p
            col = k * p + l
            for i in range(q):
                a = alpha[i]
                b = (m * (1 + a)) % p
                base = ((i + k) % q) * p
                row = i * p
                for j in range(p):
                    T[row + j][col] = base + ((a * j + b) % p)
    return T


def table_knv(p, q, u):
    """The SAME table from KNV's printed normal form, with theta_i = u_i^{-1}:
    (i,j)(k,l) = (i+k, l(1+theta_k)^{-1} + (j + l(1+theta_k)^{-1})
                        theta_i^{-1} theta_{i+k}).
    Independent transcription; Step 1 of the proof is the claim that the two
    agree."""
    n = p * q
    th = [inv_mod(x, p) for x in u]
    T = [[0] * n for _ in range(n)]
    for i in range(q):
        for j in range(p):
            for k in range(q):
                c = inv_mod((1 + th[k]) % p, p)
                for l in range(p):
                    s = (l * c) % p
                    v = (s + (j + s) * inv_mod(th[i], p) % p * th[(i + k) % q]) % p
                    T[i * p + j][k * p + l] = ((i + k) % q) * p + v
    return T


def is_loop(T, n):
    for x in range(n):
        if len(set(T[x])) != n:
            return False
        if len(set(T[y][x] for y in range(n))) != n:
            return False
    for e in range(n):
        if all(T[e][y] == y for y in range(n)) and all(T[y][e] == y for y in range(n)):
            return True
    return False


def is_assoc(T, n):
    for x in range(n):
        Tx = T[x]
        for y in range(n):
            xy = Tx[y]
            Ty = T[y]
            Txy = T[xy]
            for z in range(n):
                if Txy[z] != Tx[Ty[z]]:
                    return False
    return True


def is_abelian(T, n):
    for x in range(n):
        for y in range(x + 1, n):
            if T[x][y] != T[y][x]:
                return False
    return True


def is_right_bol(T, n):
    """((z x) y) x = z ((x y) x) for all x, y, z. Exhaustive."""
    for x in range(n):
        for y in range(n):
            xyx = T[T[x][y]][x]
            row = [T[T[T[z][x]][y]][x] for z in range(n)]
            col = [T[z][xyx] for z in range(n)]
            if row != col:
                return False
    return True


# ---------------------------------------------------------------------------
# |Mlt_r| by brute force -- no use of the proof
# ---------------------------------------------------------------------------
def mlt_r(T, n, want_elems=False, cap=4000000):
    """<R_x : x in Q> where R_x : w -> w*x, as a set of permutations encoded as
    bytes. Closes under a small seed, then FORCES closure by verifying that
    every one of the n right translations lies in the enumerated set, adding
    any that does not as a further generator. So an undercount is impossible."""
    assert n < 256
    ident = bytes(range(n))
    trans = [bytes(T[w][x] for w in range(n)) for x in range(n)]

    def comp(a, b):                     # a then b
        return bytes(b[v] for v in a)

    seed = [t for t in (trans[1], trans[min(n - 1, 2)], trans[n // 2]) if t != ident]
    elems = {ident}
    frontier = [ident]
    gens = list(dict.fromkeys(seed))
    while True:
        while frontier:
            nf = []
            for e in frontier:
                for g in gens:
                    h = comp(e, g)
                    if h not in elems:
                        elems.add(h)
                        nf.append(h)
                        if len(elems) > cap:
                            raise RuntimeError('group larger than cap')
            frontier = nf
        missing = [r for r in trans if r not in elems]
        if not missing:
            return (len(elems), elems) if want_elems else (len(elems), None)
        gens.append(missing[0])
        frontier = list(elems)


def kernel_of_shift(elems, p, q):
    """The subgroup of Mlt_r acting trivially on the Z_q coordinate: exactly the
    pure translations, by Step 2. Returned as a set of permutations."""
    return set(e for e in elems if all(e[i * p] // p == i for i in range(q)))


# ---------------------------------------------------------------------------
# dim R and the row space R itself
# ---------------------------------------------------------------------------
def _rref(rows, p, q):
    basis = []
    for r in rows:
        r = list(r)
        for b in basis:
            piv = next(i for i, v in enumerate(b) if v)
            if r[piv] % p:
                f = r[piv]
                r = [(r[c] - f * b[c]) % p for c in range(q)]
        if any(v % p for v in r):
            piv = next(i for i, v in enumerate(r) if v % p)
            iv = inv_mod(r[piv], p)
            r = [(v * iv) % p for v in r]
            basis.append(r)
            basis.sort(key=lambda b: next(i for i, v in enumerate(b) if v))
    for i in range(len(basis)):
        piv = next(c for c, v in enumerate(basis[i]) if v)
        for j in range(len(basis)):
            if i != j and basis[j][piv] % p:
                f = basis[j][piv]
                basis[j] = [(basis[j][c] - f * basis[i][c]) % p for c in range(q)]
    return tuple(tuple(b) for b in basis)


def gens_R(p, q, u):
    """x^(k)_i = u_i + u_{i+k}, and all cyclic shifts of each."""
    rows = []
    for k in range(q):
        x = [(u[i] + u[(i + k) % q]) % p for i in range(q)]
        for s in range(q):
            rows.append([x[(i + s) % q] for i in range(q)])
    return rows


def space_R(p, q, u):
    return _rref(gens_R(p, q, u), p, q)


def dim_R(p, q, u):
    return len(space_R(p, q, u))


def admissible_cells(maxn, qs=None):
    out = []
    for q in range(3, maxn):
        if not is_prime(q):
            continue
        if qs is not None and q not in qs:
            continue
        for p in range(q + 1, maxn // q + 1):
            if is_prime(p) and (p * p - 1) % q == 0:
                out.append((p, q))
    return sorted(out)


# ===========================================================================
def main():
    print('checks of the computational content of the claim: FOR ODD PRIMES p > q WITH')
    print('q DIVIDING p^2 - 1, every nonassociative right Bol loop of order pq has right')
    print('multiplication group of order p^2 q -- the two hypotheses are part of the')
    print('claim, exactly as in the title, abstract and theorem of the paper')
    print('(the paper\'s affirmative answer to KNV Problem 7.1). This program does NOT')
    print('prove that claim for all admissible p and q. It GRANTS the imported KNV')
    print('Theorem 1.1(ii)-(iii) classification -- that every such loop is Z_{pq} or a')
    print('member of the theta-family -- which is assumed here and NOT re-proved; it')
    print('brute-forces |Mlt_r| over that family only for pq <= 255; and it checks the')
    print('steps of the paper\'s hand proof as identities on finitely many cells.')
    print('See the SCOPE note at the end for what is not established.')
    print('python %s -- exact integer arithmetic in F_p and F_{p^2}, no floats,'
          % sys.version.split()[0])
    print('no randomness, standard library only')
    print('exactness guard: Fraction(1,3)*3 == 1 is %s' % (Fraction(1, 3) * 3 == 1))
    print('')

    # -------------------------------------------------------------------
    print('=== SECTION A -- the printed objects come out of the printed normal form ===')
    # A1. p=5, q=3: the paper prints F_25 = F_5[z]/(z^2-2), w = 2+3z,
    #     gamma = 3 + b z, u = (1, 2+2b, 2-2b), admissible b in {0,2,3}.
    F5 = F2(5)
    check('A1_least_nonresidue_mod_5_is_2', F5.t == 2, 't = %d' % F5.t)
    w53 = (2, 3)
    check('A1_w_equals_2_plus_3z_has_order_3_in_F25',
          F5.order(w53) == 3 and F5.add(F5.add(F5.mul(w53, w53), w53), F5.one()) == (0, 0),
          'w^2 + w + 1 = 0, order(w) = 3')
    F, w, winv, raw = setup(5, 3)
    check('A1_the_programs_own_root_of_unity_at_5_3_is_that_w', w == w53, 'w = %s' % (w,))
    printed = {}
    for b in range(5):
        g = (3, b)
        printed[b] = ((1, (2 + 2 * b) % 5, (2 - 2 * b) % 5))
    got = {}
    for b in range(5):
        got[b] = uvec(F, w, winv, (3, b), 3)
    check('A1_u_equals_1_2plus2b_2minus2b_for_every_b',
          all(got[b] == printed[b] for b in range(5)),
          'b=0..4 -> ' + ', '.join(str(got[b]) for b in range(5)))
    check('A1_gamma_3_plus_bz_has_trace_one_so_u_lies_in_F5',
          all(uvec(F, w, winv, (3, b), 3) is not None for b in range(5)),
          'Tr(3+bz) = 6 = 1 mod 5')
    adm = [b for b in range(5) if admissible(printed[b], 5)[0]]
    check('A1_admissible_b_is_exactly_0_2_3', adm == [0, 2, 3],
          'b in %s; b=1 gives u_2=0, b=4 gives u_1=0' % adm)
    F, w, winv, mem53 = members(5, 3)
    check('A1_the_two_Gamma_classes_at_5_3_are_1_2_2_and_1_1_3',
          [u for _, u in mem53] == [(1, 2, 2), (1, 1, 3)],
          'classes: %s' % ([u for _, u in mem53],))
    check('A1_class_b_0_is_symmetric_so_it_is_the_Bruck_loop_B_5_3',
          mem53[0][1] == tuple(mem53[0][1][(-i) % 3] for i in range(3)),
          'u = (1,2,2), u_i = u_{-i}, gamma = 1/2 = 3')

    # A2. p=7, q=3: KNV's Gamma = {1, 3, 4} with w = 2.
    check('A2_least_order_3_element_mod_7_is_2', root_in_Fp(7, 3) == 2, 'w = 2')
    F7, w7, wi7, raw7 = setup(7, 3)
    gs7 = sorted(g[0] for g in raw7)
    check('A2_admissible_gamma_mod_7_is_1_3_4_5', gs7 == [1, 3, 4, 5],
          'gamma in %s; gamma=2 and 6 give some u_i = 0' % gs7)
    _, _, _, mem73 = members(7, 3)
    check('A2_Gamma_classes_at_7_3_are_gamma_1_3_4_as_KNV_print',
          sorted(g[0] for g, _ in mem73) == [1, 3, 4],
          'representatives %s, since 5 = 1-3' % sorted(g[0] for g, _ in mem73))
    check('A2_their_u_vectors_are_1_2_4_and_1_5_1_and_1_3_3',
          sorted(u for _, u in mem73) == sorted([(1, 2, 4), (1, 5, 1), (1, 3, 3)]),
          '%s' % ([u for _, u in mem73],))

    # A3. the loop conditions, over the whole census range
    cells = admissible_cells(255)
    check('A3_admissible_cells_with_pq_at_most_255_number_29', len(cells) == 29,
          '%d cells, %s' % (len(cells), cells))
    check('A3_every_cell_has_p_congruent_to_plus_or_minus_1_mod_q',
          all(branch(p, q) is not None for p, q in cells),
          'q | p^2-1 <=> p = +-1 mod q')
    allmem = {}
    bad = []
    for (p, q) in cells:
        _, w, wi, mm = members(p, q)
        allmem[(p, q)] = mm
        for g, u in mm:
            ok, why = admissible(u, p)
            if not ok:
                bad.append((p, q, u, why))
    check('A3_every_member_satisfies_u_0_1_and_u_i_nonzero_and_u_i_plus_u_j_nonzero',
          not bad, '%d members over 29 cells, 0 violations'
          % sum(len(v) for v in allmem.values()))
    # A4. Niederreiter-Robinson: u obeys u_{i+2} = lambda u_{i+1} - u_i, lambda = w + w^{-1}
    badrec = []
    for (p, q) in cells:
        F, w, wi, mm = members(p, q)
        lam = F.add(w, wi)
        assert lam[1] % p == 0
        lam = lam[0] % p
        for g, u in mm:
            for i in range(q):
                if (u[(i + 2) % q] - lam * u[(i + 1) % q] + u[i]) % p:
                    badrec.append((p, q, u))
    check('A4_every_member_obeys_the_linear_recurrence_u_i2_eq_lambda_u_i1_minus_u_i',
          not badrec, 'lambda = w + w^{-1} in F_p; 0 violations over all 29 cells')

    # -------------------------------------------------------------------
    print('')
    print('=== SECTION B -- the tables are right Bol loops, by the identity ===')
    small = [(p, q) for (p, q) in cells if p * q <= 105]
    check('B0_cells_checked_by_the_full_right_Bol_identity',
          small == [(5, 3), (7, 3), (11, 3), (11, 5), (13, 3), (13, 7), (17, 3),
                    (19, 3), (19, 5), (23, 3), (29, 3), (31, 3)],
          '%s (pq <= 105)' % (small,))
    nb = na = 0
    bfail = []
    for (p, q) in small:
        n = p * q
        for g, u in allmem[(p, q)]:
            T = table_affine(p, q, u)
            if not is_loop(T, n):
                bfail.append((p, q, u, 'not a loop'))
            if not is_right_bol(T, n):
                bfail.append((p, q, u, 'not right Bol'))
            nb += 1
            if is_assoc(T, n):
                na += 1
    check('B1_every_member_of_those_cells_is_a_loop_and_satisfies_the_right_Bol_identity',
          not bfail, '%d members, all right Bol loops' % nb)
    check('B2_the_associative_members_of_those_cells_number_exactly_the_cells_with_q_dividing_p_minus_1',
          na == len([1 for (p, q) in small if (p - 1) % q == 0]),
          '%d associative among %d; cells with q | p-1: %s'
          % (na, nb, [(p, q) for (p, q) in small if (p - 1) % q == 0]))
    # the two transcriptions of the normal form agree -- this IS Step 1
    dis = []
    for (p, q) in small:
        for g, u in allmem[(p, q)]:
            if table_affine(p, q, u) != table_knv(p, q, u):
                dis.append((p, q, u))
    check('B3_the_affine_form_of_Step_1_equals_KNVs_printed_normal_form_entry_by_entry',
          not dis, 'two independent transcriptions agree on all %d tables of those cells' % nb)

    # -------------------------------------------------------------------
    print('')
    print('=== SECTION C -- the census over the ASSUMED theta-family (plus Z_{pq}), '
          'pq <= 255; the classification that this family is exhaustive is imported, '
          'not checked here ===')
    tot = 0
    assoc = 0
    nonassoc = 0
    okp2q = 0
    okdim2 = 0
    anomalies = []
    rinn_ok = 0
    for (p, q) in cells:
        mm = allmem[(p, q)]
        tot += len(mm) + 1              # + the cyclic group Z_{pq}, KNV Thm 1.1(ii)
        assoc += 1
        for g, u in mm:
            n = p * q
            T = table_affine(p, q, u)
            d = dim_R(p, q, u)
            mo, _ = mlt_r(T, n)
            if d == 1:
                assoc += 1
                if mo != p * q:
                    anomalies.append((p, q, u, d, mo))
            else:
                nonassoc += 1
                if mo == p * p * q:
                    okp2q += 1
                    if mo // (p * q) == p:
                        rinn_ok += 1
                else:
                    anomalies.append((p, q, u, d, mo))
                if d == 2:
                    okdim2 += 1
                else:
                    anomalies.append((p, q, u, d, mo))
    check('C1_loops_in_the_census_number_534_equal_to_sum_of_p_minus_q_plus_4_over_2',
          tot == 534 and tot == sum((p - q + 4) // 2 for p, q in cells),
          '%d objects generated up to isomorphism = the theta-family members plus Z_{pq} '
          'over the 29 cells with pq <= 255; that these EXHAUST the right Bol loops of '
          'those orders is granted from KNV Thm 1.1(ii)-(iii), not established here' % tot)
    check('C2_members_with_dim_R_equal_1_number_44_being_29_copies_of_Z_pq_plus_15_one_per_cell_with_q_dividing_p_minus_1',
          assoc == 44 and assoc == 29 + len([1 for p, q in cells if (p - 1) % q == 0]),
          '%d = 29 cyclic Z_{pq} (added once per cell by construction, not re-tested here) '
          '+ %d theta-family members with dim R = 1, one for each cell with q | p-1. That '
          'such a member is exactly the nonabelian group of order pq is the classification '
          'statement quoted in Section 2 of the paper (Q_gamma is associative iff gamma = 1, '
          'which needs q | p-1), GRANTED here; is_assoc and is_abelian are actually run on '
          'one of these 15 only at (7,3) in E9, and as a count over the pq <= 105 cells in B2'
          % (assoc, len([1 for p, q in cells if (p - 1) % q == 0])))
    check('C3_members_with_dim_R_not_equal_1_the_criterion_used_here_for_nonassociativity_number_490',
          nonassoc == 490,
          '%d members have dim R != 1, which is the criterion this census uses for '
          'nonassociativity; the test not is_assoc(T, pq) is NOT run member by member. The '
          'equivalence is granted from the paper\'s Section 2, and each of these members is '
          'separately found in C4 to have |Mlt_r| = p^2 q rather than the pq of a group'
          % nonassoc)
    check('C4_all_490_have_Mlt_r_of_order_p_squared_q_by_brute_force', okp2q == 490,
          '%d/%d, |Mlt_r| computed by forced closure of ALL pq right translations'
          % (okp2q, nonassoc))
    check('C5_all_490_have_dim_R_equal_2_by_the_formula', okdim2 == 490,
          '%d/%d' % (okdim2, nonassoc))
    check('C6_anomalies_are_zero', not anomalies, '0 anomalies')
    check('C7_RInn_has_order_p_throughout_since_Mlt_r_equals_pq_times_RInn',
          rinn_ok == 490, '|Mlt_r|/pq = p for all %d' % rinn_ok)
    check('C8_the_two_routines_are_independent_and_agree_on_every_census_loop',
          okp2q == okdim2 == 490,
          'mlt_r never reads the proof; dim_R never builds the group')

    # C9. the kernel P: elementary abelian of order p^2 (Step 2), on the small cells
    kok = 0
    ktot = 0
    for (p, q) in small:
        for g, u in allmem[(p, q)]:
            if dim_R(p, q, u) != 2:
                continue
            n = p * q
            T = table_affine(p, q, u)
            mo, elems = mlt_r(T, n, want_elems=True)
            P = kernel_of_shift(elems, p, q)
            ktot += 1
            # order p^2, abelian, exponent p, and every element a pure translation
            ident = bytes(range(n))

            def comp(a, b):
                return bytes(b[v] for v in a)
            expp = all(_pow_perm(x, p, n) == ident for x in P)
            ab = all(comp(x, y) == comp(y, x) for x in P for y in P)
            trans_only = all(all(x[i * p] // p == i for i in range(q)) for x in P)
            if len(P) == p * p and expp and ab and trans_only:
                kok += 1
    check('C9_the_kernel_of_the_shift_action_is_elementary_abelian_of_order_p_squared',
          kok == ktot and ktot > 0,
          '%d/%d nonassociative members of the pq<=105 cells: |P| = p^2, abelian, exponent p, '
          'pure translations only' % (kok, ktot))

    # -------------------------------------------------------------------
    print('')
    print('=== SECTION D -- each step of the proof, checked as an identity on finitely '
          'many cells only, NOT symbolically for all p and q: D2-D4 stop at pq <= 105, '
          'the rest at pq <= 255, and D4 uses one sampled family of translations ===')
    # D1 (Step 1). R_{(k,l)} is the affine map (i,j) -> (i+k, alpha_i j + beta_i).
    d1 = []
    for (p, q) in cells:
        for g, u in allmem[(p, q)]:
            n = p * q
            T = table_affine(p, q, u)
            for k in range(q):
                uk = u[k]
                mfac = (uk * inv_mod(uk + 1, p)) % p
                for l in range(p):
                    m = (l * mfac) % p
                    for i in range(q):
                        a = (u[i] * inv_mod(u[(i + k) % q], p)) % p
                        b = (m * (1 + a)) % p
                        for j in range(p):
                            if T[i * p + j][k * p + l] != ((i + k) % q) * p + (a * j + b) % p:
                                d1.append((p, q, u, i, j, k, l))
    check('D1_Step1_every_right_translation_is_the_stated_affine_map_with_alpha_i_u_i_over_u_i_plus_k',
          not d1, 'all 29 cells, all members, all (i,j,k,l): 0 mismatches')
    # D2 (Step 3). R_{(1,0)}^k = R_{(k,0)}.
    d2 = []
    for (p, q) in small:
        for g, u in allmem[(p, q)]:
            n = p * q
            T = table_affine(p, q, u)
            s = bytes(T[x][1 * p + 0] for x in range(n))
            for k in range(q):
                if _pow_perm(s, k, n) != bytes(T[x][k * p + 0] for x in range(n)):
                    d2.append((p, q, u, k))
    check('D2_Step3_R_1_0_to_the_k_equals_R_k_0', not d2,
          's = R_{(1,0)} has beta = 0; s^k = R_{(k,0)} for all k, over the pq<=105 cells')
    # D3 (Step 3). R_{(k,l)} = n_{k,l} s^k with n translating by t_i = m(1 + u_{i+k}/u_i).
    d3 = []
    for (p, q) in small:
        for g, u in allmem[(p, q)]:
            n = p * q
            T = table_affine(p, q, u)
            s = bytes(T[x][1 * p + 0] for x in range(n))
            for k in range(q):
                uk = u[k]
                mfac = (uk * inv_mod(uk + 1, p)) % p
                sk = _pow_perm(s, k, n)
                for l in range(p):
                    m = (l * mfac) % p
                    t = [(m * (1 + u[(i + k) % q] * inv_mod(u[i], p))) % p for i in range(q)]
                    nn = bytes(i * p + (j + t[i]) % p for i in range(q) for j in range(p))
                    # n_{k,l} then s^k
                    if bytes(sk[v] for v in nn) != bytes(T[x][k * p + l] for x in range(n)):
                        d3.append((p, q, u, k, l))
    check('D3_Step3_R_k_l_factors_as_the_translation_by_t_i_m_1_plus_u_i_plus_k_over_u_i_times_s_to_the_k',
          not d3, '0 mismatches over the pq<=105 cells')
    # D4 (Step 3). conjugation by s acts on t as t -> (t_{i+1} u_{i+1}/u_i)_i, i.e. as the
    #             plain cyclic shift in the coordinates r_i = t_i u_i.
    d4 = []
    for (p, q) in small:
        for g, u in allmem[(p, q)]:
            n = p * q
            T = table_affine(p, q, u)
            s = bytes(T[x][1 * p + 0] for x in range(n))
            sinv = _inv_perm(s, n)
            for t0 in range(1, p):
                t = [(t0 * (i + 1)) % p for i in range(q)]     # an arbitrary translation
                nn = bytes(i * p + (j + t[i]) % p for i in range(q) for j in range(p))
                # s^{-1} n s, as maps: apply s, then n, then s^{-1}
                conj = bytes(sinv[nn[s[x]]] for x in range(n))
                tp = [(conj[i * p] - i * p) % p for i in range(q)]
                pred = [(t[(i + 1) % q] * u[(i + 1) % q] % p * inv_mod(u[i], p)) % p
                        for i in range(q)]
                r = [(t[i] * u[i]) % p for i in range(q)]
                rp = [(tp[i] * u[i]) % p for i in range(q)]
                # the opposite conjugate, which must be the shift the other way
                conj2 = bytes(s[nn[sinv[x]]] for x in range(n))
                tp2 = [(conj2[i * p] - i * p) % p for i in range(q)]
                rp2 = [(tp2[i] * u[i]) % p for i in range(q)]
                if (tp != pred
                        or rp != [r[(i + 1) % q] for i in range(q)]
                        or rp2 != [r[(i - 1) % q] for i in range(q)]):
                    d4.append((p, q, u, t0))
    check('D4_Step3_conjugation_by_s_is_the_plain_cyclic_shift_in_the_coordinates_r_i_t_i_u_i',
          not d4, 's^{-1} n s shifts r forward and s n s^{-1} shifts it back; 0 mismatches over '
                  'the pq<=105 cells and every t of the sampled family')
    # D5 (Step 3). the generators of P become m * x^(k) with x^(k)_i = u_i + u_{i+k}.
    d5 = []
    for (p, q) in cells:
        for g, u in allmem[(p, q)]:
            for k in range(q):
                uk = u[k]
                mfac = (uk * inv_mod(uk + 1, p)) % p
                for l in range(p):
                    m = (l * mfac) % p
                    t = [(m * (1 + u[(i + k) % q] * inv_mod(u[i], p))) % p for i in range(q)]
                    r = [(t[i] * u[i]) % p for i in range(q)]
                    if r != [(m * ((u[i] + u[(i + k) % q]) % p)) % p for i in range(q)]:
                        d5.append((p, q, u, k, l))
    check('D5_Step3_in_r_coordinates_the_generators_are_m_times_x_k_with_x_k_i_u_i_plus_u_i_plus_k',
          not d5, 'all 29 cells, all members, all (k,l): 0 mismatches')
    # D6 (Step 4). x^(k) = gamma(1+w^k) v_1 + (1-gamma)(1+w^{-k}) v_{-1}.
    d6 = []
    for (p, q) in cells:
        F, w, wi, mm = members(p, q)
        for g, u in mm:
            gg = F.sub(F.one(), g)
            for k in range(q):
                c1 = F.mul(g, F.add(F.one(), F.pw(w, k)))
                c2 = F.mul(gg, F.add(F.one(), F.pw(wi, k)))
                for i in range(q):
                    lhs = F.emb((u[i] + u[(i + k) % q]) % p)
                    rhs = F.add(F.mul(c1, F.pw(w, i)), F.mul(c2, F.pw(wi, i)))
                    if lhs != rhs:
                        d6.append((p, q, u, k, i))
    check('D6_Step4_x_k_equals_gamma_1_plus_w_k_v_1_plus_1_minus_gamma_1_plus_w_minus_k_v_minus_1',
          not d6, 'identity in F_{p^2}, all 29 cells, all members, all k, all coordinates')
    # D7 (Step 4). v_1, v_{-1} are F_{p^2}-independent, so W has dimension 2 and dim R <= 2.
    d7 = []
    for (p, q) in cells:
        F, w, wi, mm = members(p, q)
        # a 2x2 minor of the (q x 2) matrix [v_1 | v_{-1}] is invertible
        det = F.sub(F.mul(F.one(), F.pw(wi, 1)), F.mul(F.pw(w, 1), F.one()))
        if det == F.zero():
            d7.append((p, q))
    check('D7_Step4_v_1_and_v_minus_1_are_independent_over_F_p2_so_dim_W_is_2',
          not d7, 'minor w^{-1} - w != 0 in all 29 cells, since w has odd order q > 1')
    check('D8_Step4_dim_R_never_exceeds_2_over_the_whole_census',
          all(dim_R(p, q, u) <= 2 for (p, q) in cells for _, u in allmem[(p, q)]),
          'F_p-independence inside F_p^{Z_q} survives extension to F_{p^2}')
    # D9 (Step 5). det[[2g, 2(1-g)], [2 g w, 2(1-g) w^{-1}]] = 4 g (1-g) (w^{-1} - w), != 0.
    d9 = []
    d9nz = 0
    for (p, q) in cells:
        F, w, wi, mm = members(p, q)
        for g, u in mm:
            gg = F.sub(F.one(), g)
            two = F.emb(2)
            M = [[F.mul(two, g), F.mul(two, gg)],
                 [F.mul(F.mul(two, g), w), F.mul(F.mul(two, gg), wi)]]
            det = F.sub(F.mul(M[0][0], M[1][1]), F.mul(M[0][1], M[1][0]))
            closed = F.mul(F.mul(F.emb(4), F.mul(g, gg)), F.sub(wi, w))
            if det != closed:
                d9.append((p, q, u))
            if det != F.zero():
                d9nz += 1
            elif g != F.one():
                d9.append(('unexpected zero', p, q, u))
    check('D9_Step5_the_2x2_determinant_of_x_0_and_its_shift_is_4_gamma_1_minus_gamma_w_inv_minus_w',
          not d9, 'closed form verified over all 29 cells and every member')
    nonass_mem = [(p, q, g, u) for (p, q) in cells for g, u in allmem[(p, q)]
                  if dim_R(p, q, u) == 2]
    check('D10_Step5_that_determinant_is_nonzero_for_every_nonassociative_member',
          d9nz == len(nonass_mem) == 490,
          '%d of %d members have nonzero determinant; the %d with gamma = 1 are the '
          'nonabelian groups' % (d9nz, sum(len(v) for v in allmem.values()), 534 - 29 - 490))
    # D11. the strengthening: R = W cap F_p^{Z_q} does not depend on gamma.
    d11 = []
    for (p, q) in cells:
        sp = set(space_R(p, q, u) for _, u in allmem[(p, q)] if dim_R(p, q, u) == 2)
        if len(sp) != 1:
            d11.append((p, q, len(sp)))
    check('D11_the_plane_R_is_the_same_for_every_nonassociative_gamma_in_each_cell',
          not d11, 'exactly one 2-dimensional R per cell, over all 29 cells')
    # D12. the eigenvalues of the cyclic shift on R are exactly w and w^{-1}.
    d12 = []
    for (p, q) in cells:
        F, w, wi, mm = members(p, q)
        lam = F.add(w, wi)[0] % p                      # w + w^{-1} lies in F_p
        for _, u in mm:
            B = space_R(p, q, u)
            if len(B) != 2:
                continue
            # write the shift of each basis vector in the basis B
            M = []
            for b in B:
                sb = [b[(i + 1) % q] for i in range(q)]
                co = _coords(sb, B, p, q)
                if co is None:
                    d12.append((p, q, u, 'not shift stable'))
                    break
                M.append(co)
            else:
                tr = (M[0][0] + M[1][1]) % p
                de = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) % p
                if tr != lam % p or de != 1 % p:
                    d12.append((p, q, u, tr, de, lam))
    check('D12_the_shift_acting_on_R_has_characteristic_polynomial_X2_minus_w_plus_w_inv_X_plus_1',
          not d12, 'trace = w + w^{-1}, determinant = 1, i.e. eigenvalues exactly w and w^{-1}, '
                   'over all 29 cells')

    # -------------------------------------------------------------------
    print('')
    print('=== SECTION E -- the controls ===')
    # E1. FORCED POSITIVE: u = (1,1,2) at p=5,q=3 is a loop, is NOT right Bol,
    #     and the decider returns the UPPER branch p^3 q = 375.
    uc = (1, 1, 2)
    n = 15
    Tc = table_affine(5, 3, uc)
    check('E1_control_u_1_1_2_at_p_5_q_3_meets_KNVs_loop_conditions',
          admissible(uc, 5)[0], 'u_0 = 1, all u_i != 0, all u_i + u_j != 0')
    check('E2_control_is_a_loop_of_order_15', is_loop(Tc, n), 'both translations bijective, '
          'two-sided identity present')
    check('E3_control_is_NOT_right_Bol', not is_right_bol(Tc, n),
          '((zx)y)x = z((xy)x) fails')
    check('E4_control_fails_the_recurrence_for_every_lambda_in_F_5',
          all(any((uc[(i + 2) % 3] - lam * uc[(i + 1) % 3] + uc[i]) % 5 for i in range(3))
              for lam in range(5)),
          'so Niederreiter-Robinson does not apply, consistent with E3')
    moc, _ = mlt_r(Tc, n)
    check('E5_control_has_Mlt_r_of_order_375_the_upper_branch_p_cubed_q', moc == 375,
          '|Mlt_r| = %d = 5^3 * 3, while p^2 q = 75' % moc)
    check('E6_control_has_dim_R_equal_3_so_the_formula_returns_the_upper_branch_too',
          dim_R(5, 3, uc) == 3, 'dim R = 3, and |Mlt_r| = p^{dim R} q')
    # E7. the hand determinant printed in the paper
    x0 = [(2 * t) % 5 for t in uc]
    check('E7_x_0_for_the_control_is_2_2_4', x0 == [2, 2, 4], 'x^(0) = 2u = (2,2,4)')
    M = [[x0[(i + s) % 3] for i in range(3)] for s in range(3)]
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
           - M[1][0] * (M[0][1] * M[2][2] - M[0][2] * M[2][1])
           + M[2][0] * (M[0][1] * M[1][2] - M[0][2] * M[1][1]))
    check('E8_the_circulant_of_2_2_4_has_determinant_minus_32_equal_3_mod_5',
          det == -32 and det % 5 == 3,
          'first-column expansion 2(8-4) - 2(4-8) + 4(4-16) = 8 + 8 - 48 = -32 = 3 mod 5')
    # E9. gamma = 1 at (7,3): the nonabelian group of order 21.
    u1 = [u for g, u in mem73 if g == (1, 0)][0]
    T1 = table_affine(7, 3, u1)
    mo1, _ = mlt_r(T1, 21)
    check('E9_gamma_equals_1_at_7_3_gives_the_nonabelian_group_of_order_21',
          is_assoc(T1, 21) and not is_abelian(T1, 21) and mo1 == 21 and dim_R(7, 3, u1) == 1,
          'u = %s, associative, nonabelian, |Mlt_r| = %d = pq, dim R = 1' % (u1, mo1))
    # E10. u == 1 (theta == 1): the cyclic group.
    for (p, q) in [(5, 3), (7, 3)]:
        uu = tuple([1] * q)
        Tz = table_affine(p, q, uu)
        moz, _ = mlt_r(Tz, p * q)
        check('E10_u_identically_1_at_%d_%d_gives_the_cyclic_group_of_order_%d' % (p, q, p * q),
              is_assoc(Tz, p * q) and is_abelian(Tz, p * q) and moz == p * q
              and dim_R(p, q, uu) == 1,
              'associative, abelian, |Mlt_r| = %d = pq, dim R = 1' % moz)
    # E11. B_{p,q} at gamma = 1/2 -- the case KNV proved (Theorem 1.1(v)).
    bok = []
    for (p, q) in [(5, 3), (7, 3), (11, 3), (11, 5), (13, 7)]:
        half = None
        for g, u in allmem[(p, q)]:
            if u == tuple(u[(-i) % q] for i in range(q)) and dim_R(p, q, u) == 2:
                half = u
        n = p * q
        T = table_affine(p, q, half)
        mo, elems = mlt_r(T, n, want_elems=True)
        P = kernel_of_shift(elems, p, q)
        ident = bytes(range(n))
        okk = (mo == p * p * q and len(P) == p * p
               and all(_pow_perm(x, p, n) == ident for x in P))
        bok.append((p, q, half, mo, okk))
        check('E11_B_%d_%d_matches_KNV_Theorem_1_1_v' % (p, q), okk,
              'u = %s symmetric (gamma = 1/2), |Mlt_r| = %d = p^2 q, kernel elementary '
              'abelian of order p^2' % (half, mo))

    # -------------------------------------------------------------------
    print('')
    print('=== SECTION F -- two cardinality counts the paper no longer states (F1, F2) '
          'and the wide dimension sweep ===')
    qs = (3, 5, 7, 11, 13)
    wide = [(p, q) for q in qs for p in range(q + 1, 200)
            if is_prime(p) and (p * p - 1) % q == 0]
    mm1 = mm2 = 0
    for (p, q) in wide:
        _, _, _, mem = members(p, q)
        if len(mem) != (p - q + 2) // 2:
            mm1 += 1
        want = (p - q) // 2 if (p - 1) % q == 0 else (p - q + 2) // 2
        if len([1 for _, u in mem if dim_R(p, q, u) == 2]) != want:
            mm2 += 1
    check('F1_Gamma_has_p_minus_q_plus_2_over_2_classes_a_count_the_paper_does_not_state',
          mm1 == 0,
          '%d cells, q in %s, p < 200: 0 mismatches. NOT A PAPER CLAIM: the paper states '
          'no cardinality for Gamma, so this check corresponds to nothing in it; it is kept '
          'as extra verification of the parameterisation and is not evidence for the theorem'
          % (len(wide), list(qs)))
    check('F2_nonassociative_members_number_p_minus_q_over_2_when_q_divides_p_minus_1_and_p_minus_q_plus_2_over_2_otherwise_a_count_the_paper_does_not_state',
          mm2 == 0,
          '0 mismatches over the same %d cells. NOT A PAPER CLAIM: like F1, the number of '
          'nonassociative members per cell is asserted nowhere in the paper; the check is '
          'retained because it exercises the census, not because the paper rests on it'
          % len(wide))
    sweep = [(p, q) for q in qs for p in range(q + 1, 401)
             if is_prime(p) and (p * p - 1) % q == 0]
    nsw = ngroup = ndim2 = 0
    dist = {}
    for (p, q) in sweep:
        _, _, _, mem = members(p, q)
        for g, u in mem:
            d = dim_R(p, q, u)
            nsw += 1
            key = ('group' if d == 1 else 'nonassoc', d)
            dist[key] = dist.get(key, 0) + 1
            if d == 1:
                ngroup += 1
            elif d == 2:
                ndim2 += 1
    check('F3_the_wide_sweep_covers_every_member_of_every_cell_with_q_at_most_13_and_p_at_most_400',
          nsw == ngroup + ndim2,
          '%d cells, %d members examined' % (len(sweep), nsw))
    check('F4_every_member_of_the_wide_sweep_with_dim_R_not_equal_1_has_dim_R_equal_2',
          ndim2 == nsw - ngroup and ndim2 > 0,
          '%d of %d members with dim R != 1 have dim R = 2; the other %d have dim R = 1 and '
          'are the gamma = 1 nonabelian groups by the classification statement quoted in '
          'Section 2 of the paper. This sweep builds no multiplication table, so both labels '
          'come from dim_R alone: neither associativity nor gamma = 1 is tested here'
          % (ndim2, nsw - ngroup, ngroup))
    check('F5_the_dimension_distribution_of_the_wide_sweep_has_no_third_value',
          set(dist) == {('group', 1), ('nonassoc', 2)},
          '%s' % sorted(dist.items()))
    check('F6_the_wide_sweep_reaches_p_400_far_beyond_the_brute_force_range_p_83',
          max(p for p, q in sweep) >= 397 and max(p for p, q in cells) == 83,
          'brute force reached p = 83 (pq <= 255); the formula sweep reaches p = %d'
          % max(p for p, q in sweep))

    # -------------------------------------------------------------------
    print('')
    note('SCOPE. What this program does NOT establish. (1) KNV Theorem 1.1(ii)-(iii), '
         'that every right Bol loop of order pq is isomorphic to Z_{pq} or to a member of the '
         'theta-family parameterised above, is NOT re-proved here: it is assumed, and every '
         'count in Sections C and F is a count over that family only. If the classification '
         'were incomplete the theorem would be too. (2) NOT RE-RUN: any loop outside the '
         'theta-family; any order other than pq -- and at order p^3 = 27 the analogous upper '
         'branch IS attained, so nothing here is a statement about right Bol loops in general. '
         '(3) The brute-force census is bounded at pq <= 255 (29 cells, p <= 83) and the '
         'dimension sweep at q <= 13, p <= 400; larger cells are covered only by the hand proof '
         'of Section 3 of the paper, which this program checks step by step but does not '
         'replace. (4) The '
         'control in Section E lies OUTSIDE the right Bol family, so it validates the decider, '
         'not the existence of a right Bol loop with |Mlt_r| = p^3 q -- the paper claims there '
         'is none. (5) Isotopism counts, KNV Conjecture 7.3 and the GAP LOOPS package are not '
         'touched by this program. (6) WORD USE, and it is load-bearing: the only checks that '
         'actually execute an associativity test are B2 (as a count over the pq <= 105 cells) '
         'and the individual samples E9 and E10. Everywhere else -- the census counts C2-C5, '
         'C7 and C9, the "nonassociative member" of D10 and D11, and the sweep counts F2, F4 '
         'and F5, including the "nonassoc" label printed in the F5 distribution -- membership '
         'in the nonassociative class is decided by dim R != 1 (equivalently, within this '
         'range, dim R = 2), and "the gamma = 1 nonabelian group" names the dim R = 1 members. '
         'That dim R = 1 holds exactly for the associative members, and exactly at gamma = 1, '
         'is imported from the classification statement quoted in Section 2 of the paper; it '
         'is NOT re-derived member by member here. (7) OVER-CHECK, corresponding to nothing '
         'in the paper: F1 and F2 verify the cardinality formulas |Gamma| = (p-q+2)/2 and '
         '(p-q)/2 or (p-q+2)/2 nonassociative members per cell, over 91 cells. The paper '
         'states neither formula anywhere, so those two checks establish more than the paper '
         'claims and support no statement in it; they are reported here for completeness and '
         'nothing in the theorem or its proof depends on them.')
    print('')
    if _FAILS:
        print('VERDICT: %d CHECKS FAILED (%s)' % (len(_FAILS), ', '.join(_FAILS[:5])))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % len(_PASSES))
    return 0


# ---------------------------------------------------------------------------
# permutation helpers (defined after use in main() only for readability)
# ---------------------------------------------------------------------------
def _pow_perm(s, k, n):
    r = bytes(range(n))
    b = s
    while k:
        if k & 1:
            r = bytes(b[v] for v in r)
        b = bytes(b[v] for v in b)
        k >>= 1
    return r


def _inv_perm(s, n):
    out = [0] * n
    for i, v in enumerate(s):
        out[v] = i
    return bytes(out)


def _coords(vec, basis, p, q):
    """Coordinates of vec in the reduced-row-echelon basis, or None if vec is
    not in the span."""
    co = []
    r = list(v % p for v in vec)
    for b in basis:
        piv = next(i for i, v in enumerate(b) if v)
        c = r[piv] % p
        co.append(c)
        if c:
            r = [(r[i] - c * b[i]) % p for i in range(q)]
    if any(v % p for v in r):
        return None
    return co


if __name__ == '__main__':
    sys.exit(main())
