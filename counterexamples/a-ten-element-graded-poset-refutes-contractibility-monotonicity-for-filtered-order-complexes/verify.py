#!/usr/bin/env python3
"""Referee verifier for "A ten-element graded poset refutes contractibility
monotonicity for filtered order complexes".

TARGET (Kitajima, arXiv:2606.15241v2, Remark 2.13, p. 8, second sentence):
    "It is an open problem whether Delta^{(k)}(P) is contractible whenever
     Delta^{(k-1)}(P) is contractible."
Here P is a finite graded (= pure) poset of rank n and Delta^{(k)}(P) is the
subcomplex of the order complex Delta(P) spanned by the chains of rank spread
at most k, 0 <= k <= n.

WHAT THIS PROGRAM DOES.  It reads the fifteen cover relations EXACTLY AS THE
PAPER PRINTS THEM (the block COVER_TABLE below is a transcription of the
displayed table in Section 2, and is parsed, not hand-expanded), and from
those fifteen pairs alone it re-derives every quantity the paper asserts:

  * the transitive closure and its rank-gap profile, and that the cover
    relation recomputed from the closure is again exactly those fifteen pairs;
  * that P has no minimum and no maximum and is pure of rank 3;
  * the f-vectors, Euler characteristics and INTEGRAL homology (Smith normal
    form, so torsion is visible) of Delta^{(1)}, Delta^{(2)}, Delta^{(3)};
  * that Delta^{(2)}(P) is collapsible, by REPLAYING THE PAPER'S OWN PRINTED
    27-row collapse table (the block COLLAPSE_TABLE below is a transcription
    of the table displayed in Section 4, and is parsed): at its own step each
    printed tau is checked to have exactly one proper coface in the complex
    then present, and that coface to be the printed sigma.  An independent
    free-face search is then required to return that same table, so the paper
    and the program cannot disagree silently;
  * every intermediate object of the two Mayer-Vietoris proofs in Sections 3
    and 4 -- the star decomposition of Delta(P), the complex Y, the pieces
    A, B, A cap B, the two 4-cycles, and the claim that the comparison map
    phi is an isomorphism (checked as an equality of integer lattices, not
    numerically);
  * the two printed chain identities of Section 4, as membership in the
    lattice of 2-boundaries of A;
  * the paper's scope claim that this P cannot be coned into a witness at a
    proper filtration stage;
  * five controls, in both polarities, including one against a published
    formula of the source paper and one against a space with torsion.

Exact integer arithmetic throughout.  No floating point, no randomness, no
external data file.  Python 3.9+, standard library only (itertools).  Runs in
well under a second on a laptop.

Contract: one `PASS <name>` line per check, closing with
`VERDICT: ALL <n> CHECKS PASS`, and exit status 0 if and only if every check
passed.  What is NOT covered is printed as `NOT RE-RUN:` lines at the end.
"""

import itertools
import sys

# ---------------------------------------------------------------------------
# 0. THE OBJECT, AS PRINTED IN THE PAPER
# ---------------------------------------------------------------------------
# Transcribed from the displayed table of Section 2.  Parsed below; nothing in
# this program hard-codes the expanded list of pairs.
COVER_TABLE = """
    a0 < b0
    a1 < b1      a1 < b2
    b0 < c0      b0 < c1      b0 < c2
    b1 < c0      b1 < c1
    b2 < c0      b2 < c2
    c0 < d0
    c1 < d0      c1 < d1
    c2 < d0      c2 < d1
"""
RANK_TABLE = """
    rank 0: a0 a1
    rank 1: b0 b1 b2
    rank 2: c0 c1 c2
    rank 3: d0 d1
"""
# Transcribed from the 27-row table PRINTED IN THE PAPER, Section 4, in the
# paper's own numbering (the printed table is read column by column: rows 1-9,
# then 10-18, then 19-27).  This is the paper's certificate, not this program's
# output: it is parsed and REPLAYED below against Delta^{(2)}(P), and the
# replay is what licenses the paper's Proposition.  The program's own
# independent greedy search is then required to reproduce this list exactly,
# so a discrepancy between paper and program cannot pass silently.
COLLAPSE_TABLE = """
     1  a0 < c0            a0 < b0 < c0
     2  a0 < c1            a0 < b0 < c1
     3  a0 < b0            a0 < b0 < c2
     4  a0                 a0 < c2
     5  a1 < c1            a1 < b1 < c1
     6  a1 < b1            a1 < b1 < c0
     7  a1 < c0            a1 < b2 < c0
     8  a1 < b2            a1 < b2 < c2
     9  a1                 a1 < c2
    10  b0 < c0            b0 < c0 < d0
    11  b1 < c0            b1 < c0 < d0
    12  b1 < d0            b1 < c1 < d0
    13  b1 < c1            b1 < c1 < d1
    14  b1                 b1 < d1
    15  b2 < c0            b2 < c0 < d0
    16  c0                 c0 < d0
    17  b2 < d0            b2 < c2 < d0
    18  b2 < c2            b2 < c2 < d1
    19  b2                 b2 < d1
    20  c1 < d0            b0 < c1 < d0
    21  b0 < c1            b0 < c1 < d1
    22  c1                 c1 < d1
    23  b0 < d0            b0 < c2 < d0
    24  d0                 c2 < d0
    25  b0 < c2            b0 < c2 < d1
    26  b0                 b0 < d1
    27  c2                 c2 < d1
"""


def parse_ranks(text):
    r = {}
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        head, rest = line.split(':', 1)
        lvl = int(head.split()[1])
        for name in rest.split():
            r[name] = lvl
    return r


def parse_covers(text):
    out = []
    for line in text.strip().splitlines():
        for tok in line.split():
            pass
        parts = line.replace('<', ' < ').split()
        i = 0
        while i + 2 < len(parts) + 0:
            if parts[i + 1] == '<':
                out.append((parts[i], parts[i + 2]))
                i += 3
            else:
                i += 1
    return out


def parse_collapse(text):
    """Each line is `<step> <tau chain> <sigma chain>`, chains written with
    '<' between element names.  Returns [(tau, sigma), ...] with each face a
    tuple of names sorted by name, which is this program's face convention.
    The printed step numbers are checked to run 1, 2, 3, ... in order."""
    out = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        head, rest = line.split(None, 1)
        step = int(head)
        if step != len(out) + 1:
            raise ValueError('collapse table step out of order: %r' % line)
        halves = rest.split('  ')
        halves = [h.strip() for h in halves if h.strip()]
        if len(halves) != 2:
            raise ValueError('collapse table row is not tau/sigma: %r' % line)
        faces = tuple(tuple(sorted(h.replace('<', ' ').split()))
                      for h in halves)
        out.append(faces)
    return out


RANKS = parse_ranks(RANK_TABLE)
COVERS = parse_covers(COVER_TABLE)
PRINTED_COLLAPSE = parse_collapse(COLLAPSE_TABLE)
ELEMS = sorted(RANKS)

# ---------------------------------------------------------------------------
# 1. THE CHECK HARNESS
# ---------------------------------------------------------------------------
_STATE = {'pass': 0, 'fail': 0}


def check(name, got, want, detail=''):
    """One check.  `got == want` decides it; both are printed on failure."""
    if got == want:
        _STATE['pass'] += 1
        print('PASS %s%s' % (name, ('  ' + detail) if detail else ''))
    else:
        _STATE['fail'] += 1
        print('FAIL %s  got %r  want %r%s'
              % (name, got, want, ('  ' + detail) if detail else ''))


# ---------------------------------------------------------------------------
# 2. POSET MACHINERY
# ---------------------------------------------------------------------------
def closure(covers):
    lt = set(covers)
    changed = True
    while changed:
        changed = False
        for (x, y) in list(lt):
            for (u, v) in list(lt):
                if y == u and (x, v) not in lt:
                    lt.add((x, v))
                    changed = True
    return lt


def cover_relation(lt):
    pts = set(x for p in lt for x in p)
    return set((x, y) for (x, y) in lt
               if not any((x, z) in lt and (z, y) in lt for z in pts))


def chains(elems, lt):
    """Every nonempty chain, as a tuple sorted by element NAME.

    Comparability is tested in BOTH directions, so nothing here assumes that
    the alphabetical order of the labels refines the order of the poset -- it
    happens to for the witness, and it does not for the coned posets of
    Section 5, where a silent assumption would have reported a disconnected
    complex.
    """
    out = []
    for r in range(1, len(elems) + 1):
        for c in itertools.combinations(sorted(elems), r):
            if all((c[i], c[j]) in lt or (c[j], c[i]) in lt
                   for i in range(len(c)) for j in range(i + 1, len(c))):
                out.append(c)
    return out


def spread(c, ranks):
    return max(ranks[x] for x in c) - min(ranks[x] for x in c)


def delta_k(elems, lt, ranks, k):
    return set(c for c in chains(elems, lt) if spread(c, ranks) <= k)


def induced(elems, lt):
    """The order complex of the sub-poset induced on `elems`."""
    sub = set((x, y) for (x, y) in lt if x in elems and y in elems)
    return set(chains(sorted(elems), sub))


# ---------------------------------------------------------------------------
# 3. SIMPLICIAL MACHINERY: f-vector, Euler characteristic, integral homology
# ---------------------------------------------------------------------------
def fvec(faces):
    by = {}
    for f in faces:
        by[len(f) - 1] = by.get(len(f) - 1, 0) + 1
    return tuple(by.get(d, 0) for d in range(max(by) + 1)) if by else ()


def chi(faces):
    return sum((-1) ** (len(f) - 1) for f in faces)


def snf(M):
    """Nonzero invariant factors of an integer matrix, by Smith normal form.
    Pure integer arithmetic; the length of the result is the rank."""
    A = [row[:] for row in M]
    m = len(A)
    n = len(A[0]) if A else 0
    inv, r, c = [], 0, 0
    while r < m and c < n:
        piv = None
        for i in range(r, m):
            for j in range(c, n):
                if A[i][j] != 0 and (piv is None
                                     or abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None:
            break
        pi, pj = piv
        A[r], A[pi] = A[pi], A[r]
        for row in A:
            row[c], row[pj] = row[pj], row[c]
        while True:
            done = True
            for i in range(r + 1, m):
                q = A[i][c] // A[r][c]
                if A[i][c] % A[r][c] == 0:
                    if q:
                        for j in range(c, n):
                            A[i][j] -= q * A[r][j]
                else:
                    for j in range(c, n):
                        A[i][j] -= q * A[r][j]
                    A[r], A[i] = A[i], A[r]
                    done = False
            for j in range(c + 1, n):
                q = A[r][j] // A[r][c]
                if A[r][j] % A[r][c] == 0:
                    if q:
                        for i in range(r, m):
                            A[i][j] -= q * A[i][c]
                else:
                    for i in range(r, m):
                        A[i][j] -= q * A[i][c]
                    for i in range(r, m):
                        A[i][c], A[i][j] = A[i][j], A[i][c]
                    done = False
            if done:
                break
        inv.append(abs(A[r][c]))
        r += 1
        c += 1
    return inv


def _by_dim(faces):
    by = {}
    for f in faces:
        by.setdefault(len(f) - 1, []).append(tuple(f))
    for d in by:
        by[d].sort()
    return by


def boundary(faces, d):
    """The matrix of partial_d : C_d -> C_{d-1}, plus the two index maps."""
    by = _by_dim(faces)
    rows = by.get(d - 1, [])
    cols = by.get(d, [])
    ri = {f: i for i, f in enumerate(rows)}
    M = [[0] * len(cols) for _ in rows]
    for j, f in enumerate(cols):
        for t in range(len(f)):
            g = f[:t] + f[t + 1:]
            M[ri[g]][j] = (-1) ** t
    return M, rows, cols


def homology(faces):
    """(betti list, torsion dict) of integral simplicial homology."""
    by = _by_dim(faces)
    top = max(by) if by else -1
    rank_d, invs = {}, {}
    for d in range(1, top + 1):
        M, rows, cols = boundary(faces, d)
        iv = snf(M) if rows and cols else []
        invs[d] = iv
        rank_d[d] = len(iv)
    betti, tor = [], {}
    for d in range(0, top + 1):
        betti.append(len(by.get(d, [])) - rank_d.get(d, 0) - rank_d.get(d + 1, 0))
        t = sorted(x for x in invs.get(d + 1, []) if x > 1)
        if t:
            tor[d] = t
    return betti, tor


def connected(faces):
    verts = [f[0] for f in faces if len(f) == 1]
    if not verts:
        return False
    adj = {v: set() for v in verts}
    for f in faces:
        if len(f) == 2:
            adj[f[0]].add(f[1])
            adj[f[1]].add(f[0])
    seen, stack = {verts[0]}, [verts[0]]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == len(verts)


# ---------------------------------------------------------------------------
# 4. COLLAPSING
# ---------------------------------------------------------------------------
def cofaces(faces, tau):
    s = set(tau)
    return [f for f in faces if len(f) == len(tau) + 1 and s.issubset(f)]


def greedy_collapse(faces, record=None):
    """Repeatedly delete a free face and its unique coface.  Deterministic:
    the candidate order is by (dimension, name)."""
    cur, steps = set(faces), 0
    while True:
        found = None
        for tau in sorted(cur, key=lambda f: (len(f), f)):
            cf = cofaces(cur, tau)
            if len(cf) == 1:
                found = (tau, cf[0])
                break
        if not found:
            break
        if record is not None:
            record.append(found)
        cur.discard(found[0])
        cur.discard(found[1])
        steps += 1
    return steps, sorted(cur)


def replay(faces, pairs):
    """Strict replay: every (tau, sigma) must have tau free with unique
    coface sigma in the complex that is actually present at that step."""
    cur = set(faces)
    for i, (tau, sig) in enumerate(pairs, 1):
        if cofaces(cur, tau) != [sig]:
            return False, i, sorted(cur)
        cur.discard(tau)
        cur.discard(sig)
    return True, len(pairs), sorted(cur)


# ---------------------------------------------------------------------------
# 5. INTEGER LATTICES (for the Mayer-Vietoris comparison map)
# ---------------------------------------------------------------------------
def as_matrix(vectors, dim):
    """Columns = the given coordinate vectors (dicts index -> coefficient)."""
    return [[v.get(i, 0) for v in vectors] for i in range(dim)]


def lat_invariants(vectors, dim):
    if not vectors:
        return (0, ())
    iv = snf(as_matrix(vectors, dim))
    return (len(iv), tuple(iv))


def same_lattice(gens_a, gens_b, dim):
    """True iff the two integer lattices coincide, GIVEN that one contains the
    other.  Equal rank and equal invariant factors force index 1."""
    return lat_invariants(gens_a, dim) == lat_invariants(gens_b, dim)


def chain_of_edges(edge_index, terms):
    """terms = [(sign, (u, v)), ...] with (u, v) already name-sorted."""
    v = {}
    for s, e in terms:
        v[edge_index[e]] = v.get(edge_index[e], 0) + s
    return {k: c for k, c in v.items() if c}


def fundamental_cycles(faces):
    """A Z-basis of the 1-cycle lattice: the fundamental cycles of a spanning
    forest.  Returned in edge coordinates."""
    verts = sorted(f[0] for f in faces if len(f) == 1)
    edges = sorted(f for f in faces if len(f) == 2)
    eidx = {e: i for i, e in enumerate(edges)}
    adj = {v: [] for v in verts}
    for e in edges:
        adj[e[0]].append(e)
        adj[e[1]].append(e)
    parent, tree, seen = {}, set(), set()
    for root in verts:
        if root in seen:
            continue
        seen.add(root)
        parent[root] = None
        stack = [root]
        while stack:
            x = stack.pop()
            for e in adj[x]:
                y = e[1] if e[0] == x else e[0]
                if y not in seen:
                    seen.add(y)
                    parent[y] = (x, e)
                    tree.add(e)
                    stack.append(y)

    def path_to_root(x):
        out = []
        while parent[x] is not None:
            p, e = parent[x]
            out.append((1 if e[0] == x else -1, e))   # oriented from x to p
            x = p
        return out

    cycles = []
    for e in edges:
        if e in tree:
            continue
        u, w = e
        # e oriented u -> w, then w -> root, then root -> u
        terms = [(1, e)]
        terms += [(s, ee) for s, ee in path_to_root(w)]
        terms += [(-s, ee) for s, ee in path_to_root(u)]
        cycles.append(chain_of_edges(eidx, terms))
    return cycles, eidx, edges


def boundary_2_columns(faces, eidx):
    """The boundaries of the 2-faces, in edge coordinates."""
    cols = []
    for f in sorted(f for f in faces if len(f) == 3):
        x, y, z = f
        cols.append(chain_of_edges(eidx, [(1, (y, z)), (-1, (x, z)), (1, (x, y))]))
    return cols


# ---------------------------------------------------------------------------
# 6. THE RUN
# ---------------------------------------------------------------------------
def main():
    print('=== Kitajima Remark 2.13 (arXiv:2606.15241v2, p. 8): '
          'a ten-element counterexample at k = 3 ===')
    print('object read from the paper: %d elements, %d printed covers'
          % (len(ELEMS), len(COVERS)))
    print()

    # ---- Section 2: the poset -------------------------------------------
    print('--- Section 2: the object is a graded poset of rank 3 ---')
    check('covers-parsed-15', len(COVERS), 15,
          'the fifteen pairs of the printed table')
    check('covers-distinct', len(set(COVERS)), 15)
    check('rank-step-one',
          sorted(set(RANKS[y] - RANKS[x] for (x, y) in COVERS)), [1],
          'every printed cover joins consecutive levels')
    check('level-sizes',
          tuple(sum(1 for e in ELEMS if RANKS[e] == r) for r in range(4)),
          (2, 3, 3, 2))

    lt = closure(COVERS)
    check('closure-31-pairs', len(lt), 31,
          'the fifteen covers plus sixteen further strict pairs')
    gaps = {}
    for (x, y) in lt:
        g = RANKS[y] - RANKS[x]
        gaps[g] = gaps.get(g, 0) + 1
    check('closure-rank-gap-profile', dict(sorted(gaps.items())),
          {1: 15, 2: 12, 3: 4})
    check('cover-relation-recovered', cover_relation(lt) == set(COVERS), True,
          'the covers of the closure are again exactly the printed fifteen')

    minimal = [x for x in ELEMS if not any(y == x for (_, y) in lt)]
    maximal = [x for x in ELEMS if not any(y == x for (y, _) in lt)]
    check('minimal-elements', minimal, ['a0', 'a1'])
    check('maximal-elements', maximal, ['d0', 'd1'])
    check('no-zero-hat-no-one-hat',
          (len(minimal) > 1, len(maximal) > 1), (True, True),
          'so Delta(P) is not a cone for the trivial reason')

    allch = chains(ELEMS, lt)
    maxch = [c for c in allch if not any(set(c) < set(d) for d in allch)]
    check('maximal-chains-count', len(maxch), 11)
    check('pure-of-rank-3', sorted(set(len(c) - 1 for c in maxch)), [3],
          'every maximal chain has length 3, so P is graded of rank n = 3')

    D = {k: delta_k(ELEMS, lt, RANKS, k) for k in (0, 1, 2, 3)}
    check('delta3-is-the-full-order-complex', D[3] == set(allch), True,
          'k = n = 3, so Delta^{(3)}(P) = Delta(P)')
    check('delta0-is-the-vertex-set', fvec(D[0]), (10,))
    check('filtration-is-increasing',
          (D[0] < D[1], D[1] < D[2], D[2] < D[3]), (True, True, True))

    # ---- the three filtration stages ------------------------------------
    print()
    print('--- the filtration stages, integral homology by Smith normal form ---')
    b1, t1 = homology(D[1])
    check('delta1-f-vector', fvec(D[1]), (10, 15))
    check('delta1-euler', chi(D[1]), -5)
    check('delta1-betti', b1, [1, 6], 'torsion %s' % (t1 or 'none'))
    check('delta1-NOT-contractible', b1 != [1], True,
          'so the witness is not swept up by the k = 2 sentence Kitajima '
          'proves one line earlier, whose hypothesis is Delta^{(1)} '
          'contractible')

    b2, t2 = homology(D[2])
    check('delta2-f-vector', fvec(D[2]), (10, 27, 18))
    check('delta2-face-count', len(D[2]), 55)
    check('delta2-euler', chi(D[2]), 1)
    check('delta2-betti-acyclic', (b2, t2), ([1, 0, 0], {}),
          'Delta^{(2)}(P) has the reduced homology of a point, and no torsion')

    b3, t3 = homology(D[3])
    check('delta3-f-vector', fvec(D[3]), (10, 31, 34, 11))
    check('delta3-face-count', len(D[3]), 86)
    check('delta3-euler', chi(D[3]), 2,
          'a contractible complex has Euler characteristic 1')
    check('delta3-betti-H2-is-Z', (b3, t3), ([1, 0, 1, 0], {}),
          'H_2(Delta(P); Z) = Z, so Delta^{(3)}(P) is NOT contractible')

    # ---- the headline ---------------------------------------------------
    print()
    print('--- the headline: contractibility is not monotone in k ---')
    # The PAPER'S OWN printed table is the primary certificate.  It is parsed
    # from COLLAPSE_TABLE above and replayed here; the program's independent
    # greedy search is then required to reproduce it, so this is a check ON THE
    # PAPER and not a replay of the program's own output.
    check('printed-collapse-table-has-27-rows', len(PRINTED_COLLAPSE), 27,
          'the table printed in Section 4 of the paper, transcribed and parsed')
    check('printed-collapse-faces-are-faces-of-delta2',
          sorted({f for pair in PRINTED_COLLAPSE for f in pair}
                 - set(D[2])), [],
          'every face named in the printed table really is a face of '
          'Delta^{(2)}(P), and there are 54 distinct ones')
    check('printed-collapse-names-54-distinct-faces',
          len({f for pair in PRINTED_COLLAPSE for f in pair}), 54)
    ok, n, rest = replay(D[2], PRINTED_COLLAPSE)
    check('printed-collapse-replays-strictly', (ok, n, len(rest)), (True, 27, 1),
          'THE PAPER\'S TABLE, replayed row by row: at its own step each '
          'printed tau has exactly one proper coface in the complex then '
          'present, and it is the printed sigma')
    check('printed-collapse-ends-at-one-vertex',
          (len(rest), rest[0] if rest else None), (1, ('d1',)),
          'after the paper\'s 27 rows exactly one face is left, the vertex d1')
    seq = []
    steps2, left2 = greedy_collapse(D[2], record=seq)
    check('delta2-collapsible', (steps2, len(left2)), (27, 1),
          '27 elementary free-face collapses remove 54 of the 55 faces, '
          'leaving the single vertex %s' % (left2[0][0] if left2 else '?'))
    check('greedy-search-reproduces-the-printed-table',
          seq == PRINTED_COLLAPSE, True,
          'an independent free-face search of Delta^{(2)}(P), run without '
          'reference to the table, returns the paper\'s 27 rows in the '
          'paper\'s order -- so paper and program cannot disagree silently')
    print('    the 27 collapse pairs AS PRINTED IN THE PAPER (tau -> sigma), '
          'each replayed above;')
    print('    each tau has exactly one proper coface at its step, and it is '
          'sigma, so a referee')
    print('    can redo the collapse by hand from this list:')
    for i, (tau, sig) in enumerate(PRINTED_COLLAPSE, 1):
        print('      %2d  %-12s ->  %s' % (i, '<'.join(tau), '<'.join(sig)))
    check('delta2-contractible-because-collapsible',
          (steps2, len(left2), b2), (27, 1, [1, 0, 0]),
          'collapsible => contractible, and the homology agrees')
    steps3, left3 = greedy_collapse(D[3])
    check('delta3-greedy-collapse-stalls', len(left3) > 1, True,
          'CONTROL: the same routine does NOT reduce Delta^{(3)} to a point '
          '(%d collapses, %d faces left), so it is not answering yes '
          'indiscriminately' % (steps3, len(left3)))
    check('THE-REFUTATION',
          (b2 == [1, 0, 0] and len(left2) == 1, b3[2]), (True, 1),
          'Delta^{(2)}(P) contractible and H_2(Delta^{(3)}(P)) = Z: the '
          'implication of Remark 2.13 fails at k = 3')

    # ---- Section 3: the Mayer-Vietoris proof for Delta^{(3)} ------------
    print()
    print('--- Section 3: Delta(P) = St(d0) cup St(d1), and the complex Y ---')
    check('d0-d1-incomparable',
          (('d0', 'd1') in lt, ('d1', 'd0') in lt), (False, False))
    check('no-face-contains-both-maxima',
          any('d0' in f and 'd1' in f for f in D[3]), False)
    below = {v: set(x for x in ELEMS if (x, v) in lt) for v in ('d0', 'd1')}
    check('elements-below-d0', sorted(below['d0']),
          ['a0', 'a1', 'b0', 'b1', 'b2', 'c0', 'c1', 'c2'])
    check('elements-below-d1', sorted(below['d1']),
          ['a0', 'a1', 'b0', 'b1', 'b2', 'c1', 'c2'])
    st = {}
    for v in ('d0', 'd1'):
        st[v] = set(f for f in D[3] if v in f or set(f) <= below[v])
    check('stars-cover-delta3', st['d0'] | st['d1'] == D[3], True)
    for v in ('d0', 'd1'):
        link = set(f for f in st[v] if v not in f)
        check('star-%s-is-a-cone' % v,
              all(tuple(sorted(f + (v,))) in st[v] for f in link), True,
              'every face of the link joins %s, so the star is contractible' % v)
    Y = st['d0'] & st['d1']
    ranks02 = set(x for x in ELEMS if RANKS[x] <= 2)
    check('Y-is-the-order-complex-of-P02-minus-c0',
          Y == induced(ranks02 - {'c0'}, lt), True)
    check('Y-f-vector', fvec(Y), (7, 11, 4))
    printed_Y_edges = ['a0 b0', 'a1 b1', 'a1 b2', 'b0 c1', 'b0 c2', 'b1 c1',
                       'b2 c2', 'a0 c1', 'a0 c2', 'a1 c1', 'a1 c2']
    check('Y-edges-match-the-printed-list',
          sorted(f for f in Y if len(f) == 2),
          sorted(tuple(s.split()) for s in printed_Y_edges))
    printed_Y_tris = ['a0 b0 c1', 'a0 b0 c2', 'a1 b1 c1', 'a1 b2 c2']
    check('Y-triangles-match-the-printed-list',
          sorted(f for f in Y if len(f) == 3),
          sorted(tuple(s.split()) for s in printed_Y_tris))
    check('Y-euler', chi(Y), 0)
    check('Y-connected', connected(Y), True)
    bY, tY = homology(Y)
    check('Y-betti', bY, [1, 1, 0], 'torsion %s' % (tY or 'none'))

    # Y itself is a union of two cones, over the two minimal elements, and
    # their intersection is S^0.  That is the second half of the suspension
    # chain S^0 -> Y -> Delta(P) of Section 3.
    above = {v: set(x for x in ELEMS if (v, x) in lt) for v in ('a0', 'a1')}
    stY = {}
    for v in ('a0', 'a1'):
        stY[v] = set(f for f in Y if v in f or set(f) <= above[v])
    check('Y-is-covered-by-two-cones', stY['a0'] | stY['a1'] == Y, True)
    for v in ('a0', 'a1'):
        link = set(f for f in stY[v] if v not in f)
        check('Y-star-%s-is-a-cone' % v,
              all(tuple(sorted(f + (v,))) in stY[v] for f in link), True)
    S0 = stY['a0'] & stY['a1']
    check('Y-cone-intersection-is-S0', sorted(S0), [('c1',), ('c2',)],
          'two points, so Y is homotopy equivalent to the suspension of S^0, '
          'that is to S^1')
    check('mayer-vietoris-shifts-Y-into-degree-2', (bY[1], b3[2]), (1, 1),
          'H_2(Delta P) = H_1(Y) = Z: Delta(P) is the union of two '
          'contractible stars along Y, hence homotopy equivalent to the '
          'suspension of Y, hence to S^2')

    # ---- Section 4: the Mayer-Vietoris proof for Delta^{(2)} ------------
    print()
    print('--- Section 4: Delta^{(2)} = A cup B with A cap B = Delta(P_[1,2]) ---')
    A = induced(set(x for x in ELEMS if RANKS[x] <= 2), lt)
    B = induced(set(x for x in ELEMS if RANKS[x] >= 1), lt)
    AB = induced(set(x for x in ELEMS if RANKS[x] in (1, 2)), lt)
    check('A-cup-B-is-delta2', A | B == D[2], True,
          'a chain of rank spread <= 2 lies in levels 0-2 or in levels 1-3')
    check('A-cap-B-is-the-middle', A & B == AB, True)
    check('A-f-vector', fvec(A), (8, 16, 7))
    check('A-euler', chi(A), -1)
    bA, tA = homology(A)
    check('A-betti-H1-is-Z2', (bA, tA), ([1, 2, 0], {}))
    stA = {}
    for v in ('a0', 'a1'):
        stA[v] = set(f for f in A if v in f or set(f) <= above[v])
    check('A-is-covered-by-two-cones', stA['a0'] | stA['a1'] == A, True)
    for v in ('a0', 'a1'):
        link = set(f for f in stA[v] if v not in f)
        check('A-star-%s-is-a-cone' % v,
              all(tuple(sorted(f + (v,))) in stA[v] for f in link), True)
    check('A-cone-intersection-is-three-points',
          sorted(stA['a0'] & stA['a1']), [('c0',), ('c1',), ('c2',)],
          'so A is homotopy equivalent to the suspension of three points, a '
          'wedge of two circles, which is why H_1(A) = Z^2')
    check('AB-f-vector', fvec(AB), (6, 7))
    check('AB-connected-H1-is-Z2', (connected(AB), homology(AB)[0]),
          (True, [1, 2]))
    check('B-f-vector', fvec(B), (8, 18, 11))
    bB, tB = homology(B)
    check('B-acyclic', (bB, tB), ([1, 0, 0], {}),
          'B is contractible; the proof in the paper is that it is the union '
          'of the two cones St_B(d0), St_B(d1) along a tree')
    stB = {}
    for v in ('d0', 'd1'):
        stB[v] = set(f for f in B if v in f or set(f) <= below[v])
    check('B-is-covered-by-its-two-cones', stB['d0'] | stB['d1'] == B, True)
    for v in ('d0', 'd1'):
        link = set(f for f in stB[v] if v not in f)
        check('B-star-%s-is-a-cone' % v,
              all(tuple(sorted(f + (v,))) in stB[v] for f in link), True)
    T = stB['d0'] & stB['d1']
    check('B-cone-intersection-is-a-tree', fvec(T), (5, 4),
          'five vertices b0,b1,b2,c1,c2 and four edges: connected and acyclic, '
          'hence contractible')
    check('B-cone-intersection-connected-acyclic',
          (connected(T), homology(T)[0]), (True, [1, 0]))

    # phi : H_1(A cap B) -> H_1(A) is an isomorphism, checked as a lattice
    # identity in the 1-chain group of A.
    cyc, eidx, edgesA = fundamental_cycles(A)
    dimA = len(edgesA)
    check('A-cycle-lattice-rank', len(cyc), 9,
          '16 edges - 8 vertices + 1 component')
    g1 = chain_of_edges(eidx, [(1, ('b1', 'c1')), (-1, ('b0', 'c1')),
                               (1, ('b0', 'c0')), (-1, ('b1', 'c0'))])
    g2 = chain_of_edges(eidx, [(1, ('b2', 'c2')), (-1, ('b0', 'c2')),
                               (1, ('b0', 'c0')), (-1, ('b2', 'c0'))])
    d2cols = boundary_2_columns(A, eidx)
    check('A-has-seven-2-faces', len(d2cols), 7)

    def is_cycle(v):
        # boundary of a 1-chain, in vertex coordinates
        acc = {}
        for i, c in v.items():
            u, w = edgesA[i]
            acc[w] = acc.get(w, 0) + c
            acc[u] = acc.get(u, 0) - c
        return not any(acc.values())

    check('gamma1-is-a-4-cycle-in-A-cap-B',
          (is_cycle(g1), len(g1),
           all(edgesA[i] in [tuple(sorted(f)) for f in AB if len(f) == 2]
               for i in g1)), (True, 4, True))
    check('gamma2-is-a-4-cycle-in-A-cap-B',
          (is_cycle(g2), len(g2),
           all(edgesA[i] in [tuple(sorted(f)) for f in AB if len(f) == 2]
               for i in g2)), (True, 4, True))

    r_tri = lat_invariants(d2cols, dimA)
    r_full = lat_invariants(d2cols + [g1, g2], dimA)
    check('boundaries-of-A-are-independent', r_tri[0], 7)
    check('phi-is-injective', r_full[0] - r_tri[0], 2,
          'gamma1, gamma2 are independent modulo the 2-boundaries of A')
    check('phi-is-surjective',
          same_lattice(d2cols + [g1, g2], cyc, dimA), True,
          'the 2-boundaries together with gamma1, gamma2 generate the whole '
          '1-cycle lattice of A -- equal rank and equal invariant factors, '
          'hence index 1')
    check('phi-is-an-isomorphism',
          (r_full[0] - r_tri[0] == 2,
           same_lattice(d2cols + [g1, g2], cyc, dimA)), (True, True),
          'so Mayer-Vietoris gives H_2(Delta^{(2)}) = ker phi = 0 and '
          'H_1(Delta^{(2)}) = coker phi = 0')

    # the two printed chain identities
    e1 = chain_of_edges(eidx, [(1, ('a1', 'c1')), (-1, ('a0', 'c1')),
                               (1, ('a0', 'c0')), (-1, ('a1', 'c0'))])
    e2 = chain_of_edges(eidx, [(1, ('a1', 'c2')), (-1, ('a0', 'c2')),
                               (1, ('a0', 'c0')), (-1, ('a1', 'c0'))])
    for nm, g, e in (('gamma1', g1, e1), ('gamma2', g2, e2)):
        diff = dict(g)
        for i, c in e.items():
            diff[i] = diff.get(i, 0) - c
        diff = {k: c for k, c in diff.items() if c}
        check('printed-identity-for-%s' % nm,
              same_lattice(d2cols, d2cols + [diff], dimA), True,
              'the difference is a 2-boundary of A, so the printed '
              'congruence modulo the a<b<c triangles holds')

    # ---- Section 5: the scope claim -------------------------------------
    print()
    print('--- Section 5: this P cannot be coned into a witness at k < n ---')
    for tag, extra_cov, extra_rank, want, why in (
            ('adjoin-a-minimum',
             [('m', x) for x in ELEMS if RANKS[x] == 0], 0, [1, 1, 0],
             'H_1 = Z'),
            ('adjoin-a-maximum',
             [(x, 'M') for x in ELEMS if RANKS[x] == 3], 4, [1, 0, 1],
             'H_2 = Z')):
        new = 'm' if tag.endswith('minimum') else 'M'
        r2 = dict((k, v + (1 if extra_rank == 0 else 0)) for k, v in RANKS.items())
        r2[new] = extra_rank if extra_rank else 0
        cov2 = [(x, y) for (x, y) in COVERS] + extra_cov
        lt2 = closure(cov2)
        el2 = sorted(r2)
        ch2 = chains(el2, lt2)
        mx2 = [c for c in ch2 if not any(set(c) < set(d) for d in ch2)]
        check('scope-%s-is-pure-of-rank-4' % tag,
              sorted(set(len(c) - 1 for c in mx2)), [4],
              'so it is a legitimate rank-4 graded poset and k = 2 < n = 4 '
              'would be a proper filtration stage')
        Dk = delta_k(el2, lt2, r2, 2)
        bb, tt = homology(Dk)
        check('scope-%s-breaks-delta2' % tag, (bb, tt), (want, {}),
              '%s: the coned poset has rank 4, but its Delta^{(2)} is no '
              'longer contractible, so the hypothesis of the implication is '
              'lost and this P cannot be lifted to a witness at k < n' % why)

    # ---- controls -------------------------------------------------------
    print()
    print('--- controls, both polarities ---')
    cr = {'x0': 0, 'x1': 1, 'x2': 2, 'x3': 3}
    clt = closure([('x0', 'x1'), ('x1', 'x2'), ('x2', 'x3')])
    got = []
    for k in (1, 2, 3):
        F = delta_k(sorted(cr), clt, cr, k)
        got.append((fvec(F), greedy_collapse(F)[0], len(greedy_collapse(F)[1])))
    check('CTRL+chain-collapsible-at-every-k', got,
          [((4, 3), 3, 1), ((4, 5, 2), 5, 1), ((4, 6, 4, 1), 7, 1)],
          'FORCED POSITIVE: for the 3-chain every stage is collapsible, so '
          'the collapser can say yes')

    kr = {'p0': 0, 'p1': 0, 'q0': 1, 'q1': 1}
    klt = closure([('p0', 'q0'), ('p0', 'q1'), ('p1', 'q0'), ('p1', 'q1')])
    F = delta_k(sorted(kr), klt, kr, 1)
    check('CTRL-crown-not-collapsible',
          (fvec(F), homology(F)[0], len(greedy_collapse(F)[1]) == 1),
          ((4, 4), [1, 1], False),
          'FORCED NEGATIVE: the 2+2 crown has Delta^{(1)} a 4-cycle')

    pr, pc = {}, []
    for i in range(4):
        for j in (0, 1):
            pr['e%d%d' % (i, j)] = i
    for i in range(3):
        for j in (0, 1):
            for j2 in (0, 1):
                pc.append(('e%d%d' % (i, j), 'e%d%d' % (i + 1, j2)))
    F = delta_k(sorted(pr), closure(pc), pr, 2)
    bP, _ = homology(F)
    check('CTRL=kitajima-published-formula', (fvec(F), bP[2]), ((8, 20, 16), 3),
          'FORCED SILENT: for Kitajima\'s Example P_n with n = 4 levels of two, '
          'Delta^{(2)} must be a wedge of 2n-2k-1 = 3 copies of S^2, and it is '
          '-- an independent engine check against a published value')

    rp2 = [(1, 2, 3), (1, 3, 4), (1, 4, 5), (1, 5, 6), (1, 2, 6), (2, 3, 5),
           (3, 4, 6), (2, 4, 5), (3, 5, 6), (2, 4, 6)]
    faces = set()
    for f in rp2:
        for r in range(1, 4):
            faces.update(itertools.combinations(sorted(f), r))
    bR, tR = homology(faces)
    check('CTRL-torsion-not-blind', (bR, tR.get(1)), ([1, 0, 0], [2]),
          'the six-vertex minimal triangulation of the real projective plane '
          'has H_1 = Z/2; the engine reports it, so "no torsion" above is a '
          'measurement and not a blind spot')

    # ---- what this program does NOT cover -------------------------------
    print()
    print('NOT RE-RUN: the exhaustive census behind the minimality remark of '
          'Section 5 (all rank-3 graded posets with every level of size at '
          'most three: 22,930,191 labelled posets, 49,843 isomorphism '
          'classes, 201 witness classes). That search ran elsewhere, its '
          'stdout was not preserved, and nothing here reproduces it. The '
          'paper therefore asserts no minimality, and the refutation does '
          'not depend on it.')
    print('NOT RE-RUN: the correction of Section 6 to Kitajima\'s Remark 4.8, '
          'p. 15 (three printed groups Z^12 that must '
          'read Z^11). That is a long-exact-sequence argument on the source\'s '
          'own printed integers, is logically independent of the refutation, '
          'and is not a computation this program performs.')
    print('NOT RE-RUN: any search of the literature. Whether the k = 3 case '
          'was settled elsewhere is not a question a program can answer, and '
          'no claim of priority is checked here.')
    print('NOT RE-RUN: the open case k < n. Nothing here bears on whether the '
          'implication of Remark 2.13 holds at a proper filtration stage; the '
          'witness has k = n = 3.')
    print()

    n = _STATE['pass']
    if _STATE['fail']:
        print('VERDICT: %d of %d CHECKS FAILED' % (_STATE['fail'], n + _STATE['fail']))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % n)
    return 0


if __name__ == '__main__':
    sys.exit(main())
