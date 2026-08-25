#!/usr/bin/env python3
"""INDEPENDENT verification of the two graph6 strings printed in the accompanying
paper (paper.tex): the 19-vertex triangulation T and the 20-vertex planar graph W.

The strings are decoded from scratch, T is rebuilt from the construction the paper
describes in words, both graphs are tested for induced-P6-freeness by two mutually
independent algorithms, and planarity is certified by exhibiting a sphere face set.

Standard library only. Integer arithmetic only. graph6 decoder written from the
format specification (McKay, "Description of graph6, sparse6 and digraph6 encodings").

Everything is printed; nothing is asserted silently.
"""
import itertools
import sys

# ----------------------------------------------------------------------------- input
G6_T = "R~vMLaWbCKM?U?R?Go@a?EG?KO?KO?"
G6_W = "S~vMLaWbCKM?U?R?Go@a?EG?KO?KOE???"

NAMES19 = ['a', 'b', 'c0', 'c1', 'c2', 'c3', 'c4',
           'u0', 'u1', 'u2', 'u3', 'v0', 'v1', 'v2', 'v3',
           'w0', 'w1', 'w2', 'w3']
NAMES20 = NAMES19 + ['x']

FAIL = []


def line(s=''):
    print(s, flush=True)


def check(label, ok, detail=''):
    tag = 'OK  ' if ok else 'FAIL'
    if not ok:
        FAIL.append(label + (' :: ' + detail if detail else ''))
    line('  [%s] %s%s' % (tag, label, ('   ' + detail if detail else '')))
    return ok


# ------------------------------------------------------------------- graph6 decoding
def pair_index(i, j):
    """graph6 bit position of the pair {i,j}, i<j.  Column-major upper triangle:
    (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),...  -> k = j(j-1)/2 + i."""
    if i > j:
        i, j = j, i
    return j * (j - 1) // 2 + i


def decode_g6(s, label):
    """Strict graph6 decode.  Returns (n, edges_set, ok_flag)."""
    line('--- graph6 STRICT decode of %s: %r  (len=%d)' % (label, s, len(s)))
    ok = True

    # 0. no header, no whitespace
    ok &= check('%s: no >>graph6<< header present (none expected)' % label,
                not s.startswith('>>graph6<<'))
    ok &= check('%s: no whitespace / newline characters' % label,
                all(not c.isspace() for c in s))

    # 1. every byte in the legal printable range 63..126
    bad = [(k, c, ord(c)) for k, c in enumerate(s) if not (63 <= ord(c) <= 126)]
    ok &= check('%s: every byte in legal range [63,126]' % label, not bad,
                'offenders=%r' % bad if bad else 'min=%d max=%d' %
                (min(ord(c) for c in s), max(ord(c) for c in s)))
    if bad:
        return None, None, False

    # 2. order field N(n)
    b0 = ord(s[0]) - 63
    if b0 == 63:                       # 126 == '~' escape
        if len(s) < 4:
            check('%s: multi-byte order field truncated' % label, False)
            return None, None, False
        n = 0
        for c in s[1:4]:
            n = (n << 6) | (ord(c) - 63)
        data = s[4:]
        line('      order field: multi-byte ("~" escape), n=%d' % n)
    else:
        n = b0
        data = s[1:]
        line('      order field: single byte %r -> n = %d - 63 = %d'
             % (s[0], ord(s[0]), n))
        ok &= check('%s: single-byte order field requires n<=62' % label, n <= 62)

    # 3. exactly the right number of data bytes
    nbits = n * (n - 1) // 2
    need = (nbits + 5) // 6
    ok &= check('%s: data byte count' % label, len(data) == need,
                'n=%d -> %d bits -> ceil/6 = %d data bytes; string has %d'
                % (n, nbits, need, len(data)))
    if len(data) != need:
        return n, None, False

    # 4. bit vector
    bits = []
    for c in data:
        v = ord(c) - 63
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    assert len(bits) == 6 * need

    # 5. trailing padding bits must be zero
    pad = bits[nbits:]
    ok &= check('%s: %d trailing padding bit(s) all zero' % (label, len(pad)),
                all(p == 0 for p in pad), 'pad=%r' % pad)

    edges = set()
    for j in range(1, n):
        for i in range(j):
            if bits[pair_index(i, j)]:
                edges.add((i, j))
    line('      decoded: n=%d  m=%d' % (n, len(edges)))
    return n, edges, bool(ok)


def encode_g6(n, edges):
    """Re-encode, so we can compare byte-for-byte against the printed string."""
    nbits = n * (n - 1) // 2
    bits = [0] * nbits
    for (i, j) in edges:
        bits[pair_index(i, j)] = 1
    while len(bits) % 6:
        bits.append(0)
    if n <= 62:
        out = chr(n + 63)
    else:
        out = chr(126) + ''.join(chr(63 + ((n >> k) & 63)) for k in (12, 6, 0))
    for t in range(0, len(bits), 6):
        v = 0
        for b in bits[t:t + 6]:
            v = (v << 1) | b
        out += chr(63 + v)
    return out


# ------------------------------------------------------------- graph helpers (stdlib)
def adj_of(n, edges):
    A = [set() for _ in range(n)]
    for (i, j) in edges:
        A[i].add(j)
        A[j].add(i)
    return A


def components(n, A, alive=None):
    alive = set(range(n)) if alive is None else set(alive)
    seen, comps = set(), []
    for s in sorted(alive):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in A[v]:
                if w in alive and w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(sorted(comp))
    return comps


def connectivity_at_least_3(n, A):
    """Exhaustive: connected, no cut vertex, no 2-cut.  Returns (bool, witness)."""
    if len(components(n, A)) != 1:
        return False, 'disconnected'
    for v in range(n):
        if len(components(n, A, set(range(n)) - {v})) != 1:
            return False, '1-cut {%d}' % v
    for v, w in itertools.combinations(range(n), 2):
        rest = set(range(n)) - {v, w}
        if len(rest) >= 2 and len(components(n, A, rest)) != 1:
            return False, '2-cut {%d,%d}' % (v, w)
    return True, 'no 1-cut and no 2-cut exists (all %d vertices and %d pairs tested)' % (
        n, n * (n - 1) // 2)


def all_triangles(n, A):
    return [t for t in itertools.combinations(range(n), 3)
            if t[1] in A[t[0]] and t[2] in A[t[0]] and t[2] in A[t[1]]]


# --------------------------------------------------------------- induced P6 detection
def induced_p6_by_subsets(n, A, names):
    """Exhaustive over all 6-subsets.  A 6-subset induces P6 iff the induced graph is
    connected, has 5 edges, and has degree multiset {1,1,2,2,2,2}."""
    found = []
    for S in itertools.combinations(range(n), 6):
        Sset = set(S)
        deg = {}
        m = 0
        for v in S:
            d = len(A[v] & Sset)
            deg[v] = d
            m += d
        m //= 2
        if m != 5:
            continue
        if sorted(deg.values()) != [1, 1, 2, 2, 2, 2]:
            continue
        # connected?
        start = S[0]
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in A[v] & Sset:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        if len(seen) != 6:
            continue
        # recover path order
        ends = [v for v in S if deg[v] == 1]
        path = [ends[0]]
        prev = None
        cur = ends[0]
        while len(path) < 6:
            nxt = [w for w in A[cur] & Sset if w != prev][0]
            path.append(nxt)
            prev, cur = cur, nxt
        found.append(tuple(path))
    return found


def induced_p6_by_extension(n, A):
    """Second, structurally different implementation: grow a path z1..zt, extending only
    to y in N(z_t) with N(y) disjoint from {z1..z_{t-1}}.  Counts each P6 twice
    (once per direction)."""
    out = []

    def grow(path):
        if len(path) == 6:
            out.append(tuple(path))
            return
        zt = path[-1]
        pref = set(path[:-1])
        for y in sorted(A[zt]):
            if y in path:
                continue
            if A[y] & pref:
                continue
            path.append(y)
            grow(path)
            path.pop()

    for v in range(n):
        grow([v])
    return out


# ------------------------------------------ sphere-embedding (planarity) verification
def find_face_sets(n, A, edges, cap=None):
    """Search for every set F of triangles with EVERY edge in EXACTLY TWO members.
    Such an F is then tested elsewhere for being a closed surface (every vertex link a
    single cycle) and for Euler characteristic 2.  A connected closed surface with chi=2
    is the sphere, so such an F is a genuine plane triangulation embedding -- a rigorous
    planarity certificate that does not depend on how the paper's wording is read.

    Solutions are returned as a list of DISTINCT face sets.  The backtracking below
    assigns one triangle at a time, so it reaches one and the same face set by many
    different orders of choice: giving a deficient edge its two faces as (t1,t2) or as
    (t2,t1) doubles the number of completing paths.  Counting completions therefore
    measures the search order, not the number of embeddings -- on T it reports 64 = 2**6
    completions of a face set that is in fact unique.  We deduplicate on the chosen
    triangle set so that len(sols) is a count of embeddings.

    `cap` bounds the number of DISTINCT solutions collected; the second return value
    says whether that bound cut the search short, in which case the list is NOT
    exhaustive and no uniqueness conclusion may be drawn from it.  cap=None searches
    exhaustively."""
    tris = all_triangles(n, A)
    tedges = []
    for t in tris:
        tedges.append([(t[0], t[1]), (t[0], t[2]), (t[1], t[2])])
    e_tris = {}
    for e in edges:
        e_tris[tuple(sorted(e))] = []
    for ti, es in enumerate(tedges):
        for e in es:
            e_tris[e].append(ti)
    cnt = {e: 0 for e in e_tris}
    chosen = []
    sols = []
    seen = set()
    truncated = []

    def bt():
        if cap is not None and len(sols) >= cap:
            truncated.append(True)      # the cap, not exhaustion, ended this branch
            return
        # any edge still short of 2 faces?
        best, bestcands = None, None
        for e, c in cnt.items():
            if c < 2:
                cands = [ti for ti in e_tris[e]
                         if ti not in chosen_set and all(cnt[x] < 2 for x in tedges[ti])]
                if len(cands) < 2 - c:
                    return                      # dead
                if best is None or len(cands) < len(bestcands):
                    best, bestcands = e, cands
        if best is None:
            key = frozenset(chosen)             # a face SET, independent of choice order
            if key not in seen:
                seen.add(key)
                sols.append(sorted(tris[ti] for ti in chosen))
            return
        for ti in bestcands:
            chosen.append(ti)
            chosen_set.add(ti)
            for x in tedges[ti]:
                cnt[x] += 1
            bt()
            for x in tedges[ti]:
                cnt[x] -= 1
            chosen_set.discard(ti)
            chosen.pop()

    chosen_set = set()
    bt()
    return tris, sols, bool(truncated)


def edges_in_exactly_two_faces(edges, faces):
    inc = {}
    for f in faces:
        for e in itertools.combinations(sorted(f), 2):
            inc[e] = inc.get(e, 0) + 1
    return all(inc.get(tuple(sorted(e)), 0) == 2 for e in edges) and len(inc) == len(edges)


def links_are_single_cycles(n, A, faces):
    """Every vertex link a single cycle == the face set is a closed surface, not just a
    pseudomanifold with a pinch point."""
    for v in range(n):
        link = {}
        for f in faces:
            if v in f:
                o = [z for z in f if z != v]
                link.setdefault(o[0], []).append(o[1])
                link.setdefault(o[1], []).append(o[0])
        if set(link) != A[v] or any(len(x) != 2 for x in link.values()):
            return False
        # single cycle?
        start = next(iter(link))
        prev, cur, cnt = None, start, 0
        while True:
            nxt = link[cur][0] if link[cur][0] != prev else link[cur][1]
            prev, cur = cur, nxt
            cnt += 1
            if cur == start:
                break
        if cnt != len(A[v]):
            return False
    return True


def is_sphere_face_set(n, A, edges, faces):
    """Quiet form of surface_report: does `faces` embed the graph in the SPHERE?"""
    return (edges_in_exactly_two_faces(edges, faces)
            and links_are_single_cycles(n, A, faces)
            and n - len(edges) + len(faces) == 2)


def surface_report(n, A, edges, faces, names):
    line('      candidate face set: %d faces' % len(faces))
    ok1 = edges_in_exactly_two_faces(edges, faces)
    check('every edge lies in exactly two faces', ok1)
    ok2 = links_are_single_cycles(n, A, faces)
    check('every vertex link is a SINGLE cycle (=> closed surface)', ok2)
    chi = n - len(edges) + len(faces)
    ok3 = check('Euler characteristic n-m+f = 2 (=> sphere, not a torus)', chi == 2,
                'n-m+f = %d-%d+%d = %d' % (n, len(edges), len(faces), chi))
    return ok1 and ok2 and ok3


# ------------------------------------------------------------------- isomorphism test
def refine_invariant(n, A):
    inv = [len(A[v]) for v in range(n)]
    for _ in range(n):
        new = [(inv[v], tuple(sorted(inv[w] for w in A[v]))) for v in range(n)]
        codes = {c: i for i, c in enumerate(sorted(set(new)))}
        nxt = [codes[c] for c in new]
        if nxt == inv:
            break
        inv = nxt
    return inv


def isomorphic(n, A, B):
    if sorted(len(x) for x in A) != sorted(len(x) for x in B):
        return None
    ia, ib = refine_invariant(n, A), refine_invariant(n, B)
    if sorted(ia) != sorted(ib):
        return None
    order = sorted(range(n), key=lambda v: (-len(A[v]), ia[v]))
    mapping = {}
    used = set()

    def bt(k):
        if k == n:
            return True
        v = order[k]
        for w in range(n):
            if w in used or len(B[w]) != len(A[v]) or ib[w] != ia[v]:
                continue
            good = True
            for u in order[:k]:
                if (u in A[v]) != (mapping[u] in B[w]):
                    good = False
                    break
            if not good:
                continue
            mapping[v] = w
            used.add(w)
            if bt(k + 1):
                return True
            used.discard(w)
            del mapping[v]
        return False

    return dict(mapping) if bt(0) else None


# ============================================================================== MAIN
line('=' * 78)
line('independent check of the two graph6 strings')
line('=' * 78)

# ---- (a) decode both -------------------------------------------------------------
nT, eT, okT = decode_g6(G6_T, 'T')
line()
nW, eW, okW = decode_g6(G6_W, 'W')
line()

AT = adj_of(nT, eT)
AW = adj_of(nW, eW)

line('DECODED T: n=%d  m=%d' % (nT, len(eT)))
line('  degrees: ' + ', '.join('%s=%d' % (NAMES19[v], len(AT[v])) for v in range(nT)))
line('  degree sequence sorted: %s  sum=%d' %
     (sorted(len(x) for x in AT), sum(len(x) for x in AT)))
line('  adjacency (names):')
for v in range(nT):
    line('    %-3s -> %s' % (NAMES19[v], ' '.join(NAMES19[w] for w in sorted(AT[v]))))
line()
line('DECODED W: n=%d  m=%d' % (nW, len(eW)))
line('  degrees: ' + ', '.join('%s=%d' % (NAMES20[v], len(AW[v])) for v in range(nW)))
line('  adjacency of x (index 19): %s' % ' '.join(NAMES20[w] for w in sorted(AW[19])))
line()

line('--- (b) claim-by-claim tests')
line('order of T: T has 19 vertices')
check('decoded T has n = 19', nT == 19, 'n=%d' % nT)
line('size of T: T has 51 edges')
check('decoded T has m = 51', len(eT) == 51, 'm=%d' % len(eT))
check('51 == 3n-6 for n=19', len(eT) == 3 * nT - 6, '3*19-6 = %d' % (3 * nT - 6))
line('order and size of W: W has 20 vertices and 53 edges')
check('decoded W has n = 20', nW == 20, 'n=%d' % nW)
check('decoded W has m = 53', len(eW) == 53, 'm=%d' % len(eW))
line('how W is obtained from T: W = T + x with edges xa, xb only')
check('W restricted to vertices 0..18 equals T exactly (labelled)',
      set((i, j) for (i, j) in eW if j <= 18) == eT)
check('x (index 19) has degree exactly 2', len(AW[19]) == 2, 'deg=%d' % len(AW[19]))
check('x is adjacent to a and b, and NOT to c0',
      AW[19] == {0, 1}, 'N(x) = %s' % sorted(NAMES20[w] for w in AW[19]))
line('the target value: Euler bound at n=20, and the conjectured value')
check('3n-6 at n=20 equals 54', 3 * 20 - 6 == 54, '3*20-6 = %d' % (3 * 20 - 6))
check('floor(50/17*(20-2))+1 == 53', (50 * 18) // 17 + 1 == 53,
      'floor(900/17)+1 = %d+1 = %d' % (900 // 17, 900 // 17 + 1))
check('m(W) equals the conjectured 53', len(eW) == 53)
line()

# ---- (c) build T from the construction described in the paper ---------------------
line('--- (c) building T INDEPENDENTLY from the construction described in the paper')
IX = {nm: k for k, nm in enumerate(NAMES19)}
pe = set()


def ae(x, y):
    i, j = sorted((IX[x], IX[y]))
    if (i, j) in pe:
        raise SystemExit('prose recipe would add a duplicate edge %s-%s' % (x, y))
    pe.add((i, j))


# "Begin with a path c0c1c2c3c4."
for i in range(4):
    ae('c%d' % i, 'c%d' % (i + 1))
# "Add adjacent vertices a,b, and join each of them to every c_i."
ae('a', 'b')
for i in range(5):
    ae('a', 'c%d' % i)
    ae('b', 'c%d' % i)
prose_faces = [('a', 'c%d' % i, 'c%d' % (i + 1)) for i in range(4)] + \
              [('b', 'c%d' % i, 'c%d' % (i + 1)) for i in range(4)] + \
              [('a', 'b', 'c0'), ('a', 'b', 'c4')]
line('      base graph a,b,c0..c4 : n=7 m=%d f=%d  chi=%d'
     % (len(pe), len(prose_faces), 7 - len(pe) + len(prose_faces)))


def insert(vertex, face):
    """Insert `vertex` into triangular face `face`: 3 new edges, face -> 3 faces."""
    if tuple(sorted(face)) not in {tuple(sorted(f)) for f in prose_faces}:
        raise SystemExit('prose recipe inserts %s into %r which is NOT a face' % (vertex, face))
    for z in face:
        ae(vertex, z)
    fs = [f for f in prose_faces if tuple(sorted(f)) == tuple(sorted(face))][0]
    prose_faces.remove(fs)
    prose_faces.append((vertex, face[0], face[1]))
    prose_faces.append((vertex, face[0], face[2]))
    prose_faces.append((vertex, face[1], face[2]))


# "For 0<=i<=3, insert u_i into the face a c_i c_{i+1} ..."
for i in range(4):
    insert('u%d' % i, ('a', 'c%d' % i, 'c%d' % (i + 1)))
# "... and v_i into the face b c_i c_{i+1};"
for i in range(4):
    insert('v%d' % i, ('b', 'c%d' % i, 'c%d' % (i + 1)))
# "then insert w_i into the face u_i c_i c_{i+1}."
for i in range(4):
    insert('w%d' % i, ('u%d' % i, 'c%d' % i, 'c%d' % (i + 1)))

AP = adj_of(19, pe)
line('      PROSE T: n=19  m=%d  faces=%d  chi=%d'
     % (len(pe), len(prose_faces), 19 - len(pe) + len(prose_faces)))
check('prose T has 51 edges', len(pe) == 51, 'm=%d' % len(pe))
check('prose T plane embedding has 34 triangular faces', len(prose_faces) == 34,
      'f=%d' % len(prose_faces))
check('prose T is a SPHERE triangulation (chi=2)',
      19 - len(pe) + len(prose_faces) == 2)
line('      degrees of prose T: ' +
     ', '.join('%s=%d' % (NAMES19[v], len(AP[v])) for v in range(19)))

line('  >>> the decisive comparison')
same_labelled = (pe == eT)
check('prose T equals decoded T as LABELLED graphs (edge sets identical)', same_labelled)
if not same_labelled:
    only_prose = sorted(pe - eT)
    only_dec = sorted(eT - pe)
    line('      edges in PROSE but not in the string (%d): %s' %
         (len(only_prose), ', '.join('%s-%s' % (NAMES19[i], NAMES19[j]) for i, j in only_prose)))
    line('      edges in the STRING but not in prose (%d): %s' %
         (len(only_dec), ', '.join('%s-%s' % (NAMES19[i], NAMES19[j]) for i, j in only_dec)))
    mp = isomorphic(19, AP, AT)
    check('prose T and decoded T are at least ISOMORPHIC', mp is not None,
          ('mapping prose->decoded: ' +
           ', '.join('%s->%s' % (NAMES19[k], NAMES19[v]) for k, v in sorted(mp.items())))
          if mp else 'no isomorphism exists')

line('  >>> byte-for-byte re-encode of the prose graph')
reT = encode_g6(19, pe)
line('      paper T : %s' % G6_T)
line('      mine  T : %s' % reT)
check('re-encoding prose T reproduces the paper T string exactly', reT == G6_T)
prose_W = set(pe) | {(0, 19), (1, 19)}
reW = encode_g6(20, prose_W)
line('      paper W : %s' % G6_W)
line('      mine  W : %s' % reW)
check('re-encoding prose W (= prose T + x~a,b) reproduces the paper W string exactly',
      reW == G6_W)
line('      round-trip of decoded strings: T %s   W %s'
     % (encode_g6(nT, eT) == G6_T, encode_g6(nW, eW) == G6_W))
line()

# ---- (d) induced P6 --------------------------------------------------------------
line('--- (d) induced P6 (INDUCED, not merely a subgraph)')
for nm, n_, A_, names in (('T', nT, AT, NAMES19), ('W', nW, AW, NAMES20)):
    sub = induced_p6_by_subsets(n_, A_, names)
    ext = induced_p6_by_extension(n_, A_)
    line('  %s: 6-subset scan over C(%d,6)=%d subsets -> %d induced P6'
         % (nm, n_, len(list(itertools.combinations(range(n_), 6))), len(sub)))
    line('  %s: independent path-extension scan -> %d directed = %d undirected'
         % (nm, len(ext), len(ext) // 2))
    check('%s: the two implementations agree' % nm,
          len(ext) == 2 * len(sub) and
          {frozenset(p) for p in ext} == {frozenset(p) for p in sub})
    check('%s is induced-P6-FREE' % nm, len(sub) == 0,
          'count=%d' % len(sub))
    for p in sub[:10]:
        line('      induced P6 witness: %s' % '-'.join(names[v] for v in p))
line()

# ---- (e) is T a maximal planar graph? -------------------------------------------
line('--- (e) is T a maximal planar graph (a triangulation)?')
check('m = 3n-6', len(eT) == 3 * nT - 6, 'm=%d, 3n-6=%d' % (len(eT), 3 * nT - 6))
mind = min(len(x) for x in AT)
check('minimum degree >= 3 (every simple triangulation has this)', mind >= 3,
      'delta=%d' % mind)
line('      min degree = %d, max degree = %d' % (mind, max(len(x) for x in AT)))
c3, wit = connectivity_at_least_3(nT, AT)
check('T is 3-CONNECTED (exhaustive 1-cut and 2-cut search)', c3, wit)
tris = all_triangles(nT, AT)
line('      triangle (3-clique) count of T = %d ; a plane triangulation on 19 '
     'vertices has 2m/3 = %d faces' % (len(tris), 2 * len(eT) // 3))
check('triangle count >= number of faces (34)', len(tris) >= 34,
      'triangles=%d, faces would be %d' % (len(tris), 34))
line('      separating/non-facial triangles = %d' % (len(tris) - 34))
line('  >>> PLANARITY: searched from scratch for a sphere face set of the DECODED graph')
tris2, sols, truncated = find_face_sets(nT, AT, eT)     # cap=None: exhaustive
line('      exhaustive search over triangle sets with every edge in exactly 2: '
     '%d DISTINCT solution(s) found (search %s)'
     % (len(sols), 'TRUNCATED by cap' if truncated else 'ran to exhaustion'))
if sols:
    good = surface_report(nT, AT, eT, sols[0], NAMES19)
    check('DECODED T is a genuine plane (sphere) triangulation -- planarity PROVED, '
          'not inferred', good)
    # Whitney asserts a unique embedding in the SPHERE, so filter the pseudomanifold
    # solutions down to genuine sphere face sets before asserting uniqueness.  As
    # unordered vertex triples a face set is reflection-invariant, so "unique up to
    # reflection and relabelling" is literally "one face set" for this fixed labelling.
    spheres = [s for s in sols if is_sphere_face_set(nT, AT, eT, s)]
    line('      of those, %d is/are genuine SPHERE face sets (links single cycles, '
         'chi=2)' % len(spheres))
    check('the sphere face set is unique (Whitney: 3-connected planar has a unique '
          'embedding)', len(spheres) == 1 and not truncated,
          '%d sphere face set(s) from a search that %s'
          % (len(spheres), 'was TRUNCATED' if truncated else 'ran to exhaustion'))
    fs = {tuple(sorted(f)) for f in sols[0]}
    ps = {tuple(sorted(IX[z] for z in f)) for f in prose_faces}
    allfs = [{tuple(sorted(f)) for f in s} for s in sols]
    check('the face set my prose embedding produced is among the solutions found',
          ps in allfs, 'sym-diff vs sol[0] = %d' % len(fs ^ ps))
    check('T has exactly 34 faces => exactly 34 facial insertions',
          len(sols[0]) == 34, 'f=%d' % len(sols[0]))
    abc0 = tuple(sorted((IX['a'], IX['b'], IX['c0'])))
    check('a-b-c0 IS a face of T, so x can legally be inserted there',
          abc0 in fs)
else:
    check('DECODED T is a plane triangulation', False,
          'NO valid sphere face set exists -- T is NOT a planar triangulation')
line()

# ---- (f) does W support the headline lower bound? -------------------------------
line('--- (f) W and the headline lower bound ex_P(20, P6^ind) >= 53')
line('      W: n = %d, m = %d' % (nW, len(eW)))
c3w, witw = connectivity_at_least_3(nW, AW)
line('      W 3-connected? %s (%s)  [not required; x has degree 2 so 2-connected '
     'at best]' % (c3w, witw))
if sols:
    wfaces = [f for f in sols[0] if tuple(sorted(f)) != abc0]
    wfaces = [tuple(NAMES19[z] for z in f) for f in wfaces]
    fW = len(wfaces) + 2      # abc0 splits into triangle abx and quadrilateral a-x-b-c0
    check('W is PLANAR: T sphere embedding with x drawn inside face a-b-c0; '
          'chi = n-m+f = 2', nW - len(eW) + fW == 2,
          '20 - 53 + %d = %d' % (fW, nW - len(eW) + fW))
subW = induced_p6_by_subsets(nW, AW, NAMES20)
check('W is induced-P6-free', len(subW) == 0, 'count=%d' % len(subW))
check('W therefore certifies ex_P(20,P6^ind) >= 53', len(eW) == 53 and len(subW) == 0)
line('      Euler bound 3n-6 at n=20 = 54, so W is exactly ONE edge below it.')
line()

# ---- the paper's further remark: repeated degree-3 deletion ----------------------
line('--- further remark in the paper: repeatedly deleting a degree-3 vertex from T '
     'reaches every smaller n')
reach = set()
seen = set()
stack = [frozenset(range(nT))]
seen.add(stack[0])
budget = 400000
while stack and budget > 0:
    S = stack.pop()
    budget -= 1
    reach.add(len(S))
    if len(S) <= 4:
        continue
    for v in S:
        if len(AT[v] & set(S)) == 3:
            T2 = S - {v}
            if T2 not in seen:
                seen.add(T2)
                stack.append(T2)
line('      orders reachable by successive degree-3 deletions: %s'
     % sorted(reach, reverse=True))
check('every order 4..19 reachable', set(range(4, 20)) <= reach,
      'missing %s' % sorted(set(range(4, 20)) - reach))
line('      (induced-P6-freeness is hereditary, so freeness of the smaller '
     'triangulations is automatic once T is free -- the only content is that a '
     'degree-3 vertex always remains.)')
line()

line('=' * 78)
if FAIL:
    line('VERDICT: %d CHECK(S) FAILED' % len(FAIL))
    for f in FAIL:
        line('   FAILED: %s' % f)
else:
    line('VERDICT: ALL CHECKS PASSED')
line('=' * 78)
line('summary numbers:')
line('  T: n=%d m=%d  induced-P6 count=%d' % (nT, len(eT), len(induced_p6_by_subsets(nT, AT, NAMES19))))
line('  W: n=%d m=%d  induced-P6 count=%d' % (nW, len(eW), len(subW)))
line('  well-formed: T=%s  W=%s' % (okT, okW))

# exit 0 if and only if every check passed, so a failing verdict cannot be read as
# success by a script that only looks at the exit status
raise SystemExit(1 if FAIL else 0)
