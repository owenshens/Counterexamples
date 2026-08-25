#!/usr/bin/env python3
"""Independently verify the two graph6 strings printed in the induced planar P6-Turan paper.

TAKEN FROM THE PAPER (inputs):
  * the two graph6 strings themselves;
  * the prose recipe for building T, transcribed into build_T() below.

DERIVED HERE (the actual checks):
  * n and m of each decoded graph;
  * whether the decoded T equals the graph the recipe builds -- tested BOTH as label-equality and,
    independently, by re-encoding the recipe graph to graph6 and comparing the string byte for byte;
  * 3-connectivity of T by exhaustive 2-cut search;
  * PLANARITY of T, proved rather than inferred: an exhaustive search for a set of triangles
    meeting every edge exactly twice, then the closed-surface test (every vertex link a single
    cycle) and the Euler characteristic n - m + f = 2.  Together with m = 3n-6 that makes T a
    maximal planar graph.  The edge count ALONE would not: K_{3,3} plus a triangle on one side
    is 3-connected with n = 6 and m = 12 = 3n-6, and is not planar;
  * planarity of W, from that embedding with the new vertex drawn inside the face a b c_0;
  * induced-P6-freeness of both graphs by exhaustive search over all 6-subsets;
  * that W is exactly T plus one vertex joined to a and b and NOT to c_0, as the paper states;
  * that each of the 34 facial insertions of T contains an induced P6, with a witness path
    printed and independently re-verified for each of the 34.

NOT DERIVED HERE: the plantri census behind the upper bound.  See the closing NOT RE-RUN line,
which states the shortfall in full.

Standard library only, integer arithmetic only.
"""
from itertools import combinations

G6_T = 'R~vMLaWbCKM?U?R?Go@a?EG?KO?KO?'
G6_W = 'S~vMLaWbCKM?U?R?Go@a?EG?KO?KOE???'

checks = []


def ck(name, ok, detail=''):
    checks.append((name, ok))
    print(('PASS ' if ok else 'FAIL ') + name + (('  [' + detail + ']') if detail else ''), flush=True)
    return ok


# ---------------------------------------------------------------- graph6 codec
def g6_decode(s):
    """graph6 -> (n, set of frozenset edges). Implemented from the format spec, not a library."""
    b = [ord(c) - 63 for c in s]
    for i, v in enumerate(b):
        if not 0 <= v <= 63:
            raise ValueError(f'byte {i} out of range: {s[i]!r}')
    if b[0] == 63:                       # '~' prefix: n in the next three 6-bit groups
        n = (b[1] << 12) | (b[2] << 6) | b[3]
        data = b[4:]
    else:
        n = b[0]
        data = b[1:]
    nbits = n * (n - 1) // 2
    need = (nbits + 5) // 6
    if len(data) != need:
        raise ValueError(f'n={n} needs {need} data bytes, got {len(data)}')
    bits = []
    for v in data:
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    # padding bits beyond nbits must be zero, else the string is malformed
    pad_ok = all(x == 0 for x in bits[nbits:])
    E = set()
    idx = 0
    for j in range(1, n):                # column-major upper triangle
        for i in range(j):
            if bits[idx]:
                E.add(frozenset((i, j)))
            idx += 1
    return n, E, pad_ok


def g6_encode(n, E):
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if frozenset((i, j)) in E else 0)
    while len(bits) % 6:
        bits.append(0)
    out = chr(n + 63) if n <= 62 else '~' + ''.join(
        chr(((n >> sh) & 63) + 63) for sh in (12, 6, 0))
    for i in range(0, len(bits), 6):
        v = 0
        for bit in bits[i:i + 6]:
            v = (v << 1) | bit
        out += chr(v + 63)
    return out


# ---------------------------------------------------------------- graph helpers
def adj(n, E):
    A = {i: set() for i in range(n)}
    for e in E:
        i, j = tuple(e)
        A[i].add(j)
        A[j].add(i)
    return A


def connected(verts, A):
    verts = set(verts)
    if not verts:
        return True
    st = [next(iter(verts))]
    seen = {st[0]}
    while st:
        for w in A[st.pop()] & verts:
            if w not in seen:
                seen.add(w)
                st.append(w)
    return seen == verts


def is_3_connected(n, E):
    """No cutvertex and no 2-cut, by exhaustive removal. n is 19, so this is cheap and exact."""
    A = adj(n, E)
    V = set(range(n))
    if not connected(V, A):
        return False, 'disconnected'
    for k in (1, 2):
        for S in combinations(range(n), k):
            rest = V - set(S)
            if len(rest) > 1 and not connected(rest, A):
                return False, f'{k}-cut {S}'
    return True, 'no cutvertex, no 2-cut'


def induced_p6_count(n, E):
    """Number of 6-subsets whose induced subgraph is a path on 6 vertices."""
    A = adj(n, E)
    cnt, witness = 0, None
    for S in combinations(range(n), 6):
        Ss = set(S)
        deg = [len(A[v] & Ss) for v in S]
        if sum(deg) != 10:                       # a path on 6 vertices has 5 edges
            continue
        if sorted(deg) != [1, 1, 2, 2, 2, 2]:
            continue
        if not connected(Ss, A):
            continue
        cnt += 1
        if witness is None:
            witness = S
    return cnt, witness


def triangles(n, E):
    A = adj(n, E)
    return sum(1 for a, b, c in combinations(range(n), 3)
               if b in A[a] and c in A[a] and c in A[b])


def all_triangles(n, E):
    A = adj(n, E)
    return [t for t in combinations(range(n), 3)
            if t[1] in A[t[0]] and t[2] in A[t[0]] and t[2] in A[t[1]]]


# ------------------------------------------------- planarity by exhibited sphere embedding
def find_face_sets(n, E):
    """Every set F of triangles of the graph in which EVERY edge lies in EXACTLY TWO members,
    by exhaustive backtracking with no cap, so the enumeration is complete.

    The backtracking assigns one triangle at a time and so reaches one and the same set by
    many orders of choice; solutions are therefore deduplicated on the chosen triangle SET,
    which is what makes len(sols) a count of face sets rather than of search paths.

    Such an F is only a candidate embedding: it is tested afterwards by face_system_verdict
    for being a closed surface of Euler characteristic 2.  That, and not the edge count, is
    what proves planarity here."""
    tris = all_triangles(n, E)
    tedges = [[frozenset((t[0], t[1])), frozenset((t[0], t[2])), frozenset((t[1], t[2]))]
              for t in tris]
    e_tris = {e: [] for e in E}
    for ti, es in enumerate(tedges):
        for e in es:
            e_tris[e].append(ti)
    cnt = {e: 0 for e in E}
    chosen, chosen_set, sols, seen = [], set(), [], set()

    def bt():
        best, bestc = None, None
        for e in cnt:                            # the edge with the fewest live candidates
            if cnt[e] < 2:
                cands = [ti for ti in e_tris[e]
                         if ti not in chosen_set and all(cnt[x] < 2 for x in tedges[ti])]
                if len(cands) < 2 - cnt[e]:
                    return                       # this edge can never reach two faces: dead
                if bestc is None or len(cands) < len(bestc):
                    best, bestc = e, cands
        if best is None:                         # every edge has exactly two faces
            key = frozenset(chosen)
            if key not in seen:
                seen.add(key)
                sols.append(sorted(tris[ti] for ti in chosen))
            return
        for ti in bestc:
            chosen.append(ti)
            chosen_set.add(ti)
            for x in tedges[ti]:
                cnt[x] += 1
            bt()
            for x in tedges[ti]:
                cnt[x] -= 1
            chosen_set.discard(ti)
            chosen.pop()

    bt()
    return tris, sols


def face_system_verdict(n, E, faces):
    """Is `faces` -- each a tuple of vertices in cyclic boundary order, of any length >= 3 --
    a SPHERE embedding of the graph (n, E)?  The conditions checked are:
      (0) the graph is connected;
      (1) every face is a closed walk on distinct vertices, using edges of the graph;
      (2) every edge lies on exactly two face boundaries, and no other pair does;
      (3) at every vertex the corners (consecutive edge pairs of the incident faces) form a
          SINGLE cycle on that vertex's edges, so the identification space is a closed
          surface with no pinch point;
      (4) n - m + f = 2.
    A connected closed surface of Euler characteristic 2 is the sphere, so (0)-(4) together
    are a planarity certificate.  Returns (ok, chi, why)."""
    A = adj(n, E)
    m, f = len(E), len(faces)
    chi = n - m + f
    if not connected(range(n), A):
        return False, chi, 'the graph is disconnected'
    for fc in faces:
        if len(fc) < 3 or len(set(fc)) != len(fc):
            return False, chi, f'malformed face {fc}'
        for k in range(len(fc)):
            if fc[(k + 1) % len(fc)] not in A[fc[k]]:
                return False, chi, f'face {fc} uses a non-edge'
    inc = {}
    for fc in faces:
        for k in range(len(fc)):
            e = frozenset((fc[k], fc[(k + 1) % len(fc)]))
            inc[e] = inc.get(e, 0) + 1
    if set(inc) != set(E):
        return False, chi, 'the face boundaries do not cover exactly the edge set'
    bad = [e for e in E if inc[e] != 2]
    if bad:
        return False, chi, f'{len(bad)} edge(s) do not lie in exactly two faces'
    for v in range(n):
        ev = {frozenset((v, w)) for w in A[v]}
        corners = []
        for fc in faces:
            L = len(fc)
            for k in range(L):
                if fc[k] == v:
                    corners.append((frozenset((v, fc[k - 1])),
                                    frozenset((v, fc[(k + 1) % L]))))
        if len(corners) != len(ev):
            return False, chi, f'vertex {v} has {len(corners)} corners for degree {len(ev)}'
        touched = {e: 0 for e in ev}
        parent = {e: e for e in ev}

        def root(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for e1, e2 in corners:
            if e1 not in touched or e2 not in touched or e1 == e2:
                return False, chi, f'degenerate corner at vertex {v}'
            touched[e1] += 1
            touched[e2] += 1
            r1, r2 = root(e1), root(e2)
            if r1 != r2:
                parent[r1] = r2
        if any(c != 2 for c in touched.values()):
            return False, chi, f'the link of vertex {v} is not 2-regular'
        if len({root(e) for e in ev}) != 1:
            return False, chi, f'the link of vertex {v} is disconnected (pinch point)'
    if chi != 2:
        return False, chi, f'n-m+f = {n}-{m}+{f} = {chi}, not 2'
    return True, chi, (f'every edge in exactly 2 faces, every link a single cycle, '
                       f'n-m+f = {n}-{m}+{f} = 2')


def induced_p6_witness(n, E):
    """The first 6-subset inducing a path on six vertices, returned in path order, else None."""
    A = adj(n, E)
    for S in combinations(range(n), 6):
        Ss = set(S)
        deg = {v: len(A[v] & Ss) for v in S}
        if sum(deg.values()) != 10 or sorted(deg.values()) != [1, 1, 2, 2, 2, 2]:
            continue
        if not connected(Ss, A):
            continue
        cur = min(v for v in S if deg[v] == 1)
        path = [cur]
        while len(path) < 6:
            nxt = [w for w in A[cur] & Ss if w not in path]
            if not nxt:
                break
            path.append(nxt[0])
            cur = path[-1]
        return tuple(path)
    return None


def verifies_as_induced_p6(E, path):
    """Independent re-check of a printed witness: six distinct vertices, the five consecutive
    pairs are edges, and NO other pair among the six is an edge."""
    if path is None or len(path) != 6 or len(set(path)) != 6:
        return False
    want = {frozenset((path[k], path[k + 1])) for k in range(5)}
    have = {frozenset(p) for p in combinations(path, 2) if frozenset(p) in E}
    return len(want) == 5 and have == want


# ---------------------------------------------------------------- the paper's recipe for T
def build_T():
    """c0..c4 a path; a,b adjacent and joined to every c_i; u_i in face a c_i c_{i+1};
    v_i in face b c_i c_{i+1}; then w_i in face u_i c_i c_{i+1}."""
    # The vertex INDEXING is the one the paper states beside the strings: "with the vertices in the
    # order a, b, c_0..c_4, u_0..u_3, v_0..v_3, w_0..w_3, and x last". A first run of this check
    # invented its own order; label-equality then failed on 60 edges while the degree fingerprint
    # matched -- the graphs were fine, the assumed labelling was not.
    order = (['a', 'b'] + ['c%d' % i for i in range(5)] + ['u%d' % i for i in range(4)]
             + ['v%d' % i for i in range(4)] + ['w%d' % i for i in range(4)])
    name = {v: i for i, v in enumerate(order)}
    E = set()

    def V(x):
        return name[x]

    def add(x, y):
        E.add(frozenset((V(x), V(y))))

    for i in range(4):
        add(f'c{i}', f'c{i+1}')
    add('a', 'b')
    for i in range(5):
        add('a', f'c{i}')
        add('b', f'c{i}')
    for i in range(4):
        for t, anchor in (('u', 'a'), ('v', 'b')):
            add(f'{t}{i}', anchor)
            add(f'{t}{i}', f'c{i}')
            add(f'{t}{i}', f'c{i+1}')
    for i in range(4):
        add(f'w{i}', f'u{i}')
        add(f'w{i}', f'c{i}')
        add(f'w{i}', f'c{i+1}')
    return len(name), E, name


def canon_key(n, E):
    """Cheap isomorphism-invariant fingerprint: sorted multiset of sorted neighbour-degree lists."""
    A = adj(n, E)
    return sorted(tuple(sorted(len(A[w]) for w in A[v])) for v in range(n))


print('=' * 78)
print('T :', G6_T)
print('W :', G6_W)
print('=' * 78)

nT, ET, padT = g6_decode(G6_T)
nW, EW, padW = g6_decode(G6_W)
print(f'decoded T: n={nT} m={len(ET)}   decoded W: n={nW} m={len(EW)}')
print()

ck('T is well-formed graph6 (padding bits zero)', padT)
ck('W is well-formed graph6 (padding bits zero)', padW)
ck('T has 19 vertices, as the paper states', nT == 19, f'n={nT}')
ck('T has 51 edges, as the paper states', len(ET) == 51, f'm={len(ET)}')
ck('T has m = 3n-6, the edge count of a maximal planar graph (planarity itself is proved '
   'below, not inferred from this count)', len(ET) == 3 * nT - 6, f'3n-6={3*nT-6}')
ck('W has 20 vertices, as the paper states', nW == 20, f'n={nW}')
ck('W has 53 edges, as the paper states', len(EW) == 53, f'm={len(EW)}')

# --- does the decoded T agree with the paper's own recipe?
nR, ER, name = build_T()
print()
print(f'recipe graph: n={nR} m={len(ER)}')
ck('the recipe builds 19 vertices', nR == 19, f'n={nR}')
ck('the recipe builds 51 edges', len(ER) == 51, f'm={len(ER)}')
ck('decoded T is LABEL-EQUAL to the recipe graph', ET == ER,
   f'symmetric difference {len(ET ^ ER)}')
ck('re-encoding the recipe graph reproduces the printed string byte for byte',
   g6_encode(nR, ER) == G6_T, g6_encode(nR, ER))
ck('decoded T and the recipe graph share a degree fingerprint',
   canon_key(nT, ET) == canon_key(nR, ER))

# --- structure of T
ok3, why3 = is_3_connected(nT, ET)
ck('T is 3-connected', ok3, why3)
AT = adj(nT, ET)
degs = sorted(len(AT[v]) for v in range(nT))
print(f'       T degree sequence: {degs}')
print(f'       T min degree = {degs[0]}, triangles = {triangles(nT, ET)}, 2n-4 = {2*nT-4}')

# --- PLANARITY of T, searched for from scratch rather than inferred from the edge count
print()
_tris, _sols = find_face_sets(nT, ET)
_spheres = [s for s in _sols if face_system_verdict(nT, ET, s)[0]]
print(f'       T has {len(_tris)} triangles as subgraphs; a 19-vertex plane triangulation has '
      f'2m/3 = {2*len(ET)//3} faces')
print(f'       exhaustive search over triangle sets meeting every edge exactly twice: '
      f'{len(_sols)} distinct solution(s), {len(_spheres)} of them a sphere face set')
ck('a set of triangles of T exists with every edge of T in exactly two of them',
   len(_sols) >= 1, f'{len(_sols)} such set(s) found')
if _spheres:
    FT = _spheres[0]
    okPT, chiT, whyPT = face_system_verdict(nT, ET, FT)
else:
    FT, okPT, chiT = None, False, None
    whyPT = 'the exhaustive search found NO sphere face set, so T is not planar'
ck('T is PLANAR: the exhibited face set is a closed surface with Euler characteristic 2, '
   'i.e. the sphere', okPT, whyPT)
ck('T is a maximal planar graph: planar AND m = 3n-6',
   okPT and len(ET) == 3 * nT - 6, f'planar={okPT}, m={len(ET)}, 3n-6={3*nT-6}')
ck('the sphere embedding of T is unique, the search having run to exhaustion',
   len(_spheres) == 1, f'{len(_spheres)} sphere face set(s) from an uncapped search')
ck('T has exactly 34 faces, hence exactly 34 facial insertions',
   FT is not None and len(FT) == 34, f'f={len(FT) if FT is not None else 0}')
_abc0 = tuple(sorted((name['a'], name['b'], name['c0'])))
ck('a, b, c_0 bound a FACE of T, so the new vertex x may legally be inserted there',
   FT is not None and _abc0 in FT, f'face {_abc0} in the embedding: '
   f'{FT is not None and _abc0 in FT}')

# --- W as T plus one vertex
AW = adj(nW, EW)
extra = [v for v in range(nW) if v not in range(nT)]
newv = nW - 1
Wminus = {e for e in EW if newv not in e}
ck('W restricted to the first 19 vertices is exactly T', Wminus == ET,
   f'symmetric difference {len(Wminus ^ ET)}')
nb = sorted(AW[newv])
inv = {v: k for k, v in name.items()}
ck('the added vertex has exactly two neighbours', len(nb) == 2, str([inv.get(v, v) for v in nb]))
ck("the added vertex is joined to a and b", {inv.get(v) for v in nb} == {'a', 'b'},
   str(sorted(inv.get(v, v) for v in nb)))
ck("the added vertex is NOT joined to c_0, as the paper says",
   name['c0'] not in AW[newv])
ck('a, b and c_0 are mutually adjacent in T, so abc_0 is a triangle of T',
   frozenset((name['a'], name['b'])) in ET
   and frozenset((name['a'], name['c0'])) in ET
   and frozenset((name['b'], name['c0'])) in ET)

# --- the load-bearing property
print()
cT, wT = induced_p6_count(nT, ET)
ck('T is induced-P6-free', cT == 0, f'induced P6 count = {cT}' + (f', e.g. {wT}' if wT else ''))
cW, wW = induced_p6_count(nW, EW)
ck('W is induced-P6-free', cW == 0, f'induced P6 count = {cW}' + (f', e.g. {wW}' if wW else ''))
# --- planarity of W: the sphere embedding of T with x drawn inside the face a b c_0.  The
#     edges xa and xb cut that triangle into the triangle a b x and the quadrilateral
#     a x b c_0, so W has 34 - 1 + 2 = 35 faces.  Planarity is part of the planar Turan
#     bound, so it belongs in the predicate of the ex_P check below, not in its prose.
if FT is not None and _abc0 in FT:
    _a, _b, _c0 = name['a'], name['b'], name['c0']
    FW = [tuple(t) for t in FT if tuple(t) != _abc0]
    FW.append((_a, _b, newv))                      # the triangle a b x
    FW.append((_a, newv, _b, _c0))                 # the quadrilateral a x b c_0
    okPW, chiW, whyPW = face_system_verdict(nW, EW, FW)
    detPW = f'f={len(FW)} (33 triangles + a b x + the quadrilateral a x b c_0), {whyPW}'
else:
    okPW, chiW, detPW = False, None, ('the sphere embedding of T that this needs is '
                                      'unavailable, so planarity of W is NOT established')
ck('W is PLANAR: the sphere embedding of T with x drawn inside the face a b c_0', okPW, detPW)

ck('W is planar, has 53 edges on 20 vertices and is induced-P6-free, so ex_P(20,P6^ind) >= 53',
   nW == 20 and len(EW) == 53 and cW == 0 and okPW,
   f'n={nW}, m={len(EW)}, induced P6 count={cW}, planar={okPW}')
ck("m(W) is one below Euler's bound 3n-6 = 54",
   3 * nW - 6 == 54 and len(EW) + 1 == 3 * nW - 6, f'm(W)={len(EW)}, 3n-6={3*nW-6}')
ck("the conjectured value floor(50(n-2)/17)+1 at n=20 equals m(W)",
   (50 * (nW - 2)) // 17 + 1 == len(EW),
   f'floor(50*18/17)+1 = floor(900/17)+1 = {(50*(nW-2))//17}+1 = {(50*(nW-2))//17 + 1}, '
   f'm(W)={len(EW)}')

# --- the 34 facial insertions of T, each of which the paper says contains an induced P6.
#     A witness path is printed for every one and then re-verified from the edge set.  If the
#     embedding above failed, no insertion is enumerated and BOTH checks below FAIL: absence
#     of the input is a failure here, never a silently skipped check.
print()
print('the 34 facial insertions of T (the inserted vertex is called z, index 19):')
ins_total = ins_with = ins_verified = 0


def _nm(z):
    return inv.get(z, 'z')


if FT is not None:
    for _fc in sorted(FT):
        ins_total += 1
        _Ein = set(ET) | {frozenset((nT, _z)) for _z in _fc}
        _w = induced_p6_witness(nT + 1, _Ein)
        if _w is not None:
            ins_with += 1
            if verifies_as_induced_p6(_Ein, _w):
                ins_verified += 1
        print('       face %-12s  n=%d m=%d  induced P6: %s'
              % ('-'.join(_nm(t) for t in _fc), nT + 1, len(_Ein),
                 ' - '.join(_nm(t) for t in _w) if _w is not None else 'NONE FOUND'))
ck('each of the 34 facial insertions of T contains an induced P6',
   ins_total == 34 and ins_with == 34, f'{ins_with} of {ins_total} insertions contain one')
ck('each printed witness re-verifies as an induced P6 (five path edges, no chord)',
   ins_total == 34 and ins_verified == 34, f'{ins_verified} of {ins_total} re-verified')

nf = sum(1 for _, o in checks if not o)
print()
print(f'VERDICT: {"ALL " + str(len(checks)) + " CHECKS PASS" if nf == 0 else str(nf) + " OF " + str(len(checks)) + " CHECKS FAILED"}')
print('NOT RE-RUN: the plantri census behind the UPPER bound -- no plantri run was made here, '
      'none of the 78,435,562 minimum-degree-four isomorphism classes over 6 <= n <= 20 was '
      'generated or tested, neither published count 11,284,042 nor 64,719,885 was reproduced, '
      'and the recursion giving (|L_n|)_{n=4..20}, in particular L_19 = {T} and L_20 = empty, '
      'was not re-run; the 34 facial insertions of T ARE checked above, but they close the '
      'degree-three branch of the upper bound only given L_19 = {T}, which is not verified here.')
raise SystemExit(1 if nf else 0)
