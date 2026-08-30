#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# verify.py -- checker for
#   "A Completely Uniform Nested Steiner Quadruple System of Order 56"
#
# Python 3.9+, STANDARD LIBRARY ONLY (itertools, base64, hashlib, sys).  No
# third-party package, no external data file, no solver, no network.  Exact
# integer arithmetic throughout; there is not a single floating-point number
# in this program and therefore no float decision.
#
# EVERY input it consumes is printed in the paper: the thirteen base blocks of
# Section 3, the two development/doubling recipes of Section 3, and the 29
# base-64 lines of Section 4.  Nothing is read from disk.  The program
# re-derives every quantity the paper asserts, prints one `PASS <name>` line
# per check, states what it did NOT cover, and exits 0 iff every check passed.
#
#   python3 verify.py            (about one second, single process)
# ---------------------------------------------------------------------------
import base64
import hashlib
import itertools
import sys

# ---------------------------------------------------------------------------
# The object, exactly as printed in the paper
# ---------------------------------------------------------------------------

# Paper, Section 3, display (3.1): the thirteen base 4-sets on {0,...,13}.
SQS14_BASE = (
    (0, 1, 10, 11), (0, 1, 4, 2), (0, 2, 9, 7), (0, 4, 13, 10), (0, 7, 4, 11),
    (0, 8, 7, 1), (0, 8, 13, 2), (0, 9, 11, 8), (0, 9, 1, 5), (0, 10, 12, 2),
    (0, 12, 9, 4), (0, 13, 12, 1), (7, 13, 12, 10),
)

# Paper, Section 4: the witness, as the 29 printed lines.  The base-64 string
# is their concatenation with no separator of any kind.
WITNESS_LINES = (
    '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',
    'CqJTpqzgP4gW0/KuKrI5mO92ezE8LFaBl3bkY3c=',
)

# Paper, Section 4: the SHA-256 digests printed there, and the redundant
# decimal prefix of the trit vector (three printed lines, concatenated).
WITNESS_B64_SHA256 = 'f802277ca20184e1b026bddc89ed1b30ed791885fa832ef17b657fcc5b2970b8'
WITNESS_TRIT_SHA256 = 'e5a8368b4a0b9f21d4e8e44be375a1cc9e00a364804e2b2366fa04e60aaf4016'
TRIT_PREFIX_DECIMAL = (
    '20110111202222101002110202021011002010221201210012021202220021000102010122211210021121200002'
    '111021121201022022101020012211221221212220212210100110101000112212100122112120002020112200211'
    '220001120110000'
)

# Paper, Section 4: the first twenty nested blocks, printed there in full.
FIRST_20_PRINTED = (
    ((0, 4), (1, 2)), ((0, 1), (3, 6)), ((0, 5), (1, 9)), ((0, 7), (1, 8)),
    ((0, 1), (10, 11)), ((0, 12), (1, 13)), ((0, 14), (1, 15)), ((0, 16), (1, 26)),
    ((0, 25), (1, 17)), ((0, 1), (18, 24)), ((0, 23), (1, 19)), ((0, 22), (1, 20)),
    ((0, 27), (1, 21)), ((0, 29), (1, 28)), ((0, 30), (1, 54)), ((0, 1), (31, 53)),
    ((0, 32), (1, 52)), ((0, 1), (33, 51)), ((0, 1), (34, 50)), ((0, 49), (1, 35)),
)

NBLOCKS = 6930          # b(56) = 56*55*54/24
MULT = 9                # (56-2)/6

# ---------------------------------------------------------------------------
# Check bookkeeping
# ---------------------------------------------------------------------------
_passed = []
_failed = []


def check(name, ok, detail=''):
    """Record one check.  `ok` must be a genuine bool from exact arithmetic."""
    if ok:
        _passed.append(name)
        print('PASS %s%s' % (name, (' | ' + detail) if detail else ''))
    else:
        _failed.append(name)
        # deliberately not the token the harness scans for a passing run
        print('CHECK-FAILED %s%s' % (name, (' | ' + detail) if detail else ''))
    return bool(ok)


# ---------------------------------------------------------------------------
# The two recipes of Section 3
# ---------------------------------------------------------------------------
def develop_z7(block):
    """x -> x+1 (mod 7) acting on Z_7 x {0,1} identified with {0,...,13} by
    (x,y) -> x + 7y.  Seven translates of a base block."""
    out = []
    for t in range(7):
        out.append(tuple(sorted(((x % 7 + t) % 7) + 7 * (x // 7) for x in block)))
    return out


def round_robin_1_factorization(n):
    """The 1-factorization of K_n printed in Section 3: vertex n-1 plays
    infinity, and for i = 0..n-2,
        F_i = {{n-1,i}} u {{(i+j) mod (n-1), (i-j) mod (n-1)} : j = 1..(n-2)/2}."""
    m = n - 1
    factors = []
    for i in range(m):
        f = [tuple(sorted((m, i)))]
        for j in range(1, (m - 1) // 2 + 1):
            f.append(tuple(sorted(((i + j) % m, (i - j) % m))))
        factors.append(f)
    return factors


def is_1_factorization(n, factors):
    """(ok, detail) -- exact: n-1 factors, each a perfect matching on n points,
    together each of the C(n,2) edges exactly once."""
    if len(factors) != n - 1:
        return False, 'got %d factors, want %d' % (len(factors), n - 1)
    seen = {}
    for f in factors:
        pts = [x for e in f for x in e]
        if len(f) != n // 2 or sorted(pts) != list(range(n)):
            return False, 'a factor is not a perfect matching'
        for e in f:
            seen[e] = seen.get(e, 0) + 1
    want = n * (n - 1) // 2
    if len(seen) != want or set(seen.values()) != {1}:
        return False, 'edges covered %d distinct, multiplicities %s' % (len(seen), sorted(set(seen.values())))
    return True, '%d factors, %d edges each exactly once' % (n - 1, want)


def hanani_double(n, blocks_n):
    """SQS(n) -> SQS(2n) on {0,...,2n-1} via (x,lev) -> x + n*lev.
    Returns (blocks, (c_a, c_b, c_c)) with the three family sizes."""
    factors = round_robin_1_factorization(n)
    L = lambda x, lev: x + n * lev
    fam_a, fam_b, fam_c = [], [], []
    for lev in (0, 1):
        for b in blocks_n:
            fam_a.append(tuple(sorted(L(x, lev) for x in b)))
    for x in range(n):
        for y in range(x + 1, n):
            fam_b.append(tuple(sorted((L(x, 0), L(y, 0), L(x, 1), L(y, 1)))))
    for f in factors:
        # ORDERED pairs of distinct edges of a common 1-factor: {x,y} at level 0
        # and {z,w} at level 1 is a DIFFERENT block from the reverse.
        for e1, e2 in itertools.permutations(f, 2):
            fam_c.append(tuple(sorted((L(e1[0], 0), L(e1[1], 0),
                                       L(e2[0], 1), L(e2[1], 1)))))
    out = fam_a + fam_b + fam_c
    return sorted(out), (len(fam_a), len(fam_b), len(fam_c))


def sqs_report(v, blocks):
    """(ok, detail) -- blocks form an SQS(v): the right number of distinct
    4-subsets of {0,...,v-1} covering every 3-subset exactly once."""
    b = v * (v - 1) * (v - 2) // 24
    if len(blocks) != b or len(set(blocks)) != b:
        return False, '%d blocks (%d distinct), want %d' % (len(blocks), len(set(blocks)), b)
    seen = set()
    for blk in blocks:
        if len(set(blk)) != 4 or not all(0 <= x < v for x in blk):
            return False, 'block %r is not a 4-subset of the point set' % (blk,)
        for t in itertools.combinations(blk, 3):
            if t in seen:
                return False, 'triple %r covered twice' % (t,)
            seen.add(t)
    t3 = v * (v - 1) * (v - 2) // 6
    if len(seen) != t3:
        return False, 'covered %d triples, want %d' % (len(seen), t3)
    return True, '%d blocks, all %d triples exactly once' % (b, t3)


def matchings(block):
    """The three perfect matchings of a 4-set (a<b<c<d), indexed by the trit:
    0 -> {a,b|c,d}, 1 -> {a,c|b,d}, 2 -> {a,d|b,c}."""
    a, b, c, d = block
    return (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c)))


def nested_pair_profile(blocks, trits):
    """The multiset of nested pairs of the nesting `trits` of `blocks`."""
    cnt = {}
    for blk, t in zip(blocks, trits):
        for e in matchings(blk)[t]:
            cnt[e] = cnt.get(e, 0) + 1
    return cnt


def is_completely_uniform(v, blocks, trits):
    """(ok, detail): all C(v,2) pairs occur as nested pairs and every
    multiplicity equals (v-2)/6."""
    mu = (v - 2) // 6
    cnt = nested_pair_profile(blocks, trits)
    allp = v * (v - 1) // 2
    if len(cnt) != allp:
        return False, '%d of %d pairs occur; %d missing' % (len(cnt), allp, allp - len(cnt))
    bad = sorted(e for e, m in cnt.items() if m != mu)
    if bad:
        return False, '%d pairs have multiplicity != %d (e.g. %r -> %d)' % (
            len(bad), mu, bad[0], cnt[bad[0]])
    return True, '%d pairs, every multiplicity == %d' % (allp, mu)


def exhaustive_nesting(v, blocks, node_cap):
    """A deterministic depth-first search for a completely uniform nesting,
    used only for the small forced-positive controls.  No randomness, no
    solver; returns (found, nodes, trits)."""
    mu = (v - 2) // 6
    remaining = {}
    for x in range(v):
        for y in range(x + 1, v):
            remaining[(x, y)] = mu
    choice = [0] * len(blocks)
    nodes = [0]

    def rec(i):
        nodes[0] += 1
        if nodes[0] > node_cap:
            raise RuntimeError('node cap')
        if i == len(blocks):
            return True
        for k, (p, q) in enumerate(matchings(blocks[i])):
            if remaining[p] > 0 and remaining[q] > 0:
                remaining[p] -= 1
                remaining[q] -= 1
                choice[i] = k
                if rec(i + 1):
                    return True
                remaining[p] += 1
                remaining[q] += 1
        return False

    try:
        found = rec(0)
    except RuntimeError:
        return False, nodes[0], None
    return found, nodes[0], list(choice)


# ===========================================================================
print('=' * 74)
print('A COMPLETELY UNIFORM NESTED STEINER QUADRUPLE SYSTEM OF ORDER 56')
print('checker for the paper in this folder; exact integer arithmetic only')
print('=' * 74)

# --- 1. the arithmetic of the cell -----------------------------------------
b56 = 56 * 55 * 54 // 24
pairs56 = 56 * 55 // 2
check('admissibility-of-56',
      56 % 6 == 2 and b56 == NBLOCKS and pairs56 == 1540
      and 2 * b56 == 1540 * MULT and (56 - 2) % 6 == 0 and (56 - 2) // 6 == MULT,
      '56 = 2 mod 6; b = 56*55*54/24 = %d; C(56,2) = %d; 2b/C(56,2) = %d = (56-2)/6'
      % (b56, pairs56, 2 * b56 // pairs56))

# --- 2. Lu's SQS(14), rebuilt from the thirteen printed base blocks --------
check('base-blocks-are-thirteen-4-sets',
      len(SQS14_BASE) == 13
      and all(len(set(b)) == 4 and all(0 <= x <= 13 for x in b) for b in SQS14_BASE)
      and len(set(tuple(sorted(b)) for b in SQS14_BASE)) == 13,
      '13 distinct 4-subsets of {0,...,13}')

sqs14 = sorted(set(e for b in SQS14_BASE for e in develop_z7(b)))
check('z7-development-gives-91-blocks',
      len(sqs14) == 91 and 13 * 7 == 91,
      '13 base blocks x 7 translates = 91 distinct blocks, no collision')

ok, detail = sqs_report(14, sqs14)
check('sqs14-is-a-steiner-quadruple-system', ok, detail)

# --- 3. the two doublings --------------------------------------------------
f14 = round_robin_1_factorization(14)
ok, detail = is_1_factorization(14, f14)
check('round-robin-1-factorization-of-k14', ok, detail)

sqs28, (a28, b28f, c28) = hanani_double(14, sqs14)
check('doubling-counts-14-to-28',
      (a28, b28f, c28) == (182, 91, 546) and a28 + b28f + c28 == 819
      and 819 == 28 * 27 * 26 // 24,
      '%d + %d + %d = 819 = 28*27*26/24' % (a28, b28f, c28))
ok, detail = sqs_report(28, sqs28)
check('sqs28-is-a-steiner-quadruple-system', ok, detail)

f28 = round_robin_1_factorization(28)
ok, detail = is_1_factorization(28, f28)
check('round-robin-1-factorization-of-k28', ok, detail)

sqs56, (a56, b56f, c56) = hanani_double(28, sqs28)
check('doubling-counts-28-to-56',
      (a56, b56f, c56) == (1638, 378, 4914) and a56 + b56f + c56 == NBLOCKS
      and NBLOCKS == 56 * 55 * 54 // 24,
      '%d + %d + %d = %d = 56*55*54/24' % (a56, b56f, c56, NBLOCKS))
ok, detail = sqs_report(56, sqs56)
check('sqs56-is-a-steiner-quadruple-system', ok, detail)

check('sqs56-blocks-are-in-canonical-order',
      sqs56 == sorted(sqs56) and all(list(b) == sorted(b) for b in sqs56),
      'increasing 4-tuples, sorted lexicographically: B_0 = %r, B_%d = %r'
      % (sqs56[0], NBLOCKS - 1, sqs56[-1]))

deg = {}
for blk in sqs56:
    for e in itertools.combinations(blk, 2):
        deg[e] = deg.get(e, 0) + 1
check('sqs56-pair-degree-is-27',
      len(deg) == 1540 and set(deg.values()) == {27} and 27 == (56 - 2) // 2,
      'every one of the 1540 pairs lies in exactly 27 = (56-2)/2 blocks, and 27/3 = %d'
      % (27 // 3))

# --- 4. the witness -------------------------------------------------------
b64 = ''.join(WITNESS_LINES)
check('witness-base64-length-and-digest',
      len(b64) == 1832 and len(WITNESS_LINES) == 29
      and hashlib.sha256(b64.encode('ascii')).hexdigest() == WITNESS_B64_SHA256,
      '29 printed lines concatenate to 1832 ASCII characters, sha256 = %s'
      % WITNESS_B64_SHA256)

num = int.from_bytes(base64.b64decode(b64), 'big')
trits = []
for _ in range(NBLOCKS):
    trits.append(num % 3)
    num //= 3
check('witness-decodes-to-exactly-6930-trits',
      len(trits) == NBLOCKS and num == 0 and set(trits) <= {0, 1, 2},
      'little-endian base 3, residue above t_%d is 0, values in {0,1,2}' % (NBLOCKS - 1))

check('witness-trit-vector-digest',
      hashlib.sha256(bytes(trits)).hexdigest() == WITNESS_TRIT_SHA256,
      'sha256 of the 6930 raw trit bytes = %s' % WITNESS_TRIT_SHA256)

check('witness-decimal-prefix-agrees-with-base64',
      len(TRIT_PREFIX_DECIMAL) == 200
      and ''.join(map(str, trits[:200])) == TRIT_PREFIX_DECIMAL,
      'the two independent encodings printed in Section 4 agree digit for digit '
      'over all 200 trits')

# --- 5. the nesting -------------------------------------------------------
nested = [matchings(blk)[t] for blk, t in zip(sqs56, trits)]

wf = True
for blk, (p, q) in zip(sqs56, nested):
    if set(p) & set(q) or tuple(sorted(p + q)) != blk or len(p) != 2 or len(q) != 2:
        wf = False
        break
check('nested-blocks-are-well-formed', wf,
      'each of the %d nested blocks is an unordered pair of DISJOINT pairs whose '
      'union is its own 4-set' % NBLOCKS)

check('nested-blocks-underlie-the-sqs56',
      sorted(tuple(sorted(p + q)) for p, q in nested) == sorted(sqs56),
      'the multiset of underlying 4-sets is exactly the SQS(56) rebuilt above')

cnt = nested_pair_profile(sqs56, trits)
allpairs = set((x, y) for x in range(56) for y in range(x + 1, 56))
check('all-1540-pairs-occur-as-nested-pairs',
      len(cnt) == 1540 and not (allpairs - set(cnt)),
      'distinct nested pairs = %d of C(56,2) = 1540, missing = 0' % len(cnt))

check('nesting-is-uniform-of-multiplicity-9',
      set(cnt.values()) == {MULT},
      'multiplicity set = %s, target (56-2)/6 = %d' % (sorted(set(cnt.values())), MULT))

check('multiplicity-sum-is-consistent',
      sum(cnt.values()) == 2 * NBLOCKS == 1540 * MULT == 13860,
      'sum = %d = 2*%d = 1540*%d' % (sum(cnt.values()), NBLOCKS, MULT))

check('first-twenty-nested-blocks-as-printed',
      tuple(nested[:20]) == FIRST_20_PRINTED,
      'B_0..B_19 decode to the twenty nested blocks displayed in Section 4')

# --- 6. anti-controls: the checker can say NO ------------------------------
bad = list(trits)
bad[0] = (bad[0] + 1) % 3
ok, detail = is_completely_uniform(56, sqs56, bad)
check('anti-control-one-flipped-trit-is-rejected', not ok,
      'trit t_0 changed from %d to %d: correctly REJECTED -- %s' % (trits[0], bad[0], detail))

ok, detail = is_completely_uniform(56, sqs56, [0] * NBLOCKS)
check('anti-control-trivial-all-zero-nesting-is-rejected', not ok,
      'nesting every block as {a,b|c,d}: correctly REJECTED -- %s' % detail)

# --- 7. forced positives: two settled orders, re-derived independently -----
sqs4 = [(0, 1, 2, 3)]
sqs8, cnt8 = hanani_double(4, sqs4)
ok8, d8 = sqs_report(8, sqs8)
found8, nodes8, tr8 = exhaustive_nesting(8, sqs8, 10 ** 6)
ok8b, d8b = (is_completely_uniform(8, sqs8, tr8) if found8 else (False, 'no nesting found'))
check('control-v8-forced-positive',
      ok8 and found8 and ok8b and (8 - 2) // 6 == 1,
      'SQS(8) rebuilt by the same doubling (%s); exhaustive search settled it in %d nodes; '
      '%s -- independently reproduces the v=8 entry of Lu\'s Table 4' % (d8, nodes8, d8b))

found14, nodes14, tr14 = exhaustive_nesting(14, sqs14, 5 * 10 ** 6)
ok14, d14 = (is_completely_uniform(14, sqs14, tr14) if found14 else (False, 'no nesting found'))
check('control-v14-forced-positive',
      found14 and ok14 and (14 - 2) // 6 == 2 and tr14 != None,
      'the SQS(14) above, nested by our own exhaustive search and NOT by Lu\'s printed '
      'pairings: settled in %d nodes; %s -- independently reproduces the v=14 entry of '
      'Lu\'s Table 4' % (nodes14, d14))

# --- 8. proved-silent controls: orders where no such nesting can exist -----
b16 = 16 * 15 * 14 // 24
check('control-v16-arithmetically-silent',
      16 % 6 == 4 and 2 * b16 == 280 and 16 * 15 // 2 == 120 and 280 % 120 != 0,
      '16 = 4 mod 6, b(16) = %d, 2b = 280, C(16,2) = 120, 280/120 is not an integer, '
      'so no completely uniform nested SQS(16) exists' % b16)

b28v = 28 * 27 * 26 // 24
check('control-v28-arithmetically-silent',
      28 % 6 == 4 and 2 * b28v == 1638 and 28 * 27 // 2 == 378 and 1638 % 378 != 0,
      'the SQS(28) built above has b = %d, 2b = 1638, C(28,2) = 378, 1638/378 is not an '
      'integer, so it admits no completely uniform nesting' % b28v)

# --- 9. the surrounding census row, recomputed ----------------------------
check('sister-cells-of-the-v56-row',
      924 * 15 == 13860 and 1260 * 11 == 13860 and 2 * NBLOCKS == 13860
      and 924 < 1540 and 1260 < 1540,
      '924*15 = 1260*11 = 13860 = 2b(56): the two cells of the v=56 row that remain open '
      'are arithmetically consistent and are NOT settled here')

check('order-62-already-in-print',
      62 * 61 // 2 == 1891 and (62 - 2) // 6 == 10 and 1891 * 10 == 2 * (62 * 61 * 60 // 24),
      'C(62,2) = 1891 and (62-2)/6 = 10, so v = 62 -- which is settled in the published '
      'literature -- is LARGER than 56: 56 is a hole, not the frontier')

adm = [v for v in range(8, 51) if v % 6 == 2]
check('admissible-orders-8-to-50-and-the-gap-to-56',
      adm == [8, 14, 20, 26, 32, 38, 44, 50]
      and not any(v % 6 == 2 for v in range(51, 56))
      and 56 % 6 == 2,
      'the admissible orders in [8,50] are exactly %s; none of 51..55 is 2 mod 6, so 56 is '
      'the next admissible order after 50' % adm)

# ---------------------------------------------------------------------------
# What this program does NOT cover
# ---------------------------------------------------------------------------
print('')
print('NOT RE-RUN: nonexistence of anything. The two remaining open cells of the v=56')
print('NOT RE-RUN: row -- 924 nested pairs at multiplicity 15, and 1260 at multiplicity 11')
print('NOT RE-RUN: -- are only checked for arithmetic consistency above; no search for or')
print('NOT RE-RUN: against them was run, and they stay open.')
print('NOT RE-RUN: whether this SQS(56) is rotational, i.e. admits an automorphism of')
print('NOT RE-RUN: order 55 fixing a point. No automorphism computation was performed on')
print('NOT RE-RUN: the 6930 blocks; the paper claims only that rotationality is not USED.')
print('NOT RE-RUN: uniqueness, isomorphism, or any enumeration of the nestings of this or')
print('NOT RE-RUN: any other SQS(56). One witness is exhibited; 0% of the space is')
print('NOT RE-RUN: exhausted, and no claim of minimality or canonicity is made.')
print('NOT RE-RUN: Table 4 of Lu at v = 20, 26, 32, 38, 44, 50. Only v = 8 and v = 14 are')
print('NOT RE-RUN: re-derived here as forced positives; the other six rows are taken on')
print('NOT RE-RUN: the published authority and bear on novelty, never on correctness.')
print('NOT RE-RUN: the search that produced the witness. This program re-checks the')
print('NOT RE-RUN: witness only; the randomised descent that found it is not reproduced,')
print('NOT RE-RUN: and none of the above depends on it.')
print('')

n = len(_passed)
if _failed:
    print('%d CHECK(S) DID NOT PASS: %s' % (len(_failed), ', '.join(_failed)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % n)
sys.exit(0)
