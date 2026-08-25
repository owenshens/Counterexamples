from itertools import combinations
from itertools import combinations_with_replacement as cwr

def parts(n, mx=None):            # partitions, parts decreasing
    mx = n if mx is None else mx
    if n == 0:
        yield ()
        return
    for k in range(min(n, mx), 0, -1):
        for t in parts(n - k, k):
            yield (k,) + t

def evens(m):                     # even partitions of m
    return [] if m % 2 else [tuple(2*p for p in t)
                             for t in parts(m//2)]

def key(ls):
    return tuple(sorted(ls, reverse=True))

def sd(a, b):                     # a superdominates b
    A, B, sa, sb = sorted(a), sorted(b), 0, 0
    for i in range(min(len(A), len(B))):
        sa, sb = sa + A[i], sb + B[i]
        if sa > sb:
            return False
    return True

def alphas(r, s, m):              # ordered: r odd, then s even
    if r:
        for a in range(1, m + 1, 2):
            for t in alphas(r - 1, s, m - a):
                yield (a,) + t
    elif s:
        for a in range(2, m + 1, 2):
            for t in alphas(0, s - 1, m - a):
                yield (a,) + t
    else:
        yield ()

Y = {(18,):0,(16,2):9,(14,4):6,(14,2,2):20,(12,6):3,
 (12,4,2):14,(12,2,2,2):32,(10,8):0,(10,6,2):12,(10,4,4):12,
 (10,4,2,2):24,(10,2,2,2,2):44,(8,8,2):12,(8,6,4):10,
 (8,6,2,2):24,(8,4,4,2):21,(8,4,2,2,2):36,(8,2,2,2,2,2):56,
 (6,6,6):10,(6,6,4,2):20,(6,6,2,2,2):36,(6,4,4,4):20,
 (6,4,4,2,2):30,(6,4,2,2,2,2):48,(6,2,2,2,2,2,2):68,
 (4,4,4,4,2):30,(4,4,4,2,2,2):40,(4,4,2,2,2,2,2):60,
 (4,2,2,2,2,2,2,2):80,(2,2,2,2,2,2,2,2,2):100}

d, L = 9, evens(18)
assert len(L) == 30 and sorted(L) == sorted(Y)
one = [(a, b) for a in L for b in L if a != b and sd(a, b)]
print(len(one), min(Y[a] - Y[b] for a, b in one))    # 417 0

blk = {}
for r in range(d + 1):
    for s in range(d//2 + 1):
        V = [(al, la) for al in alphas(r, s, d)
             for la in evens(d - sum(al))] if r + s else []
        if V:
            blk[(r, s)] = V
print(len(blk), sum(len(V) for V in blk.values()))   # 15 131

A = lambda u: key([2*a for a in u[0]] + 2*list(u[1]))
def Bs(u):
    t = sorted([2*a for a in u[0]], reverse=True)
    return key([t[0] + t[1]] + t[2:] + 2*list(u[1]))
C = lambda u, v: key([a + b for a, b in zip(u[0], v[0])]
                     + list(u[1]) + list(v[1]))

dg = [(A(u), Bs(u)) for (r, s), V in blk.items()
      if r + s >= 2 for u in V]
of = [(A(u), A(v), C(u, v)) for V in blk.values()
      for u, v in combinations(V, 2)]
print(len(dg), len(of), len(dg) + len(of))       # 119 937 1056
print(min(Y[a] - Y[b] for a, b in dg),
      min(Y[a] + Y[b] - 2*Y[c] for a, b, c in of))   # 0 0
print(all(sd(a, b) for a, b in dg),
      all(sd(a + b, c + c) for a, b, c in of))       # True True

two = [(a, b, m) for a, b in cwr(L, 2) for m in L
       if sd(a + b, m + m)]
bad = [t for t in two if Y[t[0]] + Y[t[1]] - 2*Y[t[2]] < 0]
print(len(two), bad)
# 7025 [((10, 8), (4, 4, 4, 2, 2, 2), (8, 4, 4, 2))]
print(Y[(10,8)] + Y[(4,4,4,2,2,2)] - 2*Y[(8,4,4,2)])      # -2
