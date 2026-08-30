#!/usr/bin/env python3
"""verify.py -- checks the computational claims of

    "An Explicit Degree-m Minimal Polynomial: a Proof of Conjecture 11 of Lasjaunias"

Python 3.9+, STANDARD LIBRARY ONLY (sys, time, math for comb only via a local integer routine).
All arithmetic is exact in F_p, F_p[T] and F_p[T][X] on integer tuples.  There is no floating-point
decision anywhere, no solver, no search, no external data file.

WHAT IT DOES
    PART A  reads the objects PRINTED IN THE PAPER for the witness cell (p,k,m) = (17,4,4) -- the word
            W, omega_4, eps_A, the four continuants, H(V), the conjectured P and the full quotient Q --
            re-derives each of them from the definitions transcribed from the source, and checks the
            factorisation H = P*Q by multiplying the printed quotient out.
    PART B  runs the same derivations over EVERY admissible triple (p, k, m) = (p, (p-1)/m, m) with p an
            odd prime <= BOUND (default 300), and checks each lemma of the paper's proof on every one of
            them, aggregated one PASS line per lemma.
    PART C  negative controls: objects that MUST NOT divide H(V) (a perturbed P; the source's T-less
            Proposition 5.1 polynomial; the conjectured P against eps_B instead of eps_A).
    PART D  published pins -- the source author's own printed integers and polynomials, reproduced.
    PART E  RECONSTRUCTS the root alpha of Proposition 9 as a Laurent series in F_p((1/T)) (Newton in
            u = 1/T, started at alpha_0 = w_{k+1}T), on every admissible triple with p <= 60, and
            certifies to order u^20: H(V)(alpha) = 0; the Hensel hypothesis v(H) > 2v(H'), which turns
            that into an exact root; that the first partial quotient is w_{k+1}T, so alpha has a pole
            at T = infinity; the diagonalised Frobenius relation alpha^p = (A alpha + B)/(B alpha + A);
            beta = (alpha+1)/(alpha-1) = 1 mod 1/T and beta^p = a^k beta; beta^m = a EXACTLY, which is
            the selection zeta = 1; P(alpha) = 0; and that NONE of the k-1 sibling factors, the ones
            beta^m = zeta a with zeta^k = 1, zeta != 1, vanishes at alpha -- so the census now
            distinguishes the conjectured P from the other degree-m factors of H(V), which mere exact
            division cannot do.  Three controls -- a perturbed alpha, a wrong initial partial quotient,
            the eps_B form of H -- and the paper's PRINTED P evaluated at the reconstructed alpha.

WHAT IT DOES NOT DO -- printed as `NOT RE-RUN:` lines at the end, and repeated in REVIEW_NOTE.md's
`## Scope`.  The program still does not PROVE the theorem: it re-runs it on a finite census, it quotes
Hensel's lemma for the passage from order u^20 to an exact root, and it quotes Capelli's criterion for
the irreducibility of P (only the hypothesis v_{T-1}(a) = +-1 is re-derived).

    python3 verify.py            # BOUND = 300
    python3 verify.py 100        # smaller census
"""
import sys
import time

BOUND = int(sys.argv[1]) if len(sys.argv) > 1 else 300

# ===========================================================================
# 0.  THE OBJECTS AS PRINTED IN THE PAPER.  Transcribed character-for-character
#     from paper.tex, section "The witness cell".  Nothing below is recomputed
#     into these strings; they are the program's INPUT.
# ===========================================================================
PAPER_p, PAPER_k, PAPER_m = 17, 4, 4
PAPER_omega = 15
PAPER_W = [7, 8, 6, 4, 13, 11, 9, 10]          # w_1, ..., w_8 in F_17
PAPER_epsA = 1
PAPER_epsB = 16
PAPER_K13 = "13T^3+13T"
PAPER_K14 = "T^4+6T^2+1"
PAPER_K58 = "T^4+6T^2+1"
PAPER_K68 = "4T^3+4T"
PAPER_H = "(4T^3+4T)X^18 + (16T^4+11T^2+16)X^17 + (T^4+6T^2+1)X + (13T^3+13T)"
PAPER_P = "X^4 + 4TX^3 + 6X^2 + 4TX + 1"
PAPER_Q = ("(4T^3+4T)X^14 + (12T^2+16)X^13 + (13T^3+14T)X^12 + (9T^2+6)X^11 + (4T^3+11T)X^10"
           " + (4T^2+16)X^9 + (13T^3+2T)X^8 + (4T^3+15T)X^6 + (13T^2+1)X^5 + (13T^3+6T)X^4"
           " + (8T^2+11)X^3 + (4T^3+3T)X^2 + (5T^2+1)X + (13T^3+13T)")
# the two hand specialisations printed in the paper
PAPER_P_AT_T0_ROOTS = [3, -3, 6, -6]           # P(X) mod (T) = (X-3)(X+3)(X-6)(X+6) over F_17
PAPER_H_AT_T1_SCALAR = 8                       # H(V)|_{T=1} = 8 (X-1)(X+1)^17 over F_17

# ===========================================================================
# 1.  F_p[T]:  tuple of ints, low degree first, no trailing zeros
# ===========================================================================
def tnorm(a):
    i = len(a)
    while i > 0 and a[i - 1] == 0:
        i -= 1
    return tuple(a[:i])


def tadd(a, b, p):
    n = max(len(a), len(b))
    return tnorm([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)])


def tsub(a, b, p):
    n = max(len(a), len(b))
    return tnorm([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(n)])


def tmul(a, b, p):
    if not a or not b:
        return ()
    r = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    r[i + j] = (r[i + j] + ai * bj) % p
    return tnorm(r)


def tscal(c, a, p):
    c %= p
    return () if c == 0 else tnorm([(c * x) % p for x in a])


def teval(a, t, p):
    v = 0
    for c in reversed(a):
        v = (v * t + c) % p
    return v


def tdivmod(a, b, p):
    """exact division with remainder in F_p[T]; b != 0."""
    assert b
    q = [0] * max(0, len(a) - len(b) + 1)
    r = list(a)
    li = pow(b[-1], p - 2, p)
    for i in range(len(r) - 1, len(b) - 2, -1):
        c = r[i] % p
        if c == 0:
            continue
        c = c * li % p
        q[i - len(b) + 1] = c
        for j in range(len(b)):
            r[i - len(b) + 1 + j] = (r[i - len(b) + 1 + j] - c * b[j]) % p
    return tnorm(q), tnorm(r)


def tstr(a):
    if not a:
        return "0"
    out = []
    for i in range(len(a) - 1, -1, -1):
        c = a[i]
        if c == 0:
            continue
        if i == 0:
            out.append(str(c))
        elif i == 1:
            out.append(("" if c == 1 else str(c)) + "T")
        else:
            out.append(("" if c == 1 else str(c)) + "T^" + str(i))
    return "+".join(out)


# ===========================================================================
# 2.  F_p[T][X]:  list of F_p[T] elements, index = degree in X
# ===========================================================================
def xnorm(A):
    i = len(A)
    while i > 0 and A[i - 1] == ():
        i -= 1
    return A[:i]


def xadd(A, B, p):
    n = max(len(A), len(B))
    return xnorm([tadd(A[i] if i < len(A) else (), B[i] if i < len(B) else (), p) for i in range(n)])


def xsub(A, B, p):
    n = max(len(A), len(B))
    return xnorm([tsub(A[i] if i < len(A) else (), B[i] if i < len(B) else (), p) for i in range(n)])


def xmul(A, B, p):
    if not A or not B:
        return []
    R = [()] * (len(A) + len(B) - 1)
    for i, ai in enumerate(A):
        if ai:
            for j, bj in enumerate(B):
                if bj:
                    R[i + j] = tadd(R[i + j], tmul(ai, bj, p), p)
    return xnorm(R)


def xdivmod_monic(A, P, p):
    """A = P*Q + R in F_p[T][X], P MONIC in X.  Exact: no inverses are taken."""
    assert P and P[-1] == (1,)
    Q = [()] * max(0, len(A) - len(P) + 1)
    R = list(A)
    d = len(P) - 1
    for i in range(len(R) - 1, d - 1, -1):
        c = R[i]
        if c == ():
            continue
        Q[i - d] = c
        for j in range(d + 1):
            if P[j]:
                R[i - d + j] = tsub(R[i - d + j], tmul(c, P[j], p), p)
    return xnorm(Q), xnorm(R)


def xstr(A):
    if not A:
        return "0"
    out = []
    for i in range(len(A) - 1, -1, -1):
        if A[i] == ():
            continue
        cs = tstr(A[i])
        out.append("(%s)X^%d" % (cs, i) if i else "(%s)" % cs)
    return " + ".join(out)


# ===========================================================================
# 3.  THE PARSER for the strings printed in the paper.  Strict: it accepts
#     exactly the normal form the paper uses and raises on anything else, so a
#     mis-transcription is an error and not a silently different object.
# ===========================================================================
def parse_T(s, p):
    """'13T^3+13T' -> F_p[T] element."""
    s = s.replace(" ", "").replace("*", "")
    if not s:
        raise ValueError("empty T-polynomial")
    out = ()
    for term in s.split("+"):
        if not term:
            raise ValueError("empty term in %r" % s)
        if "T" in term:
            co, _, ex = term.partition("T")
            c = 1 if co == "" else int(co)
            if ex == "":
                e = 1
            elif ex.startswith("^"):
                e = int(ex[1:])
            else:
                raise ValueError("bad exponent %r" % ex)
        else:
            c, e = int(term), 0
        mono = [0] * e + [c % p]
        out = tadd(out, tnorm(mono), p)
    return out


def parse_X(s, p):
    """'(4T^3+4T)X^18 + ... + (13T^3+13T)' -> F_p[T][X] element."""
    s = s.replace(" ", "").replace("*", "")
    terms, depth, cur = [], 0, ""
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "+" and depth == 0:
            terms.append(cur)
            cur = ""
        else:
            cur += ch
    terms.append(cur)
    if depth != 0:
        raise ValueError("unbalanced parentheses in %r" % s)
    A = []
    for term in terms:
        if not term:
            raise ValueError("empty X-term in %r" % s)
        if term.startswith("("):
            close = term.index(")")
            coef = parse_T(term[1:close], p)
            rest = term[close + 1:]
        else:
            i = term.find("X")
            if i < 0:
                coef, rest = parse_T(term, p), ""
            else:
                head = term[:i]
                coef = (1,) if head == "" else parse_T(head, p)
                rest = term[i:]
        if rest == "":
            e = 0
        elif rest == "X":
            e = 1
        elif rest.startswith("X^"):
            e = int(rest[2:])
        else:
            raise ValueError("bad X part %r" % rest)
        while len(A) <= e:
            A.append(())
        A[e] = tadd(A[e], coef, p)
    return xnorm(A)


# ===========================================================================
# 4.  THE SOURCE'S DEFINITIONS, transcribed from Lasjaunias, J. Integer Seq. 21
#     (2018) Art. 18.8.3 (= arXiv:1803.01739), lines 43-47, 60-73, 182-192, 202.
# ===========================================================================
def finv(a, p):
    return pow(a % p, p - 2, p)


def omega_and_word(p, k):
    """omega_k = (-1)^(k-1) prod_{1<=i<=k}(1 - 1/(2i));  w_1 = (2k-1)(2k omega_k)^{-1};
    w_{i+1} w_i = (2k-2i-1)(2k-2i+1)(i(2k-i))^{-1}.  Returns (omega_k, [w_1..w_{2k}])."""
    om = 1 if (k - 1) % 2 == 0 else p - 1
    for i in range(1, k + 1):
        om = om * ((1 - finv(2 * i, p)) % p) % p
    w = [None] * (2 * k + 1)
    w[1] = (2 * k - 1) * finv(2 * k * om % p, p) % p
    for i in range(1, 2 * k):
        rhs = (2 * k - 2 * i - 1) * (2 * k - 2 * i + 1) % p * finv(i * (2 * k - i) % p, p) % p
        if w[i] == 0:
            raise ValueError("w_%d vanishes at p=%d k=%d" % (i, p, k))
        w[i + 1] = rhs * finv(w[i], p) % p
    return om, w[1:]


def K(w, p, n, r):
    """K_{n,r} = <w_n T, ..., w_r T> in F_p[T], with K_{n,n-1} = 1 and K_{n,n-2} = 0."""
    if r <= n - 2:
        return ()
    A, B = (), (1,)
    for i in range(n, r + 1):
        A, B = B, tadd(tmul((0, w[i - 1]), B, p), A, p)
    return B


def const_K(w, p, lo, hi):
    """<w_lo, ..., w_hi> in F_p (bare entries, no T).  = K_{lo,hi}(1)."""
    A, B = 0, 1
    for i in range(lo, hi + 1):
        A, B = B, (w[i - 1] * B + A) % p
    return B


def cf_value(seq, p):
    """[a_1,...,a_n] = a_1 + 1/(a_2 + 1/(...)) in F_p; None if undefined."""
    v = seq[-1] % p
    for a in reversed(seq[:-1]):
        if v == 0:
            return None
        v = (a + finv(v, p)) % p
    return v


def eps_AB(p, k, j):
    """eps_A = (-1)^(k+j+1) [w_{j+1},...,w_{2k},omega_k];  eps_B drops the omega_k tail."""
    om, w = omega_and_word(p, k)
    a = cf_value(list(w[j:2 * k]) + [om], p)
    b = cf_value(list(w[j:2 * k]), p)
    s = 1 if (k + j + 1) % 2 == 0 else -1
    return (None if a is None else s * a % p), (None if b is None else s * b % p)


def build_H(p, k, j, eps):
    """H([p,k,j,eps]) = K_{j+2,2k} X^{p+1} - K_{j+1,2k} X^p + eps(K_{1,j} X + K_{1,j-1})."""
    om, w = omega_and_word(p, k)
    H = [()] * (p + 2)
    H[p + 1] = K(w, p, j + 2, 2 * k)
    H[p] = tsub((), K(w, p, j + 1, 2 * k), p)
    H[1] = tscal(eps, K(w, p, 1, j), p)
    H[0] = tscal(eps, K(w, p, 1, j - 1), p)
    return xnorm(H)


def ibinom(n, i):
    """exact integer binomial coefficient."""
    if i < 0 or i > n:
        return 0
    r = 1
    for t in range(i):
        r = r * (n - t) // (t + 1)
    return r


def P_conj(p, k, m):
    """P(X) = sum_{0<=i<=m} (-1)^(ki) binom(m,i) T^((1-(-1)^i)/2) X^(m-i)  -- the conjectured factor."""
    A = [()] * (m + 1)
    for i in range(m + 1):
        c = ibinom(m, i) % p
        if (k * i) % 2:
            c = (-c) % p
        if c == 0:
            continue
        A[m - i] = tadd(A[m - i], ((0, c) if i % 2 else (c,)), p)
    return xnorm(A)


def P_split(p, k, m):
    """((1+c)(X+1)^m + (1-c)(X-1)^m)/2 with c = (-1)^k T -- the parity-split form."""
    c = (0, (-1) ** k % p)
    Xp1, Xm1 = [(1,), (1,)], [tsub((), (1,), p), (1,)]
    A, B = [(1,)], [(1,)]
    for _ in range(m):
        A = xmul(A, Xp1, p)
        B = xmul(B, Xm1, p)
    lhs = xadd(xmul([tadd((1,), c, p)], A, p), xmul([tsub((1,), c, p)], B, p), p)
    return lhs           # equals 2P


def R_of(p, k):
    """R_k = T^p mod (T^2-1)^k in F_p[T], and the quotient Q_k with T^p = (T^2-1)^k Q_k + R_k."""
    Pk = (1,)
    for _ in range(k):
        Pk = tmul(Pk, (p - 1, 0, 1), p)
    Tp = [0] * (p + 1)
    Tp[p] = 1
    q, r = tdivmod(tnorm(Tp), Pk, p)
    return Pk, q, r


def primes_upto(n):
    s = [True] * (n + 1)
    s[0] = s[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            for j in range(i * i, n + 1, i):
                s[j] = False
    return [i for i in range(3, n + 1) if s[i]]


def valuation_at(a, root, p):
    """multiplicity of (T - root) in a nonzero a in F_p[T]."""
    assert a
    v = 0
    cur = a
    lin = tnorm([(-root) % p, 1])
    while True:
        q, r = tdivmod(cur, lin, p)
        if r != ():
            return v
        v += 1
        cur = q


# ===========================================================================
# 4b.  F_p((1/T)):  truncated Laurent series in u = 1/T, with TRACKED precision
#
#      A series is the triple (lo, co, prec):  sum_{i} co[i] u^(lo+i), every
#      coefficient with exponent < prec is EXACT, nothing is known at exponents
#      >= prec.  co carries no leading and no trailing zero, so for a nonzero
#      series lo is exactly its valuation v (v(u) = 1, v(T) = -1); the zero
#      series is (prec, (), prec), i.e. "0 mod u^prec".  Every operation
#      propagates prec honestly, so a claim "X = 0 mod u^N" is only ever made
#      when the computed object actually knows its coefficients that far.
#      SCAP is the storage cap: no coefficient at exponent >= SCAP is retained.
# ===========================================================================
SCAP = 64


def snorm(lo, co, prec, p):
    co = [c % p for c in co]
    if lo + len(co) > prec:
        co = co[:max(0, prec - lo)]
    i = 0
    while i < len(co) and co[i] == 0:
        i += 1
    lo += i
    co = co[i:]
    while co and co[-1] == 0:
        co.pop()
    if not co:
        return (prec, (), prec)
    return (lo, tuple(co), prec)


def sconst(c, p):
    return snorm(0, [c], SCAP, p)


def sfromT(a, p):
    """F_p[T] (low degree first) -> Laurent series in u = 1/T, exact to SCAP."""
    if not a:
        return (SCAP, (), SCAP)
    d = len(a) - 1
    return snorm(-d, [a[d - i] for i in range(d + 1)], SCAP, p)


def sadd(A, B, p):
    prec = min(A[2], B[2])
    lo = min(A[0], B[0])
    n = prec - lo
    if n <= 0:
        return (prec, (), prec)
    co = [0] * n
    for (l, c, _) in (A, B):
        for i, v in enumerate(c):
            e = l + i - lo
            if 0 <= e < n:
                co[e] += v
    return snorm(lo, co, prec, p)


def sneg(A, p):
    return snorm(A[0], [(-c) % p for c in A[1]], A[2], p)


def ssub(A, B, p):
    return sadd(A, sneg(B, p), p)


def sscal(c, A, p):
    return snorm(A[0], [c * x for x in A[1]], A[2], p)


def smul(A, B, p):
    prec = min(A[2] + B[0], B[2] + A[0])
    lo = A[0] + B[0]
    n = prec - lo
    if n <= 0 or not A[1] or not B[1]:
        return (prec, (), prec)
    co = [0] * n
    for i, x in enumerate(A[1]):
        if not x:
            continue
        for j, y in enumerate(B[1]):
            e = i + j
            if e >= n:
                break
            co[e] += x * y
    return snorm(lo, co, prec, p)


def sinv(A, p):
    """1/A.  A must be nonzero; relative precision is preserved."""
    lo, c, prec = A
    assert c, "cannot invert the zero series"
    r = prec - lo
    cc = list(c) + [0] * (r - len(c))
    i0 = finv(cc[0], p)
    h = [i0]
    for n in range(1, r):
        s = 0
        for i in range(1, n + 1):
            if cc[i]:
                s += cc[i] * h[n - i]
        h.append((-i0 * s) % p)
    return snorm(-lo, h, prec - 2 * lo, p)


def sdiv(A, B, p):
    return smul(A, sinv(B, p), p)


def sfrob(A, p):
    """A -> A^p.  In characteristic p this only spreads the exponents, because
    (sum c_i u^i)^p = sum c_i^p u^(pi) and c^p = c for c in F_p."""
    lo, c, prec = A
    nlo, nprec = lo * p, min(SCAP, prec * p)
    n = nprec - nlo
    if n <= 0 or not c:
        return (nprec, (), nprec)
    co = [0] * n
    for i, v in enumerate(c):
        e = i * p
        if e < n:
            co[e] = v
    return snorm(nlo, co, nprec, p)


def spow(A, n, p):
    r = sconst(1, p)
    B = A
    while n:
        if n & 1:
            r = smul(r, B, p)
        n >>= 1
        if n:
            B = smul(B, B, p)
    return r


def spolypart(A):
    """the polynomial part of A: the terms with exponent <= 0, as F_p[T]."""
    lo, co, prec = A
    d = -lo
    if d < 0:
        return ()
    out = [0] * (d + 1)
    for i, v in enumerate(co):
        e = lo + i
        if e <= 0:
            out[-e] = v
    return tnorm(out)


def sval(A):
    """the valuation of A, or its precision when A is 0 to that precision."""
    return A[0]


def svanishes(A, N):
    """certified: every coefficient of A below u^N is zero AND A knows that far."""
    return A[0] >= N and A[2] >= N


def sstr(A, n=4):
    lo, co, prec = A
    if not co:
        return "0 mod u^%d" % prec
    out = []
    for i, v in enumerate(co[:n]):
        e = lo + i
        if v == 0:
            continue
        if e < 0:
            out.append("%dT^%d" % (v, -e) if e < -1 else "%dT" % v)
        elif e == 0:
            out.append("%d" % v)
        else:
            out.append("%d/T^%d" % (v, e) if e > 1 else "%d/T" % v)
    return " + ".join(out) + " + O(1/T^%d)" % (lo + min(n, len(co)))


# ===========================================================================
# 5.  THE CHECK LEDGER
# ===========================================================================
CHECKS = 0
FAILED = []


def ck(name, cond, detail=""):
    global CHECKS
    CHECKS += 1
    if cond:
        print("PASS %s%s" % (name, (" " + detail) if detail else ""))
    else:
        FAILED.append(name)
        print("FAIL %s%s" % (name, (" " + detail) if detail else ""))
    return bool(cond)


t0 = time.time()
p, k, m = PAPER_p, PAPER_k, PAPER_m
print("=== PART A -- the witness cell (p,k,m) = (%d,%d,%d), objects READ FROM THE PAPER ===" % (p, k, m))

H_paper = parse_X(PAPER_H, p)
P_paper = parse_X(PAPER_P, p)
Q_paper = parse_X(PAPER_Q, p)
K13_paper = parse_T(PAPER_K13, p)
K14_paper = parse_T(PAPER_K14, p)
K58_paper = parse_T(PAPER_K58, p)
K68_paper = parse_T(PAPER_K68, p)

om, w = omega_and_word(p, k)
ck("A01-omega-from-definition", om == PAPER_omega,
   "omega_4 = %d; paper prints %d" % (om, PAPER_omega))
ck("A02-word-from-recurrence", list(w) == PAPER_W,
   "w_1..w_8 = %s; paper prints %s" % (list(w), PAPER_W))
ck("A03-antisymmetry-eq2", all((w[2 * k - i] + w[i - 1]) % p == 0 for i in range(1, 2 * k + 1)),
   "w_{9-i} = -w_i for i = 1..8")
ck("A04-K13", K(w, p, 1, 3) == K13_paper, "K_{1,3} = %s" % tstr(K13_paper))
ck("A05-K14", K(w, p, 1, 4) == K14_paper, "K_{1,4} = %s" % tstr(K14_paper))
ck("A06-K58", K(w, p, 5, 8) == K58_paper, "K_{5,8} = %s = (-1)^4 K_{1,4}" % tstr(K58_paper))
ck("A07-K68", K(w, p, 6, 8) == K68_paper, "K_{6,8} = %s = (-1)^3 K_{1,3}" % tstr(K68_paper))
EA, EB = eps_AB(p, k, k)
ck("A08-epsA-from-cf-definition", EA == PAPER_epsA and EB == PAPER_epsB,
   "eps_A = %s, eps_B = %s; paper prints %d and %d" % (EA, EB, PAPER_epsA, PAPER_epsB))
ck("A09-epsA-equals-minus-one-to-k", EA == (-1) ** k % p, "(-1)^4 = 1")
H = build_H(p, k, k, EA)
ck("A10-H-from-definition", H == H_paper, "deg_X H = %d = p+1" % (len(H) - 1))
Pc = P_conj(p, k, m)
ck("A11-P-from-closed-form", Pc == P_paper, "P = %s" % xstr(P_paper))
ck("A12-degP-equals-m", len(P_paper) - 1 == m, "deg_X P = %d" % (len(P_paper) - 1))
ck("A13-parity-split-identity", P_split(p, k, m) == xmul([(2,)], P_paper, p),
   "2P = (1+T)(X+1)^4 + (1-T)(X-1)^4")
ck("A14-printed-quotient-multiplies-out", xmul(P_paper, Q_paper, p) == H_paper,
   "P*Q = H exactly, deg_X Q = %d" % (len(Q_paper) - 1))
Qd, Rd = xdivmod_monic(H_paper, P_paper, p)
ck("A15-division-remainder-zero", Rd == [] and Qd == Q_paper, "remainder 0 and quotient = printed Q")
K14, K13 = K(w, p, 1, 4), K(w, p, 1, 3)
U, V = tadd(K14, K13, p), tsub(K14, K13, p)
Tm1, Tp1 = (1,), (1,)
for _ in range(k):
    Tm1 = tmul(Tm1, (p - 1, 1), p)
    Tp1 = tmul(Tp1, (1, 1), p)
ck("A16-U-is-T-minus-1-to-k", U == Tm1, "K_{1,4}+K_{1,3} = (T-1)^4")
ck("A17-V-is-T-plus-1-to-k", V == Tp1, "K_{1,4}-K_{1,3} = (T+1)^4")
ck("A18-difference-of-squares", tsub(tmul(K14, K14, p), tmul(K13, K13, p), p) ==
   tmul(Tm1, Tp1, p), "K_{1,4}^2 - K_{1,3}^2 = (T^2-1)^4")
# a = (c-1)/(c+1) with c = (-1)^k T = T ; lambda = U/V must equal a^k
ck("A19-multiplier-is-a-to-the-k", tmul(U, Tp1, p) == tmul(V, Tm1, p),
   "(K_{1,4}+K_{1,3})(T+1)^4 = (K_{1,4}-K_{1,3})(T-1)^4, i.e. lambda = a^4")
ck("A20-valuation-of-a-at-T-1", valuation_at((p - 1, 1), 1, p) == 1 and valuation_at((1, 1), 1, p) == 0,
   "v_{T-1}(T-1) = 1 and v_{T-1}(T+1) = 0, so v_{T-1}(a) = 1, coprime to m = 4")
# hand specialisations printed in the paper
H_T0 = xnorm([(teval(c, 0, p),) if teval(c, 0, p) else () for c in H_paper])
want_T0 = [()] * (p + 1)
want_T0[p] = (p - 1,)
want_T0[1] = (1,)
ck("A21-H-at-T-0", H_T0 == xnorm(want_T0), "H|_{T=0} = -(X^17 - X)")
P_T0 = xnorm([(teval(c, 0, p),) if teval(c, 0, p) else () for c in P_paper])
prod = [(1,)]
for r in PAPER_P_AT_T0_ROOTS:
    prod = xmul(prod, [((-r) % p,), (1,)], p)
ck("A22-P-at-T-0-splits", P_T0 == prod, "P|_{T=0} = (X-3)(X+3)(X-6)(X+6) over F_17")
H_T1 = xnorm([(teval(c, 1, p),) if teval(c, 1, p) else () for c in H_paper])
fac = [((-1) % p,), (1,)]
for _ in range(p):
    fac = xmul(fac, [(1,), (1,)], p)
ck("A23-H-at-T-1", H_T1 == xmul([(PAPER_H_AT_T1_SCALAR,)], fac, p),
   "H|_{T=1} = 8 (X-1)(X+1)^17 over F_17")
P_T1 = xnorm([(teval(c, 1, p),) if teval(c, 1, p) else () for c in P_paper])
q4 = [(1,)]
for _ in range(k):
    q4 = xmul(q4, [(1,), (1,)], p)
ck("A24-P-at-T-1", P_T1 == q4, "P|_{T=1} = (X+1)^4")
print()

# ===========================================================================
print("=== PART B -- the same derivations on every admissible triple with p <= %d ===" % BOUND)
cells = []
for q in primes_upto(BOUND):
    for mm_ in range(2, q):
        if (q - 1) % mm_:
            continue
        kk = (q - 1) // mm_
        if not (1 <= kk < q / 2):
            continue
        cells.append((q, kk, mm_))
N = len(cells)
tally = {}


def bump(name, ok):
    a, b = tally.get(name, (0, 0))
    tally[name] = (a + (1 if ok else 0), b + 1)


for (q, kk, mm_) in cells:
    om_, w_ = omega_and_word(q, kk)
    Pk, Qk, Rk = R_of(q, kk)
    K1k, K1k1 = K(w_, q, 1, kk), K(w_, q, 1, kk - 1)
    K2k, K2k1 = K(w_, q, 2, kk), K(w_, q, 2, kk - 1)
    bump("B01-admissible-1-le-k-lt-p-over-2", 1 <= kk < q / 2)
    bump("B02-eq2-antisymmetry",
         all((w_[2 * kk - i] + w_[i - 1]) % q == 0 for i in range(1, 2 * kk + 1)))
    # the reversal lemma, on the eight index pairs the proof uses
    okrev = True
    for (n_, r_) in ((1, kk), (2, kk), (1, kk - 1), (2, kk - 1),
                     (kk + 1, 2 * kk), (kk + 2, 2 * kk), (kk + 1, 2 * kk - 1), (kk + 2, 2 * kk - 1)):
        if n_ > r_ + 1:
            continue
        if K(w_, q, n_, r_) != tscal((-1) ** (r_ - n_ + 1), K(w_, q, 2 * kk + 1 - r_, 2 * kk + 1 - n_), q):
            okrev = False
    bump("B03-reversal-lemma", okrev)
    bump("B04-eq6-K-k+1-2k", K(w_, q, kk + 1, 2 * kk) == tscal((-1) ** kk, K1k, q))
    bump("B05-eq6-K-k+2-2k", K(w_, q, kk + 2, 2 * kk) == tscal((-1) ** (kk - 1), K1k1, q))
    bump("B06-eq4a-K-1-2k", K(w_, q, 1, 2 * kk) == tscal((-1) ** kk, Pk, q))
    bump("B07-eq4b-K-2-2k", K(w_, q, 2, 2 * kk) == tscal((-1) ** kk, Rk, q))
    bump("B08-lemma-a-K1k-squares", tsub(tmul(K1k, K1k, q), tmul(K1k1, K1k1, q), q) == Pk)
    bump("B09-lemma-b-mixed", tsub(tmul(K1k, K2k, q), tmul(K1k1, K2k1, q), q) == Rk)
    num = tsub(tmul(Rk, Rk, q), (1,), q)
    qq, rr = tdivmod(num, Pk, q) if num else ((), ())
    bump("B10-lemma-c-exact-division", rr == ())
    bump("B11-lemma-c-K2k-squares", tsub(tmul(K2k, K2k, q), tmul(K2k1, K2k1, q), q) == qq)
    gk, gk1 = teval(K1k, 1, q), teval(K1k1, 1, q)
    hk, hk1 = teval(K2k, 1, q), teval(K2k1, 1, q)
    bump("B12-at-T-1-g-squares-agree", (gk * gk - gk1 * gk1) % q == 0)
    bump("B13-at-T-1-mixed-is-one", (gk * hk - gk1 * hk1) % q == 1)
    bump("B14-at-T-1-h-squares-give-2omega", (hk * hk - hk1 * hk1) % q == 2 * om_ % q)
    bump("B15-h-squares-give-minus-2-binom-hockey-stick",
         (hk * hk - hk1 * hk1) % q == (-2 * ibinom((q - 1) // 2, kk)) % q)
    bump("B16-omega-closed-form-central-binomial",
         (2 * (-1) ** (kk + 1) * om_) % q == 2 * ibinom(2 * kk, kk) * pow(pow(4, kk, q), q - 2, q) % q)
    bump("B17-Q-at-1-is-binom", teval(Qk, 1, q) == ibinom((q - 1) // 2, kk) % q)
    if kk >= 2:
        bump("B18-eq5-determinant-at-T-1", (gk * hk1 - gk1 * hk) % q == (-1) ** kk % q)
    bump("B19-sign-g-k-is-minus-one-to-k+1-times-g-k-1", gk == ((-1) ** (kk + 1) * gk1) % q)
    bump("B20-g-k-1-nonzero", gk1 != 0)
    EA_, EB_ = eps_AB(q, kk, kk)
    bump("B21-epsA-equals-minus-one-to-k", EA_ == (-1) ** kk % q)
    U_, V_ = tadd(K1k, K1k1, q), tsub(K1k, K1k1, q)
    tm, tp = (1,), (1,)
    for _ in range(kk):
        tm = tmul(tm, (q - 1, 1), q)
        tp = tmul(tp, (1, 1), q)
    lead = 1
    for i in range(kk):
        lead = lead * w_[i] % q
    bump("B22-degrees-of-U-and-V", len(U_) - 1 == kk and len(V_) - 1 == kk)
    bump("B23-lead-squared-is-one", lead * lead % q == 1)
    if kk % 2 == 0:
        bump("B24-U-V-are-unit-times-T-minus-plus-1-to-k",
             U_ == tscal(lead, tm, q) and V_ == tscal(lead, tp, q))
    else:
        bump("B24-U-V-are-unit-times-T-minus-plus-1-to-k",
             V_ == tscal(lead, tm, q) and U_ == tscal(lead, tp, q))
    # lambda = U/V = a^k with a = (c-1)/(c+1), c = (-1)^k T
    c_ = (0, (-1) ** kk % q)
    am, ap = (1,), (1,)
    for _ in range(kk):
        am = tmul(am, tsub(c_, (1,), q), q)
        ap = tmul(ap, tadd(c_, (1,), q), q)
    bump("B25-multiplier-is-a-to-the-k", tmul(U_, ap, q) == tmul(V_, am, q))
    bump("B26-valuation-of-a-at-T-1-is-plus-minus-one",
         valuation_at(tsub(c_, (1,), q), 1, q) - valuation_at(tadd(c_, (1,), q), 1, q) in (1, -1))
    Pc_ = P_conj(q, kk, mm_)
    bump("B27-degP-equals-m", len(Pc_) - 1 == mm_)
    bump("B28-parity-split-identity", P_split(q, kk, mm_) == xmul([(2,)], Pc_, q))
    H_ = build_H(q, kk, kk, EA_)
    Qx, Rx = xdivmod_monic(H_, Pc_, q)
    bump("B29-P-divides-H", Rx == [])
    bump("B30-quotient-degree-is-p+1-minus-m", len(Qx) - 1 == q + 1 - mm_)
    # the initial partial quotient of alpha is forced by beta^m = a: matching the
    # coefficient of 1/T in 1 + 2m/(w_{k+1}T) + ... = a = 1 - 2(-1)^k/T + ... gives
    # w_{k+1} = (-1)^(k+1) m, which the word recurrence must reproduce on its own
    bump("B31-w-k+1-is-minus-one-to-k+1-times-m", w_[kk] == (-1) ** (kk + 1) * mm_ % q)

print("PASS B00-census-region %d admissible triples (p,(p-1)/m,m), p odd prime <= %d" % (N, BOUND))
CHECKS += 1
for name in sorted(tally):
    good, tot = tally[name]
    extra = ""
    if name.startswith("B18"):
        extra = (" (k >= 2 only: eq. (5) at (n,t,r) = (1,2,k) needs t <= r; the %d cells with k = 1 have "
                 "K_{1,0} = 1 and are covered by B19)" % (N - tot))
    ck(name, good == tot, "%d/%d cells%s" % (good, tot, extra))
print()

# ===========================================================================
print("=== PART C -- negative controls: these must NOT divide ===")
Pperturb = parse_X("X^4 + 4TX^3 + 6X^2 + 4TX + 2", 17)
_, Rp = xdivmod_monic(H_paper, Pperturb, 17)
ck("C01-perturbed-constant-term-fails", Rp != [],
   "X^4+4TX^3+6X^2+4TX+2 leaves remainder %s" % xstr(Rp))
H312 = build_H(3, 1, 1, eps_AB(3, 1, 1)[0])
Pbad = parse_X("X^2 + X + 1", 3)
_, Rb = xdivmod_monic(H312, Pbad, 3)
ck("C02-T-less-Prop-5-1-poly-fails", Rb != [],
   "over F_3 the T-less X^2+X+1 leaves remainder %s" % xstr(Rb))
H17b = build_H(17, 4, 4, PAPER_epsB)
_, Rb2 = xdivmod_monic(H17b, P_paper, 17)
ck("C03-epsB-instead-of-epsA-fails", Rb2 != [],
   "H([17,4,4,eps_B]) leaves remainder %s" % xstr(Rb2))
print()

# ===========================================================================
print("=== PART D -- published pins: the source author's own printed numbers ===")
om3, w3 = omega_and_word(3, 1)
ck("D01-omega-1-and-word-at-p-3", om3 == 2 and list(w3) == [1, 2],
   "omega_1 = 1/2 = 2 and W = T, -T over F_3")
ck("D02-printed-H1", build_H(3, 1, 1, 1) == parse_X("X^4 + TX^3 + TX + 1", 3),
   "H([3,1,1,+1]) = X^4+TX^3+TX+1")
ck("D03-printed-H2", build_H(3, 1, 1, 2) == parse_X("X^4 + TX^3 + 2TX + 2", 3),
   "H([3,1,1,-1]) = X^4+TX^3-TX-1")
ck("D04-printed-eps-at-p-3", eps_AB(3, 1, 1) == (2, 1), "eps_A = -1, eps_B = +1 at (3,1,1)")
Qa, Ra = xdivmod_monic(H312, P_conj(3, 1, 2), 3)
ck("D05-printed-factorisation-at-p-3", Ra == [] and Qa == parse_X("X^2 + 2", 3),
   "H([3,1,1,-1]) = (X^2+TX+1)(X^2-1)")
ck("D06-Conj-5-4-epsA", eps_AB(13, 4, 2)[0] == 12, "eps_A(13,4,2) = -1")
_, R54 = xdivmod_monic(build_H(13, 4, 2, 12), parse_X("X^4 + 6TX^3 + 2X^2 + 4", 13), 13)
ck("D07-Conj-5-4-printed-quartic-divides", R54 == [],
   "X^4+6TX^3+2X^2+4 divides H([13,4,2,-1]) with remainder 0")
ext_ok = True
for q in (7, 11, 13, 17, 19, 23, 31, 43):
    kk = (q - 1) // 2
    om_, w_ = omega_and_word(q, kk)
    if om_ != q - 1:
        ext_ok = False
    if not all(w_[i] == 2 * (-1) ** (i + 1) % q for i in range(2 * kk)):
        ext_ok = False
    if not all(eps_AB(q, kk, jj)[0] == (-1) ** kk % q for jj in range(1, 2 * kk)):
        ext_ok = False
ck("D08-extremal-case-k-equals-p-minus-1-over-2", ext_ok,
   "omega_k = -1, W = -2T,2T,...,2T and eps_A = (-1)^k for every j, at p = 7..43")
print()

# ===========================================================================
EBOUND = min(BOUND, 60)
ETARGET = 20
ecells = [c for c in cells if c[0] <= EBOUND]
print("=== PART E -- the Proposition 9 root alpha, RECONSTRUCTED in F_p((1/T)) on the %d admissible"
      " triples with p <= %d ===" % (len(ecells), EBOUND))
print("    (Newton in u = 1/T from alpha_0 = w_{k+1}T; every identity below is certified to order"
      " u^%d, and every series knows the exponent below which it is exact)" % ETARGET)
etally = {}


def ebump(name, ok):
    a, b = etally.get(name, (0, 0))
    etally[name] = (a + (1 if ok else 0), b + 1)


def enewton(x, Hs, dHs, q):
    """Newton iteration for H in F_p((1/T)).  With e = x - alpha one has the exact
    identity H(x) = e*H'(alpha) + e^p*(K2*alpha - K1) + K2*e^(p+1), whence the next
    error is -e^p (K2 alpha - K1)/H'(x): the valuation of the error is multiplied by
    p and shifted by p+k-1 > 0, so v(e) >= 0 suffices for convergence."""
    for _ in range(40):
        hx = Hs(x)
        if hx[1] == ():
            break
        corr = smul(hx, sinv(dHs(x), q), q)
        if corr[1] == ():
            break
        nx = ssub(x, corr, q)
        if nx == x:
            break
        x = nx
    return x


ctrl_perturb = ctrl_wrongpq = ctrl_paperP = ctrl_epsB = None
alpha_witness = None
for (q, kk, mm_) in ecells:
    SCAP = q + 2 * kk + 64
    om_, w_ = omega_and_word(q, kk)
    eps, epsB = eps_AB(q, kk, kk)              # eps_A from its continued-fraction definition
    # the four coefficient continuants of H(V), V = [p,k,k,eps_A], from the definition
    K2, K1 = K(w_, q, kk + 2, 2 * kk), K(w_, q, kk + 1, 2 * kk)
    Kk, Kk1 = K(w_, q, 1, kk), K(w_, q, 1, kk - 1)
    sK2, sK1, sKk, sKk1 = (sfromT(K2, q), sfromT(K1, q), sfromT(Kk, q), sfromT(Kk1, q))

    def Hs(x, q=q, sK2=sK2, sK1=sK1, sKk=sKk, sKk1=sKk1, eps=eps):
        """H(V)(x) = K_{k+2,2k} x^(p+1) - K_{k+1,2k} x^p + eps(K_{1,k} x + K_{1,k-1})."""
        return sadd(smul(sfrob(x, q), ssub(smul(sK2, x, q), sK1, q), q),
                    sscal(eps, sadd(smul(sKk, x, q), sKk1, q), q), q)

    def dHs(x, q=q, sK2=sK2, sKk=sKk, eps=eps):
        """H'(V)(x) = K_{k+2,2k} x^p + eps K_{1,k}, since p+1 = 1 and p = 0 in F_p."""
        return sadd(smul(sK2, sfrob(x, q), q), sscal(eps, sKk, q), q)

    x0 = snorm(-1, [w_[kk]], SCAP, q)          # alpha_0 = w_{k+1} T
    al = enewton(x0, Hs, dHs, q)
    hx, dhx = Hs(al), dHs(al)
    ebump("E01-alpha-is-a-root-of-H", svanishes(hx, ETARGET))
    ebump("E02-hensel-criterion-gives-an-exact-root",
          sval(hx) - 2 * sval(dhx) > 0 and sval(hx) - sval(dhx) >= 1)
    ebump("E03-first-partial-quotient-is-w-k+1-T",
          spolypart(al) == tnorm([0, w_[kk] % q]) and sval(al) == -1)
    # the diagonalised Frobenius relation of the proof: alpha^p = (A alpha + B)/(B alpha + A)
    sA, sB = sfromT(tscal(-eps, Kk, q), q), sfromT(tscal(-eps, Kk1, q), q)
    ebump("E04-frobenius-relation-is-diagonal",
          svanishes(ssub(smul(sfrob(al, q), sadd(smul(sB, al, q), sA, q), q),
                         sadd(smul(sA, al, q), sB, q), q), ETARGET))
    one = sconst(1, q)
    sc = sfromT(tnorm([0, (-1) ** kk % q]), q)          # c = (-1)^k T
    beta = sdiv(sadd(al, one, q), ssub(al, one, q), q)
    aa = sdiv(ssub(sc, one, q), sadd(sc, one, q), q)    # a = (c-1)/(c+1)
    ebump("E05-beta-and-a-are-1-mod-1-over-T",
          sval(ssub(beta, one, q)) >= 1 and sval(ssub(aa, one, q)) >= 1)
    bm = spow(beta, mm_, q)
    ebump("E06-beta-to-the-m-equals-a-so-zeta-is-1", svanishes(ssub(bm, aa, q), ETARGET))
    ebump("E07-beta-to-the-p-equals-a-to-the-k-beta",
          svanishes(ssub(sfrob(beta, q), smul(spow(aa, kk, q), beta, q), q), ETARGET))
    Pc_ = P_conj(q, kk, mm_)
    acc = sfromT(Pc_[-1], q)
    for i in range(len(Pc_) - 2, -1, -1):
        acc = sadd(smul(acc, al, q), sfromT(Pc_[i], q), q)
    ebump("E08-P-vanishes-at-alpha", svanishes(acc, ETARGET))
    roots = [z for z in range(1, q) if pow(z, kk, q) == 1]
    ebump("E09-k-th-roots-of-unity-split-in-F-p", len(roots) == kk)
    if kk >= 2:
        # the k-1 sibling degree-m factors of H(V), one per zeta != 1 with zeta^k = 1,
        # are the ones exact division cannot exclude.  P_zeta(alpha) is a nonzero
        # multiple of (alpha+1)^m - zeta a (alpha-1)^m = (alpha-1)^m (beta^m - zeta a),
        # of valuation exactly -m: none of them vanishes at alpha.
        xp1m, xm1m = spow(sadd(al, one, q), mm_, q), spow(ssub(al, one, q), mm_, q)
        okz = True
        for z in roots:
            if z == 1:
                continue
            if sval(ssub(xp1m, smul(sscal(z, aa, q), xm1m, q), q)) != -mm_:
                okz = False
        ebump("E10-sibling-factors-do-not-vanish-at-alpha", okz)
    if (q, kk, mm_) == (PAPER_p, PAPER_k, PAPER_m):
        alpha_witness = sstr(al, 5)
        # the P PRINTED IN THE PAPER, not the recomputed one, evaluated at alpha
        accp = sfromT(P_paper[-1], q)
        for i in range(len(P_paper) - 2, -1, -1):
            accp = sadd(smul(accp, al, q), sfromT(P_paper[i], q), q)
        ctrl_paperP = svanishes(accp, ETARGET)
        # control 1: move alpha by 1/T^3 -- it must stop being a root
        bad1 = sadd(al, snorm(3, [1], SCAP, q), q)
        hbad1 = Hs(bad1)
        ctrl_perturb = (not svanishes(hbad1, ETARGET), sval(hbad1))
        # control 2: start from the WRONG initial partial quotient (w_{k+1}+1)T -- no
        # root of H(V) has that polynomial part, so E01+E03 must not both hold
        bad2 = enewton(snorm(-1, [(w_[kk] + 1) % q], SCAP, q), Hs, dHs, q)
        ctrl_wrongpq = not (svanishes(Hs(bad2), ETARGET)
                            and spolypart(bad2) == tnorm([0, (w_[kk] + 1) % q]))
        # control 3: alpha belongs to eps_A.  H([p,k,k,eps_B])(alpha) = (eps_B - eps_A)
        # (K_{1,k} alpha + K_{1,k-1}), of valuation -k-1: alpha is not its root either
        hB = sadd(smul(sfrob(al, q), ssub(smul(sK2, al, q), sK1, q), q),
                  sscal(epsB, sadd(smul(sKk, al, q), sKk1, q), q), q)
        ctrl_epsB = (not svanishes(hB, ETARGET), sval(hB))
for name in sorted(etally):
    good, tot = etally[name]
    extra = ""
    if name.startswith("E10"):
        extra = " (k >= 2 only: at k = 1 there is no sibling factor)"
    ck(name, good == tot, "%d/%d cells%s" % (good, tot, extra))
ck("E11-paper-printed-P-vanishes-at-witness-alpha", bool(ctrl_paperP),
   "P as printed in the paper kills the reconstructed alpha = %s at (17,4,4)" % alpha_witness)
ck("E12-control-perturbed-alpha-is-not-a-root", bool(ctrl_perturb and ctrl_perturb[0]),
   "alpha + 1/T^3 gives H(V) of valuation %s, not = 0 mod u^%d"
   % ((ctrl_perturb[1] if ctrl_perturb else "?"), ETARGET))
ck("E13-control-wrong-first-partial-quotient-fails", bool(ctrl_wrongpq),
   "Newton from (w_{k+1}+1)T does not deliver a root of H(V) with that polynomial part")
ck("E14-control-alpha-is-not-a-root-of-the-epsB-form", bool(ctrl_epsB and ctrl_epsB[0]),
   "H([17,4,4,eps_B])(alpha) has valuation %s, not = 0 mod u^%d"
   % ((ctrl_epsB[1] if ctrl_epsB else "?"), ETARGET))
print()

# ===========================================================================
print("NOT RE-RUN: the continued-fraction CONSTRUCTION of alpha in Proposition 5.1 of the source.")
print("NOT RE-RUN:   PART E does not follow that construction; it solves H(V)(alpha) = 0 in F_p((1/T))")
print("NOT RE-RUN:   by Newton from alpha_0 = w_{k+1}T and then verifies the two properties the paper")
print("NOT RE-RUN:   actually uses -- H(V)(alpha) = 0 and first partial quotient w_{k+1}T, hence a pole")
print("NOT RE-RUN:   at T = infinity.  Uniqueness of the root with that polynomial part is not proved.")
print("NOT RE-RUN: exactness of alpha. PART E certifies its identities to order u^%d only; the step" % ETARGET)
print("NOT RE-RUN:   from that to an exact root is Hensel's lemma in the complete field F_p((1/T)),")
print("NOT RE-RUN:   whose hypothesis v(H) > 2 v(H') is re-derived (E02) and whose conclusion is cited.")
print("NOT RE-RUN: PART E for p > %d: the alpha census stops there, PARTS A-D run to p <= %d."
      % (EBOUND, BOUND))
print("NOT RE-RUN: irreducibility of P over F_p(T). Only the hypothesis of Capelli's criterion is")
print("NOT RE-RUN:   re-derived here (A20, B26: v_{T-1}(a) = +-1, coprime to m); the criterion is cited.")
print("NOT RE-RUN: nu(alpha) = m, the third clause of the conjecture. That is Corollary 4.2 of the")
print("NOT RE-RUN:   source (published Corollary 5, from Ayadi-Lasjaunias 2016) and is not our result.")
print("NOT RE-RUN: the full factorisation multisets of H(V) over F_p(T) -- no factorisation algorithm")
print("NOT RE-RUN:   is implemented; only exact division by a given monic P.")
print("NOT RE-RUN: primes p > %d. The theorem is proved for the whole infinite family; the census" % BOUND)
print("NOT RE-RUN:   above exhausts 0% of it and is corroboration, not evidence for the theorem.")
print()
print("elapsed %.1f s" % (time.time() - t0))
if FAILED:
    print("FAILURES: %s" % ", ".join(FAILED))
    print("VERDICT: %d OF %d CHECKS FAILED" % (len(FAILED), CHECKS))
    sys.exit(1)
print("VERDICT: ALL %d CHECKS PASS" % CHECKS)
sys.exit(0)
