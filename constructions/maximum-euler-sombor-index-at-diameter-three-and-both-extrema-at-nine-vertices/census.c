/* census.c -- the two enumeration kernels behind the nine-vertex diameter-three census
 * of "The Maximum Euler-Sombor Index at Diameter Three, and Both Extrema at Nine Vertices".
 *
 * It computes no square-root comparison that any conclusion rests on: the kernels only
 * ENUMERATE and BUCKET.  All exact arithmetic, all published-value controls and the verdict
 * are done by the driver census.py, which calls this binary twice.
 *
 *   mode "classes":  read graph6 lines for n=9 on stdin (produced by `geng -q -c 9`),
 *                    compute for each graph its exact diameter (bitmask ball growth, no
 *                    floating point at all) and its Euler-Sombor EDGE PROFILE, i.e. the
 *                    multiset of unordered degree pairs {d(x),d(y)} over the edges.  The
 *                    profile determines EU(G) exactly, so the driver can evaluate the whole
 *                    census in 80-digit Decimal arithmetic by evaluating one representative
 *                    per distinct profile.  Emits one line per (diameter, profile) bucket
 *                    with the number of ISOMORPHISM CLASSES in it (geng emits one graph6
 *                    per class) and up to four representative graph6 strings.
 *                    => this pass is where the 148,229 count and the uniqueness/tie counts
 *                       come from.
 *
 *   mode "labelled": an independent, isomorphism-engine-free sweep of ALL 2^36 =
 *                    68,719,476,736 labelled graphs on the vertex set {0,...,8}.  Reports
 *                    the labelled count per diameter (the connected total is checked by the
 *                    driver against the exponential-formula value for A001187(9)), the
 *                    double-precision extrema of EU over the diameter-3 cell, and -- for
 *                    windows of width `win` around two thresholds supplied on the command
 *                    line -- the number of labelled graphs and the distinct profiles inside
 *                    them.  Nothing outside the windows may beat the thresholds; if anything
 *                    does, the driver fails the run.  This pass uses no geng, no canonical
 *                    form and no isomorphism test whatsoever.
 *
 * build: gcc -O3 -march=native -fopenmp -o census census.c -lm
 * usage: geng -q -c 9 | ./census classes  > classes.txt
 *        ./census labelled <emin> <emax> <win> > labelled.txt
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>

#define N 9
#define FULL ((uint16_t)0x1FFu)
#define NEDGE 36
#define NSLOT 81 /* degree pair (a,b), 1<=a<=b<=8, stored at a*9+b */

/* ------------------------------------------------------------------ *
 * exact diameter of a graph on 9 vertices; -1 if disconnected.
 * cur[v] is the ball of radius k about v; integer arithmetic only.
 * ------------------------------------------------------------------ */
static inline int diameter9(const uint16_t *adj)
{
    uint16_t Nb[N], cur[N], nxt[N];
    int v, k, u;
    for (v = 0; v < N; v++) Nb[v] = (uint16_t)(adj[v] | (uint16_t)(1u << v));
    for (v = 0; v < N; v++) cur[v] = Nb[v];
    for (k = 1; k <= N; k++) {
        int all = 1, changed = 0;
        for (v = 0; v < N; v++) if (cur[v] != FULL) { all = 0; break; }
        if (all) return k;
        for (v = 0; v < N; v++) {
            uint16_t s = cur[v], acc = cur[v];
            while (s) { u = __builtin_ctz(s); s &= (uint16_t)(s - 1); acc |= Nb[u]; }
            nxt[v] = acc;
            if (acc != cur[v]) changed = 1;
        }
        if (!changed) return -1;               /* balls stopped growing before covering V */
        memcpy(cur, nxt, sizeof(nxt));
    }
    return -1;
}

/* the 36 vertex pairs, in graph6 bit order */
static int PI[NEDGE], PJ[NEDGE];
static void init_pairs(void)
{
    int k = 0, i, j;
    for (j = 1; j < N; j++) for (i = 0; i < j; i++) { PI[k] = i; PJ[k] = j; k++; }
}

static double EUT[9][9];
static void init_eu(void)
{
    int a, b;
    for (a = 1; a <= 8; a++) for (b = 1; b <= 8; b++)
        EUT[a][b] = sqrt((double)(a * a + a * b + b * b));
}

/* ================================================================== *
 * mode "classes"
 * ================================================================== */

typedef struct {
    uint8_t  diam;
    uint8_t  prof[NSLOT];
    uint32_t count;
    uint8_t  nrep;
    char     rep[4][16];
    uint8_t  used;
} Ent;

#define TBITS 21
#define TSIZE (1u << TBITS)

static Ent *tab;

static uint64_t fnv(const uint8_t *p, int n)
{
    uint64_t h = 1469598103934665603ULL;
    int i;
    for (i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}

static void bump(uint8_t diam, const uint8_t *prof, const char *g6)
{
    uint8_t key[NSLOT + 1];
    uint64_t h;
    uint32_t idx;
    key[0] = diam;
    memcpy(key + 1, prof, NSLOT);
    h = fnv(key, NSLOT + 1);
    idx = (uint32_t)(h & (TSIZE - 1));
    for (;;) {
        Ent *e = &tab[idx];
        if (!e->used) {
            e->used = 1; e->diam = diam; memcpy(e->prof, prof, NSLOT);
            e->count = 1; e->nrep = 1;
            snprintf(e->rep[0], 16, "%s", g6);
            return;
        }
        if (e->diam == diam && memcmp(e->prof, prof, NSLOT) == 0) {
            e->count++;
            if (e->nrep < 4) { snprintf(e->rep[e->nrep], 16, "%s", g6); e->nrep++; }
            return;
        }
        idx = (idx + 1) & (TSIZE - 1);
    }
}

static int run_classes(void)
{
    char line[64];
    uint64_t nread = 0, nbad = 0, ndisc = 0;
    tab = (Ent *)calloc(TSIZE, sizeof(Ent));
    if (!tab) { fprintf(stderr, "oom\n"); return 2; }

    while (fgets(line, sizeof(line), stdin)) {
        uint16_t adj[N];
        uint8_t deg[N], prof[NSLOT];
        int k, v, d;
        size_t L = strlen(line);
        while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
        if (L == 0) continue;
        if (L != 7 || line[0] != 'H') { nbad++; continue; }   /* n=9 graph6 is 'H' + 6 bytes */
        for (v = 0; v < N; v++) adj[v] = 0;
        for (k = 0; k < NEDGE; k++) {
            int byte = 1 + k / 6, bit = 5 - (k % 6);
            if (((line[byte] - 63) >> bit) & 1) {
                adj[PI[k]] |= (uint16_t)(1u << PJ[k]);
                adj[PJ[k]] |= (uint16_t)(1u << PI[k]);
            }
        }
        nread++;
        d = diameter9(adj);
        if (d < 0) { ndisc++; continue; }
        for (v = 0; v < N; v++) deg[v] = (uint8_t)__builtin_popcount(adj[v]);
        memset(prof, 0, sizeof(prof));
        for (k = 0; k < NEDGE; k++)
            if ((adj[PI[k]] >> PJ[k]) & 1) {
                int a = deg[PI[k]], b = deg[PJ[k]], lo = a < b ? a : b, hi = a < b ? b : a;
                prof[lo * 9 + hi]++;
            }
        bump((uint8_t)d, prof, line);
    }

    printf("READ %llu\n", (unsigned long long)nread);
    printf("MALFORMED %llu\n", (unsigned long long)nbad);
    printf("DISCONNECTED %llu\n", (unsigned long long)ndisc);
    {
        uint32_t i;
        for (i = 0; i < TSIZE; i++) {
            Ent *e = &tab[i];
            int a, b, r, m = 0;
            if (!e->used) continue;
            for (a = 1; a <= 8; a++) for (b = a; b <= 8; b++) m += e->prof[a * 9 + b];
            printf("BUCKET %u %d %u %u", e->diam, m, e->count, e->nrep);
            for (r = 0; r < e->nrep; r++) printf(" %s", e->rep[r]);
            printf(" |");
            for (a = 1; a <= 8; a++) for (b = a; b <= 8; b++)
                if (e->prof[a * 9 + b]) printf(" %d,%d:%u", a, b, e->prof[a * 9 + b]);
            printf("\n");
        }
    }
    printf("END\n");
    free(tab);
    return 0;
}

/* ================================================================== *
 * mode "labelled": all 2^36 labelled graphs on {0..8}
 * ================================================================== */

/* chunk tables: 6 chunks of 6 bits; T[c][w][v] = adjacency contribution of chunk value w */
static uint16_t T[6][64][N];
static void init_chunks(void)
{
    int c, w, b, v;
    for (c = 0; c < 6; c++) for (w = 0; w < 64; w++) {
        for (v = 0; v < N; v++) T[c][w][v] = 0;
        for (b = 0; b < 6; b++) if ((w >> b) & 1) {
            int k = c * 6 + b;
            T[c][w][PI[k]] |= (uint16_t)(1u << PJ[k]);
            T[c][w][PJ[k]] |= (uint16_t)(1u << PI[k]);
        }
    }
}

#define MAXWIN 64
typedef struct {
    uint64_t percut[11];        /* index 0 = disconnected, d = diameter */
    double   emin, emax;
    uint64_t lomask, himask;    /* a witness mask at each extreme */
    uint64_t nlo, nhi;          /* labelled graphs inside the two windows */
    uint64_t below, above;      /* labelled graphs strictly beating the thresholds */
    int      nplo, nphi;
    uint8_t  plo[MAXWIN][NSLOT], phi[MAXWIN][NSLOT];
    uint64_t clo[MAXWIN], chi[MAXWIN];
} Acc;

static void addprof(uint8_t P[][NSLOT], uint64_t *C, int *n, const uint8_t *prof)
{
    int i;
    for (i = 0; i < *n; i++) if (memcmp(P[i], prof, NSLOT) == 0) { C[i]++; return; }
    if (*n < MAXWIN) { memcpy(P[*n], prof, NSLOT); C[*n] = 1; (*n)++; }
}

static int run_labelled(double emin, double emax, double win)
{
    int nth, t, i;
    Acc *acc;
    long long hi;
    const long long NHI = 1LL << 12;     /* bits 24..35 */
    const long long NLO = 1LL << 24;     /* bits 0..23  */

    nth = omp_get_max_threads();
    acc = (Acc *)calloc((size_t)nth, sizeof(Acc));
    for (t = 0; t < nth; t++) { acc[t].emin = 1e300; acc[t].emax = -1e300; }

    #pragma omp parallel for schedule(dynamic, 1)
    for (hi = 0; hi < NHI; hi++) {
        Acc *A = &acc[omp_get_thread_num()];
        uint16_t base[N];
        int v;
        uint16_t *T4 = T[4][(hi) & 63], *T5 = T[5][(hi >> 6) & 63];
        for (v = 0; v < N; v++) base[v] = (uint16_t)(T4[v] | T5[v]);
        for (long long lo = 0; lo < NLO; lo++) {
            uint16_t adj[N];
            uint16_t *t0 = T[0][lo & 63], *t1 = T[1][(lo >> 6) & 63],
                     *t2 = T[2][(lo >> 12) & 63], *t3 = T[3][(lo >> 18) & 63];
            int d, k;
            for (v = 0; v < N; v++) adj[v] = (uint16_t)(base[v] | t0[v] | t1[v] | t2[v] | t3[v]);
            d = diameter9(adj);
            A->percut[d < 0 ? 0 : d]++;
            if (d != 3) continue;
            {
                int deg[N];
                double s = 0.0;
                for (v = 0; v < N; v++) deg[v] = __builtin_popcount(adj[v]);
                for (k = 0; k < NEDGE; k++)
                    if ((adj[PI[k]] >> PJ[k]) & 1) s += EUT[deg[PI[k]]][deg[PJ[k]]];
                if (s < A->emin) { A->emin = s; A->lomask = (uint64_t)((hi << 24) | lo); }
                if (s > A->emax) { A->emax = s; A->himask = (uint64_t)((hi << 24) | lo); }
                if (s < emin - win) A->below++;
                if (s > emax + win) A->above++;
                if (s <= emin + win || s >= emax - win) {
                    uint8_t prof[NSLOT];
                    memset(prof, 0, sizeof(prof));
                    for (k = 0; k < NEDGE; k++)
                        if ((adj[PI[k]] >> PJ[k]) & 1) {
                            int a = deg[PI[k]], b = deg[PJ[k]];
                            int l = a < b ? a : b, h2 = a < b ? b : a;
                            prof[l * 9 + h2]++;
                        }
                    if (s <= emin + win) { A->nlo++; addprof(A->plo, A->clo, &A->nplo, prof); }
                    else                 { A->nhi++; addprof(A->phi, A->chi, &A->nphi, prof); }
                }
            }
        }
    }

    {
        Acc R;
        memset(&R, 0, sizeof(R));
        R.emin = 1e300; R.emax = -1e300;
        for (t = 0; t < nth; t++) {
            for (i = 0; i <= 10; i++) R.percut[i] += acc[t].percut[i];
            if (acc[t].emin < R.emin) { R.emin = acc[t].emin; R.lomask = acc[t].lomask; }
            if (acc[t].emax > R.emax) { R.emax = acc[t].emax; R.himask = acc[t].himask; }
            R.nlo += acc[t].nlo; R.nhi += acc[t].nhi;
            R.below += acc[t].below; R.above += acc[t].above;
            for (i = 0; i < acc[t].nplo; i++) {
                int j, f = -1;
                for (j = 0; j < R.nplo; j++) if (memcmp(R.plo[j], acc[t].plo[i], NSLOT) == 0) f = j;
                if (f < 0 && R.nplo < MAXWIN) { f = R.nplo++; memcpy(R.plo[f], acc[t].plo[i], NSLOT); R.clo[f] = 0; }
                if (f >= 0) R.clo[f] += acc[t].clo[i];
            }
            for (i = 0; i < acc[t].nphi; i++) {
                int j, f = -1;
                for (j = 0; j < R.nphi; j++) if (memcmp(R.phi[j], acc[t].phi[i], NSLOT) == 0) f = j;
                if (f < 0 && R.nphi < MAXWIN) { f = R.nphi++; memcpy(R.phi[f], acc[t].phi[i], NSLOT); R.chi[f] = 0; }
                if (f >= 0) R.chi[f] += acc[t].chi[i];
            }
        }
        printf("THREADS %d\n", nth);
        printf("SWEPT %llu\n", 1ULL << 36);
        printf("DISCONNECTED %llu\n", (unsigned long long)R.percut[0]);
        for (i = 1; i <= 8; i++)
            printf("LABELLED_DIAM %d %llu\n", i, (unsigned long long)R.percut[i]);
        printf("DOUBLE_MIN %.17g mask %llu\n", R.emin, (unsigned long long)R.lomask);
        printf("DOUBLE_MAX %.17g mask %llu\n", R.emax, (unsigned long long)R.himask);
        printf("BELOW_THRESHOLD %llu\n", (unsigned long long)R.below);
        printf("ABOVE_THRESHOLD %llu\n", (unsigned long long)R.above);
        printf("WINDOW_LO_COUNT %llu\n", (unsigned long long)R.nlo);
        printf("WINDOW_HI_COUNT %llu\n", (unsigned long long)R.nhi);
        for (i = 0; i < R.nplo; i++) {
            int a, b;
            printf("WINDOW_LO_PROFILE %llu |", (unsigned long long)R.clo[i]);
            for (a = 1; a <= 8; a++) for (b = a; b <= 8; b++)
                if (R.plo[i][a * 9 + b]) printf(" %d,%d:%u", a, b, R.plo[i][a * 9 + b]);
            printf("\n");
        }
        for (i = 0; i < R.nphi; i++) {
            int a, b;
            printf("WINDOW_HI_PROFILE %llu |", (unsigned long long)R.chi[i]);
            for (a = 1; a <= 8; a++) for (b = a; b <= 8; b++)
                if (R.phi[i][a * 9 + b]) printf(" %d,%d:%u", a, b, R.phi[i][a * 9 + b]);
            printf("\n");
        }
        printf("END\n");
    }
    free(acc);
    return 0;
}

int main(int argc, char **argv)
{
    init_pairs();
    init_eu();
    if (argc >= 2 && strcmp(argv[1], "classes") == 0) return run_classes();
    if (argc >= 5 && strcmp(argv[1], "labelled") == 0) {
        init_chunks();
        return run_labelled(atof(argv[2]), atof(argv[3]), atof(argv[4]));
    }
    fprintf(stderr, "usage: %s classes   (graph6 on stdin)\n"
                    "       %s labelled <emin> <emax> <win>\n", argv[0], argv[0]);
    return 2;
}
