/*
 * D01_bfs_entropy.c — BFS for 64GB RAM machines
 *
 * Hash table sized at 2^(D+2) bits, capped at 2^29 (~45GB total).
 * Tracks S_0(k) reliably through the peak for D up to ~28.
 *
 * Memory usage per D (two hash tables, 21 bytes/entry each):
 *   D=22: 2^24 = 16M entries  →   672 MB
 *   D=24: 2^26 = 64M entries  →   2.7 GB
 *   D=26: 2^28 = 268M entries →  11.2 GB
 *   D=28: 2^29 = 537M entries →  22.5 GB  (capped)
 *   D=30: 2^29 = 537M entries →  22.5 GB  (may abort for hard semiprimes)
 *
 * Compile: cc -O3 -o D01_bfs_entropy D01_bfs_entropy.c -lm
 * Run:     ./D01_bfs_entropy [max_D] [min_D]
 *          default: min_D=10, max_D=30
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

typedef struct {
    uint64_t pp;
    uint64_t qq;
    uint32_t carry;
} State;

static State *ht[2] = {NULL, NULL};
static uint8_t *occ[2] = {NULL, NULL};
static int cnt[2] = {0, 0};
static uint64_t ht_mask = 0;
static int64_t ht_size = 0;
static int current_hbits = 0;

static void ht_alloc(int bits) {
    if (bits == current_hbits) return;
    ht_size = 1LL << bits;
    ht_mask = ht_size - 1;
    for (int i = 0; i < 2; i++) {
        free(ht[i]); free(occ[i]);
        ht[i] = (State*)malloc((size_t)ht_size * sizeof(State));
        occ[i] = (uint8_t*)malloc((size_t)ht_size);
        if (!ht[i] || !occ[i]) {
            fprintf(stderr, "OOM allocating 2^%d = %lld entries (%.1f GB per table)\n",
                    bits, (long long)ht_size,
                    (double)ht_size * sizeof(State) / 1e9);
            exit(1);
        }
    }
    current_hbits = bits;
    fprintf(stderr, "  Allocated 2^%d = %lld entries (%.1f GB total)\n",
            bits, (long long)ht_size,
            2.0 * (double)ht_size * (sizeof(State) + 1) / 1e9);
}

static inline void ht_clear(int idx) {
    memset(occ[idx], 0, ht_size);
    cnt[idx] = 0;
}

static inline uint64_t hash3(uint64_t p, uint64_t q, uint32_t c) {
    uint64_t h = p * 2654435761ULL ^ q * 2246822519ULL ^ (uint64_t)c * 3266489917ULL;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    return h & ht_mask;
}

static inline int ht_insert(int idx, uint64_t pp, uint64_t qq, uint32_t carry) {
    uint64_t i = hash3(pp, qq, carry);
    for (int probe = 0; probe < 2048; probe++) {
        uint64_t j = (i + probe) & ht_mask;
        if (!occ[idx][j]) {
            ht[idx][j] = (State){pp, qq, carry};
            occ[idx][j] = 1;
            cnt[idx]++;
            return 1;
        }
        if (ht[idx][j].pp == pp && ht[idx][j].qq == qq && ht[idx][j].carry == carry)
            return 0;
    }
    return -1;
}

static int is_prime(uint64_t n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (uint64_t i = 5; i * i <= n; i += 6)
        if (n % i == 0 || n % (i+2) == 0) return 0;
    return 1;
}

static uint64_t random_prime_range(uint64_t lo, uint64_t hi) {
    if (lo < 3) lo = 3;
    if (lo % 2 == 0) lo++;
    if (hi <= lo) { if (is_prime(lo)) return lo; return 0; }
    for (int a = 0; a < 100000; a++) {
        uint64_t r = (hi - lo) / 2;
        if (r == 0) r = 1;
        uint64_t c = lo + 2 * ((uint64_t)rand() % r);
        if (c % 2 == 0) c++;
        if (c > hi) c = hi;
        if (c < lo) c = lo;
        if (is_prime(c)) return c;
    }
    for (uint64_t c = lo; c <= hi; c += 2)
        if (is_prime(c)) return c;
    return 0;
}

typedef struct {
    double h_peak;
    int max_c0;
    int peak_k;
    int last_k;
    int peak_reliable;
    int aborted;
} BfsResult;

static BfsResult bfs_entropy(uint64_t N, int hash_bits) {
    int D = 64 - __builtin_clzll(N);
    int max_pos = D + 3;
    uint8_t Nb[64];
    for (int i = 0; i < 64; i++) Nb[i] = (N >> i) & 1;

    ht_alloc(hash_bits);
    int cur = 0, nxt = 1;
    ht_clear(cur);

    uint64_t seed = hash3(0, 0, 0);
    ht[cur][seed] = (State){0, 0, 0};
    occ[cur][seed] = 1;
    cnt[cur] = 1;

    BfsResult res = {0.0, 1, 0, 0, 0, 0};
    int steps_past_peak = 0;

    for (int k = 0; k < max_pos; k++) {
        uint8_t target = (k < 64) ? Nb[k] : 0;
        ht_clear(nxt);

        for (int64_t idx = 0; idx < ht_size; idx++) {
            if (!occ[cur][idx]) continue;
            uint64_t pp = ht[cur][idx].pp;
            uint64_t qq = ht[cur][idx].qq;
            uint32_t cin = ht[cur][idx].carry;

            /* Optimized convolution: use popcount on AND of p with bit-reversed q */
            int base_conv = 0;
            if (k >= 2) {
                uint64_t rev_qq = 0;
                for (int i = 1; i < k; i++)
                    if ((qq >> (k - i)) & 1) rev_qq |= (1ULL << i);
                base_conv = __builtin_popcountll(pp & rev_qq);
            }
            uint32_t qq0 = qq & 1, pp0 = pp & 1;

            for (int a = 0; a <= 1; a++) {
                for (int b = 0; b <= 1; b++) {
                    int conv = (k == 0) ? (a * b) : (base_conv + a * qq0 + pp0 * b);
                    int tot = conv + cin;
                    if ((tot & 1) == target) {
                        int rc = ht_insert(nxt,
                            pp | ((uint64_t)a << k),
                            qq | ((uint64_t)b << k),
                            (uint32_t)(tot >> 1));
                        if (rc < 0) {
                            res.aborted = 1;
                            goto done;
                        }
                    }
                }
            }
        }

        /* Count carry-zero states */
        int c0 = 0;
        for (int64_t i = 0; i < ht_size; i++)
            if (occ[nxt][i] && ht[nxt][i].carry == 0) c0++;

        double h = (c0 > 1) ? log2((double)c0) : 0.0;

        if (c0 > res.max_c0) {
            res.max_c0 = c0;
            res.h_peak = h;
            res.peak_k = k;
            steps_past_peak = 0;
        } else {
            steps_past_peak++;
        }
        res.last_k = k;

        /* Swap */
        int tmp = cur; cur = nxt; nxt = tmp;

        if (cnt[cur] == 0) break;
        /* Abort if table > 75% full */
        if (cnt[cur] > ht_size * 3 / 4) {
            res.aborted = 1;
            goto done;
        }
        /* Early stop: well past peak and past midpoint */
        if (steps_past_peak >= 5 && k > D / 2 + 2) break;
    }
done:
    res.peak_reliable = (res.last_k >= res.peak_k + 2) ? 1 : 0;
    return res;
}

int main(int argc, char **argv) {
    srand(time(NULL));

    int max_D = 30, min_D = 10;
    if (argc > 1) max_D = atoi(argv[1]);
    if (argc > 2) min_D = atoi(argv[2]);

    fprintf(stderr, "D01 BFS Entropy (64GB edition): D=%d..%d\n", min_D, max_D);
    printf("D,N,p,q,H_peak,max_c0,peak_k,last_k,reliable,aborted\n");
    fflush(stdout);

    double sum_h[64] = {0}, sum_h2[64] = {0};
    int count_all[64] = {0}, count_rel[64] = {0};
    double sum_h_rel[64] = {0};

    for (int D = min_D; D <= max_D; D++) {
        int hbits = D + 2;
        if (hbits > 29) hbits = 29;  /* cap at 2^29 = 537M entries ≈ 22.5GB */
        if (hbits < 16) hbits = 16;

        int samples;
        if (D <= 16) samples = 30;
        else if (D <= 20) samples = 20;
        else if (D <= 24) samples = 15;
        else if (D <= 28) samples = 10;
        else samples = 5;

        uint64_t N_lo = 1ULL << (D - 1);
        uint64_t N_hi = (1ULL << D) - 1;
        uint64_t sq_lo = (uint64_t)sqrt((double)N_lo);
        uint64_t sq_hi = (uint64_t)sqrt((double)N_hi);
        if (sq_lo < 3) sq_lo = 3;

        int good = 0;
        fprintf(stderr, "\nD=%d (hbits=%d, target=%d samples):\n", D, hbits, samples);

        for (int att = 0; att < samples * 10 && good < samples; att++) {
            uint64_t p, q, N;
            /* Mix balanced and unbalanced factors */
            if (good % 3 < 2) {
                p = random_prime_range(sq_lo * 2 / 3, sq_hi);
            } else {
                uint64_t pm = sq_hi / 3;
                if (pm < 5) pm = 5;
                p = random_prime_range(3, pm);
            }
            if (!p) continue;
            uint64_t ql = (N_lo + p - 1) / p;
            uint64_t qh = N_hi / p;
            if (ql < p) ql = p;
            if (qh < ql) continue;
            q = random_prime_range(ql, qh);
            if (!q) continue;
            N = p * q;
            if ((64 - __builtin_clzll(N)) != D) continue;

            BfsResult r = bfs_entropy(N, hbits);

            printf("%d,%llu,%llu,%llu,%.4f,%d,%d,%d,%d,%d\n",
                D, (unsigned long long)N,
                (unsigned long long)(p < q ? p : q),
                (unsigned long long)(p < q ? q : p),
                r.h_peak, r.max_c0, r.peak_k, r.last_k,
                r.peak_reliable, r.aborted);
            fflush(stdout);

            sum_h[D] += r.h_peak;
            sum_h2[D] += r.h_peak * r.h_peak;
            count_all[D]++;
            if (r.peak_reliable && !r.aborted) {
                sum_h_rel[D] += r.h_peak;
                count_rel[D]++;
            }
            good++;

            fprintf(stderr, "  [%d/%d] N=%llu (%llu×%llu) H=%.2f c0=%d %s%s\n",
                    good, samples, (unsigned long long)N,
                    (unsigned long long)(p < q ? p : q),
                    (unsigned long long)(p < q ? q : p),
                    r.h_peak, r.max_c0,
                    r.peak_reliable ? "reliable" : "UNRELIABLE",
                    r.aborted ? " ABORTED" : "");
        }
    }

    /* Final summary and regression */
    printf("\n# SUMMARY\n");
    printf("# %3s | %5s(%4s) | %10s %7s | %6s %7s\n",
           "D", "n", "rel", "avg_Hpeak", "std", "D/2", "delta");

    double sx = 0, sy = 0, sxy = 0, sxx = 0;
    int nt = 0;
    /* Reliable-only regression */
    double rx = 0, ry = 0, rxy = 0, rxx = 0;
    int rn = 0;

    for (int D = min_D; D <= max_D; D++) {
        if (!count_all[D]) continue;
        double avg = sum_h[D] / count_all[D];
        double var = sum_h2[D] / count_all[D] - avg * avg;
        double std = var > 0 ? sqrt(var) : 0;
        double delta = avg - D / 2.0;

        printf("# %3d | %5d(%4d) | %10.4f %7.4f | %6.1f %+7.3f\n",
            D, count_all[D], count_rel[D], avg, std, D / 2.0, delta);

        sx += D; sy += avg; sxy += D * avg; sxx += D * (double)D;
        nt++;

        if (count_rel[D] >= 3) {
            double avg_r = sum_h_rel[D] / count_rel[D];
            rx += D; ry += avg_r; rxy += D * avg_r; rxx += D * (double)D;
            rn++;
        }
    }

    if (nt > 2) {
        double mx = sx / nt, my = sy / nt;
        double a = (sxy - nt * mx * my) / (sxx - nt * mx * mx);
        double b = my - a * mx;
        double ss_r = 0, ss_t = 0;
        for (int D = min_D; D <= max_D; D++) {
            if (!count_all[D]) continue;
            double avg = sum_h[D] / count_all[D];
            ss_r += (avg - (a * D + b)) * (avg - (a * D + b));
            ss_t += (avg - my) * (avg - my);
        }
        printf("\n# REGRESSION (all): H_peak = %.4f * D + (%.4f), R^2=%.4f\n",
               a, b, 1 - ss_r / ss_t);
        printf("# SLOPE = %.4f (0.5000 = sqrt(N) barrier)\n", a);
    }
    if (rn > 2) {
        double rmx = rx / rn, rmy = ry / rn;
        double ra = (rxy - rn * rmx * rmy) / (rxx - rn * rmx * rmx);
        double rb = rmy - ra * rmx;
        printf("# REGRESSION (reliable only): H_peak = %.4f * D + (%.4f)\n", ra, rb);
    }

    for (int i = 0; i < 2; i++) { free(ht[i]); free(occ[i]); }
    return 0;
}
