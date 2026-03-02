#!/usr/bin/env python3
"""
D04_phantom_stats.py — Detailed per-depth phantom statistics.

For specific semiprimes, traces S_0(k), tau_odd(M_k), and Phi(k) at
every bit position k, showing the depth profile of the phantom halo.
"""
import math

import numpy as np


def analyze_phantoms(N):
    D = N.bit_length()
    N_bits = [(N >> i) & 1 for i in range(D + 2)]

    states = {(1, 1, 0)}

    for k in range(1, D + 2):
        next_states = set()
        target = N_bits[k]

        M_k = N % (1 << (k + 1))

        tau = 0
        if M_k > 0:
            d = 1
            while d * d <= M_k:
                if M_k % d == 0:
                    e = M_k // d
                    if d % 2 == 1 and e % 2 == 1 and d < (1 << (k + 1)) and e < (1 << (k + 1)):
                        tau += 1
                        if d != e:
                            tau += 1
                d += 1

        c0_states = []
        for pp, qq, cin in states:
            for a in (0, 1):
                for b in (0, 1):
                    np_ = pp | (a << k)
                    nq_ = qq | (b << k)

                    conv = 0
                    for i in range(k + 1):
                        conv += ((np_ >> i) & 1) * ((nq_ >> (k - i)) & 1)

                    total = conv + cin
                    if (total & 1) == target:
                        cout = total >> 1
                        next_states.add((np_, nq_, cout))
                        if cout == 0:
                            c0_states.append((np_, nq_))

        S0 = len(c0_states)
        phantoms = S0 - tau

        print(
            f"  k={k:2d} | S_0={S0:5d}  tau={tau:5d}  "
            f"Phi={phantoms:5d} | ratio={S0 / tau if tau > 0 else 0:5.1f}x"
        )

        states = next_states
        if len(states) > 500_000:
            break


if __name__ == "__main__":
    print("=" * 70)
    print("  Per-depth phantom statistics")
    print("=" * 70)

    print("\nN = 10403 = 101 x 103 (D = 14)")
    analyze_phantoms(10403)

    print("\nN = 525274 = 601 x 874 ... checking")
    N = 525274
    print(f"  bit-length = {N.bit_length()}")
    analyze_phantoms(N)
