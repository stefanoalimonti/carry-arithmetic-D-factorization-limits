#!/usr/bin/env python3
"""
D03_phantom_decomposition.py — Decomposition S_0(k) = tau_odd(M_k) + Phi(k).

Computes the carry-zero envelope S_0(k) via BFS, then separates the
sub-polynomial divisor core tau_odd from the exponential phantom halo Phi.
Measures the growth exponent alpha_Phi = log2(Phi_peak) / D across bit-lengths.
"""
import math
import random

import numpy as np

random.seed(42)
np.random.seed(42)


def get_bits(n, D):
    return [(n >> i) & 1 for i in range(D)]


def fast_bfs_entropy(N, max_pos=None):
    D = N.bit_length()
    if max_pos is None:
        max_pos = D + 2
    N_bits = [(N >> i) & 1 for i in range(max_pos + 1)]

    states = [(1, 1, 1)]
    curve = []

    for k in range(1, max_pos):
        next_states = []
        target = N_bits[k]

        c0 = 0
        Mk = N % (1 << (k + 1))

        for pp, qq, P_prev in states:
            P_old_k_bit = ((pp * qq) >> k) & 1
            xor_ab = target ^ P_old_k_bit

            valid_pairs = [(0, 0), (1, 1)] if xor_ab == 0 else [(0, 1), (1, 0)]

            for a, b in valid_pairs:
                np_ = pp | (a << k)
                nq_ = qq | (b << k)

                ck = 0
                for i in range(k + 1):
                    ai = (np_ >> i) & 1
                    bki = (nq_ >> (k - i)) & 1
                    ck += ai * bki

                Pk = P_prev + (ck << k)
                next_states.append((np_, nq_, Pk))

                # Pk is the carry-free convolution sum Σ conv_j·2^j.
                # For BFS states, (np_·nq_) mod 2^{k+1} = Mk always holds.
                # Pk == Mk iff there is no carry propagation at columns 0..k,
                # i.e., this is equivalent to the carry-zero condition c_{k+1}=0
                # of Definition 2 in the paper.
                if Pk == Mk:
                    c0 += 1

        curve.append((k, len(next_states), c0))
        states = next_states

        if len(states) > 1_500_000:
            break

    return curve


def generate_prime(bits):
    while True:
        p = random.getrandbits(bits)
        p |= (1 << (bits - 1)) | 1
        d = 3
        is_p = True
        while d * d <= p:
            if p % d == 0:
                is_p = False
                break
            d += 2
        if is_p:
            return p


if __name__ == "__main__":
    print("=" * 80)
    print("  S_0(k) = tau_odd(M_k) + Phi(k): Phantom Growth Exponent")
    print("=" * 80)
    print()

    for D in range(10, 24, 2):
        phantoms_peaks = []
        tau_peaks = []

        for _ in range(5):
            dp = D // 2
            dq = D - dp
            p = generate_prime(dp)
            q = generate_prime(dq)
            N = p * q

            curve = fast_bfs_entropy(N, max_pos=D + 2)

            max_p = 0
            max_tau = 0

            for k, stot, c0 in curve:
                Mk = N % (1 << (k + 1))
                bound = 1 << (k + 1)

                tau = 0
                if Mk > 0:
                    d = 1
                    while d * d <= Mk:
                        if Mk % d == 0:
                            e = Mk // d
                            if d % 2 == 1 and e % 2 == 1 and d < bound and e < bound:
                                tau += 1
                                if d != e:
                                    tau += 1
                        d += 1

                phantoms = c0 - tau
                if phantoms > max_p:
                    max_p = phantoms
                if tau > max_tau:
                    max_tau = tau

            phantoms_peaks.append(max_p)
            tau_peaks.append(max_tau)

        avg_p = sum(phantoms_peaks) / len(phantoms_peaks)
        avg_t = sum(tau_peaks) / len(tau_peaks)

        alpha_p = math.log2(avg_p) / D if avg_p > 0 else 0
        alpha_t = math.log2(avg_t) / D if avg_t > 0 else 0
        pct = avg_p / (avg_p + avg_t) * 100 if (avg_p + avg_t) > 0 else 0

        print(
            f"D={D:2d} | <Phi> = {avg_p:8.1f} (alpha = {alpha_p:.4f}) | "
            f"<tau> = {avg_t:6.1f} (alpha = {alpha_t:.4f}) | %Phantom = {pct:.1f}%"
        )
