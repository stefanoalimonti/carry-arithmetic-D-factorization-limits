#!/usr/bin/env python3
"""
D02_hensel_lifting.py — Hensel lifting demonstration for integer factorization.

Shows that p-adic lifting from a small modulus base m produces a single
deterministic path mod m, m^2, m^4, ..., but hitting the true factor
requires m ~ sqrt(N), giving no speedup over trial division.

Demonstrates the core mechanism behind the phantom halo: in base 2,
Hensel lifting is degenerate (every bit choice produces a valid lift),
generating 2^(D/2) modular solutions that mask the true factors.
"""
import math
import random


def ext_gcd(a, b):
    if a == 0:
        return b, 0, 1
    g, y, x = ext_gcd(b % a, a)
    return g, x - (b // a) * y, y


def mod_inv(a, m):
    g, x, _ = ext_gcd(a, m)
    if g != 1:
        raise ValueError(f"modular inverse does not exist for {a} mod {m}")
    return x % m


def is_prime(n):
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def hensel_lift(N, m, max_depth=10):
    """
    Perform Hensel lifting from base modulus m.

    For each pair (p0, q0) with p0*q0 = N mod m, lifts the solution
    to mod m^2, m^4, m^8, ... checking for exact factorization.
    """
    base_sols = []
    for p0 in range(1, m):
        if math.gcd(p0, m) == 1:
            q0 = (N * mod_inv(p0, m)) % m
            base_sols.append((p0, q0))

    results = []

    for p0, q0 in base_sols:
        pk, qk = p0, q0
        current_mod = m
        path = [(current_mod, pk, qk)]
        success = False

        for step in range(1, max_depth + 1):
            error = (N - pk * qk) // current_mod

            if error == 0 and pk > 1 and qk > 1:
                success = True
                break

            try:
                qk_inv = mod_inv(qk, current_mod)
                t = (error * qk_inv) % current_mod

                pk_next = pk + t * current_mod
                if pk_next > 0 and N % pk_next == 0:
                    qk_next = N // pk_next
                else:
                    qk_next = (N * mod_inv(pk_next, current_mod**2)) % (current_mod**2)

                current_mod = current_mod ** 2
                pk, qk = pk_next, qk_next
                path.append((current_mod, pk, qk))

            except (ValueError, ZeroDivisionError):
                break

        results.append({
            'base': (p0, q0),
            'path': path,
            'success': success,
            'final_p': pk if success else None,
            'final_q': qk if success else None,
        })

    return results


if __name__ == "__main__":
    print("=" * 70)
    print("  Hensel Lifting Demonstration")
    print("=" * 70)

    test_cases = [
        (11, 13, 143),
        (13, 17, 221),
        (23, 29, 667),
        (41, 43, 1763),
    ]

    for p, q, N in test_cases:
        print(f"\nN = {N} = {p} x {q}")

        for m in [3, 5, 7]:
            print(f"  Base modulus m = {m}")
            results = hensel_lift(N, m, max_depth=5)

            for r in results:
                p0, q0 = r['base']
                is_correct = (p % m == p0) or (q % m == p0)
                marker = " <-- correct seed" if is_correct else ""

                print(f"    Seed (mod {m}): p0={p0}, q0={q0}{marker}")

                for mod_k, pk, qk in r['path']:
                    match_p = pk % mod_k == p % mod_k
                    match_q = pk % mod_k == q % mod_k
                    target = "p" if match_p else ("q" if match_q else "?")
                    print(f"      mod {mod_k:<6d}: pk={pk:<5d} (tracks {target})")

                if r['success']:
                    print(f"      FACTORED: {r['final_p']} x {r['final_q']} = {N}")
                else:
                    print(f"      Lift stalled (no exact factorization found)")
