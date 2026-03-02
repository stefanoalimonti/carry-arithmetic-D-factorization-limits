# Experiments — Paper D

| Script | Description | Language | Referenced in |
|--------|-------------|----------|---------------|
| D01_bfs_entropy.c | Hash-table BFS computing S₀(k) for D ≤ 28 | C | §5.1 (Table 2) |
| D02_hensel_lifting.py | Hensel lifting from small moduli, barrier demonstration | Python | §6 |
| D03_phantom_decomposition.py | S₀ = τ + Φ decomposition and growth exponent α | Python | Thm. 1 (§3), §5.1 (Table 1) |
| D04_phantom_stats.py | Per-depth phantom profile for specific semiprimes | Python | §5.1 |

## Requirements

- Python >= 3.8, NumPy
- C compiler (gcc/clang) for D01
