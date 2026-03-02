# carry-arithmetic-D-factorization-limits

**The Carry-Zero Entropy Bound: Structural Limits of Bitwise Factorization**

*Author: Stefano Alimonti* · [ORCID 0009-0009-1183-1698](https://orcid.org/0009-0009-1183-1698)

## Main Result

The carry-zero BFS tree of binary multiplication decomposes cleanly:

$$S_0(k) = \tau_{\text{odd}}(M_k) + \Phi(k)$$

- **Core** $\tau_{\text{odd}}$: real odd-divisor pairs, sub-polynomial (Wigert-Ramanujan).
- **Halo** $\Phi$: phantom modular aliases, exponential ($\alpha \in [0.38, 0.50]$).

Phantoms dominate at > 93% for $D \geq 12$, arising from degenerate Hensel lifting in $\mathbf{Z}_2$.

## Status

Ready for submission. Two theorems proved; phantom scaling empirical ($D \leq 28$).

## Repository Structure

```
paper/carry_entropy_bound.md        The paper
experiments/
  D01_bfs_entropy.c                 C hash-table BFS for D <= 28
  D02_hensel_lifting.py             Hensel lifting demonstration
  D03_phantom_decomposition.py      Phantom vs divisor decomposition
  D04_phantom_stats.py              Per-depth phantom statistics
```

## Reproduction

```bash
# C experiment (requires gcc/clang)
cc -O3 -o bfs_entropy experiments/D01_bfs_entropy.c -lm
./bfs_entropy 28 10

# Python experiments
pip install numpy
python experiments/D02_hensel_lifting.py
python experiments/D03_phantom_decomposition.py
python experiments/D04_phantom_stats.py
```

## Dependencies

- Python >= 3.8, NumPy
- C compiler (gcc or clang) for `D01_bfs_entropy.c`

## Companion Papers

| Label | Title | Repository |
|-------|-------|------------|
| [A] | Spectral Theory of Carries in Positional Multiplication | [`carry-arithmetic-A-spectral-theory`](https://github.com/stefanoalimonti/carry-arithmetic-A-spectral-theory) |
| [C] | Eigenvalue Statistics of Carry Companion Matrices: Markov-Driven GOE↔GUE Transition | [`carry-arithmetic-C-matrix-statistics`](https://github.com/stefanoalimonti/carry-arithmetic-C-matrix-statistics) |
| [E] | The Trace Anomaly of Binary Multiplication | [`carry-arithmetic-E-trace-anomaly`](https://github.com/stefanoalimonti/carry-arithmetic-E-trace-anomaly) |
| [G] | The Angular Uniqueness of Base 2 | [`carry-arithmetic-G-angular-uniqueness`](https://github.com/stefanoalimonti/carry-arithmetic-G-angular-uniqueness) |
| [H] | Carry Polynomials and the Partial Euler Product | [`carry-arithmetic-H-euler-control`](https://github.com/stefanoalimonti/carry-arithmetic-H-euler-control) |
| [P1] | Pi from Pure Arithmetic | [`carry-arithmetic-P1-pi-spectral`](https://github.com/stefanoalimonti/carry-arithmetic-P1-pi-spectral) |
| [P2] | The Sector Ratio in Binary Multiplication | [`carry-arithmetic-P2-sector-ratio`](https://github.com/stefanoalimonti/carry-arithmetic-P2-sector-ratio) |

### Citation

```bibtex
@article{alimonti2026entropy_bound,
  author  = {Alimonti, Stefano},
  title   = {The Carry-Zero Entropy Bound: Structural Limits of Bitwise Factorization},
  year    = {2026},
  note    = {Preprint},
  url     = {https://github.com/stefanoalimonti/carry-arithmetic-D-factorization-limits}
}
```

## License

Paper: CC BY 4.0. Code: MIT License.
