# mvndist 0.1.0

## Initial Release

### Features

* **Wasserstein-2 distance**: Fast closed-form implementation using matrix square roots
* **Fisher-Rao distance**: Optimized implementation using the exact method from Skovgaard (1984)
  - Solves optimization problem f(y) = tr(C(y)²) = 0
  - Exact closed-form solution for univariate case
  - Numerical optimization for multivariate case
* **Hellinger distance**: Bounded similarity measure [0,1]
* **Bhattacharyya distance**: Related to classification error bounds
* **KL divergence**: Information-theoretic divergence (asymmetric)
* **Symmetrized KL**: Symmetric version of KL divergence
* **Affine-invariant distance**: Riemannian metric on SPD matrices
* **compare_distances()**: Utility function to compute all distances at once

### Documentation

* Comprehensive README with examples and guidance on which distance to use
* Full roxygen2 documentation for all exported functions
* Unit tests for all distance measures

### Performance

* All distances except Fisher-Rao have closed-form solutions (O(p³) complexity)
* Fisher-Rao uses efficient optimization with good convergence properties
* Suitable for embeddings up to 100+ dimensions

### Dependencies

* expm: For matrix exponential and logarithm operations
* MASS: For multivariate normal utilities

### Notes

* Fisher-Rao and Wasserstein-2 measure fundamentally different geometric properties
* Wasserstein-2 recommended for most ML/LLM applications due to speed and stability
* Fisher-Rao recommended for theoretical statistical work requiring information geometry
