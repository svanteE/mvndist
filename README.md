# mvndist

Distance Measures Between Multivariate Normal Distributions

## Overview

`mvndist` provides a comprehensive set of distance and divergence measures between multivariate normal distributions, including:

- **Wasserstein-2 (Bures)** - Optimal transport distance with closed form (FAST, recommended for most uses)
- **Fisher-Rao** - Information geometry metric via optimization (exact but slower)
- **Hellinger** - Bounded distance [0,1]
- **Bhattacharyya** - Related to classification error
- **Kullback-Leibler** - Information-theoretic divergence
- **Affine-invariant** - Riemannian metric on SPD matrices

## Key Distinction: Fisher-Rao vs Wasserstein-2

Both measure "distance" between distributions, but they capture different concepts:

### Fisher-Rao Distance (Information Geometry)
- **What it measures**: Geodesic distance on the statistical manifold
- **Computation**: Uses Eriksen (1987) method solving the geodesic equations
- **Speed**: Slower, especially for high dimensions
- **When to use**: Theoretical statistical work, information-theoretic analysis
- **Implementation**: Optimization of f(y) = tr(C(y)²) = 0 from Eriksen (1987)

### Wasserstein-2 Distance (Optimal Transport)
- **What it measures**: Cost of optimally transporting probability mass
- **Computation**: Closed-form solution via matrix square roots
- **Speed**: Very fast, even in high dimensions
- **When to use**: Machine learning, LLM applications, practical work
- **Implementation**: Direct formula, no iteration needed

**Rule of thumb**: Use Wasserstein-2 unless you specifically need the Fisher information metric for theoretical reasons.

## Installation

You can install the development version from source:

```r
# Install dependencies
install.packages(c("expm", "MASS"))

# Install from source
install.packages("path/to/mvndist", repos = NULL, type = "source")

# Or using devtools
devtools::install_local("path/to/mvndist")
```

Note: The package no longer requires `deSolve` - Fisher-Rao distance now uses the optimized Eriksen (1987) method based on Skovgaard's (1984) Riemannian geometry framework.

## Quick Start

```r
library(mvndist)

# Define two bivariate normal distributions
mu1 <- c(0, 0)
Sigma1 <- diag(2)

mu2 <- c(2, 1)
Sigma2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)

# Compute individual distances
wasserstein2_distance(mu1, Sigma1, mu2, Sigma2)
hellinger_distance(mu1, Sigma1, mu2, Sigma2)
bhattacharyya_distance(mu1, Sigma1, mu2, Sigma2)

# Compare all distances at once
compare_distances(mu1, Sigma1, mu2, Sigma2)
```

## Key Functions

### Fast, Closed-Form Distances

- `wasserstein2_distance()` - Optimal transport (recommended for most applications)
- `hellinger_distance()` - Bounded similarity measure
- `bhattacharyya_distance()` - Pattern recognition, classification
- `kl_divergence()` - Information theory (asymmetric)
- `symmetrized_kl()` - Symmetric version of KL
- `affine_invariant_distance()` - Riemannian metric on covariances

### Numerical Methods

- `fisher_rao_distance()` - Information geometry metric (slow, may be unstable)
- `fisher_rao_geodesic_numerical()` - Compute geodesic path

### Utilities

- `compare_distances()` - Compute all distances for easy comparison

## Which Distance Should I Use?

| Application | Recommended Distance | Why |
|-------------|---------------------|-----|
| **Machine Learning** | Wasserstein-2 | Fast closed form, meaningful interpolation |
| **LLM embeddings** | Wasserstein-2 | Efficient for high-dimensional spaces |
| **Model merging/interpolation** | Wasserstein-2 | Optimal transport interpretation |
| **Classification** | Bhattacharyya | Related to Bayes error |
| **Bounded similarity** | Hellinger | Normalized [0,1] |
| **Information theory** | Fisher-Rao | Theoretically optimal for statistical inference |
| **Fast approximation** | Fisher-Rao | Only when exact information metric is needed |
| **Variational inference** | KL divergence | Natural for VAEs |
| **Covariance smoothing** | Affine-invariant | SPD matrix operations |

### Performance Comparison

```r
# Example: Compare computation time
system.time(wasserstein2_distance(mu1, Sigma1, mu2, Sigma2))
# Typical: < 0.01 seconds

system.time(fisher_rao_distance(mu1, Sigma1, mu2, Sigma2))
# Typical: 0.1-1 seconds (requires optimization)
```

**Bottom line**: Wasserstein-2 is 10-100x faster than Fisher-Rao and sufficient for most practical applications. Use Fisher-Rao only when you specifically need the information geometry perspective.

## Examples

### LLM Applications

```r
# Compare embedding distributions
embedding1_mu <- colMeans(embeddings_model1)
embedding1_Sigma <- cov(embeddings_model1)

embedding2_mu <- colMeans(embeddings_model2)
embedding2_Sigma <- cov(embeddings_model2)

# Wasserstein distance for model comparison
wasserstein2_distance(embedding1_mu, embedding1_Sigma, 
                      embedding2_mu, embedding2_Sigma)
```

### Out-of-Distribution Detection

```r
# Train distribution
train_mu <- colMeans(train_features)
train_Sigma <- cov(train_features)

# Test sample
test_mu <- colMeans(test_features)
test_Sigma <- cov(test_features)

# Detect distribution shift
d <- bhattacharyya_distance(train_mu, train_Sigma, test_mu, test_Sigma)
is_ood <- d > threshold
```

## References

### Key Papers

- **Wasserstein-2**: Takatsu, A. (2011). Wasserstein geometry of Gaussian measures. *Osaka Journal of Mathematics*, 48(4), 1005-1026.

- **Fisher-Rao**: Skovgaard, L. T. (1984). A Riemannian geometry of the multivariate normal model. *Scandinavian Journal of Statistics*, 11(4), 211-223.

- **Affine-invariant**: Pennec, X., Fillard, P., & Ayache, N. (2006). A Riemannian framework for tensor computing. *International Journal of Computer Vision*, 66(1), 41-66.

- **Information Geometry**: Amari, S. (2016). *Information Geometry and Its Applications*. Springer.

### Included Documentation

A detailed paper on the Fisher-Rao geodesic implementation is included with the package:
```r
# View the paper location
system.file("fisher_rao_geodesic.pdf", package = "mvndist")
```

## License

MIT License - see LICENSE file for details

## Author

Poul Svante Eriksen - Aalborg University

## Contributing

Contributions are welcome! Please open an issue or pull request on GitHub.
