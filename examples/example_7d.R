# Nontrivial example: 7-dimensional multivariate normals
# Scenario: Compare two neural network embedding distributions

library(mvndist)
set.seed(42)

# ============================================================
# Simulate two 7-dimensional embedding distributions
# (e.g., from different fine-tuned language models)
# ============================================================

# Distribution 1: Base model embeddings
n1 <- 500
mu1 <- c(0, 0.5, -0.3, 0.2, 0, -0.1, 0.4)
Sigma1 <- matrix(c(
  1.0, 0.3, 0.1, -0.2, 0.05, 0.1, 0.0,
  0.3, 0.9, 0.2, 0.1, -0.1, 0.05, 0.15,
  0.1, 0.2, 1.2, 0.0, 0.2, 0.1, -0.05,
  -0.2, 0.1, 0.0, 0.8, 0.3, 0.0, 0.2,
  0.05, -0.1, 0.2, 0.3, 1.1, 0.15, 0.1,
  0.1, 0.05, 0.1, 0.0, 0.15, 0.95, 0.05,
  0.0, 0.15, -0.05, 0.2, 0.1, 0.05, 1.3
), nrow = 7, byrow = TRUE)

# Distribution 2: Fine-tuned model embeddings (shifted and rotated)
mu2 <- c(0.5, 0.8, -0.5, 0.5, 0.3, 0.2, 0.6)
Sigma2 <- matrix(c(
  1.1, 0.25, 0.15, -0.1, 0.1, 0.15, 0.05,
  0.25, 1.0, 0.25, 0.2, -0.05, 0.1, 0.2,
  0.15, 0.25, 1.3, 0.05, 0.25, 0.15, 0.0,
  -0.1, 0.2, 0.05, 0.9, 0.35, 0.05, 0.25,
  0.1, -0.05, 0.25, 0.35, 1.2, 0.2, 0.15,
  0.15, 0.1, 0.15, 0.05, 0.2, 1.05, 0.1,
  0.05, 0.2, 0.0, 0.25, 0.15, 0.1, 1.4
), nrow = 7, byrow = TRUE)

# ============================================================
# Compare all distances
# ============================================================

cat("7-Dimensional Multivariate Normal Distributions\n")
cat("==============================================\n")
cat("Model 1 (Base): μ1 =", mu1, "\n")
cat("Model 2 (Fine-tuned): μ2 =", mu2, "\n")
cat("\nΣ1 (covariance matrix 1):\n")
print(Sigma1)
cat("\nΣ2 (covariance matrix 2):\n")
print(Sigma2)
cat("\n")

# Compute all distances
distances <- compare_distances(mu1, Sigma1, mu2, Sigma2)

cat("Distance Comparison:\n")
print(distances)
cat("\n")

# ============================================================
# Individual detailed results
# ============================================================

# Helper function: D_0 from equation (7) of the paper
D0_fisher_rao <- function(mu1, Sigma1, mu2, Sigma2) {
  p <- length(mu1)
  
  # Transform to canonical parameters relative to (0, I)
  U <- chol(solve(Sigma1))
  U_inv <- solve(U)
  
  mu_transformed <- U %*% (mu2 - mu1)
  Sigma_transformed <- U %*% Sigma2 %*% t(U)
  
  De <- solve(Sigma_transformed)
  de <- De %*% mu_transformed
  
  # Equation (7): D_0 formula
  Dei <- solve(De)
  mu <- Dei %*% de
  tau <- sum(mu * de)
  phi <- -de %*% t(mu) / 2
  Ga <- Dei + (1 + tau/4) * mu %*% t(mu)
  A <- (De + Ga - phi - t(phi)) / 2
  
  # Matrix logarithm
  G <- expm::logm(A)
  sqrt(sum(G^2) / 2)
}

cat("Individual Distance Measures:\n")
cat(strrep("-", 40), "\n", sep = "")

# D_0 from equation (7)
d_0 <- D0_fisher_rao(mu1, Sigma1, mu2, Sigma2)
cat("D_0 (Equation 7, Fisher-Rao theory):", round(d_0, 4), "\n\n")

w2 <- wasserstein2_distance(mu1, Sigma1, mu2, Sigma2)
cat("Wasserstein-2 (optimal transport):", round(w2, 4), "\n")

hell <- hellinger_distance(mu1, Sigma1, mu2, Sigma2)
cat("Hellinger distance [0-1]:", round(hell, 4), "\n")

bhat <- bhattacharyya_distance(mu1, Sigma1, mu2, Sigma2)
cat("Bhattacharyya distance:", round(bhat, 4), "\n")

kl <- kl_divergence(mu1, Sigma1, mu2, Sigma2)
cat("KL divergence (asymmetric):", round(kl, 4), "\n")

skl <- symmetrized_kl(mu1, Sigma1, mu2, Sigma2)
cat("Symmetrized KL:", round(skl, 4), "\n")

ai <- affine_invariant_distance(mu1, Sigma1, mu2, Sigma2)
cat("Affine-invariant (Riemannian):", round(ai, 4), "\n")

fr <- fisher_rao_distance(mu1, Sigma1, mu2, Sigma2)
cat("Fisher-Rao (information geometry):", round(fr, 4), "\n")

# Comparison
cat("\n")
cat("Comparison to Fisher-Rao:\n")
cat(strrep("-", 40), "\n", sep = "")
cat("D_0 / Fisher-Rao ratio:", round(d_0 / fr, 4), "\n")
cat("Wasserstein-2 / Fisher-Rao ratio:", round(w2 / fr, 4), "\n")
cat("Hellinger / Fisher-Rao ratio:", round(hell / fr, 4), "\n")

# Get the A matrix used in Fisher-Rao computations
cat("\nFisher-Rao A Matrix:\n")
A <- fisher_rao_A_matrix(mu1, Sigma1, mu2, Sigma2)

if (!is.null(A)) {
  cat("A matrix dimensions:", nrow(A), "×", ncol(A), "\n")
  cat("Frobenius norm of A:", round(sqrt(sum(A^2)), 4), "\n")
  cat("Distance from A: sqrt(tr(A²))/2 =", round(sqrt(sum(A^2)) / 2, 4), "\n\n")
  
  cat("A matrix (full):\n")
  print(round(A, 4))
} else {
  cat("Failed to compute A matrix\n")
}

cat("\n")
cat("Interpretation:\n")
cat("- Wasserstein-2 ≈", round(w2, 4), ": Fast metric, optimal transport cost\n")
cat("- Hellinger ≈", round(hell, 4), ": Bounded [0,1], statistical distance\n")
cat("- Fisher-Rao ≈", round(fr, 4), ": Information geometry, slower computation\n")
