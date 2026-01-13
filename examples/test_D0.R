# Test D_0 with specific De and de values

library(mvndist)
library(expm)

# Define the helper function D_0 from equation (7)
D0_fisher_rao_direct <- function(de, De) {
  # Equation (7): D_0 formula
  Dei <- solve(De)
  mu <- Dei %*% de
  tau <- sum(mu * de)
  phi <- -de %*% t(mu) / 2
  Ga <- Dei + (1 + tau/4) * mu %*% t(mu)
  A <- (De + Ga - phi - t(phi)) / 2
  
  # Matrix logarithm
  G <- logm(A)
  sqrt(sum(G^2) / 2)
}

# Fisher-Rao distance from canonical (De, de) space
fisher_rao_from_canonical <- function(de, De) {
  # The A matrix used in Fisher-Rao
  # Construct T and compute A = log(T'T)
  p <- length(de)
  
  # For this simple test, use the standard formula
  # A is from the optimization problem, but for canonical space we can approximate
  A_FR <- matrix(0, 2*p + 1, 2*p + 1)
  
  # This is approximate - the real computation requires optimization
  # Just show the distance formula
  Dei <- solve(De)
  mu <- Dei %*% de
  A_simple <- cbind(rbind(-diag(p), -t(mu)), rep(0, p+1))
  
  sqrt(sum(A_simple^2) / 2)
}

# Test values
De <- matrix(c(
  1.0, 0.1, 0.1,
  0.1, 1.0, 0.1,
  0.1, 0.1, 1.0
), nrow = 3, byrow = TRUE)

de <- matrix(c(1, 1, 1), ncol = 1)

cat("Test De matrix:\n")
print(De)
cat("\nTest de vector:", as.vector(de), "\n\n")

# Compute D_0
d_0 <- D0_fisher_rao_direct(de, De)
cat("D_0 (from equation 7):", round(d_0, 6), "\n")

# Show intermediate computations
cat("\n--- Intermediate computations ---\n")
Dei <- solve(De)
mu <- Dei %*% de
tau <- sum(mu * de)

cat("Dei (inverse of De):\n")
print(round(Dei, 4))
cat("\nmu = Dei %*% de:", as.vector(round(mu, 4)), "\n")
cat("tau = sum(mu * de):", round(tau, 6), "\n")

phi <- -de %*% t(mu) / 2
cat("\nphi matrix:\n")
print(round(phi, 4))

Ga <- Dei + (1 + tau/4) * mu %*% t(mu)
cat("\nGa matrix:\n")
print(round(Ga, 4))

A <- (De + Ga - phi - t(phi)) / 2
cat("\nA matrix:\n")
print(round(A, 4))

G <- logm(A)
cat("\nG = logm(A):\n")
print(round(G, 4))

cat("\nFrobenius norm squared of G:", round(sum(G^2), 6), "\n")
cat("D_0 = sqrt(sum(G^2)/2) =", round(sqrt(sum(G^2)/2), 6), "\n")

# Also compute Fisher-Rao distance formula
cat("\n--- Fisher-Rao distance formula ---\n")
cat("Fisher-Rao distance = sqrt(tr(A^2)/2) where A = log(T'T)\n")
cat("For D_0, distance =", round(d_0, 6), "\n")

# ------------------------------------------------------------
# Eigen-aligned case: de is an eigenvector of De
# In this case, D0 equals the Fisher–Rao distance
# Construction: set Sigma1 = I, mu1 = 0, choose symmetric PD De and
# set de = alpha * v where v is an eigenvector of De. Then
# Sigma2 = solve(De), mu2 = solve(De) %*% de realizes (De, de) as the
# canonical parameters of (mu2, Sigma2) relative to (mu1, Sigma1).
# ------------------------------------------------------------

cat("\n=== Eigen-aligned verification (D0 == Fisher–Rao) ===\n")

# Pick a simple diagonal De with known eigenvectors (standard basis)
De_eig <- diag(c(2, 3, 5))
v1 <- matrix(c(1, 0, 0), ncol = 1)     # eigenvector for eigenvalue 2
alpha <- 1
de_eig <- alpha * v1                    # de aligned with eigenvector of De

# Build (mu1, Sigma1) and (mu2, Sigma2) that yield (De, de) canonically
mu1_eig <- c(0, 0, 0)
Sigma1_eig <- diag(3)
Sigma2_eig <- solve(De_eig)
mu2_eig <- as.vector(solve(De_eig) %*% de_eig)

# Compute D0 directly from (De, de)
d0_eig <- D0_fisher_rao_direct(de_eig, De_eig)

# Compute Fisher–Rao via mvndist
fr_eig <- fisher_rao_distance(mu1_eig, Sigma1_eig, mu2_eig, Sigma2_eig)

cat("De (eigen-aligned):\n"); print(De_eig)
cat("de (aligned with eigenvector):", as.vector(de_eig), "\n")
cat("mu2 (constructed):", round(mu2_eig, 6), "\n")
cat("Sigma2 (constructed):\n"); print(Sigma2_eig)

cat("\nD0 (eq. 7):", round(d0_eig, 6), "\n")
cat("Fisher–Rao:", round(fr_eig, 6), "\n")
cat("Abs. diff:", format(abs(d0_eig - fr_eig), scientific = TRUE), "\n")
cat("Matches within 1e-8:", abs(d0_eig - fr_eig) < 1e-8, "\n")
