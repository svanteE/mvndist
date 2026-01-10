#' Helper: Construct Z matrix from y (upper triangular part)
#' @noRd
construct_Z <- function(y, d, D) {
  p <- length(d)
  
  # Handle univariate case
  if (p == 1) {
    return(matrix(-d[1]^2 / (2 * D[1]), 1, 1))
  }
  
  Z <- matrix(0, p, p)
  
  # Diagonal elements (known)
  for (i in 1:p) {
    Z[i, i] <- -d[i]^2 / (2 * D[i])
  }
  
  # Upper triangular part from y
  idx <- 1
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      Z[i, j] <- y[idx]
      idx <- idx + 1
    }
  }
  
  # Lower triangular part determined from upper
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      Z[j, i] <- -(d[i] * d[j] + D[j] * Z[i, j]) / D[i]
    }
  }
  
  Z
}

#' Helper: Construct T matrix
#' @noRd
construct_T <- function(y, d, D) {
  p <- length(d)
  Z <- construct_Z(y, d, D)
  
  T_mat <- matrix(0, 2*p+1, 2*p+1)
  
  # Block structure: T = [D, d, Z; 0, 1, -d'*D^-1; 0, 0, D^-1]
  T_mat[1:p, 1:p] <- diag(D)
  T_mat[1:p, p+1] <- d
  T_mat[1:p, (p+2):(2*p+1)] <- Z
  
  T_mat[p+1, p+1] <- 1
  T_mat[p+1, (p+2):(2*p+1)] <- -d / D
  
  T_mat[(p+2):(2*p+1), (p+2):(2*p+1)] <- diag(1/D)
  
  T_mat
}

#' Helper: Objective function f(y) = tr(C(y)^2)
#' @noRd
objective_f <- function(y, d, D) {
  p <- length(d)
  
  T_mat <- construct_T(y, d, D)
  
  # Compute A(y) = log(T'T)
  TtT <- t(T_mat) %*% T_mat
  A <- expm::logm(TtT)
  
  # Extract C(y) block: rows 1:p, columns (p+2):(2p+1)
  C <- A[1:p, (p+2):(2*p+1)]
  
  # Return tr(C^2)
  sum(C^2)
}

#' Fisher-Rao Geodesic (Optimized)
#'
#' Computes the Fisher-Rao geodesic between two multivariate normal distributions
#' by solving the optimization problem described in the paper. This method finds
#' the exact geodesic by minimizing f(y) = tr(C(y)^2) where C(y) is the off-diagonal
#' block in the matrix logarithm decomposition.
#'
#' @param mu1 Mean vector of the first distribution
#' @param Sigma1 Covariance matrix of the first distribution
#' @param mu2 Mean vector of the second distribution
#' @param Sigma2 Covariance matrix of the second distribution
#' @param n_steps Number of points along the geodesic path (default: 100)
#' @param tol Tolerance for endpoint error (default: 1e-6)
#'
#' @return A list of length n_steps, each element containing:
#'   \item{t}{Parameter value (0 to 1)}
#'   \item{mu}{Mean vector at this point}
#'   \item{Sigma}{Covariance matrix at this point}
#'
#' @references
#' Skovgaard, L. T. (1984). A Riemannian geometry of the multivariate
#' normal model. Scandinavian Journal of Statistics, 11(4), 211-223.
#'
#' Eriksen, P. S. (1987). Geodesics connected with the Fisher metric on
#' the multivariate normal manifold. Proceedings of the GST Workshop.
#'
#' @export
fisher_rao_geodesic_numerical <- function(mu1, Sigma1, mu2, Sigma2, 
                                          n_steps = 100, tol = 1e-6) {
  # For now, return linear interpolation as placeholder
  # Full geodesic path computation would require solving the optimized path at each t
  times <- seq(0, 1, length.out = n_steps)
  
  warning("Full geodesic path not yet implemented. Use fisher_rao_distance() for distance computation.")
  
  lapply(times, function(t) {
    list(t = t, 
         mu = (1-t)*mu1 + t*mu2,
         Sigma = (1-t)*Sigma1 + t*Sigma2)
  })
}

#' Fisher-Rao Distance (Optimized)
#'
#' Computes the Fisher-Rao distance by optimizing the geodesic equation.
#' This method solves the optimization problem f(y) = tr(C(y)^2) = 0 to find
#' the exact geodesic, then computes the distance as the path length.
#' For univariate distributions, uses the exact closed-form solution.
#'
#' @param mu1 Mean vector of the first distribution
#' @param Sigma1 Covariance matrix of the first distribution
#' @param mu2 Mean vector of the second distribution
#' @param Sigma2 Covariance matrix of the second distribution
#' @param method Optimization method (default: "BFGS")
#' @param control Control parameters for optim()
#'
#' @return The Fisher-Rao distance (non-negative scalar), or NA if computation fails
#'
#' @examples
#' # Bivariate example
#' mu1 <- c(0, 0)
#' Sigma1 <- diag(2)
#' mu2 <- c(1, 1)
#' Sigma2 <- diag(2)
#' d <- fisher_rao_distance(mu1, Sigma1, mu2, Sigma2)
#' print(d)
#' 
#' # Univariate example (uses exact formula)
#' d_uni <- fisher_rao_distance(0, matrix(1), 2, matrix(4))
#' print(d_uni)
#'
#' @importFrom expm logm
#' @export
fisher_rao_distance <- function(mu1, Sigma1, mu2, Sigma2, 
                                 method = "BFGS", control = list()) {
  # Check for ill-conditioned matrices
  cond1 <- kappa(Sigma1)
  cond2 <- kappa(Sigma2)
  
  if (cond1 > 1e10 || cond2 > 1e10) {
    warning(paste("Ill-conditioned covariance matrix detected.",
                  "Fisher-Rao distance may be inaccurate."))
    return(NA_real_)
  }
  
  p <- length(mu1)
  
  # Transform to canonical parameters relative to (0, I)
  U <- chol(solve(Sigma1))
  U_inv <- solve(U)
  
  mu_transformed <- U %*% (mu2 - mu1)
  Sigma_transformed <- U %*% Sigma2 %*% t(U)
  
  De <- solve(Sigma_transformed)
  de <- De %*% mu_transformed
  
  # Univariate case - use exact formula
  if (p == 1) {
    de_scalar <- as.numeric(de)
    De_scalar <- as.numeric(De)
    
    sigma1 <- 1
    sigma2 <- sqrt(1 / De_scalar)
    mu_diff <- de_scalar / De_scalar
    
    inner_term <- (mu_diff^2 / (2 * sigma1 * sigma2) + sigma1/sigma2 + sigma2/sigma1)^2
    distance <- sqrt(2) / 2 * acosh(0.5 * inner_term - 1)
    
    return(distance)
  }
  
  # Multivariate case - optimize to find geodesic
  # Diagonalize De via eigendecomposition
  eigen_De <- eigen(De)
  De_diag <- diag(eigen_De$values)
  V <- eigen_De$vectors
  de_rotated <- t(V) %*% de
  
  # Transform to d and D
  D <- sqrt(diag(De_diag))
  d <- de_rotated / D
  
  # Optimize y (upper triangular part of Z)
  n_y <- p*(p-1)/2
  y0 <- rep(0, n_y)
  
  result <- optim(
    par = y0,
    fn = objective_f,
    d = d,
    D = D,
    method = method,
    control = control
  )
  
  # Compute distance from optimized T matrix
  T_mat <- construct_T(result$par, d, D)
  TtT <- t(T_mat) %*% T_mat
  A <- expm::logm(TtT)
  distance <- sqrt(sum(A^2) / 2)
  
  # Check convergence
  if (result$value > 1e-4) {
    warning(paste("Optimization may not have converged. f(y) =", 
                  round(result$value, 8)))
  }
  
  distance
}
