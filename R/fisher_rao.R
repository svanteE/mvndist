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

#' Helper: Gradient of objective function (optional, for speed)
#' @noRd
gradient_f <- function(y, d, D) {
  # Numerical gradient via finite differences
  p <- length(d)
  n_y <- length(y)
  grad <- numeric(n_y)
  h <- 1e-8
  
  f0 <- objective_f(y, d, D)
  
  for (i in 1:n_y) {
    y_plus <- y
    y_plus[i] <- y[i] + h
    grad[i] <- (objective_f(y_plus, d, D) - f0) / h
  }
  
  grad
}

#' Fisher-Rao Geodesic (Optimized)
#'
#' Computes the Fisher-Rao geodesic path between two multivariate normal distributions
#' using the method of Eriksen (1987), which solves the geodesic equations for the
#' Fisher information metric formulated by Skovgaard (1984). The method minimizes
#' f(y) = tr(C(y)^2) where C(y) is the off-diagonal block in the matrix logarithm
#' decomposition. Once optimal y is found, it constructs A(y) and computes Λ(t) = exp(t·A)
#' to extract the geodesic path μ(t) and Σ(t).
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
#' (Established the Riemannian geometry framework and differential equations)
#'
#' Eriksen, P. S. (1987). Geodesics connected with the Fisher metric on
#' the multivariate normal manifold. Proceedings of the GST Workshop.
#' (Solved the geodesic equations, enabling practical computation)
#'
#' @importFrom expm logm expm
#' @export
fisher_rao_geodesic_numerical <- function(mu1, Sigma1, mu2, Sigma2, 
                                          n_steps = 100, tol = 1e-6) {
  p <- length(mu1)
  
  # Check for ill-conditioned matrices
  cond1 <- kappa(Sigma1)
  cond2 <- kappa(Sigma2)
  
  if (cond1 > 1e10 || cond2 > 1e10) {
    warning("Ill-conditioned covariance matrix detected. Geodesic may be inaccurate.")
    return(NULL)
  }
  
  # Transform to canonical parameters relative to (0, I)
  U <- chol(solve(Sigma1))
  U_inv <- solve(U)
  
  mu_transformed <- U %*% (mu2 - mu1)
  Sigma_transformed <- U %*% Sigma2 %*% t(U)
  
  De <- solve(Sigma_transformed)
  de <- De %*% mu_transformed
  
  # Univariate case
  if (p == 1) {
    de_scalar <- as.numeric(de)
    De_scalar <- as.numeric(De)
    
    sigma1 <- 1
    sigma2 <- sqrt(1 / De_scalar)
    mu_diff <- de_scalar / De_scalar
    
    # Use exact geodesic formula for univariate case
    times <- seq(0, 1, length.out = n_steps)
    
    geodesic <- lapply(times, function(t) {
      # Interpolate mean linearly
      mu_t <- (1 - t) * 0 + t * mu_diff
      
      # Interpolate sigma on log scale (geometric mean)
      sigma_t <- sigma1^(1-t) * sigma2^t
      Sigma_t <- matrix(sigma_t^2, 1, 1)
      
      # Transform back to original coordinates
      mu_original <- mu1 + U_inv %*% matrix(mu_t, ncol = 1)
      Sigma_original <- U_inv %*% Sigma_t %*% t(U_inv)
      
      list(t = t, 
           mu = as.vector(mu_original),
           Sigma = Sigma_original)
    })
    
    return(geodesic)
  }
  
  # Multivariate case - use optimized method
  # Diagonalize De via eigendecomposition
  eigen_De <- eigen(De)
  De_diag <- diag(eigen_De$values)
  V <- eigen_De$vectors
  de_rotated <- t(V) %*% de
  
  # Transform to d and D
  D <- sqrt(diag(De_diag))
  d <- de_rotated / D
  
  # Optimize y to minimize tr(C(y)^2)
  n_y <- p*(p-1)/2
  y0 <- rnorm(n_y, mean = 0, sd = 0.01)
  
  # Multi-restart optimization
  best_result <- list(value = Inf)
  
  for (attempt in 1:3) {
    result <- optim(
      par = if(attempt == 1) y0 else rnorm(n_y, sd = 0.05),
      fn = objective_f,
      gr = gradient_f,
      d = d,
      D = D,
      method = "BFGS",
      control = list()
    )
    
    if (result$value < best_result$value) {
      best_result <- result
    }
    
    if (result$value < 1e-6) break
  }
  
  y_opt <- best_result$par
  
  # Check convergence
  if (best_result$value > 1e-4) {
    warning(paste("Optimization may not have converged. f(y) =", 
                  round(best_result$value, 8)))
  }
  
  # Construct optimal T matrix and compute A = log(T'T)
  T_mat <- construct_T(y_opt, d, D)
  TtT <- t(T_mat) %*% T_mat
  A <- expm::logm(TtT)
  
  # At optimum, C(y) = 0, so A has the block structure:
  # A = [A11, a12,  0   ]
  #     [a21,  0,  -a12^t]
  #     [0,  -a12, -A11]
  # where A11 is p×p, a12 is p×1
  
  # Generate geodesic path: Λ(t) = exp(t·A)
  times <- seq(0, 1, length.out = n_steps)
  
  geodesic <- lapply(times, function(t) {
    # Compute Λ(t) = exp(t·A)
    Lambda_t <- expm::expm(t * A)
    
    # Extract De(t) and de(t) from Λ(t) - these are already in the rotated (d,D) space
    # Λ[1:p, 1:p] = De(t) in rotated coordinates
    # Λ[1:p, p+1] = de(t) in rotated coordinates
    De_t_rotated <- Lambda_t[1:p, 1:p]
    de_t_rotated <- Lambda_t[1:p, p+1, drop=FALSE]
    
    # Rotate back using V to get to canonical space (before the d,D rotation)
    De_t_canonical <- V %*% De_t_rotated %*% t(V)
    de_t_canonical <- V %*% de_t_rotated
    
    # Convert (De, de) to (mu, Sigma) in canonical space
    Sigma_t_canonical <- solve(De_t_canonical)
    mu_t_canonical <- Sigma_t_canonical %*% de_t_canonical
    
    # Transform back to original coordinates
    mu_original <- mu1 + U_inv %*% mu_t_canonical
    Sigma_original <- U_inv %*% Sigma_t_canonical %*% t(U_inv)
    
    list(t = t, 
         mu = as.vector(mu_original),
         Sigma = Sigma_original)
  })
  
  geodesic
}

#' Fisher-Rao Distance (Optimized)
#'
#' Computes the Fisher-Rao distance using the geodesic solution method of Eriksen (1987).
#' This method solves the optimization problem f(y) = tr(C(y)^2) = 0 to find the exact
#' geodesic for the Fisher information metric (Skovgaard 1984), then computes the distance
#' as the path length. For univariate distributions, uses the exact closed-form solution.
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
  
  # Better initial guess: use small random perturbations
  # This helps escape poor local minima
  y0 <- rnorm(n_y, mean = 0, sd = 0.01)
  
  # Try optimization with multiple restarts if it doesn't converge
  best_result <- list(value = Inf)
  
  for (attempt in 1:3) {
    result <- optim(
      par = if(attempt == 1) y0 else rnorm(n_y, sd = 0.05),
      fn = objective_f,
      gr = gradient_f,  # Use gradient for faster convergence
      d = d,
      D = D,
      method = method,
      control = control
    )
    
    if (result$value < best_result$value) {
      best_result <- result
    }
    
    # Break early if converged well
    if (result$value < 1e-6) break
  }
  
  result <- best_result
  
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
