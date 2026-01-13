#' Wasserstein-2 (Bures) Distance Between Multivariate Normals
#'
#' Computes the Wasserstein-2 distance (also known as Bures or Frechet distance)
#' between two multivariate normal distributions. This metric has a closed-form
#' solution for Gaussian distributions and is based on optimal transport theory.
#'
#' @param mu1 Mean vector of the first distribution (numeric vector)
#' @param Sigma1 Covariance matrix of the first distribution (positive definite matrix)
#' @param mu2 Mean vector of the second distribution (numeric vector)
#' @param Sigma2 Covariance matrix of the second distribution (positive definite matrix)
#'
#' @return The Wasserstein-2 distance (non-negative scalar)
#'
#' @references
#' Takatsu, A. (2011). Wasserstein geometry of Gaussian measures.
#' Osaka Journal of Mathematics, 48(4), 1005-1026.
#'
#' @examples
#' mu1 <- c(0, 0)
#' Sigma1 <- diag(2)
#' mu2 <- c(1, 1)
#' Sigma2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
#' wasserstein2_distance(mu1, Sigma1, mu2, Sigma2)
#'
#' @export
wasserstein2_distance <- function(mu1, Sigma1, mu2, Sigma2) {
  # Check for identical distributions (avoid numerical precision issues)
  if (isTRUE(all.equal(mu1, mu2)) && isTRUE(all.equal(Sigma1, Sigma2))) {
    return(0)
  }

  # Mean part
  mu_dist_sq <- sum((mu1 - mu2)^2)
  
  # Covariance part (Bures metric)
  Sigma1_sqrt <- expm::sqrtm(Sigma1)
  M <- Sigma1_sqrt %*% Sigma2 %*% Sigma1_sqrt
  M_sqrt <- expm::sqrtm(M)
  
  Sigma_dist_sq <- sum(diag(Sigma1)) + sum(diag(Sigma2)) - 2 * sum(diag(M_sqrt))
  
  sqrt(mu_dist_sq + Sigma_dist_sq)
}

#' Hellinger Distance Between Multivariate Normals
#'
#' Computes the Hellinger distance between two multivariate normal distributions.
#' This distance is bounded between 0 and 1, making it interpretable as a similarity measure.
#'
#' @param mu1 Mean vector of the first distribution
#' @param Sigma1 Covariance matrix of the first distribution
#' @param mu2 Mean vector of the second distribution
#' @param Sigma2 Covariance matrix of the second distribution
#'
#' @return The Hellinger distance (scalar in [0, 1])
#'
#' @examples
#' mu1 <- c(0, 0)
#' Sigma1 <- diag(2)
#' mu2 <- c(1, 1)
#' Sigma2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
#' hellinger_distance(mu1, Sigma1, mu2, Sigma2)
#'
#' @export
hellinger_distance <- function(mu1, Sigma1, mu2, Sigma2) {
  Sigma_sum <- (Sigma1 + Sigma2) / 2
  det_term <- sqrt(det(Sigma1) * det(Sigma2)) / det(Sigma_sum)
  
  exp_term <- exp(-0.125 * t(mu1 - mu2) %*% solve(Sigma_sum) %*% (mu1 - mu2))
  
  sqrt(1 - det_term^0.25 * as.numeric(exp_term))
}

#' Bhattacharyya Distance Between Multivariate Normals
#'
#' Computes the Bhattacharyya distance between two multivariate normal distributions.
#' This distance is related to the Bhattacharyya coefficient and Bayes classification error.
#'
#' @param mu1 Mean vector of the first distribution
#' @param Sigma1 Covariance matrix of the first distribution
#' @param mu2 Mean vector of the second distribution
#' @param Sigma2 Covariance matrix of the second distribution
#'
#' @return The Bhattacharyya distance (non-negative scalar)
#'
#' @examples
#' mu1 <- c(0, 0)
#' Sigma1 <- diag(2)
#' mu2 <- c(1, 1)
#' Sigma2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
#' bhattacharyya_distance(mu1, Sigma1, mu2, Sigma2)
#'
#' @export
bhattacharyya_distance <- function(mu1, Sigma1, mu2, Sigma2) {
  Sigma_avg <- (Sigma1 + Sigma2) / 2
  
  as.numeric(0.125 * t(mu1 - mu2) %*% solve(Sigma_avg) %*% (mu1 - mu2) +
    0.5 * log(det(Sigma_avg) / sqrt(det(Sigma1) * det(Sigma2))))
}

#' KL Divergence Between Multivariate Normals
#'
#' Computes the Kullback-Leibler divergence from distribution 1 to distribution 2.
#' Note: KL divergence is asymmetric (not a true distance metric).
#'
#' @param mu1 Mean vector of the first distribution
#' @param Sigma1 Covariance matrix of the first distribution
#' @param mu2 Mean vector of the second distribution
#' @param Sigma2 Covariance matrix of the second distribution
#'
#' @return The KL divergence (non-negative scalar)
#'
#' @examples
#' mu1 <- c(0, 0)
#' Sigma1 <- diag(2)
#' mu2 <- c(1, 1)
#' Sigma2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
#' kl_divergence(mu1, Sigma1, mu2, Sigma2)
#'
#' @export
kl_divergence <- function(mu1, Sigma1, mu2, Sigma2) {
  d <- length(mu1)
  Sigma2_inv <- solve(Sigma2)
  
  as.numeric(0.5 * (sum(diag(Sigma2_inv %*% Sigma1)) + 
         t(mu2 - mu1) %*% Sigma2_inv %*% (mu2 - mu1) - 
         d + log(det(Sigma2) / det(Sigma1))))
}

#' Symmetrized KL Divergence
#'
#' Computes the symmetrized Kullback-Leibler divergence, which is the average
#' of KL(P||Q) and KL(Q||P). This provides a symmetric measure.
#'
#' @param mu1 Mean vector of the first distribution
#' @param Sigma1 Covariance matrix of the first distribution
#' @param mu2 Mean vector of the second distribution
#' @param Sigma2 Covariance matrix of the second distribution
#'
#' @return The symmetrized KL divergence (non-negative scalar)
#'
#' @export
symmetrized_kl <- function(mu1, Sigma1, mu2, Sigma2) {
  0.5 * (kl_divergence(mu1, Sigma1, mu2, Sigma2) + 
         kl_divergence(mu2, Sigma2, mu1, Sigma1))
}

#' Affine-Invariant Distance
#'
#' Computes the affine-invariant Riemannian distance between two multivariate
#' normal distributions. Uses the log-Euclidean metric on covariance matrices.
#'
#' @param mu1 Mean vector of the first distribution
#' @param Sigma1 Covariance matrix of the first distribution
#' @param mu2 Mean vector of the second distribution
#' @param Sigma2 Covariance matrix of the second distribution
#'
#' @return The affine-invariant distance (non-negative scalar)
#'
#' @references
#' Pennec, X., Fillard, P., & Ayache, N. (2006). A Riemannian framework
#' for tensor computing. International Journal of Computer Vision, 66(1), 41-66.
#'
#' @export
affine_invariant_distance <- function(mu1, Sigma1, mu2, Sigma2) {
  # Mean part
  Sigma1_inv <- solve(Sigma1)
  mu_dist_sq <- as.numeric(t(mu1 - mu2) %*% Sigma1_inv %*% (mu1 - mu2))
  
  # Covariance part (log-Euclidean metric)
  Sigma1_inv_sqrt <- expm::sqrtm(Sigma1_inv)
  M <- Sigma1_inv_sqrt %*% Sigma2 %*% Sigma1_inv_sqrt
  Sigma_dist_sq <- sum(expm::logm(M)^2)
  
  sqrt(mu_dist_sq + Sigma_dist_sq)
}
