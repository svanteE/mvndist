#' Compare Multiple Distance Measures
#'
#' Computes all available distance measures between two multivariate normal
#' distributions and returns them in a named vector for easy comparison.
#'
#' @param mu1 Mean vector of the first distribution
#' @param Sigma1 Covariance matrix of the first distribution
#' @param mu2 Mean vector of the second distribution
#' @param Sigma2 Covariance matrix of the second distribution
#' @param include_fisher Logical; include Fisher-Rao distance (slow, may fail for
#'   ill-conditioned matrices). Default: FALSE
#'
#' @return A named numeric vector with all computed distances
#'
#' @examples
#' mu1 <- c(0, 0)
#' Sigma1 <- diag(2)
#' mu2 <- c(1, 1)
#' Sigma2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
#' compare_distances(mu1, Sigma1, mu2, Sigma2)
#'
#' @export
compare_distances <- function(mu1, Sigma1, mu2, Sigma2, include_fisher = FALSE) {
  distances <- c(
    Wasserstein2 = wasserstein2_distance(mu1, Sigma1, mu2, Sigma2),
    Hellinger = hellinger_distance(mu1, Sigma1, mu2, Sigma2),
    Bhattacharyya = bhattacharyya_distance(mu1, Sigma1, mu2, Sigma2),
    KL_symmetrized = symmetrized_kl(mu1, Sigma1, mu2, Sigma2),
    Affine_Invariant = affine_invariant_distance(mu1, Sigma1, mu2, Sigma2)
  )
  
  if (include_fisher) {
    distances <- c(distances, 
                  Fisher_Rao = fisher_rao_distance(mu1, Sigma1, mu2, Sigma2))
  }
  
  distances
}
