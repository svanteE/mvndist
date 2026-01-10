library(testthat)
library(mvndist)

test_that("Wasserstein-2 distance works", {
  mu1 <- c(0, 0)
  Sigma1 <- diag(2)
  mu2 <- c(1, 1)
  Sigma2 <- diag(2)
  
  d <- wasserstein2_distance(mu1, Sigma1, mu2, Sigma2)
  
  expect_true(d >= 0)
  expect_true(is.finite(d))
  expect_equal(d, sqrt(2), tolerance = 1e-10)  # Known closed form
})

test_that("Hellinger distance is bounded [0,1]", {
  mu1 <- c(0, 0)
  Sigma1 <- diag(2)
  mu2 <- c(1, 1)
  Sigma2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  
  d <- hellinger_distance(mu1, Sigma1, mu2, Sigma2)
  
  expect_true(d >= 0)
  expect_true(d <= 1)
})

test_that("Distance is zero for identical distributions", {
  mu <- c(1, 2, 3)
  Sigma <- matrix(c(1, 0.5, 0.3, 0.5, 2, 0.4, 0.3, 0.4, 1.5), 3, 3)
  
  expect_equal(wasserstein2_distance(mu, Sigma, mu, Sigma), 0, tolerance = 1e-10)
  expect_equal(hellinger_distance(mu, Sigma, mu, Sigma), 0, tolerance = 1e-10)
  expect_equal(bhattacharyya_distance(mu, Sigma, mu, Sigma), 0, tolerance = 1e-10)
})

test_that("KL divergence is asymmetric", {
  mu1 <- c(0, 0)
  Sigma1 <- diag(2)
  mu2 <- c(1, 1)
  Sigma2 <- diag(2) * 2
  
  kl1 <- kl_divergence(mu1, Sigma1, mu2, Sigma2)
  kl2 <- kl_divergence(mu2, Sigma2, mu1, Sigma1)
  
  expect_false(isTRUE(all.equal(kl1, kl2)))
})

test_that("Symmetrized KL is symmetric", {
  mu1 <- c(0, 0)
  Sigma1 <- diag(2)
  mu2 <- c(1, 1)
  Sigma2 <- diag(2) * 2
  
  skl1 <- symmetrized_kl(mu1, Sigma1, mu2, Sigma2)
  skl2 <- symmetrized_kl(mu2, Sigma2, mu1, Sigma1)
  
  expect_equal(skl1, skl2, tolerance = 1e-10)
})

test_that("compare_distances returns all metrics", {
  mu1 <- c(0, 0)
  Sigma1 <- diag(2)
  mu2 <- c(1, 1)
  Sigma2 <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  
  dists <- compare_distances(mu1, Sigma1, mu2, Sigma2)
  
  expect_length(dists, 5)
  expect_true(all(dists >= 0))
  expect_true(all(is.finite(dists)))
})
