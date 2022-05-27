setPriors <- function(distances, X_est, class_mat, p, G, n, m) {
  # sigma_sq (inverse-gamma prior)
  ########### Not sure why I'm using these values as priors ###############
  SSR_init <- sum((as.matrix(dist(X_est)) - distances)^2)
  prior_shape <- 5
  prior_scale <- (prior_shape - 1) * SSR_init / m
  ###################################################################

  # mu_j (normal prior)

  ## Use the means of the estimated object configuration
  prior_mean <- matrix(NA, nrow = p, ncol = G)
  for (i in 1:G) {
    index_vec <- class_mat == i
    if (sum(index_vec) == 1) {
      prior_mean[, i] <- X_est[which(index_vec), ]
    } else if (sum(index_vec) >= 1) {
      prior_mean[, i] <- colMeans(X_est[which(index_vec), , drop = FALSE])
    } else if (sum(index_vec) == 0) {
      prior_mean[, i] <- colMeans(X_est)  ############################################# ASK
    }
  }

  # T_j (Inverse-Wishart hyperprior for mu_j)
  prior_alpha <- p+4
  prior_Bj <- array(NA, c(p,p,G))
  for (i in 1:G) {
    index_vec <- class_mat == i
    if (sum(index_vec) > 1) {
      prior_Bj[,,i] <- (prior_alpha - p - 1) * cov(X_est[which(index_vec), , drop = FALSE])
    } else {
      prior_Bj[,,i] <- (prior_alpha - p - 1) * cov(X_est) ############################################# ASK
    }
  }
  list(prior_shape = prior_shape, prior_scale = prior_scale,
       prior_mean = prior_mean,
       prior_alpha = prior_alpha, prior_Bj = prior_Bj)
}


