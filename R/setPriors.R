setPriors <- function(distances, X_est, mclust_result, p, G, n, m, model_type) {
  # sigma_sq (inverse-gamma prior)
  ########### Not sure why I'm using these values as priors ###############
  SSR_init <- sum((as.matrix(dist(X_est)) - distances)^2) / 2
  prior_shape <- 5
  prior_scale <- (prior_shape - 1) * SSR_init / m
  ###################################################################

  # mu_j (normal prior)

  if (model_type == "Unequal Unrestricted" || model_type == "Equal Unrestricted") {
    ## Use the means of the estimated object configuration
    prior_alpha <- p+2
    #prior_mean <- mclust_result$parameters$mean
    #prior_Bj <- (prior_alpha - p - 1) * mclust_result$parameters$variance$sigma

    prior_mean <- matrix(NA, nrow = p, ncol = G)
    prior_Bj <- array(NA, c(p,p,G))
    for (i in 1:G) {
      prior_mean[,i] <- colMeans(X_est)
      prior_Bj[,,i] <- ((prior_alpha - p - 1) * cov(X_est)) #/ G #Changed the prior here
    }
    out <- list(prior_shape = prior_shape, prior_scale = prior_scale,
                prior_mean = prior_mean,
                prior_alpha = prior_alpha, prior_Bj = prior_Bj)
  } else if (model_type == "Unequal Diagonal" || model_type == "Equal Diagonal") {
    prior_mean <- matrix(NA, nrow = p, ncol = G)
    prior_IG_alpha <- (p+2) / 2
    prior_IG_beta <- eigen(cov(X_est))$values[1] / 2
    for (i in 1:G) {
      prior_mean[,i] <- colMeans(X_est)
    }
    out <- list(prior_shape = prior_shape, prior_scale = prior_scale,
                prior_mean = prior_mean,
                prior_IG_alpha = prior_IG_alpha, prior_IG_beta = prior_IG_beta)
  } else if (model_type == "Equal Spherical" || model_type == "Unequal Spherical") {
    prior_mean <- matrix(NA, nrow = p, ncol = G)
    prior_IG_alpha <- (p+2) / 2
    prior_IG_beta <- (sum(diag(cov(X_est))) / p) / (G^(2/p)) #eigen(cov(X_est))$values[1] / 2 #
    for (i in 1:G) {
      prior_mean[,i] <- colMeans(X_est)
    }
    out <- list(prior_shape = prior_shape, prior_scale = prior_scale,
                prior_mean = prior_mean,
                prior_IG_alpha = prior_IG_alpha, prior_IG_beta = prior_IG_beta)
  }

  return(out)

  # for (i in 1:G) {
  #   index_vec <- class_mat == i
  #   if (sum(index_vec) == 1) {
  #     prior_mean[, i] <- X_est[which(index_vec), ]
  #   } else if (sum(index_vec) >= 1) {
  #     prior_mean[, i] <- colMeans(X_est[which(index_vec), , drop = FALSE])
  #   } else if (sum(index_vec) == 0) {
  #     prior_mean[, i] <- colMeans(X_est)  ############################################# ASK
  #   }
  # }

  # # T_j (Inverse-Wishart hyperprior for mu_j)
  # prior_alpha <- p+4
  # prior_Bj <- array(NA, c(p,p,G))
  # for (i in 1:G) {
  #   index_vec <- class_mat == i
  #   if (sum(index_vec) > 1) {
  #     prior_Bj[,,i] <- (prior_alpha - p - 1) * cov(X_est[which(index_vec), , drop = FALSE])
  #   } else {
  #     prior_Bj[,,i] <- (prior_alpha - p - 1) * cov(X_est) ############################################# ASK
  #   }
  # }
}


