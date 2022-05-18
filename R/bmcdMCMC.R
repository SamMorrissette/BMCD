#' @importFrom gtools rdirichlet
#' @importFrom LaplacesDemon rinvwishart
#' @importFrom mvtnorm rmvnorm


bmcdMCMC <- function(distances, mcmc_list, priors, p, G, n, m, bmcd_iter, labelswitch_iter) {

  # Unpack list
  X_list = mcmc_list$X_list
  sigma_sq_list <- mcmc_list$sigma_sq_list
  eps_list <- mcmc_list$eps_list
  mu_list <- mcmc_list$mu_list
  n_list <- mcmc_list$n_list
  T_list <- mcmc_list$T_list
  z_list <- mcmc_list$z_list
  class_list <- mcmc_list$class_list

  tri_ind <- lower.tri(matrix(data=NA, nrow = n, ncol = n))
  X_star = X_list[[1]] # For transformation of X
  SSR_init <- sum((as.matrix(dist(X_list[[1]])) - distances)^2)

  x_prop_var <- (2.38^2) * sigma_sq_list[[1]] / (n-1) ############ how do I choose a good proposal variance?
  sigma_prop_variance <- (2.38^2) * 2 * (SSR_init / 2)^2 / ((m-2)*(m-2)*(m-4))

  accept_x <- 0
  accept_sigma <- 0

  for (t in 2:bmcd_iter) {


    # Generating X using random walk M-H algorithm ----------------------------

    for (i in 1:n) {
      j <- class_list[t-1, i]


      x_old <- X_list[[t-1]][i,]
      x_new <- rnorm(p, mean = x_old, sd = sqrt(x_prop_var))

      X_new <- X_list[[t-1]]
      X_new[i,] = x_new

      ## First term
      A_old <- as.matrix(x_old - mu_list[[t-1]][ ,j])
      Q1_old <- t(A_old) %*% solve(T_list[[t-1]][,,j]) %*% A_old

      A_new <- as.matrix(x_new - mu_list[[t-1]][ ,j])
      Q1_new <- t(A_new) %*% solve(T_list[[t-1]][,,j]) %*% A_new

      ## Second term
      delta_new <- as.matrix(dist(X_new))
      delta_old <- as.matrix(dist(X_list[[t-1]]))

      diff_new <- delta_new - distances
      diff_old <- delta_old - distances

      Q2_new <- (1/sigma_sq_list[t-1]) * sum(diff_new[tri_ind]^2)
      Q2_old <- (1/sigma_sq_list[t-1]) * sum(diff_old[tri_ind]^2)

      ## Third term

      norm_new <- delta_new / sqrt(sigma_sq_list[t-1])
      norm_old <- delta_old / sqrt(sigma_sq_list[t-1])

      Q3_new <- prod(pnorm(norm_new[tri_ind]))
      Q3_old <- prod(pnorm(norm_old[tri_ind]))

      ## Calculate the ratio (using log instead to simplify things)
      ratio <- -0.5*(Q1_new - Q1_old) - (Q2_new - Q2_old) + (log(Q3_new) - log(Q3_old))

      ## Have to take log of random uniform variable since we are working on log scale
      if (log(runif(1)) < ratio) {
        X_list[[t]][i, ] <- x_new
        accept_x <- accept_x + 1
      } else {
        X_list[[t]][i, ] <- x_old
      }
    }

    # Transform X using Procrustean Similarity Transformation -----------------

    ## The following code uses the notation used in the main paper

    vec_1 <- matrix(1, nrow = n, ncol = 1)
    J = diag(1, nrow = n, ncol = n) - ((1/n) * (vec_1 %*% t(vec_1)))
    C = t(X_star) %*% J %*% X_list[[t]]
    svd_C = svd(C)
    mat_T = svd_C$v %*% t(svd_C$u)
    little_t = (1/n) * (t(X_star - X_list[[t]] %*% mat_T) %*% vec_1)
    X_list[[t]] = (X_list[[t]] %*% mat_T) + (vec_1 %*% t(little_t))


    # Generate sigma_sq -------------------------------------------------------

    delta <- as.matrix(dist(X_list[[t]]))
    SSR <- sum((delta - distances)^2)

    sigma_sq_old <- sigma_sq_list[t-1]
    sigma_sq_new <- rnorm(1, mean = sigma_sq_old, sd = sqrt(sigma_prop_variance))
    while (sigma_sq_new < 0) {
      sigma_sq_new <- rnorm(1, mean = sigma_sq_old, sd = sqrt(sigma_prop_variance)) # Restricted to be positive
    }


    ## Third term
    norm_old <- delta / sqrt(sigma_sq_old)
    norm_new <- delta / sqrt(sigma_sq_new)


    ratio <- -sum(log(pnorm(norm_new[tri_ind])/pnorm(norm_old[tri_ind]))) -
      ((0.5*SSR + priors$prior_scale) * ((1/sigma_sq_new) - (1/sigma_sq_old))) -
      (0.5*m + priors$prior_shape + 1) * log(sigma_sq_new/sigma_sq_old)

    if (log(runif(1)) < ratio) {
      sigma_sq_list[t] <- sigma_sq_new
      accept_sigma <- accept_sigma + 1
    } else {
      sigma_sq_list[t] <- sigma_sq_old
    }


    # Generate epsilon --------------------------------------------------------

    eps_list[t, ] <- gtools::rdirichlet(1, n_list[t-1,] + 1)

    # Generate mu and T for each component ------------------------------------

    for (k in 1:G) {
      if (n_list[t-1, k] > 0) {
        x_j <- X_list[[t]][which(class_list[t-1, ] == k), ]
        x_j <- matrix(x_j, nrow = n_list[t-1, k])
        centered_x <- sweep(x_j, 2, mu_list[[t-1]][, k])
        S_j <- 0
        for (q in 1:nrow(centered_x)) {
          S_j = S_j + (centered_x[q, ] %*% t(centered_x[q,]))
        }
      } else {
        S_j <- 0
      }

      T_list_pst <- priors$prior_Bj[,,k] + (S_j/2)
      tryCatch({
        T_list[[t]][,,k] <<- LaplacesDemon::rinvwishart(priors$prior_alpha + (n_list[t-1, k] / 2), T_list_pst)
      }, error = function(e) {
        diag(T_list_pst) <- diag(T_list_pst) + 1e-05
        T_list[[t]][,,k] <<- LaplacesDemon::rinvwishart(priors$prior_alpha + (n_list[t-1, k] / 2), T_list_pst)
      }
      )

      if (n_list[t-1, k] > 1) {
        x_bar_j <- colMeans(x_j)
      } else if (n_list[t-1, k] == 1) {
        x_bar_j <- x_j
      } else if (n_list[t-1, k] == 0) {
        x_bar_j <- 0
      }

      pst_mean = (n_list[t-1, k]  * x_bar_j + priors$prior_mean[, k]) / (n_list[t-1, k] + 1)

      pst_var =  T_list[[t]][,,k] / (n_list[t-1, k]  + 1)

      mu_list[[t]][,k] = mvtnorm::rmvnorm(1, mean = pst_mean, sigma = pst_var)
    }


    # Relabeling procedure ----------------------------------------------------



  }
}
