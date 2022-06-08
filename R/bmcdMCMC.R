#' @importFrom gtools rdirichlet
#' @importFrom LaplacesDemon rinvwishart
#' @importFrom mvtnorm rmvnorm
#' @importFrom RcppHungarian HungarianSolver
#' @importFrom combinat permn


bmcdMCMC <- function(distances, mcmc_list, priors, p, G, n, m, bmcd_iter, bmcd_burn, labelswitch_iter) {

  # Unpack list
  X_list = mcmc_list$X_list
  sigma_sq_list <- mcmc_list$sigma_sq_list
  eps_list <- mcmc_list$eps_list
  mu_list <- mcmc_list$mu_list
  n_list <- mcmc_list$n_list
  T_list <- mcmc_list$T_list
  z_list <- mcmc_list$z_list
  class_list <- mcmc_list$class_list

  # Initialize various MCMC variables (proposal variances, acceptance rates, transformation parameters, label-switching parameters)
  tri_ind <- lower.tri(matrix(data=NA, nrow = n, ncol = n))
  X_star = X_list[[1]] # For transformation of X
  SSR_init <- sum((as.matrix(dist(X_list[[1]])) - distances)^2)

  x_prop_var <- (2.38^2) * sigma_sq_list[[1]] / (n-1) ############ how do I choose a good proposal variance?
  sigma_prop_variance <- (2.38^2) * 2 * (SSR_init / 2)^2 / ((m-2)*(m-2)*(m-4))

  accept_x <- 0
  accept_sigma <- 0

  init_theta <-rep(list(list()), G)
  init_s <- init_theta
  qwe <- c()

  for (t in 2:bmcd_iter) {
    # Iteration printing  -----------------------------------------------------
    if (t %% 100 == 0) {
      print(t)
    }

    # Generating X using random walk M-H algorithm ----------------------------

    for (i in 1:n) {
      j <- class_list[t-1, i]


      x_old <- X_list[[t-1]][i,]
      x_new <- rnorm(p, mean = x_old, sd = sqrt(x_prop_var))

      X_new <- X_list[[t-1]]
      X_new[i,] = x_new

      ## First term
      A_old <- as.matrix(x_old - mu_list[[t-1]][ ,j])
      Q1_old <- t(A_old) %*% solve(T_list[[t-1]][,,j], tol = 1e-17) %*% A_old

      A_new <- as.matrix(x_new - mu_list[[t-1]][ ,j])
      Q1_new <- t(A_new) %*% solve(T_list[[t-1]][,,j], tol = 1e-17) %*% A_new

      ## Second term
      delta_new <- as.matrix(dist(X_new))
      delta_old <- as.matrix(dist(X_list[[t-1]]))

      diff_new <- delta_new - distances
      diff_old <- delta_old - distances

      Q2_new <- (1/(2*sigma_sq_list[t-1])) * sum(diff_new[tri_ind]^2)
      Q2_old <- (1/(2*sigma_sq_list[t-1])) * sum(diff_old[tri_ind]^2)

      ## Third term

      norm_new <- delta_new / sqrt(sigma_sq_list[t-1])
      norm_old <- delta_old / sqrt(sigma_sq_list[t-1])

      Q3_new <- sum((pnorm(norm_new[tri_ind], log.p = TRUE)))
      Q3_old <- sum((pnorm(norm_old[tri_ind], log.p = TRUE)))

      ## Calculate the ratio (using log instead to simplify things)
      ratio <- -0.5*(Q1_new - Q1_old) - (Q2_new - Q2_old) + (Q3_new - Q3_old)

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
    SSR <- sum((delta - distances)^2) / 2

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

      T_list_pst <- as.matrix(priors$prior_Bj[,,k] + (S_j/2))
      T_list[[t]][,,k] <- LaplacesDemon::rinvwishart(priors$prior_alpha + (n_list[t-1, k] / 2), T_list_pst)
      # tryCatch({
      #   T_list[[t]][,,k] <<- LaplacesDemon::rinvwishart(priors$prior_alpha + (n_list[t-1, k] / 2), T_list_pst)
      # }, error = function(e) {
      #   diag(T_list_pst) <- diag(T_list_pst) + 1e-05
      #   T_list[[t]][,,k] <<- LaplacesDemon::rinvwishart(priors$prior_alpha + (n_list[t-1, k] / 2), T_list_pst)
      # }
      # )

      if (n_list[t-1, k] > 1) {
        x_bar_j <- colMeans(x_j)
      } else if (n_list[t-1, k] == 1) {
        x_bar_j <- x_j
      } else if (n_list[t-1, k] == 0) {
        x_bar_j <- 0
      }

      pst_mean = (n_list[t-1, k]  * x_bar_j + priors$prior_mean[, k]) / (n_list[t-1, k] + 1)

      pst_var =  T_list[[t]][,,k] / (n_list[t-1, k]  + 1)

      if (p > 1) {
        mu_list[[t]][,k] = mvtnorm::rmvnorm(1, mean = pst_mean, sigma = pst_var)
      } else {
        mu_list[[t]][,k] = rnorm(1, mean = pst_mean, sd = sqrt(pst_var))
      }
    }


    # Relabeling procedure ----------------------------------------------------

    if (t == labelswitch_iter) {
      for (comp in 1:G) {

        init_theta[[comp]][[1]] <- (1 / labelswitch_iter) * sum(eps_list[1:labelswitch_iter, comp])
        init_theta[[comp]][[2]] <- (1 / labelswitch_iter) * Reduce(`+`,
                                                          rapply(mu_list[1:labelswitch_iter],
                                                                 classes = 'matrix', how = 'list',
                                                                 f = function(x) x[, comp, drop = FALSE]))
        init_theta[[comp]][[3]] <- (1 / labelswitch_iter) *
          Reduce(`+`, lapply(T_list[1:labelswitch_iter], function(x) x[,,comp]))


        init_s[[comp]][[1]] <- (1 / labelswitch_iter) * sum((eps_list[1:labelswitch_iter,comp] - init_theta[[comp]][[1]])^2)

        mu_temp <- rapply(mu_list[1:labelswitch_iter],
                          classes = 'matrix', how = 'list',
                          f = function(x) x[, comp, drop = FALSE])

        init_s[[comp]][[2]] <- (1 / labelswitch_iter) *
          Reduce(`+`,
                 lapply(mu_temp,
                        function(x) (x[, 1] - init_theta[[comp]][[2]])^2))

        T_temp <- lapply(T_list[1:labelswitch_iter], function(x) x[,,comp])

        init_s[[comp]][[3]] <- (1 / labelswitch_iter) *
          Reduce(`+`, lapply(T_temp, function(x) (x - init_theta[[comp]][[3]])^2))
      }

      theta_perms <- permn(init_theta)
      perm_labels <- permn(1:G)

    } else if (t > labelswitch_iter) {
      r <- t - labelswitch_iter

      ## Testing
      eps_theta <- rep(NA, G)

      for (comp in 1:G) {
        eps_theta[comp] <- init_theta[[comp]][[1]]

      }
      eps_list[labelswitch_iter+r,]
      test_mat <- matrix(NA, nrow = G,  ncol = G)
      for (comp in 1:G) {
        test_mat[comp, ] <- (eps_list[labelswitch_iter+r, comp] - eps_theta)^2 / init_s[[comp]][[1]]
      }

      ## mu
      mu_theta <- matrix(NA, nrow = p, ncol = G)

      for (comp in 1:G) {
        mu_theta[,comp] <- init_theta[[comp]][[2]]
      }

      test_mat2 <- matrix(NA, nrow = G,  ncol = G)
      for (comp in 1:G) {
        test_mat2[comp, ] <- colSums(
          sweep((mu_list[[labelswitch_iter+r]][,comp] - mu_theta)^2, 1, init_s[[comp]][[2]], FUN = '/'))
      }


      ## T
      T_theta <- array(NA, c(p,p,G))
      for (comp in 1:G) {
        T_theta[,,comp] <- init_theta[[comp]][[3]]
      }

      test_mat3 <- matrix(NA, nrow = G, ncol = G)
      for (comp in 1:G) {
        test_mat3[comp, ] <- apply(T_theta, 3, function(x) {
          sum((T_list[[labelswitch_iter+r]][,,comp] - x)^2 / init_s[[comp]][[3]])
        }, simplify = TRUE)
      }


      ## Assignment problem from here
      assignment_solution <- HungarianSolver(test_mat+test_mat2+test_mat3)
      new_ind <- assignment_solution$pairs[,2]

      if(sum(new_ind != c(1:G)) != 0) qwe <- c(qwe, t)
      eps_list[labelswitch_iter+r,] <- eps_list[labelswitch_iter+r, new_ind]
      mu_list[[labelswitch_iter+r]] <- mu_list[[labelswitch_iter+r]][ , new_ind, drop = FALSE]
      T_list[[labelswitch_iter+r]] <- T_list[[labelswitch_iter+r]][ ,, new_ind, drop = FALSE]

      ## Step 2
      old_theta <- init_theta
      for (comp in 1:G) {
        init_theta[[comp]][[1]] <- (((labelswitch_iter + r - 1) / (labelswitch_iter + r)) * init_theta[[comp]][[1]]) +
          ((1 / (labelswitch_iter + r)) * eps_list[labelswitch_iter+r, comp])
        init_theta[[comp]][[2]] <- (((labelswitch_iter + r - 1) / (labelswitch_iter + r)) * init_theta[[comp]][[2]]) +
          ((1 / (labelswitch_iter + r)) * mu_list[[labelswitch_iter+r]][, comp])
        init_theta[[comp]][[3]] <- (((labelswitch_iter + r - 1) / (labelswitch_iter + r)) * init_theta[[comp]][[3]]) +
          ((1 / (labelswitch_iter + r)) * T_list[[labelswitch_iter+r]][,, comp])

        init_s[[comp]][[1]] <- (((labelswitch_iter + r - 1) / (labelswitch_iter + r)) * init_s[[comp]][[1]]) +
          (((labelswitch_iter + r - 1) / (labelswitch_iter + r)) * ((old_theta[[comp]][[1]] - init_theta[[comp]][[1]])^2)) +
          ((1 / (labelswitch_iter+r)) * ((eps_list[labelswitch_iter+r, comp] - init_theta[[comp]][[1]])^2))
        init_s[[comp]][[2]] <- (((labelswitch_iter + r - 1) / (labelswitch_iter + r)) * init_s[[comp]][[2]]) +
          (((labelswitch_iter + r - 1) / (labelswitch_iter + r)) * (old_theta[[comp]][[2]] - init_theta[[comp]][[2]])^2) +
          ((1 / (labelswitch_iter+r)) * ((mu_list[[labelswitch_iter+r]][, comp] - init_theta[[comp]][[2]])^2))
        init_s[[comp]][[3]] <- (((labelswitch_iter + r - 1) / (labelswitch_iter + r)) * init_s[[comp]][[3]]) +
          (((labelswitch_iter + r - 1) / (labelswitch_iter + r)) * (old_theta[[comp]][[3]] - init_theta[[comp]][[3]])^2) +
          ((1 / (labelswitch_iter+r)) * ((T_list[[labelswitch_iter+r]][,, comp] - init_theta[[comp]][[3]])^2))
      }
    }



    # Calculate cluster probabilities -----------------------------------------

    for (a in 1:n) {
      denom = 0
      for (w in 1:G) {
        denom = denom + #(eps_list[t, w] *
                           mvtnorm::dmvnorm(X_list[[t]][a, , drop = FALSE],
                                            mean = mu_list[[t]][, w, drop = FALSE],
                                            sigma = matrix(T_list[[t]][,,w, drop = FALSE], ncol = p, nrow = p))#)
      }
      for (k in 1:G) {
        z_list[[t]][a, k] = #eps_list[t, k] *
          mvtnorm::dmvnorm(X_list[[t]][a, , drop = FALSE],
                           mean = mu_list[[t]][, k, drop = FALSE],
                           sigma = matrix(T_list[[t]][,,k, drop = FALSE], ncol = p, nrow = p)) / denom
      }
    }


    # Cluster assignment ------------------------------------------------------

    clust <- apply(z_list[[t]], 1, which.max)
    for (k in 1:G) {
      n_list[t, k] <- sum(clust == k)
    }
    class_list[t, ] <- clust
  }

  # Discard burn-in ---------------------------------------------------------

  eps_list <- eps_list[(bmcd_burn+1):bmcd_iter, ]
  mu_list <- mu_list[(bmcd_burn+1):bmcd_iter]
  T_list <- T_list[((bmcd_burn+1):bmcd_iter)]
  sigma_sq_list <- sigma_sq_list[(bmcd_burn+1):bmcd_iter]
  X_list <- X_list[(bmcd_burn+1):bmcd_iter]
  z_list <- z_list[(bmcd_burn+1):bmcd_iter]
  n_list <- n_list[(bmcd_burn+1):bmcd_iter, ]
  class_list <- class_list[(bmcd_burn+1):bmcd_iter, ]


  # Return list ------------------------------------------------------------

  list(X_list = X_list,
       sigma_sq_list = sigma_sq_list,
       eps_list = eps_list,
       mu_list = mu_list,
       T_list = T_list,
       z_list = z_list,
       n_list = n_list,
       class_list = class_list)

}
