#' @importFrom mclust Mclust mclustBIC
#' @importFrom foreach %dopar% foreach
#' @importFrom iterators icount
#' @importFrom doParallel registerDoParallel
#' @export

bmcd <- function(distances, bmds_object, p, min_G, max_G, bmcd_iter, bmcd_burn, labelswitch_iter,
                 parallel = FALSE, num_cores = 0, model_type = "Unequal Unrestricted") {

  n <- nrow(distances)
  m <- n*(n-1)/2

  if(parallel == TRUE & num_cores > 0) {
    doParallel::registerDoParallel(cores=num_cores)
  }

  X_out <- bmds_object$X_out
  sigma_sq_out <- bmds_object$sigma_sq_out
  X_est <- X_out[[p]]

  # Set initial parameter values using GMM  ---------------------------------

  out <- priors <- vector("list", length = max_G - min_G + 1)



  if (parallel == FALSE & num_cores == 0) {
    ind <- 1
    for (G in min_G:max_G) {
      mclust_out <- Mclust(X_est, G = G) # Fit GMM
      print(mclust_out)

      sigma_sq <- sigma_sq_out[[p]] # Posterior mean of sigma squared (measurement error) from BMDS algorithm
      eps_init <- mclust_out$parameters$pro # Vector of component mixing proportions
      mu_init <- mclust_out$parameters$mean # Vector of component means
      T_init <- mclust_out$parameters$variance$sigma # List of component covariance matrices
      z_init <- mclust_out$z # Matrix of cluster assignment probabilities

      ## Accounting for one dimension
      mu_init <- matrix(mu_init, nrow = p, ncol = G)
      T_init <- array(T_init, c(p,p,G))

      # Set priors --------------------------------------------------------------

      priors[[ind]] <- setPriors(distances, X_est, mclust_out, p, G, n, m, model_type)

      # Initialize MCMC lists ---------------------------------------------------
      init_params <- list(sigma_sq, eps_init, mu_init, T_init, z_init)
      mcmc_list <- initLists(init_params, X_est, mclust_out$classification, p, G, n, bmcd_iter)

      # Run MCMC ----------------------------------------------------------------
      out[[ind]] <- bmcdMCMC(distances, mcmc_list, priors[[ind]], p, G, n, m, bmcd_iter, bmcd_burn, labelswitch_iter, model_type)
      out[[ind]]$G <- G
      ind <- ind + 1
    }
  } else if (parallel == TRUE & num_cores == 0) {
    stop("num_cores = 0")
  } else if (parallel == TRUE & num_cores > 0) {
    comb <- function(x, ...) {
      lapply(seq_along(x),
             function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    out_list <- foreach::foreach(G=min_G:max_G, .packages = c("BMCD", "mclust"), .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
      mclust_out <- Mclust(X_est, G = G) # Fit GMM

      sigma_sq <- sigma_sq_out[[p]] # Posterior mean of sigma squared (measurement error) from BMDS algorithm
      eps_init <- mclust_out$parameters$pro # Vector of component mixing proportions
      mu_init <- mclust_out$parameters$mean # Vector of component means
      T_init <- mclust_out$parameters$variance$sigma # List of component covariance matrices
      z_init <- mclust_out$z # Matrix of cluster assignment probabilities

      ## Accounting for one dimension
      mu_init <- matrix(mu_init, nrow = p, ncol = G)
      T_init <- array(T_init, c(p,p,G))

      # Set priors --------------------------------------------------------------

      priors <- setPriors(distances, X_est, mclust_out, p, G, n, m, model_type)

      # Initialize MCMC lists ---------------------------------------------------
      init_params <- list(sigma_sq, eps_init, mu_init, T_init, z_init)
      mcmc_list <- initLists(init_params, X_est, mclust_out$classification, p, G, n, bmcd_iter)

      # Run MCMC ----------------------------------------------------------------
      out <- bmcdMCMC(distances, mcmc_list, priors, p, G, n, m, bmcd_iter, bmcd_burn, labelswitch_iter, model_type)
      out$G <- G
      list(out, priors)
    }
    out <- out_list[[1]]
    priors <- out_list[[2]]
  }

  #out
  list(X_out, out, priors, p)
}
