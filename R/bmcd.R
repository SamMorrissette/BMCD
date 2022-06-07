#' @importFrom mclust Mclust mclustBIC
#' @importFrom foreach %dopar% foreach
#' @importFrom iterators icount
#' @importFrom doParallel registerDoParallel
#' @export

bmcd <- function(distances, max_p, min_G, max_G, bmds_iter, bmds_burn, bmcd_iter, bmcd_burn, labelswitch_iter,
                 parallel = FALSE, num_cores = NULL) {

  if(parallel == TRUE & num_cores > 0) {
    doParallel::registerDoParallel(cores=num_cores)
  }
  n <- nrow(distances)
  m <- n*(n-1) / 2

  # Obtain an initial guess for X using BMDS algorithm ----------------------

  X_out <- vector("list", length = max_p)
  sigma_sq_out <- vector("list", length = max_p)

  if (parallel == FALSE) {
    for (i in 1:max_p) {
      # Doesn't include delta matrix (reduce memory required)
      temp_bmds <- edited_bmdsMCMC(DIST = distances, p = i, nwarm = bmds_burn, niter = bmds_iter)
      X_out[[i]] <- temp_bmds$x_bmds
      sigma_sq_out[[i]] <- temp_bmds$e_sigma
      print(i)
    }
  } else if (parallel == TRUE & num_cores == 0) {
    stop("num_cores = 0")
  } else if (parallel == TRUE & num_cores > 0) {
    comb <- function(x, ...) {
      lapply(seq_along(x),
             function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    out_list <- foreach::foreach(j=1:max_p, .packages ="BMCD", .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
      output <- edited_bmdsMCMC(distances, j, nwarm = 1000, niter = 5000)
      list(output$x_bmds, output$e_sigma)
    }
    X_out <- out_list[[1]]
    sigma_sq_out <- out_list[[2]]
  }


  # Calculate MDSIC to choose best p ----------------------------------------
  mdsics <- edited_MDSIC(distances, X_out)
  print(mdsics)
  plot(mdsics)
  p <- which.min(mdsics)
  #p <- max_p
  print(p)

  X_est <- X_out[[p]]

  # Set initial parameter values using GMM  ---------------------------------

  out <- priors <- vector("list", length = max_G - min_G + 1)



  if (parallel == FALSE) {
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

      priors[[ind]] <- setPriors(distances, X_est, mclust_out, p, G, n, m)

      # Initialize MCMC lists ---------------------------------------------------
      init_params <- list(sigma_sq, eps_init, mu_init, T_init, z_init)
      mcmc_list <- initLists(init_params, X_est, mclust_out$classification, p, G, n, bmcd_iter)

      # Run MCMC ----------------------------------------------------------------
      out[[ind]] <- bmcdMCMC(distances, mcmc_list, priors[[ind]], p, G, n, m, bmcd_iter, bmcd_burn, labelswitch_iter)
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

      priors <- setPriors(distances, X_est, mclust_out, p, G, n, m)

      # Initialize MCMC lists ---------------------------------------------------
      init_params <- list(sigma_sq, eps_init, mu_init, T_init, z_init)
      mcmc_list <- initLists(init_params, X_est, mclust_out$classification, p, G, n, bmcd_iter)

      # Run MCMC ----------------------------------------------------------------
      out <- bmcdMCMC(distances, mcmc_list, priors, p, G, n, m, bmcd_iter, bmcd_burn, labelswitch_iter)
      out$G <- G
      list(out, priors)
    }
    out <- out_list[[1]]
    priors <- out_list[[2]]
  }

  #out
  list(X_out, out, priors, p)
}
