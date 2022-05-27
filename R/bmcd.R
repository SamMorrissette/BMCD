#' @importFrom mclust Mclust mclustBIC
#' @export

bmcd <- function(distances, p, G, bmds_iter, bmds_burn, bmcd_iter, bmcd_burn, labelswitch_iter) {
  n <- nrow(distances)
  m <- n*(n-1) / 2

  # Obtain an initial guess for X using BMDS algorithm ----------------------

  bmds_out <- edited_bmdsMCMC(DIST = distances, p = p, nwarm = bmds_burn, niter = bmds_iter) # Doesn't include delta matrix (reduce memory required)
  X_est <- bmds_out$x_bmds

  # Set initial parameter values using GMM  ---------------------------------

  mclust_out <- Mclust(X_est, G = G) # Fit GMM

  sigma_sq <- bmds_out$e_sigma # Posterior mean of sigma squared (measurement error) from BMDS algorithm
  eps_init <- mclust_out$parameters$pro # Vector of component mixing proportions
  mu_init <- mclust_out$parameters$mean # Vector of component means
  T_init <- mclust_out$parameters$variance$sigma # List of component covariance matrices
  z_init <- mclust_out$z # Matrix of cluster assignment probabilities

  ## Accounting for one dimension
  mu_init <- matrix(mu_init, nrow = p, ncol = G)
  T_init <- array(T_init, c(p,p,G))

  # Set priors --------------------------------------------------------------

  priors <- setPriors(distances, X_est, mclust_out$classification, p, G, n, m)


  # Initialize MCMC lists ---------------------------------------------------
  init_params <- list(sigma_sq, eps_init, mu_init, T_init, z_init)
  mcmc_list <- initLists(init_params, X_est, mclust_out$classification, p, G, n, bmcd_iter)


  # Run MCMC ----------------------------------------------------------------
  out <- bmcdMCMC(distances, mcmc_list, priors, p, G, n, m, bmcd_iter, bmcd_burn, labelswitch_iter)


  # Calculate MIC -----------------------------------------------------------
  #calc_MIC <- MIC(distances, out, priors)

}
