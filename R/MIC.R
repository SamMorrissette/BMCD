#' @importFrom gtools ddirichlet
#' @importFrom mvtnorm dmvnorm
#' @importFrom LaplacesDemon dinvwishart
#'

MIC <- function(distances, bmcd_MCMC_list, priors) {
  mcmc <- bmcd_MCMC_list

  eps_star <- colMeans(mcmc$eps_list)
  mu_star <- Reduce("+", mcmc$mu_list) / length(mcmc$mu_list)
  T_star <- Reduce("+", mcmc$T_list) / length(mcmc$T_list)

  # Calculate pi(Lambda*) ----------------------------------------------------

  ## Epsilon
  pr_eps <- gtools::ddirichlet(eps_star, rep(1, length(eps_star)))

  ## Mu_j
  pr_mu <- rep(NA, G)
  for (i in 1:G) {
    pr_mu[i] <- mvtnorm::dmvnorm(mu_star[,i], mean = priors$prior_mean[,i], sigma = T_star[,,i])
  }

  ## T_j
  pr_T <- rep(NA, G)
  for (i in 1:G) {
    pr_T[i] <- LaplacesDemon::dinvwishart(T_star[,,i], nu = priors$prior_alpha, S = priors$prior_Bj[,,i], log = FALSE) # Change to log = TRUE?
  }

  pi_Lambda <- pr_eps * pr_mu * pr_T



  # Calculate pi(X_pg | Lambda*) --------------------------------------------

  X <- Reduce("+", mcmc$X_list) / length(mcmc$X_list)
  n <- nrow(X)
  m <- n*(n-1) / 2

  x_densities <- rep(NA, n)
  for (i in 1:n) {
    total = 0
    for (comp in 1:G) {
      total = total + (eps_star[comp] * mvtnorm::dmvnorm(X[i,], mean = mu_star[,comp], sigma = T_star[,,comp]))
    }
    x_densities[i] <- total
  }


  # Calculate pi(Lambda* | X_pg) --------------------------------------------

  #n_star <- colMeans(mcmc$n_list)
  z_star <- Reduce("+", mcmc$z_list) / length(mcmc$z_list)
  K_star <- apply(z_star, 1, which.max)
  n_star <- table(K_star)

  ## Calculate pi(eps* | K)

  pi_eps <- gtools::ddirichlet(eps_star, n_star + 1)

  ## Calculate the product
  pi_mu <- rep(NA, G)
  pi_T <- rep(NA, G)

  for (i in 1:G) {
    ## Calculate pi(mu_j* | others)
    x_j = X[which(K_star==i), ]
    x_j_bar = colMeans(x_j)
    comp_mean = ((n_star[i] * x_j_bar) + priors$prior_mean[,i]) / (n_star[i] + 1)
    comp_sigma = T_star[,,i] / (n_star[i] + 1)
    pi_mu[i] <- mvtnorm::dmvnorm(mu_star[,i], mean = comp_mean, sigma = comp_sigma)

    ## Calculate S_j
    centered_x <- sweep(x_j, 2, mu_star[,i])
    S_j <- 0
    for (q in 1:nrow(centered_x)) {
      S_j = S_j + (centered_x[q, ] %*% t(centered_x[q,]))
    }

    ## Calculate pi(T_j* | others)
    pi_T[i] <- LaplacesDemon::dinvwishart(T_star[,,i], nu = priors$prior_alpha + (n_star[i]/2), priors$prior_Bj[,,i] + (S_j / 2))
  }
  expectation <- mean(pi_eps * prod(pi_mu * pi_T))

  SSR <- sum((as.matrix(dist(X)) - distances)^2) / 2

  MIC <- ((m-2) * SSR) - (2 * log(prod(x_densities) * prod(pi_Lambda) / expectation))
}
