#' @importFrom gtools ddirichlet
#' @importFrom mvtnorm dmvnorm
#' @importFrom LaplacesDemon dinvwishart
#'

MIC <- function(distances, X_out, bmcd_MCMC_list, priors, p, min_G, max_G) {
  mics <-rep(NA, max_G - min_G + 1)
  n <- nrow(distances)
  m <- n*(n-1) / 2

  # Calculate last term -----------------------------------------------------
  A <- 0
  for (q in 1:p) {
    X <- X_out[[q]]

    # Calculate Hn ------------------------------------------------------------
    Hn = -((n + 1) * log(pi)) + (2 * log(gamma((n+1) / 2)))

    # Second term -------------------------------------------------------------
    centered_X <- sweep(X, 2, colMeans(X))
    S = 0
    for (j in 1:n) {
      S = S + centered_X[j,] %*% t(centered_X[j,])
    }

    # Third term --------------------------------------------------------------

    R <- vector("list", length = q)
    for (i in 1:q) {
      X_temp <- X_out[[i]]
      R[[i]] <- t(X_temp) %*% X_temp
    }

    if (q >= 2) {
      r <- rep(NA, q-1)
      total <- 0
      for (j in 1:(q-1)) {
        r[j] <- diag(R[[q]])[j] / diag(R[[q-1]])[j]
        total <- total + log(r[j]) + ((n+1) * log(n+1) / (n+r[j]))
      }
    } else {
      total <- 0
    }
    A = Hn - (n * log((diag(S)[q]) / n)) + total
  }

  for (index in 1:(max_G - min_G + 1)) {
    mcmc <- bmcd_MCMC_list[[index]]
    G <- mcmc$G
    prior_G <- priors[[index]]

    eps_star <- colMeans(mcmc$eps_list)
    mu_star <- Reduce("+", mcmc$mu_list) / length(mcmc$mu_list)
    T_star <- Reduce("+", mcmc$T_list) / length(mcmc$T_list)

    # Calculate pi(Lambda*) ----------------------------------------------------

    ## Epsilon
    pr_eps <- gtools::ddirichlet(eps_star, rep(1, length(eps_star)))

    ## Mu_j
    pr_mu <- rep(NA, G)
    for (i in 1:G) {
      pr_mu[i] <- mvtnorm::dmvnorm(mu_star[,i], mean = prior_G$prior_mean[,i], sigma = matrix(T_star[,,i], nrow = p, ncol = p))
    }

    ## T_j
    pr_T <- rep(NA, G)
    for (i in 1:G) {
      tryCatch({
        pr_T[i] <<- LaplacesDemon::dinvwishart(T_star[,,i],
                                              nu = prior_G$prior_alpha,
                                              S = prior_G$prior_Bj[,,i]) # Change to log = TRUE?
      }, error = function(e) {
        diag(prior_G$prior_Bj[,,i]) <- diag(prior_G$prior_Bj[,,i]) + 1e-05
        pr_T[i] <<- LaplacesDemon::dinvwishart(T_star[,,i],
                                              nu = prior_G$prior_alpha,
                                              S = prior_G$prior_Bj[,,i]) # Change to log = TRUE?
      })
    }

    pi_Lambda <- pr_eps * pr_mu * pr_T

    # Calculate pi(X_pG | Lambda*) --------------------------------------------

    X <- Reduce("+", mcmc$X_list) / length(mcmc$X_list)


    x_densities <- rep(NA, n)
    for (i in 1:n) {
      total = 0
      for (comp in 1:G) {
        total = total + (eps_star[comp] * mvtnorm::dmvnorm(X[i,],
                                                           mean = mu_star[,comp],
                                                           sigma = matrix(T_star[,,comp], nrow = p, ncol = p)))
      }
      x_densities[i] <- total
    }


    # Calculate pi(Lambda* | X_pG) --------------------------------------------

    #n_star <- colMeans(mcmc$n_list)
    z_star <- Reduce("+", mcmc$z_list) / length(mcmc$z_list)
    K_star <- apply(z_star, 1, which.max)
    n_star <- table(factor(K_star, levels = 1:G))

    ## Calculate pi(eps* | K)

    pi_eps <- gtools::ddirichlet(eps_star, n_star + 1)

    ## Calculate the product
    pi_mu <- rep(NA, G)
    pi_T <- rep(NA, G)

    for (i in 1:G) {
      ## Calculate pi(mu_j* | others)
      x_j = X[which(K_star==i), , drop = FALSE]
      if (nrow(x_j) >= 1) {
        x_j_bar = colMeans(x_j)
        comp_mean = ((n_star[i] * x_j_bar) + prior_G$prior_mean[,i]) / (n_star[i] + 1)
        comp_sigma = matrix(T_star[,,i] / (n_star[i] + 1), ncol = p, nrow = p)
        pi_mu[i] <- mvtnorm::dmvnorm(mu_star[,i], mean = comp_mean, sigma = comp_sigma)

        ## Calculate S_j
        centered_x <- sweep(x_j, 2, mu_star[,i])
        S_j <- 0
        for (q in 1:nrow(centered_x)) {
          S_j = S_j + (centered_x[q, ] %*% t(centered_x[q,]))
        }

        ## Calculate pi(T_j* | others)
        Sigma <- prior_G$prior_Bj[,,i] + (S_j / 2)

        tryCatch({
          pi_T[i] <<- LaplacesDemon::dinvwishart(T_star[,,i],
                                                nu = prior_G$prior_alpha + (n_star[i]/2),
                                                S = Sigma)
        }, error = function(e) {
          diag(Sigma) <- diag(Sigma) + 1e-05
          pi_T[i] <<- LaplacesDemon::dinvwishart(T_star[,,i],
                                                 nu = prior_G$prior_alpha + (n_star[i]/2),
                                                 S = Sigma) # Change to log = TRUE?
        })
      } else {
        pi_mu[i] <- mvtnorm::dmvnorm(mu_star[,i],
                                     mean = prior_G$prior_mean[,i],
                                     sigma = matrix(T_star[,,i], ncol = p, nrow = p))
        tryCatch({
          pi_T[i] <<- LaplacesDemon::dinvwishart(T_star[,,i],
                                                nu = prior_G$prior_alpha,
                                                S = prior_G$prior_Bj[,,i]) # Change to log = TRUE?
        }, error = function(e) {
          diag(prior_G$prior_Bj[,,i]) <- diag(prior_G$prior_Bj[,,i]) + 1e-05
          pi_T[i] <<- LaplacesDemon::dinvwishart(T_star[,,i],
                                                nu = prior_G$prior_alpha,
                                                S = prior_G$prior_Bj[,,i]) # Change to log = TRUE?
        })
      }
    }
    expectation <- mean(pi_eps * prod(pi_mu * pi_T))

    SSR <- sum((as.matrix(dist(X)) - distances)^2) / 2

    print(paste("SSR:", SSR))
    print(paste("X densities:", sum(log(x_densities))))
    print(paste("pi_Lambda densities:", sum(log(pi_Lambda))))
    print(paste("expectation:", log(expectation)))
    print((2 * (sum(log(x_densities)) + sum(log(pi_Lambda)) - log(expectation))))

    if (p == 1) {
      mics[index] <- ((m-2) * SSR) - (2 * (sum(log(x_densities)) + sum(log(pi_Lambda)) - log(expectation)))
    } else if (p > 1) {
      mics[index] <- ((m-2) * SSR) - (2 * (sum(log(x_densities)) + sum(log(pi_Lambda)) - log(expectation))) + A
    }
  }
  return(mics)
}
