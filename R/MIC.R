#' @importFrom gtools ddirichlet
#' @importFrom mvtnorm dmvnorm
#' @importFrom LaplacesDemon dinvwishart
#'

MIC <- function(distances, X_out, bmcd_MCMC_list, priors, min_G, max_G) {
  mics <-rep(NA, max_G - min_G + 1)
  bics <- rep(NA, max_G - min_G + 1)
  aics <- rep(NA, max_G - min_G + 1)
  mod_aics <- rep(NA, max_G - min_G + 1)
  optim_params <- rep(list(list()), (max_G - min_G + 1))
  n <- nrow(distances)
  m <- n*(n-1) / 2
  p <- dim(bmcd_MCMC_list[[1]]$T_list[[1]][,,1])[1]
  print(paste("p:",p))

  # Calculate last term -----------------------------------------------------
  # A <- 0
  # for (q in 1:p) {
  #   X <- X_out[[q]]
  #
  #   # Calculate Hn ------------------------------------------------------------
  #   Hn = -((n + 1) * log(pi)) + (2 * log(gamma((n+1) / 2)))
  #
  #   # Second term -------------------------------------------------------------
  #   centered_X <- sweep(X, 2, colMeans(X))
  #   S = 0
  #   for (j in 1:n) {
  #     S = S + centered_X[j,] %*% t(centered_X[j,])
  #   }
  #
  #   # Third term --------------------------------------------------------------
  #
  #   R <- vector("list", length = q)
  #   for (i in 1:q) {
  #     X_temp <- X_out[[i]]
  #     R[[i]] <- t(X_temp) %*% X_temp
  #   }
  #
  #   if (q >= 2) {
  #     r <- rep(NA, q-1)
  #     total <- 0
  #     for (j in 1:(q-1)) {
  #       r[j] <- diag(R[[q]])[j] / diag(R[[q-1]])[j]
  #       total <- total + log(r[j]) + ((n+1) * log(n+1) / (n+r[j]))
  #     }
  #   } else {
  #     total <- 0
  #   }
  #   A = A + Hn - (n * log((diag(S)[q]) / n)) + total
  # }

  for (index in 1:(max_G - min_G + 1)) {
    mcmc <- bmcd_MCMC_list[[index]]
    G <- mcmc$G
    prior_G <- priors[[index]]

    eps_star <- colMeans(mcmc$eps_list)
    mu_star <- Reduce("+", mcmc$mu_list) / length(mcmc$mu_list)
    T_star <- Reduce("+", mcmc$T_list) / length(mcmc$T_list)
    X <- Reduce("+", mcmc$X_list) / length(mcmc$X_list)
    SSR <- sum((as.matrix(dist(X)) - distances)^2) / 2


    num_params <- G*(p+p+1)#G*((p*p - p) / 2  + (2*p) + 1)

    aics_iter <- rep(NA, nrow(mcmc$eps_list))
    for (iter in 1:nrow(mcmc$eps_list)) {
      xs <- c()
      for (i in 1:n) {
        total <- 0
        for (j in 1:G) {
          total <- total + mcmc$eps_list[iter, j] * dmvnorm(X[i,], mcmc$mu_list[[iter]][,j], mcmc$T_list[[iter]][,,j])
        }
        xs <- c(xs, total)
      }
      aics_iter[iter] <- 2*num_params - (2*sum(log(xs)))
    }
    aics[index] <- min(aics_iter)
    optim_params[[index]] <- list(eps = mcmc$eps_list[which.min(aics_iter), ],
                                mu = mcmc$mu_list[[which.min(aics_iter)]],
                                S = mcmc$T_list[[which.min(aics_iter)]],
                                z = mcmc$z_list[[which.min(aics_iter)]])


    print(index)



    # # Calculate pi(Lambda*) ----------------------------------------------------
    #
    # ## Epsilon
    # pr_eps <- gtools::ddirichlet(eps_star, rep(1, length(eps_star)))
    #
    # ## Mu_j
    # pr_mu <- rep(NA, G)
    # for (i in 1:G) {
    #   pr_mu[i] <- mvtnorm::dmvnorm(mu_star[,i], mean = prior_G$prior_mean[,i], sigma = matrix(T_star[,,i], nrow = p, ncol = p))
    # }
    #
    # ## T_j
    # pr_T <- rep(NA, G)
    # for (i in 1:G) {
    #   pr_T[i] <- LaplacesDemon::dinvwishart(T_star[,,i],
    #                                          nu = prior_G$prior_alpha,
    #                                          S = prior_G$prior_Bj[,,i])
    #   # tryCatch({
    #   #   pr_T[i] <<- LaplacesDemon::dinvwishart(T_star[,,i],
    #   #                                         nu = prior_G$prior_alpha,
    #   #                                         S = prior_G$prior_Bj[,,i]) # Change to log = TRUE?
    #   # }, error = function(e) {
    #   #   diag(prior_G$prior_Bj[,,i]) <- diag(prior_G$prior_Bj[,,i]) + 1e-05
    #   #   pr_T[i] <<- LaplacesDemon::dinvwishart(T_star[,,i],
    #   #                                         nu = prior_G$prior_alpha,
    #   #                                         S = prior_G$prior_Bj[,,i]) # Change to log = TRUE?
    #   # }, warning = function(w) {
    #   #   diag(prior_G$prior_Bj[,,i]) <- diag(prior_G$prior_Bj[,,i]) + 1e-05
    #   #   pr_T[i] <<- LaplacesDemon::dinvwishart(T_star[,,i],
    #   #                                          nu = prior_G$prior_alpha,
    #   #                                          S = prior_G$prior_Bj[,,i]) # Change to log = TRUE?
    #   # })
    # }
    #
    # pi_Lambda <- pr_eps * pr_mu * pr_T
    #
    # # Calculate pi(X_pG | Lambda*) --------------------------------------------
    #
    #
    #
    # x_densities <- rep(NA, n)
    # for (i in 1:n) {
    #   total = 0
    #   for (comp in 1:G) {
    #     total = total + (eps_star[comp] * mvtnorm::dmvnorm(X[i,],
    #                                                        mean = mu_star[,comp],
    #                                                        sigma = matrix(T_star[,,comp], nrow = p, ncol = p),
    #                                                        log = FALSE))
    #   }
    #   x_densities[i] <- total
    # }
    #
    #
    # # Calculate pi(Lambda* | X_pG) --------------------------------------------
    #
    # ## Calculate pi(eps* | K)
    #
    # product <- rep(NA, nrow(mcmc$n_list))
    # for (i in 1:nrow(mcmc$n_list)) {
    #   pi_eps <- gtools::ddirichlet(eps_star, mcmc$n_list[i,] + 1)
    #   pi_mu <- pi_T <- rep(NA, G)
    #   for (j in 1:G) {
    #     ## Calculate pi(mu_j* | others)
    #     x_j <- X[which(mcmc$class_list[i,] == j), , drop = FALSE]
    #     if (nrow(x_j) >= 1) {
    #       x_j_bar = colMeans(x_j)
    #       comp_mean = ((mcmc$n_list[i,j] * x_j_bar) + prior_G$prior_mean[,j]) / (mcmc$n_list[i,j] + 1)
    #       comp_sigma = matrix(T_star[,,j] / (mcmc$n_list[i,j] + 1), ncol = p, nrow = p)
    #       pi_mu[j] <- mvtnorm::dmvnorm(mu_star[,j],
    #                                    mean = comp_mean,
    #                                    sigma = comp_sigma)
    #       ## Calculate pi(T_j* | others)
    #
    #       ### Calculate S_j
    #       centered_x <- sweep(x_j, 2, mcmc$mu_list[[i]][,j])
    #       S_j <- 0
    #       for (q in 1:nrow(centered_x)) {
    #         S_j = S_j + (centered_x[q, ] %*% t(centered_x[q,]))
    #       }
    #
    #       Sigma <- prior_G$prior_Bj[,,j] + (S_j / 2)
    #       pi_T[j] <- LaplacesDemon::dinvwishart(T_star[,,j],
    #                                              nu = prior_G$prior_alpha + (mcmc$n_list[i,j]/2),
    #                                              S = Sigma)
    #       # tryCatch({
    #       #   pi_T[j] <<- LaplacesDemon::dinvwishart(T_star[,,j],
    #       #                                          nu = prior_G$prior_alpha + (mcmc$n_list[i,j]/2),
    #       #                                          S = Sigma)
    #       # }, error = function(e) {
    #       #   diag(Sigma) <- diag(Sigma) + 1e-05
    #       #   pi_T[j] <<- LaplacesDemon::dinvwishart(T_star[,,j],
    #       #                                          nu = prior_G$prior_alpha + (mcmc$n_list[i,j]/2),
    #       #                                          S = Sigma) # Change to log = TRUE?
    #       # }, warning = function(w) {
    #       #   diag(Sigma) <- diag(Sigma) + 1e-05
    #       #   pi_T[j] <<- LaplacesDemon::dinvwishart(T_star[,,j],
    #       #                                          nu = prior_G$prior_alpha + (mcmc$n_list[i,j]/2),
    #       #                                          S = Sigma) # Change to log = TRUE?
    #       # })
    #     } else {
    #       pi_mu[j] <- mvtnorm::dmvnorm(mu_star[,j],
    #                                    mean = prior_G$prior_mean[,j],
    #                                    sigma = matrix(T_star[,,j], ncol = p, nrow = p))
    #       pi_T[j] <- LaplacesDemon::dinvwishart(T_star[,,j],
    #                                              nu = prior_G$prior_alpha,
    #                                              S = prior_G$prior_Bj[,,j])
    #
    #       # tryCatch({
    #       #   pi_T[j] <<- LaplacesDemon::dinvwishart(T_star[,,j],
    #       #                                          nu = prior_G$prior_alpha,
    #       #                                          S = prior_G$prior_Bj[,,j]) # Change to log = TRUE?
    #       # }, error = function(e) {
    #       #   diag(prior_G$prior_Bj[,,j]) <- diag(prior_G$prior_Bj[,,j]) + 1e-05
    #       #   pi_T[j] <<- LaplacesDemon::dinvwishart(T_star[,,j],
    #       #                                          nu = prior_G$prior_alpha,
    #       #                                          S = prior_G$prior_Bj[,,j]) # Change to log = TRUE?
    #       # }, warning = function(w) {
    #       #   diag(prior_G$prior_Bj[,,j]) <- diag(prior_G$prior_Bj[,,j]) + 1e-05
    #       #   pi_T[j] <<- LaplacesDemon::dinvwishart(T_star[,,j],
    #       #                                          nu = prior_G$prior_alpha,
    #       #                                          S = prior_G$prior_Bj[,,j]) # Change to log = TRUE?
    #       # })
    #     }
    #   }
    #   product[i] <- pi_eps * prod(pi_mu * pi_T)
    # }
    # expectation <- mean(product)
    #
    #
    # # print(sum(log(x_densities)))
    # # print(sum(log(pi_Lambda)))
    # # print(log(expectation))
    # # print(paste("SSR:", (m-2) * log(SSR)))
    # # print(paste("X_densities:", prod((x_densities))))
    # # print(SSR)
    # #  print(prod(x_densities))
    # #  print(prod(pi_Lambda))
    # # print(paste("Numerator:", prod(x_densities) * prod(pi_Lambda)))
    # # print(paste("Denominator:", expectation))
    # # print(paste("X_densities:", sum(log(x_densities))))
    # # print(paste("pi_Lambda:", sum(log(pi_Lambda))))
    # # print(paste("Expectation:", log(expectation)))
    # # print(-(2 * (sum(log(x_densities)) + sum(log(pi_Lambda)) - log(expectation))))
    # #print(A)
    #
    # if (p == 1) {
    #   mics[index] <- ((m-2) * log(SSR)) - (2 * (sum(log(x_densities)) + sum(log(pi_Lambda)) - log(expectation)))
    # } else if (p > 1) {
    #   mics[index] <- ((m-2) * log(SSR)) - (2 * (sum(log(x_densities)) + sum(log(pi_Lambda)) - log(expectation))) + A
    # }

    # ##BIC TESTING
    #
    # num_params <- G*(p+p+1)#G*((p*p - p) / 2  + (2*p) + 1)
    # xs <- c()
    # for (i in 1:n) {
    #   total <- 0
    #   for (j in 1:G) {
    #     total <- total + eps_star[j]*dmvnorm(X[i,], mu_star[,j], T_star[,,j])
    #   }
    #   xs <- c(xs, total)
    # }
    # bics[index] <- -2 * (sum(log(xs))) + (log(n) * num_params)
    # #bics[index] <- ((m-2) * log(SSR) - 2 * (sum(log(xs))) + (log(n) * num_params))
    # aics[index] <- 2*num_params - (2*sum(log(xs)))
    # mod_aics[index] <- (m-2) * log(SSR) -(2*num_params - (2*sum(log(xs))))
  }
  # print(paste("Bics:", bics))
  # print(which.min(bics))
  # print(paste("MICS:", MICS))
  # print(which.min(mics))
  # print(paste("mod_aics:", mod_aics))
  # print(which.min(mod_aics))
  return(list(aics = aics,
              optim_params = optim_params))

}
