initLists <- function(paramList, X_est, class_mat, p, G, n, iter) {
  X_list <- lapply(1:iter, matrix, data = NA, nrow = n, ncol = p)
  X_list[[1]] <- X_est

  sigma_sq_list <- rep(NA, iter)
  sigma_sq_list[1] <- paramList[[1]]

  eps_list <- matrix(NA, nrow = iter, ncol = G)
  eps_list[1,] <- paramList[[2]]

  mu_list <- lapply(1:iter, matrix, data = NA, nrow = p, ncol = G)
  mu_list[[1]] <- paramList[[3]]


  n_list <- matrix(NA, nrow = iter, ncol = G)

  for (j in 1:G) {
    n_list[1, j] = sum(class_mat == j)
  }


  T_list <- rep(list(list()), iter)
  for (i in 1:iter) {
    T_list[[i]] <- array(NA, c(p,p,G))
  }
  T_list[[1]] <- paramList[[4]]

  z_list <- lapply(1:iter, matrix, data = NA, nrow = n, ncol = G)
  z_list[[1]] <- paramList[[5]]

  class_list <- matrix(NA, nrow = iter, ncol = n)
  class_list[1, ] <- class_mat

  list(X_list = X_list,
       sigma_sq_list = sigma_sq_list,
       eps_list = eps_list,
       mu_list = mu_list,
       n_list = n_list,
       T_list = T_list,
       z_list = z_list,
       class_list = class_list)
}
