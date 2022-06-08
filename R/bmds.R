#' @importFrom foreach %dopar% foreach
#' @importFrom doParallel registerDoParallel
#' @export

bmds <- function(distances, max_p, bmds_iter, bmds_burn, parallel = FALSE, num_cores = NULL) {
  if(parallel == TRUE & num_cores > 0) {
    doParallel::registerDoParallel(cores=num_cores)
  }
  n <- nrow(distances)

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
      output <- edited_bmdsMCMC(distances, j, nwarm = bmds_burn, niter = bmds_iter)
      list(output$x_bmds, output$e_sigma)
    }
    X_out <- out_list[[1]]
    sigma_sq_out <- out_list[[2]]
  }
  mdsics <- edited_MDSIC(distances, X_out)

  return(list(X_out = X_out,
              sigma_sq_out = sigma_sq_out,
              mdsics = mdsics))
}
