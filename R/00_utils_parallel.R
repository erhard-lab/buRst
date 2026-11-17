#' Simple lapply/mclapply switch (internal)
#'
#' @keywords internal
#' @noRd
NULL

# internal helper: choose lapply vs mclapply
.lapply_par <- function(X, FUN, ..., parallel = FALSE, ncores = 2L) {
  if (!isTRUE(parallel)) return(lapply(X, FUN, ...))
  if (.Platform$OS.type == "windows") {
    warning("parallel=TRUE requested but mclapply not supported on Windows; using lapply().")
    return(lapply(X, FUN, ...))
  }
  parallel::mclapply(X, FUN, ..., mc.cores = ncores)
}
