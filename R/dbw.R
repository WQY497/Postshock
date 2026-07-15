#' Donor Balancing Weights (DBW)
#'
#' Compute donor weights on pre-shock data so that a convex combination of
#' donors matches the target under L1/L2 distance. Supports PCA, sum-to-one,
#' box constraints, and L1/L2 regularization on weights.
#'
#' @param X List of data.frame/matrix. First is target (`X[[1L]]`);
#'   the remaining are donors (`X[[2L]]` ... `X[[n+1L]]`). Rows = time,
#'   columns = features.
#' @param dbw_indices Integer vector of feature columns to use.
#' @param shock_time_vec Integer vector, same length as `X`; each series
#'   is truncated to rows `1:shock_i` before matching.
#' @param scale,center Logical; z-score/center features before matching.
#' @param sum_to_1 Logical; if `TRUE`, enforce `sum(W)=1`.
#' @param bounded_below_by,bounded_above_by Numeric scalars; lower/upper bounds
#'   for each weight (broadcast to all weights).
#' @param princ_comp_count `NULL` or positive integer; if not `NULL`,
#'   apply PCA and keep that many principal components.
#' @param normchoice Character; `"l1"` or `"l2"` distance.
#' @param penalty_normchoice Character; `"l1"` (LASSO) or `"l2"` (Ridge)
#'   on weights.
#' @param penalty_lambda Numeric >= 0; regularization strength (0 = off).
#' @param Y Optional list (same length as `X`) of auxiliary features to merge.
#' @param Y_lookback_indices Optional integer vector/list; backward row indices
#'   for Y (e.g., `c(1,2,5)` means last, 2nd last, 5th last).
#' @param X_lookback_indices Optional integer vector/list; backward row indices
#'   for X. If `NULL`, use only the last pre-shock row.
#' @param inputted_transformation Function applied to each `Y[[i]]` before
#'   merging; defaults to `base::identity`.
#'#' @param verbose Logical. If \code{TRUE}, print optimization iteration
#'   information from \code{Rsolnp::solnp()}. Defaults to \code{FALSE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{opt_params}: numeric vector of donor weights (length = #donors)
#'   \item \code{convergence}: \code{"convergence"} or \code{"failed_convergence"}
#'   \item \code{loss}: final matching loss (L1 or L2)
#' }
#'
#' @details
#' With `sum_to_1=TRUE` and nonnegative bounds, \eqn{||W||_1 = 1}; L1
#' regularization on \eqn{||W||_1} is then a constant and does not change the
#' optimum. Use L2 if you want an effective penalty under simplex constraints.
#'
#' @examples
#' \donttest{
#'   set.seed(1)
#'   T <- 50; K <- 3
#'   target <- data.frame(x1=rnorm(T), x2=rnorm(T, 0.2), x3=rnorm(T,-0.1))
#'   donor1 <- data.frame(x1=rnorm(T, 0.1), x2=rnorm(T,-0.1), x3=rnorm(T, 0.3))
#'   donor2 <- data.frame(x1=rnorm(T,-0.2), x2=rnorm(T, 0.4), x3=rnorm(T, 0.0))
#'   X_list <- list(target, donor1, donor2)
#'   shock  <- c(40, 40, 40)
#'
#'   out <- dbw(
#'     X = X_list, dbw_indices = 1:K, shock_time_vec = shock,
#'     center = TRUE, scale = TRUE,
#'     sum_to_1 = TRUE, bounded_below_by = 0, bounded_above_by = 1,
#'     normchoice = "l2",
#'     penalty_normchoice = "l2", penalty_lambda = 1e-2
#'   )
#'   out$opt_params; out$loss; out$convergence
#' }
#'
#' @importFrom stats sd
#' @importFrom utils tail
#' @export

dbw <- function(
    X,
    dbw_indices,
    shock_time_vec,
    scale = FALSE,
    center = FALSE,
    sum_to_1 = TRUE,
    bounded_below_by = 0,
    bounded_above_by = 1,
    princ_comp_count = NULL,
    normchoice = "l2",
    penalty_normchoice = "l1",
    penalty_lambda = 0,
    Y = NULL,
    Y_lookback_indices = NULL,
    X_lookback_indices = NULL,
    inputted_transformation = base::identity,
    verbose = FALSE
) {
  # ---- basic checks ----
  if (!normchoice %in% c("l1","l2"))
    stop("normchoice must be 'l1' or 'l2'.")
  if (!penalty_normchoice %in% c("l1","l2"))
    stop("penalty_normchoice must be 'l1' or 'l2'.")
  if (length(shock_time_vec) != length(X))
    stop("length(shock_time_vec) must equal length(X).")
  if (any(is.na(shock_time_vec)) || any(shock_time_vec < 1))
    stop("shock_time_vec must be integers >= 1 (no NA).")
  for (i in seq_along(X)) {
    if (ncol(X[[i]]) < max(dbw_indices)) {
      stop(sprintf("X[[%d]] has only %d columns but dbw_indices=%s.",
                   i, ncol(X[[i]]), paste(dbw_indices, collapse=",")))
    }
  }
  
  # 1) select columns
  X_subset1 <- lapply(X, function(df) df[, dbw_indices, drop = FALSE])
  
  # 2) optional: merge transformed Y lookbacks with X
  if (!is.null(Y) && !is.null(Y_lookback_indices)) {
    X_Y_combiner <- function(y, x) {
      y_tr <- inputted_transformation(y)
      len_y <- nrow(y_tr)
      idx_y <- len_y - as.numeric(unlist(Y_lookback_indices)) + 1
      idx_y <- idx_y[idx_y >= 1 & idx_y <= len_y]
      if (length(idx_y) == 0) stop("Y_lookback_indices out of range.")
      y_vec <- as.numeric(t(as.matrix(y_tr[idx_y, , drop = FALSE])))
      y_mat <- matrix(rep(y_vec, each = nrow(x)), nrow = nrow(x))
      cbind(y_mat, x)
    }
    combined_X <- mapply(X_Y_combiner, y = Y, x = X_subset1, SIMPLIFY = FALSE)
  } else {
    combined_X <- X_subset1
  }
  
  # 3) truncate to pre-shock rows
  row_returner <- function(df, stv) {
    stv <- min(nrow(df), as.integer(stv))
    if (stv < 1L) stop("A shock_time < 1 was found.")
    df[1:stv, , drop = FALSE]
  }
  X_subset2 <- mapply(row_returner, df = combined_X, stv = shock_time_vec, SIMPLIFY = FALSE)
  
  # 4) extract lookbacks -> one row per series
  cov_extractor <- function(X_df) {
    if (is.null(X_lookback_indices)) {
      vals <- as.numeric(tail(X_df, 1))
      matrix(vals, nrow = 1)
    } else {
      len <- nrow(X_df)
      idx <- len - as.numeric(unlist(X_lookback_indices)) + 1
      idx <- idx[idx >= 1 & idx <= len]
      if (length(idx) == 0) stop("X_lookback_indices out of range after truncation.")
      vals <- as.numeric(t(as.matrix(X_df[idx, , drop = FALSE])))
      matrix(vals, nrow = 1)
    }
  }
  X_rows <- lapply(X_subset2, cov_extractor)
  
  # 5) stack target + donors
  dat <- do.call('rbind', X_rows)           # 1st row = target; others = donors
  if (nrow(dat) <= 1 || ncol(dat) == 0)
    stop("Covariate matrix is empty or has no donors.")
  
  # 6) drop constant columns
  const_cols <- which(apply(dat, 2, function(x) sd(x, na.rm = TRUE) == 0))
  if (length(const_cols) > 0) {
    warning(sprintf("Removed constant columns: %s", paste(const_cols, collapse=",")))
    dat <- dat[, -const_cols, drop = FALSE]
    if (ncol(dat) == 0) stop("All columns were constant; nothing left to match.")
  }
  
  # 7) optional centering/scaling
  if (scale || center) {
    dat <- base::scale(dat, center = center, scale = scale)
    dat <- as.matrix(dat)
  }
  
  # 8) sanity checks
  if (any(is.na(dat)) || any(is.infinite(dat))) {
    stop("[dbw] Found NA/Inf after preprocessing.")
  }
  
  # 9) optional PCA
  if (!is.null(princ_comp_count)) {
    if (princ_comp_count <= 0) stop("princ_comp_count must be positive.")
    if (ncol(dat) < princ_comp_count) stop("princ_comp_count exceeds #features.")
    if (ncol(dat) == 1) stop("PCA is not meaningful with 1 column.")
    sv <- svd(dat)
    dat <- dat %*% sv$v[, 1:princ_comp_count, drop = FALSE]
  }
  
  # 10) split target vs donors
  X1 <- dat[1, , drop = FALSE]
  donors <- dat[-1, , drop = FALSE]
  n <- nrow(donors)
  
  if (n <= 0L) {
    stop("No donors available after preprocessing.")
  }
  
  # With one donor under simplex constraints, the weight is fixed at one.
  # Rank deficiency is therefore irrelevant.
  if (n > 1L && qr(dat)$rank < min(dim(dat))) {
    warning(
      "Design matrix is rank-deficient (collinearity). Consider PCA."
    )
  }
  
  X0 <- lapply(
    seq_len(n),
    function(i) donors[i, , drop = FALSE]
  )
  
  # 11) objective: distance + optional penalty on W
  weightedX0 <- function(W) {
    XW <- Reduce(`+`, Map(function(w, x) w * as.matrix(x), W, X0))  # 1 x k
    diff_vec <- as.numeric(X1 - XW)
    loss <- if (normchoice == "l1") sum(abs(diff_vec)) else sqrt(sum(diff_vec^2))
    if (penalty_lambda > 0) {
      if (penalty_normchoice == 'l1') {
        # With simplex + nonnegativity, sum(abs(W)) == 1 (constant).
        loss <- loss + penalty_lambda * sum(abs(W))
      } else {
        loss <- loss + penalty_lambda * sum(W^2)
      }
    }
    loss
  }
  
  # 12) constraints
  eqfun_obj <- NULL; eqB_obj <- NULL
  if (isTRUE(sum_to_1)) { eqfun_obj <- function(W) sum(W) - 1; eqB_obj <- 0 }
  LB <- if (!is.na(bounded_below_by)) rep(bounded_below_by, n) else NULL
  UB <- if (!is.na(bounded_above_by)) rep(bounded_above_by, n) else NULL
  if (isTRUE(sum_to_1) && !is.null(LB) && !is.null(UB)) {
    if (sum(LB) > 1 + 1e-12 || sum(UB) < 1 - 1e-12)
      stop("Infeasible constraints: cannot satisfy sum(W)=1 under given bounds.")
  }
  
  # Single-donor fast path
  if (n == 1L && isTRUE(sum_to_1)) {
    
    # Under sum-to-one, the only feasible weight is 1.
    if (!is.null(LB) && LB[1L] > 1 + 1e-12) {
      stop(
        "Infeasible constraints: the single donor must have weight 1."
      )
    }
    
    if (!is.null(UB) && UB[1L] < 1 - 1e-12) {
      stop(
        "Infeasible constraints: the single donor must have weight 1."
      )
    }
    
    diff_vec <- as.numeric(
      X1 - donors[1L, , drop = FALSE]
    )
    
    final_loss <- if (normchoice == "l1") {
      sum(abs(diff_vec))
    } else {
      sqrt(sum(diff_vec^2))
    }
    
    return(list(
      opt_params = 1,
      convergence = "convergence",
      loss = as.numeric(round(final_loss, 6))
    ))
  }
  
  # 13) solve
  opt <- Rsolnp::solnp(
    par = rep(1 / n, n),
    fun = weightedX0,
    eqfun = eqfun_obj,
    eqB = eqB_obj,
    LB = LB,
    UB = UB,
    control = list(
      trace = if (isTRUE(verbose)) 1 else 0,
      tol = 1e-8
    )
  )
  
  # 14) final loss & return
  syn <- opt$pars %*% donors
  diff_vec <- as.numeric(X1 - syn)
  final_loss <- if (normchoice == "l1") sum(abs(diff_vec)) else sqrt(sum(diff_vec^2))
  list(
    opt_params  = as.numeric(opt$pars),
    convergence = if (opt$convergence == 0) "convergence" else "failed_convergence",
    loss        = as.numeric(round(final_loss, 6))
  )
}
