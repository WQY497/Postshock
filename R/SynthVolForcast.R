#' Synthetic Volatility Forecast with Donor Shock Adjustment
#'
#' @description
#' Produces k-step-ahead variance forecasts for a target series using a
#' pre-shock GARCH model, then adjusts that forecast by a donor-based
#' post-shock effect. Each donor is fitted with a GARCH-X model that includes
#' a post-shock indicator. The donor-specific post-shock indicator coefficients
#' are combined using donor balancing weights (DBW) to obtain the final
#' volatility adjustment.
#'
#' @param Y_series_list List of numeric series. The first element is the target
#'   series, and the remaining elements are donor series.
#' @param covariates_series_list List of data frames or matrices containing
#'   donor covariates. Its length must equal the number of donors.
#' @param target_covariates Data frame or matrix containing target covariates.
#'   These covariates are used for DBW matching. By default, the target-side
#'   GARCH model is fitted without exogenous regressors.
#' @param shock_time_vec Integer vector giving the shock index for each series.
#'   Its length must equal \code{length(Y_series_list)}.
#' @param shock_length_vec Integer vector giving the number of post-shock rows
#'   for which the donor post-shock indicator equals one.
#' @param k Integer forecast horizon.
#' @param dbw_scale Logical. Whether to scale DBW features.
#' @param dbw_center Logical. Whether to center DBW features.
#' @param dbw_indices Integer vector specifying which covariate columns are used
#'   by \code{\link{dbw}}. If \code{NULL}, all columns of
#'   \code{target_covariates} are used.
#' @param princ_comp_input Optional integer specifying the number of principal
#'   components used inside \code{\link{dbw}}. If \code{NULL}, PCA is not used.
#' @param covariate_indices Integer vector specifying which donor covariate
#'   columns are included in donor-side GARCH-X models. If \code{NULL}, donors
#'   use only the post-shock indicator.
#' @param penalty_lambda Numeric penalty weight forwarded to \code{\link{dbw}}.
#' @param penalty_normchoice Character string, either \code{"l1"} or
#'   \code{"l2"}, forwarded to \code{\link{dbw}}.
#' @param plots Logical. If \code{TRUE} and a function named
#'   \code{plot_maker_garch} exists, it is called to visualize the result.
#' @param shock_time_labels Optional labels passed to the plotting helper.
#' @param return_fits Logical. If \code{TRUE}, the fitted donor and target
#'   GARCH objects are returned.
#' @param garch_order Optional fixed GARCH-X order \code{c(q, p, o)}. If
#'   \code{NULL}, the order is selected by \code{\link{auto_garchx}}.
#' @param max.p Nonnegative integer. Maximum GARCH lag order passed to
#'   \code{\link{auto_garchx}}.
#' @param max.q Nonnegative integer. Maximum ARCH lag order passed to
#'   \code{\link{auto_garchx}}.
#' @param backcast.initial Optional positive numeric scalar used as
#'   \code{backcast.values} when \code{garch_order} is fixed.
#' @param vcov.type Character string. Either \code{"ordinary"} or
#'   \code{"robust"}, passed to \code{garchx()} and
#'   \code{lmtest::coeftest()}.
#'
#' @details
#' The target-side model is fitted only on the pre-shock target series and does
#' not use exogenous regressors. This avoids requiring future target-side
#' covariates at forecast time. Donor-side models are GARCH-X specifications
#' that include a post-shock indicator and, optionally, selected donor
#' covariates. The coefficient on the post-shock indicator is interpreted as the
#' donor-specific short-horizon volatility effect.
#'
#' The final adjusted forecast is computed as
#' \deqn{
#'   \hat{\sigma}^{2, adj}_{0,T+h}
#'   =
#'   \max\{\hat{\sigma}^{2, raw}_{0,T+h} + \hat{\omega}, \epsilon\},
#' }
#' where \eqn{\hat{\omega}} is the DBW-weighted donor shock effect and
#' \eqn{\epsilon = 10^{-8}} enforces a strictly positive variance forecast.
#'
#' Forecast accuracy is evaluated outside this function using realized variance
#' proxies, such as next-period squared returns.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{linear_combinations}: numeric vector of donor weights.
#'   \item \code{predictions}: list with \code{unadjusted}, \code{adjusted},
#'         and \code{arithmetic_mean} variance forecasts.
#'   \item \code{meta}: diagnostic information, including donor weights,
#'         donor-specific post-shock effects, standard errors, selected orders,
#'         and the combined shock effect.
#'   \item \code{donor_fits}: returned only if \code{return_fits = TRUE}.
#'   \item \code{target_fit}: returned only if \code{return_fits = TRUE}, or in
#'         the no-donor fast path.
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y_tgt <- rnorm(200)
#' y_d1  <- rnorm(200)
#' y_d2  <- rnorm(200)
#'
#' Y <- list(y_tgt, y_d1, y_d2)
#'
#' Xt  <- cbind(x1 = rnorm(200), x2 = rnorm(200))
#' Xd1 <- cbind(x1 = rnorm(200), x2 = rnorm(200))
#' Xd2 <- cbind(x1 = rnorm(200), x2 = rnorm(200))
#'
#' out <- SynthVolForecast(
#'   Y_series_list = Y,
#'   covariates_series_list = list(Xd1, Xd2),
#'   target_covariates = Xt,
#'   shock_time_vec = c(150L, 150L, 150L),
#'   shock_length_vec = c(10L, 10L, 10L),
#'   k = 2,
#'   covariate_indices = 1:2,
#'   dbw_indices = 1:2,
#'   plots = FALSE
#' )
#'
#' str(out$predictions)
#' }
#'
#' @importFrom stats coef complete.cases fitted predict vcov
#' @importFrom lmtest coeftest
#' @importFrom garchx garchx
#' @export
SynthVolForecast <- function(
    Y_series_list,
    covariates_series_list,
    target_covariates,
    shock_time_vec,
    shock_length_vec,
    k = 1,
    dbw_scale = TRUE,
    dbw_center = TRUE,
    dbw_indices = NULL,
    princ_comp_input = NULL,
    covariate_indices = NULL,
    penalty_lambda = 0,
    penalty_normchoice = "l1",
    plots = TRUE,
    shock_time_labels = NULL,
    return_fits = FALSE,
    garch_order = NULL,
    max.p = 3,
    max.q = 3,
    backcast.initial = NULL,
    vcov.type = "robust"
) {
  # ---- sanitize args ---------------------------------------------------------
  k <- max(1L, as.integer(k))
  vc <- match.arg(tolower(as.character(vcov.type)), c("ordinary", "robust"))
  
  n_total <- length(Y_series_list)
  n_donors <- n_total - 1L
  
  if (n_total < 1L) {
    stop("SynthVolForecast: 'Y_series_list' must contain at least the target series.")
  }
  
  if (length(shock_time_vec) != n_total) {
    stop("SynthVolForecast: length(shock_time_vec) must equal length(Y_series_list).")
  }
  
  if (length(shock_length_vec) != n_total) {
    stop("SynthVolForecast: length(shock_length_vec) must equal length(Y_series_list).")
  }
  
  if (!is.null(backcast.initial)) {
    stopifnot(
      is.numeric(backcast.initial),
      length(backcast.initial) == 1L,
      is.finite(backcast.initial),
      backcast.initial > 0
    )
  }
  
  if (!is.null(princ_comp_input)) {
    if (!is.numeric(princ_comp_input) || length(princ_comp_input) != 1L ||
        !is.finite(princ_comp_input) || princ_comp_input < 1) {
      stop("SynthVolForecast: 'princ_comp_input' must be NULL or a positive integer.")
    }
    princ_comp_input <- as.integer(princ_comp_input)
  }
  
  # ---- helper to coerce covariates ------------------------------------------
  as_df <- function(x) {
    if (is.numeric(x) && is.null(dim(x))) {
      data.frame(V1 = x)
    } else if (is.matrix(x)) {
      as.data.frame(x)
    } else if (inherits(x, "data.frame")) {
      x
    } else {
      stop("SynthVolForecast: covariates must be numeric, matrix, or data.frame.")
    }
  }
  
  # ---- FAST PATH: no donors --------------------------------------------------
  if (n_donors == 0L) {
    target_Y <- as.numeric(Y_series_list[[1L]])
    target_end <- as.integer(shock_time_vec[1L])
    
    if (target_end < 1L || target_end > length(target_Y)) {
      stop("SynthVolForecast: invalid target shock index.")
    }
    
    if (is.null(garch_order)) {
      sel_t <- auto_garchx(
        y = target_Y[seq_len(target_end)],
        xreg = NULL,
        max.p = max.p,
        max.q = max.q,
        vcov.type = vc
      )
      
      target_order <- c(sel_t$q, sel_t$p, sel_t$o)
      fit_target <- sel_t$fit
    } else {
      target_order <- garch_order
      
      fit_target <- garchx(
        y = target_Y[seq_len(target_end)],
        order = target_order,
        xreg = NULL,
        vcov.type = vc,
        backcast.values = if (is.null(backcast.initial)) NULL else backcast.initial
      )
    }
    
    raw_pred <- predict(fit_target, n.ahead = k)
    raw_vec <- as.numeric(if (is.list(raw_pred)) raw_pred$pred else raw_pred)
    
    eps <- 1e-8
    raw_vec <- pmax(raw_vec, eps)
    
    return(list(
      linear_combinations = numeric(0),
      predictions = list(
        unadjusted = raw_vec,
        adjusted = raw_vec,
        arithmetic_mean = raw_vec
      ),
      meta = list(
        n_donors = 0L,
        weights = numeric(0),
        omega_vec = numeric(0),
        omega_se = numeric(0),
        combined_omega = 0,
        shock_time = shock_time_vec,
        shock_length = shock_length_vec,
        target_order = target_order
      ),
      target_fit = fit_target
    ))
  }
  
  # ---- donor path checks -----------------------------------------------------
  if (length(covariates_series_list) != n_donors) {
    stop("SynthVolForecast: length(covariates_series_list) must equal the number of donors.")
  }
  
  target_covariates <- as_df(target_covariates)
  
  for (i in seq_len(n_donors)) {
    covariates_series_list[[i]] <- as_df(covariates_series_list[[i]])
  }
  
  if (is.null(dbw_indices)) {
    dbw_indices <- seq_len(ncol(target_covariates))
  }
  
  has_lmtest <- requireNamespace("lmtest", quietly = TRUE)
  
  # ---- per-donor fit and omega extraction -----------------------------------
  omega_star_hat_vec <- numeric(n_donors)
  omega_star_std_err_vec <- rep(NA_real_, n_donors)
  donor_fits <- vector("list", n_donors)
  donor_selected_orders <- vector("list", n_donors)
  
  for (i in seq_len(n_donors)) {
    donor_Y <- as.numeric(Y_series_list[[i + 1L]])
    donor_X <- covariates_series_list[[i]]
    
    len_i <- as.integer(shock_length_vec[i + 1L])
    start_i <- as.integer(shock_time_vec[i + 1L])
    end_i <- min(start_i + len_i, length(donor_Y))
    
    if (start_i < 1L || start_i > length(donor_Y)) {
      stop(sprintf("SynthVolForecast: invalid shock index for donor %d.", i))
    }
    
    post_shock_indicator <- rep(0L, length(donor_Y))
    
    if (len_i > 0L && end_i >= start_i + 1L) {
      post_shock_indicator[(start_i + 1L):end_i] <- 1L
    }
    
    rows_to_use <- length(donor_Y)
    
    if (is.null(covariate_indices)) {
      X_i_final <- matrix(
        post_shock_indicator[seq_len(rows_to_use)],
        ncol = 1,
        dimnames = list(NULL, "post_shock_indicator")
      )
      donor_y_final <- donor_Y[seq_len(rows_to_use)]
    } else {
      X_cov <- as.matrix(donor_X[seq_len(rows_to_use), covariate_indices, drop = FALSE])
      X_i_final <- cbind(X_cov, post_shock_indicator[seq_len(rows_to_use)])
      colnames(X_i_final)[ncol(X_i_final)] <- "post_shock_indicator"
      
      ok <- complete.cases(X_i_final)
      X_i_final <- X_i_final[ok, , drop = FALSE]
      donor_y_final <- donor_Y[seq_len(rows_to_use)][ok]
    }
    
    if (is.null(garch_order)) {
      sel <- auto_garchx(
        y = donor_y_final,
        xreg = X_i_final,
        max.p = max.p,
        max.q = max.q,
        vcov.type = vc
      )
      
      order_i <- c(sel$q, sel$p, sel$o)
      donor_fit <- sel$fit
    } else {
      order_i <- garch_order
      
      donor_fit <- garchx(
        y = donor_y_final,
        order = order_i,
        xreg = X_i_final,
        vcov.type = vc,
        backcast.values = if (is.null(backcast.initial)) NULL else backcast.initial
      )
    }
    
    donor_selected_orders[[i]] <- order_i
    donor_fits[[i]] <- donor_fit
    
    # omega = coefficient on the post-shock indicator
    if (has_lmtest) {
      vcov_fun <- function(obj) vcov(obj, vcov.type = vc)
      ct <- coeftest(donor_fit, vcov. = vcov_fun)
      
      hit <- grep("post_shock_indicator", rownames(ct), fixed = TRUE)
      
      if (!length(hit)) {
        stop(sprintf(
          "SynthVolForecast: post_shock_indicator coefficient not found for donor %d.",
          i
        ))
      }
      
      omega_star_hat_vec[i] <- ct[hit[1L], "Estimate"]
      
      if ("Std. Error" %in% colnames(ct)) {
        omega_star_std_err_vec[i] <- ct[hit[1L], "Std. Error"]
      }
    } else {
      cf <- coef(donor_fit)
      hit <- grep("post_shock_indicator", names(cf), fixed = TRUE)
      
      if (!length(hit)) {
        stop(sprintf(
          "SynthVolForecast: post_shock_indicator coefficient not found for donor %d.",
          i
        ))
      }
      
      omega_star_hat_vec[i] <- cf[hit[1L]]
    }
  }
  
  # ---- DBW weights (target + donors) ----------------------------------------
  X_for_dbw <- c(list(target_covariates), covariates_series_list)
  
  dbw_output <- dbw(
    X = X_for_dbw,
    dbw_indices = dbw_indices,
    shock_time_vec = shock_time_vec,
    scale = dbw_scale,
    center = dbw_center,
    sum_to_1 = TRUE,
    bounded_below_by = 0,
    bounded_above_by = 1,
    princ_comp_count = princ_comp_input,
    normchoice = "l2",
    penalty_normchoice = penalty_normchoice,
    penalty_lambda = penalty_lambda
  )
  
  w_hat <- as.numeric(dbw_output$opt_params)
  
  if (length(w_hat) != length(omega_star_hat_vec)) {
    stop("SynthVolForecast: DBW weight length does not match number of donor effects.")
  }
  
  omega_star_hat <- sum(w_hat * omega_star_hat_vec)
  
  # ---- target pre-shock model and forecast ----------------------------------
  target_Y <- as.numeric(Y_series_list[[1L]])
  target_end <- as.integer(shock_time_vec[1L])
  
  if (target_end < 1L || target_end > length(target_Y)) {
    stop("SynthVolForecast: invalid target shock index.")
  }
  
  if (is.null(garch_order)) {
    sel_t <- auto_garchx(
      y = target_Y[seq_len(target_end)],
      xreg = NULL,
      max.p = max.p,
      max.q = max.q,
      vcov.type = vc
    )
    
    target_order <- c(sel_t$q, sel_t$p, sel_t$o)
    fit_target <- sel_t$fit
  } else {
    target_order <- garch_order
    
    fit_target <- garchx(
      y = target_Y[seq_len(target_end)],
      order = target_order,
      xreg = NULL,
      vcov.type = vc,
      backcast.values = if (is.null(backcast.initial)) NULL else backcast.initial
    )
  }
  
  raw_pred <- predict(fit_target, n.ahead = k)
  raw_vec <- as.numeric(if (is.list(raw_pred)) raw_pred$pred else raw_pred)
  
  eps <- 1e-8
  raw_vec <- pmax(raw_vec, eps)
  adjusted_vec <- pmax(raw_vec + omega_star_hat, eps)
  arithmetic_vec <- pmax(raw_vec + mean(omega_star_hat_vec), eps)
  
  predictions <- list(
    unadjusted = raw_vec,
    adjusted = adjusted_vec,
    arithmetic_mean = arithmetic_vec
  )
  
  meta <- list(
    n_donors = n_donors,
    shock_time = shock_time_vec,
    shock_length = shock_length_vec,
    dbw_status = dbw_output$convergence,
    weights = w_hat,
    omega_vec = omega_star_hat_vec,
    omega_se = omega_star_std_err_vec,
    combined_omega = omega_star_hat,
    donor_orders = donor_selected_orders,
    target_order = target_order
  )
  
  out <- list(
    linear_combinations = w_hat,
    predictions = predictions,
    meta = meta
  )
  
  if (isTRUE(return_fits)) {
    out$donor_fits <- donor_fits
    out$target_fit <- fit_target
  }
  
  # ---- optional plot hook ----------------------------------------------------
  if (isTRUE(plots)) {
    f <- get0("plot_maker_garch", mode = "function", inherits = TRUE)
    
    if (!is.null(f)) {
      try(
        f(
          fitted(fit_target),
          shock_time_labels,
          shock_time_vec,
          shock_length_vec,
          raw_vec,
          w_hat,
          omega_star_hat_vec,
          omega_star_std_err_vec,
          adjusted_vec,
          arithmetic_vec
        ),
        silent = TRUE
      )
    }
  }
  
  out
}