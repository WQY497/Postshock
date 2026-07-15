#' Auto-select a GARCHX specification by BIC
#'
#' Grid-searches over GARCHX order specifications and selects the model with
#' the lowest BIC. The function optionally refits the selected model using a
#' backcast value initialized at the estimated unconditional variance, clamped
#' to a reasonable range relative to the sample variance.
#'
#' @param y Numeric vector of returns or residuals. Missing, infinite, and NaN
#'   values are not allowed.
#' @param xreg Optional matrix or data frame of exogenous regressors. If
#'   provided, it must have the same number of rows as \code{length(y)}.
#' @param max.p Nonnegative integer. Maximum GARCH lag order to search.
#' @param max.q Nonnegative integer. Maximum ARCH lag order to search.
#' @param search_o Integer vector specifying asymmetry orders to search.
#'   Supported values are \code{0} and \code{1}, where \code{0} denotes a
#'   symmetric GARCH model and \code{1} denotes a GJR-type asymmetric model.
#' @param final_refit Logical. If \code{TRUE}, the BIC-selected model is
#'   optionally refit using a backcast value based on the estimated
#'   unconditional variance.
#' @param clamp_factor Length-two numeric vector \code{c(low, high)}. If the
#'   estimated unconditional variance falls outside
#'   \code{[low, high] * var(y)}, the sample variance is used as the backcast.
#' @param vcov.type Character string. Either \code{"ordinary"} or
#'   \code{"robust"}, passed to \code{garchx::garchx()}.
#' @param verbose Logical. If \code{TRUE}, warnings and refit messages are
#'   printed.
#'
#' @details
#' The \code{garchx} package uses order \code{c(q, p, o)}, where \code{q} is the
#' ARCH order, \code{p} is the GARCH order, and \code{o} is the asymmetry order.
#' The trivial specification \code{c(0, 0, 0)} is skipped. Among all converged
#' models, the function selects the specification with the lowest BIC.
#'
#' When \code{final_refit = TRUE}, the selected model is refit using a scalar
#' \code{backcast.values} initialized from the estimated unconditional variance.
#' For a GJR specification, the unconditional variance is approximated as
#' \deqn{
#'   \omega / (1 - \alpha - \beta - 0.5\gamma).
#' }
#' If this value is invalid or outside the clamped range, the sample variance is
#' used instead.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{fit}: the selected \code{garchx} fit object.
#'   \item \code{q}, \code{p}, \code{o}: selected GARCHX orders.
#'   \item \code{bic}, \code{aic}: information criteria for the selected fit.
#'   \item \code{refit_used}: whether the final backcast refit replaced the
#'   original grid-search fit.
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- scale(rnorm(500))[, 1]
#'
#' out <- auto_garchx(
#'   y = y,
#'   max.p = 2,
#'   max.q = 2,
#'   search_o = c(0, 1),
#'   final_refit = TRUE,
#'   vcov.type = "robust",
#'   verbose = FALSE
#' )
#'
#' c(out$q, out$p, out$o)
#' out$bic
#' }
#'
#' @importFrom stats logLik coef var
#' @export
auto_garchx <- function(
    y,
    xreg         = NULL,
    max.p        = 3,
    max.q        = 3,
    search_o     = c(0L, 1L),
    final_refit  = TRUE,
    clamp_factor = c(0.1, 10),
    vcov.type    = c("ordinary", "robust"),
    verbose      = TRUE
) {
  vcov.type <- match.arg(vcov.type)
  
  # ---- Basic input checks ----------------------------------------------------
  if (!is.numeric(y) || !is.vector(y)) {
    stop("auto_garchx: 'y' must be a numeric vector.")
  }
  if (any(!is.finite(y))) {
    stop("auto_garchx: 'y' contains NA/Inf/NaN.")
  }
  
  n <- length(y)
  if (n < 20 && isTRUE(verbose)) {
    warning("auto_garchx: very short series; estimation may be unstable.")
  }
  
  max.p <- as.integer(max.p)
  max.q <- as.integer(max.q)
  if (length(max.p) != 1L || length(max.q) != 1L ||
      !is.finite(max.p) || !is.finite(max.q) ||
      max.p < 0L || max.q < 0L) {
    stop("auto_garchx: 'max.p' and 'max.q' must be nonnegative integers.")
  }
  
  if (!is.numeric(clamp_factor) || length(clamp_factor) != 2L ||
      any(!is.finite(clamp_factor)) || any(clamp_factor <= 0)) {
    stop("auto_garchx: 'clamp_factor' must be a positive numeric vector of length 2.")
  }
  
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    
    if (nrow(xreg) != n) {
      stop("auto_garchx: nrow(xreg) must equal length(y).")
    }
    if (any(!is.finite(xreg))) {
      stop("auto_garchx: 'xreg' contains NA/Inf/NaN.")
    }
  }
  
  # Normalize search list for o: 0 = symmetric, 1 = GJR
  search_o <- as.integer(unique(search_o))
  search_o <- search_o[search_o %in% c(0L, 1L)]
  if (!length(search_o)) {
    search_o <- 0L
  }
  
  # Holder for the best model by BIC
  best <- list(
    bic        = Inf,
    aic        = Inf,
    fit        = NULL,
    q          = NA_integer_,
    p          = NA_integer_,
    o          = NA_integer_,
    refit_used = FALSE
  )
  
  # ---- 1) Grid search over (q, p, o) -----------------------------------------
  # garchx::garchx(order = c(q, p, o))
  for (q in 0:max.q) {
    for (p in 0:max.p) {
      for (o in search_o) {
        
        # Skip the trivial model
        if (q == 0L && p == 0L && o == 0L) {
          next
        }
        
        fit0 <- try(
          garchx::garchx(
            y = y,
            order = c(q, p, o),
            xreg = xreg,
            vcov.type = vcov.type
          ),
          silent = TRUE
        )
        
        if (inherits(fit0, "try-error")) {
          next
        }
        
        ll <- as.numeric(logLik(fit0))
        k  <- length(coef(fit0))
        
        bic0 <- -2 * ll + k * log(n)
        aic0 <- -2 * ll + 2 * k
        
        if (is.finite(bic0) && bic0 < best$bic) {
          best <- list(
            bic        = bic0,
            aic        = aic0,
            fit        = fit0,
            q          = q,
            p          = p,
            o          = o,
            refit_used = FALSE
          )
        }
      }
    }
  }
  
  if (is.null(best$fit)) {
    stop("auto_garchx: no converged model in grid search.")
  }
  
  # If no final refit is requested, return the grid-search winner
  if (!isTRUE(final_refit)) {
    return(best)
  }
  
  # ---- 2) Optional refit using unconditional-variance backcast ---------------
  cf <- coef(best$fit)
  
  omega_idx <- grep("^(intercept|omega)$", names(cf), ignore.case = TRUE)
  if (!length(omega_idx)) {
    return(best)
  }
  
  omega <- as.numeric(cf[omega_idx[1L]])
  
  alpha <- sum(cf[grep("^arch", names(cf), ignore.case = TRUE)], na.rm = TRUE)
  beta  <- sum(cf[grep("^garch", names(cf), ignore.case = TRUE)], na.rm = TRUE)
  
  gamma <- if (best$o == 1L) {
    sum(cf[grep("asym|gamma|gjr", names(cf), ignore.case = TRUE)], na.rm = TRUE)
  } else {
    0
  }
  
  # Approximate unconditional variance for GJR:
  # omega / (1 - alpha - beta - 0.5 * gamma)
  denom <- 1 - alpha - beta - 0.5 * gamma
  
  if (!is.finite(omega) || !is.finite(denom) || denom <= 0) {
    return(best)
  }
  
  uncond_var <- as.numeric(omega / denom)
  
  # Clamp the backcast relative to sample variance
  sample_var <- stats::var(y, na.rm = TRUE)
  rng <- range(clamp_factor)
  
  if (!is.finite(uncond_var) || uncond_var <= 0 ||
      uncond_var < rng[1L] * sample_var ||
      uncond_var > rng[2L] * sample_var) {
    if (isTRUE(verbose)) {
      warning("[auto_garchx] backcast out of range; fallback to sample variance.")
    }
    uncond_var <- sample_var
  }
  
  if (isTRUE(verbose)) {
    message(sprintf(
      "[auto_garchx] refit with backcast = %.6f (q = %d, p = %d, o = %d)",
      uncond_var, best$q, best$p, best$o
    ))
  }
  
  fit2 <- try(
    garchx::garchx(
      y = y,
      order = c(best$q, best$p, best$o),
      xreg = xreg,
      backcast.values = uncond_var,
      vcov.type = vcov.type
    ),
    silent = TRUE
  )
  
  if (inherits(fit2, "try-error")) {
    if (isTRUE(verbose)) {
      warning("[auto_garchx] backcast refit failed; keeping grid-search model.")
    }
    return(best)
  }
  
  ll2 <- as.numeric(logLik(fit2))
  k2  <- length(coef(fit2))
  
  bic2 <- -2 * ll2 + k2 * log(n)
  aic2 <- -2 * ll2 + 2 * k2
  
  if (is.finite(bic2) && bic2 < best$bic) {
    best$fit        <- fit2
    best$bic        <- bic2
    best$aic        <- aic2
    best$refit_used <- TRUE
  } else {
    best$refit_used <- FALSE
    
    if (isTRUE(verbose)) {
      warning("[auto_garchx] refit did not improve BIC; keeping grid-search model.")
    }
  }
  
  best
}

#' @rdname auto_garchx
#' @export
auto_garch <- auto_garchx