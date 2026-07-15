#' Multi-pool shock aggregation for SynthPrediction
#'
#' @description
#' Extends \code{SynthPrediction()} to multiple donor pools representing
#' distinct shock channels (e.g., supply-side, demand-side).
#'
#' @inheritParams SynthPrediction
#' @param pools Named list. Each element is a character vector of donor series
#'   names (matching \code{names(Y_series_list)}). Pool names define channels.
#' @param pool_agg Aggregation rule across pools. One of "sum" or "weighted".
#' @param pool_weights Optional named numeric vector of pool weights.
#' @param base_use_dbw Logical. Whether DBW is used in the baseline forecast.
#' @param pool_use_dbw Logical. Whether DBW is used within each pool.
#'
#' @return
#' Same structure as \code{SynthPrediction()}, with extra:
#' \code{meta$multi_pool} and \code{predictions$adjusted_multipool}.
#'
#' @export
SynthPredictionMultiPool <- function(
    Y_series_list,
    covariates_series_list,
    shock_time_vec,
    shock_length_vec,
    k = 1,
    pools,
    pool_agg = c("sum", "weighted"),
    pool_weights = NULL,
    base_use_dbw = FALSE,
    pool_use_dbw = TRUE,
    ...
) {
  
  # ---- minimal checks ----
  if (!is.list(pools) || length(pools) < 1) stop("`pools` must be a non-empty named list.")
  if (is.null(names(pools)) || any(names(pools) == "")) stop("`pools` must be a *named* list.")
  if (!is.list(Y_series_list) || length(Y_series_list) < 2) stop("Y_series_list must be a list with target + at least 1 donor.")
  if (length(shock_time_vec) != length(Y_series_list)) stop("shock_time_vec length must match Y_series_list.")
  if (length(shock_length_vec) != length(Y_series_list)) stop("shock_length_vec length must match Y_series_list.")
  if (!is.list(covariates_series_list)) stop("covariates_series_list must be a list (can be Y_series_list as baseline).")
  if (length(covariates_series_list) != length(Y_series_list)) {
    stop("For multi-pool, covariates_series_list must have the same length/order as Y_series_list.")
  }
  
  # We require Y to be named because pools refer to donor series by name.
  if (is.null(names(Y_series_list)) || any(names(Y_series_list) == "")) {
    stop("For multi-pool, Y_series_list must be a *named* list so pools can reference donors by name.")
  }
  
  pool_agg <- match.arg(pool_agg)
  
  # ---- capture dots; DROP multipool-only args to avoid leaking/recursion ----
  dots <- list(...)
  drop_names <- c("pools", "pool_agg", "pool_weights", "base_use_dbw", "pool_use_dbw")
  dots <- dots[setdiff(names(dots), drop_names)]
  
  # helper: call SynthPrediction safely (force pools=NULL)
  call_SP <- function(args) {
    do.call(SynthPrediction, c(args, list(pools = NULL), dots))
  }
  
  # target is always the first element (matches SynthPrediction contract)
  target_name <- names(Y_series_list)[1]
  
  # ---- 1) baseline: common unadjusted forecast ----
  base_out <- call_SP(list(
    Y_series_list          = Y_series_list,
    covariates_series_list = covariates_series_list,
    shock_time_vec         = shock_time_vec,
    shock_length_vec       = shock_length_vec,
    k                      = k,
    use_dbw                = base_use_dbw
  ))
  yhat_unadj <- as.numeric(base_out$predictions$unadjusted)
  
  # ---- 2) pool effects ----
  pool_effects <- setNames(numeric(length(pools)), names(pools))
  
  for (nm in names(pools)) {
    
    donors <- unique(stats::na.omit(pools[[nm]]))
    donors <- as.character(donors)
    if (length(donors) < 1) stop("Pool '", nm, "' must be a non-empty character vector of donor names.")
    
    # Prevent target accidentally included as donor
    donors <- setdiff(donors, target_name)
    if (length(donors) < 1) stop("Pool '", nm, "' must contain at least one donor (not the target).")
    
    y_keep <- c(target_name, donors)
    
    if (any(!y_keep %in% names(Y_series_list))) {
      bad <- y_keep[!y_keep %in% names(Y_series_list)]
      stop("Pool '", nm, "' has unknown Y series: ", paste(bad, collapse = ", "))
    }
    
    # Subset by positions aligned with Y (NOT by covariate names)
    idx <- match(y_keep, names(Y_series_list))
    if (anyNA(idx)) stop("Internal error: match() produced NA in pool '", nm, "'.")
    
    out_nm <- call_SP(list(
      Y_series_list          = Y_series_list[idx],
      covariates_series_list = covariates_series_list[idx],
      shock_time_vec         = shock_time_vec[idx],
      shock_length_vec       = shock_length_vec[idx],
      k                      = k,
      use_dbw                = pool_use_dbw
    ))
    
    pool_effects[[nm]] <- as.numeric(out_nm$meta$combined_omega)
  }
  
  # ---- 3) aggregate ----
  agg_effect <- if (pool_agg == "sum") {
    sum(pool_effects)
  } else {
    if (is.null(pool_weights)) stop("pool_weights is required when pool_agg='weighted'.")
    if (is.null(names(pool_weights))) stop("pool_weights must be a named numeric vector.")
    w <- pool_weights[names(pool_effects)]
    if (anyNA(w)) stop("pool_weights must include all pool names: ", paste(names(pool_effects), collapse = ", "))
    sum(w * pool_effects)
  }
  
  # ---- 4) final adjusted ----
  out <- base_out
  out$predictions$adjusted_multipool <- yhat_unadj + as.numeric(agg_effect)
  
  out$meta$multi_pool <- list(
    pools             = pools,
    target_name       = target_name,
    pool_effects      = pool_effects,
    pool_agg          = pool_agg,
    pool_weights      = pool_weights,
    aggregated_effect = as.numeric(agg_effect),
    base_use_dbw      = base_use_dbw,
    pool_use_dbw      = pool_use_dbw
  )
  
  out
}
