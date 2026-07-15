#' SynthPrediction — Synthetic Prediction with Donor Balancing & Shock Adjustment
#'
#' @description
#' Build k-step-ahead forecasts for a **target** time series by:
#' (1) estimating each **donor** series' post-shock fixed effect via (S)ARIMA
#' with a post-shock indicator; (2) combining donor effects by weights
#' (uniform or \code{\link{dbw}}); and (3) adjusting the target ARIMA forecast
#' by the weighted donor effect. Indices across lists are assumed aligned.
#'
#' **Seasonality:** \code{seasonal=TRUE} by default. For *each* series, seasonality
#' is enabled only if \code{forecast::findfrequency()} detects a frequency > 1
#' (NA or <= 1 disables seasonality for that series).
#'
#' @param Y_series_list List of time series (e.g., \code{xts}). \code{[[1]]} is the
#' target; \code{[[2]]}..\code{[[n+1]]} are donors.
#' @param covariates_series_list List of time series aligned with \code{Y_series_list}.
#' Used for DBW features and (optionally) as \code{xreg} in ARIMA. You may pass
#' \code{Y_series_list} here as a simple baseline.
#' @param shock_time_vec Vector of shock start times, each either a time index present
#' in the corresponding series or a positive integer row number. Length must equal
#' \code{length(Y_series_list)}.
#' @param shock_length_vec Integer vector (same length) giving the post-shock window
#' length (in rows) used to estimate donor fixed effects.
#' @param k Integer forecast horizon for the target.
#' @param dbw_scale,dbw_center Logical flags forwarded to \code{dbw} when computing donor weights.
#' @param dbw_indices Integer vector of covariate columns used by \code{dbw}.
#' If \code{NULL}, informative columns are auto-selected at pre-shock rows.
#' @param use_dbw Logical; if \code{FALSE}, use equal donor weights; otherwise use \code{dbw}.
#' @param princ_comp_input Optional PCA dimension for \code{dbw} (forwarded only).
#' @param covariate_indices Integer vector of columns from
#' \code{covariates_series_list[[i]]} to include as \code{xreg} in ARIMA
#' (lagged for the target; contemporaneous for donors). If \code{NULL}, donors use only
#' the post-shock indicator and the target uses no \code{xreg}.
#' @param geometric_sets,days_before_shocktime_vec Reserved (ignored).
#' @param arima_order Optional fixed ARIMA order \code{c(p,d,q)} for all series;
#' if \code{NULL}, \code{forecast::auto.arima()} is used.
#' @param seasonal Logical: whether seasonal models are allowed/requested. Default \code{TRUE}.
#' @param seasonal_order Optional seasonal order \code{c(P,D,Q)} used only when
#' \code{arima_order} is supplied and seasonal detection is \code{TRUE}.
#' @param seasonal_period Optional seasonal period. If \code{NULL} and \code{seasonal=TRUE},
#' it is auto-detected per series via \code{findfrequency()}.
#' @param user_ic_choice Information criterion for \code{auto.arima}
#' (\code{"aicc"}, \code{"aic"}, \code{"bic"}).
#' @param stepwise,approximation Passed through to \code{forecast::auto.arima()}.
#' @param plots Logical; if \code{TRUE} and a function named
#' \code{plot_maker_synthprediction} exists, it is called to visualize results.
#' @param display_ground_truth_choice Logical; forwarded to the plotting helper if it exists.
#' @param penalty_lambda Numeric L1/L2 penalty forwarded to \code{dbw}.
#' @param penalty_normchoice Character \code{"l1"} or \code{"l2"}; forwarded to \code{dbw}.
#' @param control_shocks Optional list of shock specifications to control for (remove) in the
#' target series model. Typically \code{list(list(time=..., length=..., shape=...))}.
#' @param pools Optional named list of donor pools. If provided, the call is
#' dispatched to \code{SynthPredictionMultiPool()}.
#' @param pool_agg Aggregation rule across pools ("sum" or "weighted").
#' @param pool_weights Optional named numeric vector of pool weights.
#' @param base_use_dbw Logical. Whether DBW is used in the baseline forecast.
#' @param pool_use_dbw Logical. Whether DBW is used within each pool.


#'
#' @return A list with components:
#' \itemize{
#'   \item \code{linear_combinations} — numeric donor weights \eqn{W}.
#'   \item \code{predictions} — list with numeric vectors \code{unadjusted}, \code{adjusted}, \code{arithmetic_mean}.
#'   \item \code{meta} — metadata (weights, shock timing, ARIMA orders, IC used, etc.).
#' }
#'
#' If \code{pools} is provided, the function dispatches to
#' \code{\link{SynthPredictionMultiPool}} and returns its output (see that function
#' for the multi-pool return structure).

#'
#' @seealso \code{\link{dbw}}, \code{\link[forecast]{auto.arima}}, \code{\link[forecast]{Arima}}
#'
#' @examples
#' \donttest{
#' library(xts)
#' set.seed(1)
#' T  <- 80
#' ix <- seq.Date(as.Date("2022-01-01"), by = "day", length.out = T)
#' y1 <- xts(rnorm(T), ix)  # target
#' y2 <- xts(rnorm(T), ix)  # donor 1
#' y3 <- xts(rnorm(T), ix)  # donor 2
#' Y  <- list(y1, y2, y3)
#' X  <- list(y1, y2, y3)   # simple baseline for DBW/xreg
#' shock_time  <- rep(ix[60], length(Y))
#' shock_len   <- rep(5L,   length(Y))
#' out <- SynthPrediction(
#'    Y_series_list = Y,
#'    covariates_series_list = X,
#'    shock_time_vec = shock_time,
#'    shock_length_vec = shock_len,
#'    k = 2,
#'    plots = FALSE
#' )
#' str(out$predictions)
#' }
#' @importFrom forecast Arima auto.arima forecast findfrequency
#' @importFrom lmtest   coeftest
#' @importFrom stats    predict na.omit var
#' @importFrom zoo      index
#' @importFrom xts      lag.xts
#' @export
SynthPrediction <- function(
    Y_series_list,
    covariates_series_list,
    shock_time_vec,
    shock_length_vec,
    k = 1,
    dbw_scale = TRUE,
    dbw_center = TRUE,
    dbw_indices = NULL,
    use_dbw = TRUE,
    princ_comp_input = NULL,
    covariate_indices = NULL,
    geometric_sets = NULL,
    days_before_shocktime_vec = NULL,
    arima_order = NULL,
    seasonal = TRUE,                 # default TRUE; gated per series
    seasonal_order = NULL,
    seasonal_period = NULL,          # if NULL, detect via forecast::findfrequency()
    user_ic_choice = c("aicc","aic","bic")[1],
    stepwise = TRUE,
    approximation = TRUE,
    plots = TRUE,
    display_ground_truth_choice = FALSE,
    penalty_lambda = 0,
    penalty_normchoice = "l1",
    control_shocks = NULL,
    pools = NULL,
    pool_agg = "sum",
    pool_weights = NULL,
    base_use_dbw = FALSE,
    pool_use_dbw = TRUE,
    ...
){
  # ---- Multi-pool dispatch (EARLY EXIT) ----
  if (!is.null(pools)) {
    return(SynthPredictionMultiPool(
      Y_series_list          = Y_series_list,
      covariates_series_list = covariates_series_list,
      shock_time_vec         = shock_time_vec,
      shock_length_vec       = shock_length_vec,
      k                      = k,
      pools                  = pools,
      pool_agg               = pool_agg,
      pool_weights           = pool_weights,
      base_use_dbw           = base_use_dbw,
      pool_use_dbw           = pool_use_dbw,
      ...
    ))
  }
  
  # Ensure an identity transform exists for DBW transforms
  if (!exists("id")) id <- function(x) x
  
  n <- as.integer(length(Y_series_list) - 1L)
  dbw_output <- NULL
  
  ## ---------- 0) Convert shock times to integer indices ----------
  if (length(shock_time_vec) != length(Y_series_list)) {
    stop(
      "SynthPrediction: length(shock_time_vec) must equal ",
      "length(Y_series_list)."
    )
  }
  
  integer_shock_time_vec <- vapply(
    seq_len(n + 1L),
    function(i) {
      ti <- shock_time_vec[i]
      series_index <- zoo::index(Y_series_list[[i]])
      
      # Date, POSIX time, or character labels:
      # find the corresponding row in the series index.
      if (
        inherits(ti, "Date") ||
        inherits(ti, "POSIXt") ||
        is.character(ti)
      ) {
        pos <- match(
          as.character(ti),
          as.character(series_index)
        )
        
        if (is.na(pos)) {
          stop(sprintf(
            paste0(
              "SynthPrediction: shock_time_vec[%d] = '%s' ",
              "is not present in the index of series %d."
            ),
            i,
            as.character(ti),
            i
          ))
        }
        
        return(as.integer(pos))
      }
      
      # Numeric values are interpreted as row numbers.
      if (
        !is.numeric(ti) ||
        length(ti) != 1L ||
        !is.finite(ti) ||
        ti != as.integer(ti)
      ) {
        stop(sprintf(
          paste0(
            "SynthPrediction: shock_time_vec[%d] must be ",
            "a valid time label or a positive integer row number."
          ),
          i
        ))
      }
      
      pos <- as.integer(ti)
      
      if (pos < 1L || pos > NROW(Y_series_list[[i]])) {
        stop(sprintf(
          paste0(
            "SynthPrediction: shock_time_vec[%d] = %d ",
            "is outside the available rows of series %d."
          ),
          i,
          pos,
          i
        ))
      }
      
      pos
    },
    integer(1L)
  )
  
  ## ---------- helper: per-series seasonal decision ----------
  decide_seasonal <- function(y_vec, seasonal_flag, seasonal_period_user = NULL){
    if (!isTRUE(seasonal_flag)) return(list(flag = FALSE, period = NULL))
    if (!is.null(seasonal_period_user)) {
      return(list(flag = TRUE, period = as.integer(round(seasonal_period_user))))
    }
    freq <- tryCatch(forecast::findfrequency(as.numeric(y_vec)),
                     error = function(e) NA_real_)
    if (is.na(freq) || freq <= 1) list(flag = FALSE, period = NULL)
    else                          list(flag = TRUE,  period = as.integer(round(freq)))
  }
  
  ## ---------- 1) Donor (S)ARIMA with post-shock indicator ----------
  omega_star_hat_vec <- numeric(n)
  order_of_arima     <- vector("list", n+1)
  
  for (i in 2:(n+1)) {
    Tpre  <- integer_shock_time_vec[i]
    Ls    <- shock_length_vec[i]
    last  <- Tpre + Ls
    
    # Post-shock indicator: zeros through preT, ones during the shock window
    vec0  <- rep(0, Tpre)
    vec1  <- rep(1, Ls)
    post_shock_indicator <- c(vec0, vec1)    # equals 1 from (preT+1) onward
    
    # Donor design matrix; match y length to rows of X_i_final if NA removal occurs
    if (is.null(covariate_indices)) {
      X_i_final <- post_shock_indicator
      y_fit     <- Y_series_list[[i]][1:last]
    } else {
      X_i_subset <- covariates_series_list[[i]][1:last, covariate_indices, drop = FALSE]
      X_i_final  <- cbind(X_i_subset, post_shock_indicator)
      X_i_final  <- stats::na.omit(X_i_final)
      y_fit      <- Y_series_list[[i]][1:NROW(X_i_final)]
    }
    
    # Per-donor seasonal decision
    seas_dec <- decide_seasonal(y_vec = y_fit, seasonal_flag = seasonal, seasonal_period_user = seasonal_period)
    seasonal_flag_i <- seas_dec$flag
    seas_period_i   <- seas_dec$period
    
    # Fit donor model (fixed order if provided; otherwise auto.arima)
    donor_model <- if (!is.null(arima_order)) {
      if (seasonal_flag_i && !is.null(seasonal_order)) {
        forecast::Arima(y_fit, order = arima_order,
                        seasonal = list(order = seasonal_order, period = seas_period_i),
                        xreg = X_i_final)
      } else {
        forecast::Arima(y_fit, order = arima_order, seasonal = FALSE, xreg = X_i_final)
      }
    } else {
      tryCatch(
        forecast::auto.arima(y_fit, xreg = X_i_final, ic = user_ic_choice,
                             seasonal = seasonal_flag_i, stepwise = stepwise,
                             approximation = approximation),
        error = function(e){
          warning(sprintf("Donor %d: auto.arima() failed → fallback to Arima(1,1,1)", i))
          forecast::Arima(y_fit, order = c(1,1,1), seasonal = FALSE, xreg = X_i_final)
        }
      )
    }
    
    # Store ARIMA structure and extract post-shock coefficient (last xreg)
    order_of_arima[[i]] <- donor_model$arma
    ct <- lmtest::coeftest(donor_model)
    omega_star_hat_vec[i-1] <- ct[nrow(ct), "Estimate"]
  }
  
  ## ---------- 2) Donor weights (equal or DBW) ----------
  # ---- DBW lookback handling ----
  if (is.null(days_before_shocktime_vec)) {
    X_lookback_indices_dbw <- 1
  } else {
    stopifnot(length(days_before_shocktime_vec) == 1,
              days_before_shocktime_vec >= 1)
    X_lookback_indices_dbw <- seq_len(days_before_shocktime_vec)
  }
  
  n_donors <- n
  equal_w  <- rep(1/n_donors, n_donors)
  
  if (!isTRUE(use_dbw)) {
    w_hat <- equal_w
  } else {
    # Prefer covariates for DBW features; fallback to Y
    X_for_dbw <- if (!is.null(covariates_series_list)) covariates_series_list else Y_series_list
    X_for_dbw <- lapply(X_for_dbw, function(obj){
      m <- as.matrix(obj)
      if (is.null(dim(m))) m <- matrix(m, ncol = 1)
      if (is.null(colnames(m))) colnames(m) <- paste0("x", seq_len(ncol(m)))
      m[, , drop = FALSE]
    })
    
    # Auto-select informative dbw_indices using cross-series variance at preT rows
    if (is.null(dbw_indices)) {
      preT_int <- as.integer(integer_shock_time_vec)
      rows <- Map(function(m, t) m[t, , drop = FALSE], X_for_dbw, preT_int)
      M <- do.call(rbind, rows)
      keep <- apply(M, 2, function(col) stats::var(col, na.rm = TRUE) > 1e-8)
      dbw_indices <- if (any(keep)) which(keep) else integer(0)
    }
    
    if (length(dbw_indices) == 0) {
      warning("[DBW] no informative columns pre-shock → fallback to equal weights.")
      w_hat <- equal_w
    } else {
      dbw_output <- tryCatch(
        dbw(
          X = X_for_dbw,
          dbw_indices = dbw_indices,
          shock_time_vec = as.integer(integer_shock_time_vec),
          scale = dbw_scale, center = dbw_center,
          sum_to_1 = TRUE, bounded_below_by = 0, bounded_above_by = 1,
          Y = Y_series_list,
          Y_lookback_indices = NULL,
          X_lookback_indices = X_lookback_indices_dbw,
          inputted_transformation = id,
          penalty_lambda = penalty_lambda,
          penalty_normchoice = penalty_normchoice
        ),
        error = function(e){
          warning(paste0("[DBW] failed: ", e$message))
          NULL
        }
      )
      w_hat <- if (is.null(dbw_output)) equal_w else dbw_output$opt_params
    }
  }
  
  # Weighted donor shock effect
  omega_star_hat <- as.numeric(w_hat %*% omega_star_hat_vec)
  
  ## ---------- 3) Target (S)ARIMA: train through preT, forecast from preT+1 ----------
  pt1 <- integer_shock_time_vec[1]
  stopifnot(pt1 >= 1L)
  
  target_series     <- Y_series_list[[1]][ 1:pt1 ]  # includes preT → 1-step = preT+1
  forecast_horizon  <- k
  
  # --- Build Control Shocks Dummy Matrix ---
  # Note: build_control_xreg should be available from R/utils.R
  C_train <- build_control_xreg(
    last = pt1,
    control_spec = if (is.null(control_shocks)) NULL else control_shocks[[1]]
  )
  if (!is.null(C_train)) C_train <- as.matrix(C_train)
  
  C_future <- if (is.null(C_train)) NULL else
    matrix(0, nrow = forecast_horizon, ncol = ncol(C_train),
           dimnames = list(NULL, colnames(C_train)))
  
  # --- Combine with Covariates ---
  if (is.null(covariate_indices)) {
    xreg_input  <- C_train
    xreg_future <- C_future
  } else {
    X_lagged <- xts::lag.xts(
      covariates_series_list[[1]][1:pt1, covariate_indices, drop = FALSE],
      k = 1,
      na.pad = TRUE
    )
    X_mat <- as.matrix(X_lagged)
    
    # Fill first row NA from lagging
    if (NROW(X_mat) >= 1) X_mat[1, ] <- 0
    
    # Clean cbind logic
    xreg_input <- if (is.null(C_train)) X_mat else cbind(X_mat, C_train)
    
    # Future covariates
    X_last <- covariates_series_list[[1]][pt1, covariate_indices, drop = FALSE]
    X_future_vals <- matrix(rep(as.numeric(X_last), forecast_horizon),
                            nrow = forecast_horizon, byrow = TRUE)
    colnames(X_future_vals) <- colnames(X_mat)
    
    xreg_future <- if (is.null(C_future)) X_future_vals else cbind(X_future_vals, C_future)
  }
  
  # --- Safety Check: Force Matrix ---
  if (!is.null(xreg_input)) {
    xreg_input  <- as.matrix(xreg_input)
    xreg_future <- as.matrix(xreg_future)
    stopifnot(NROW(xreg_input) == length(target_series))
    if(anyNA(xreg_input)) xreg_input[is.na(xreg_input)] <- 0
  }
  
  seas_dec_t <- decide_seasonal(y_vec = target_series, seasonal_flag = seasonal, seasonal_period_user = seasonal_period)
  seasonal_flag_target   <- seas_dec_t$flag
  target_seasonal_period <- seas_dec_t$period
  
  target_model <- if (!is.null(arima_order)) {
    if (seasonal_flag_target && !is.null(seasonal_order)) {
      forecast::Arima(target_series, order = arima_order,
                      seasonal = list(order = seasonal_order, period = target_seasonal_period),
                      xreg = xreg_input)
    } else {
      forecast::Arima(target_series, order = arima_order, seasonal = FALSE, xreg = xreg_input)
    }
  } else {
    forecast::auto.arima(
      target_series, xreg = xreg_input, ic = user_ic_choice,
      seasonal = seasonal_flag_target, stepwise = stepwise, approximation = approximation
    )
  }
  
  # Forecast k steps ahead; use predict(), fallback to forecast() on error
  prediction_result <- tryCatch({
    if (is.null(xreg_future)) {
      stats::predict(target_model, n.ahead = forecast_horizon)
    } else {
      stats::predict(target_model, n.ahead = forecast_horizon, newxreg = xreg_future)
    }
  }, error = function(e){
    warning(paste("Target predict() failed:", e$message))
    fc_try <- tryCatch({
      if (is.null(xreg_future)) forecast::forecast(target_model, h = forecast_horizon)
      else                      forecast::forecast(target_model, h = forecast_horizon, xreg = xreg_future)
    }, error = function(e2) e2)
    if (inherits(fc_try, "error")) target_model else fc_try
  })
  
  if (inherits(prediction_result, "forecast")) {
    pred_vec <- as.numeric(prediction_result$mean)
  } else if (!is.null(prediction_result$pred)) {
    pred_vec <- as.numeric(prediction_result$pred)
  } else if (inherits(prediction_result, c("Arima","ARIMA"))) {
    stop("predict() returned a model; forecasting failed.")
  } else {
    stop("Unknown prediction object returned by predict/forecast().")
  }
  
  # Length guarantee
  forecast_horizon <- as.integer(forecast_horizon)
  if (length(pred_vec) < forecast_horizon) {
    stop(sprintf("Forecast length %d < k=%d. Something went wrong in predict/forecast.",
                 length(pred_vec), forecast_horizon))
  }
  if (length(pred_vec) > forecast_horizon) {
    warning(sprintf("Forecast length %d > k=%d; truncating to first k steps.",
                    length(pred_vec), forecast_horizon))
    pred_vec <- pred_vec[seq_len(forecast_horizon)]
  }
  
  ## ---------- 4) Shock adjustment & outputs ----------
  raw_vec        <- pred_vec
  adjusted_vec   <- pred_vec + omega_star_hat
  arithmetic_vec <- pred_vec + mean(omega_star_hat_vec)
  
  predictions <- list(
    unadjusted      = as.numeric(raw_vec),
    adjusted        = as.numeric(adjusted_vec),
    arithmetic_mean = as.numeric(arithmetic_vec)
  )
  
  meta <- list(
    n_donors       = n,
    shock_time     = shock_time_vec,
    shock_length   = shock_length_vec,
    dbw_status     = if (is.null(dbw_output)) NA else dbw_output[[2]],
    weights        = w_hat,
    omega_vec      = omega_star_hat_vec,
    combined_omega = as.numeric(omega_star_hat),
    arima_order    = order_of_arima,
    ic_used        = user_ic_choice
  )
  
  output_list <- list(
    linear_combinations = w_hat,
    predictions         = predictions,
    meta                = meta
  )
  
  if (isTRUE(plots)) {
    if(exists("plot_maker_synthprediction")) {
      plot_maker_synthprediction(
        Y_series_list,
        shock_time_vec,
        integer_shock_time_vec,
        shock_length_vec,
        raw_vec,
        w_hat,
        omega_star_hat_vec,
        adjusted_vec,
        arithmetic_vec,
        display_ground_truth = display_ground_truth_choice
      )
    }
  }
  
  return(output_list)
}