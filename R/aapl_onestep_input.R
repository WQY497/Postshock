#' Apple one-step control-shock experiment input
#'
#' Prepared input object for the Apple one-step post-shock forecasting
#' experiment used in the technical report.
#'
#' The object supports reproduction of the experiment combining:
#' \itemize{
#'   \item a one-step local regression baseline for the target series,
#'   \item donor-side shock effect estimation,
#'   \item donor balancing weights (DBW),
#'   \item control-shock adjustment via \code{build_control_xreg()}.
#' }
#'
#' @format A named list with components:
#' \describe{
#'   \item{Y_series_list}{List of target and donor response series.}
#'   \item{X_for_dbw}{List of target and donor covariate series used for DBW.}
#'   \item{shock_time_vec}{Local shock-time indices for the target and donor windows.}
#'   \item{shock_length_vec}{Shock lengths for the target and donor windows.}
#'   \item{Z_tgt_base}{Target baseline window used for one-step regression forecasting.}
#'   \item{control_shock_day}{Date of the control shock included in the baseline model.}
#'   \item{target_day}{Forecast target date.}
#'   \item{last_obs_day}{Last observed date used for one-step prediction.}
#'   \item{donor_effect_days}{Donor shock-effect dates used in the experiment.}
#'   \item{release_day}{MacBook release date associated with the target event.}
#' }
#'
#' @usage data(aapl_onestep_input)
"aapl_onestep_input"