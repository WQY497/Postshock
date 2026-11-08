#' Boundary-aligned Apple S1 article inputs
#'
#' Frozen inputs for the retrospective one-step Apple S1 illustration.
#' Derived without downloading or changing observations from
#' \code{aapl_onestep_input} and \code{aapl_control_input}. This object contains
#' neither fitted models nor precomputed forecast values.
#'
#' @format A named list:
#' \describe{
#' \item{Y_series_list}{Target AAPL close observations from 7 October through
#' 10 November 2020 (25 trading rows), followed by three unchanged 25-row donor
#' series. Prices retain the source objects' 2020-dollar processing.}
#' \item{covariates_series_list}{Aligned GSPC and IRX matching panels. GSPC
#' retains the source's 2020-base index scaling; IRX retains its quoted units.}
#' \item{shock_time_vec}{Last pre-effect boundaries: target row 25 and donor
#' rows obtained by locating each effect date and subtracting one (row 12).}
#' \item{shock_length_vec}{One observation for each effect window. Donor
#' response models therefore use rows 1 through 13, not all supplied rows.}
#' \item{control_shocks}{Target-only nested point specification for the
#' 9 October 2020 response. Forecast-time controls are zero.}
#' \item{settings}{Explicit S1 parameters: one-step automatic ARIMA with AICc,
#' stepwise and approximation enabled, no seasonality, no response-model
#' market regressors, ten-row centered/scaled DBW, no PCA and zero penalty.}
#' \item{dates}{Supplied and fitted window sizes, start/end dates, boundaries,
#' effect dates and matching starts for every series. Donor effect dates are
#' 7 June 2017, 31 October 2018 and 4 June 2019.}
#' \item{control_date}{Historical control response date.}
#' \item{evaluation}{11 November 2020 and its realized AAPL close; separate
#' from all target fitting and matching inputs and not passed to the model.}
#' \item{provenance}{Source names, schema and retrospective-data limitations.}
#' }
#'
#' @details
#' Overlapping dated rows of the two old objects are checked for exact equality.
#' Their union supplies the target window and separate evaluation observation.
#' Donor panels are unchanged; API boundaries immediately precede the located
#' donor effect observations. The old objects are retained unchanged.
#'
#' Market prices originate from Yahoo Finance, with annual CPI processing
#' based on FRED CPIAUCSL as documented in the article supplement. Retrieval
#' timestamps, raw provider snapshots and some original object builders are
#' unavailable. This is not a real-time-vintage backtest or an exact reproduction
#' of the original delayed-response experiment. Provider revisions need not
#' reconstruct the same historical values.
#'
#' The development-checkout builder is \file{data-raw/aapl_s1_input.R}
#' (excluded from the built tarball).
#' Its implementation is installed as \file{replication/apple_data_builder.R}.
#' The shared Apple diagnostic script is installed as
#' \file{replication/apple_diagnostic.R} and retains the article's numerical
#' assertions. Donor-effect and weight reference checks retain their original
#' absolute tolerances; their error bounds are propagated to the combined
#' adjustment, with a separate weighted-sum identity check. Its implementation
#' manifest distinguishes same-version builds.
#'
#' @examples
#' data(aapl_s1_input)
#' input_fields <- c("Y_series_list", "covariates_series_list", "shock_time_vec",
#'                   "shock_length_vec", "control_shocks")
#' arguments <- c(aapl_s1_input[input_fields], aapl_s1_input$settings)
#' fit <- do.call(SynthPrediction, arguments)
#' fit$predictions
#' aapl_s1_input$evaluation
#'
#' @usage data(aapl_s1_input)
"aapl_s1_input"
