#' One-step volatility forecasting input for IYG
#'
#' @description
#' A ready-to-use input object for the IYG one-step post-shock volatility
#' forecasting experiment. The target episode is the 2016 U.S. election, and
#' the candidate donor episodes include earlier U.S. elections and
#' Brexit-related political uncertainty events.
#'
#' The object stores three pre-assembled donor-pool specifications:
#' `all`, `restricted`, and `elections_only`. Each specification contains the
#' `Y_series_list`, `covariates_series_list`, `target_covariates`,
#' `shock_time_vec`, and `shock_length_vec` arguments required by
#' `SynthVolForecast()`.
#'
#' @format A named list with the following components:
#' \describe{
#'   \item{description}{Short description of the data object.}
#'   \item{target_symbol}{Character scalar giving the target asset symbol.}
#'   \item{target_event}{Character scalar giving the target event name.}
#'   \item{inputs}{Named list of ready-to-use `SynthVolForecast()` input objects. The current object includes `all`, `restricted`, and `elections_only`.}
#'   \item{donor_pool_specs}{Named list describing the donor events included in each donor-pool specification.}
#'   \item{evaluation}{List describing the realized-variance proxy, forecast horizon, and proxy index.}
#' }
#'
#' Each element of `inputs` is itself a named list with the following components:
#' \describe{
#'   \item{target_event}{Character scalar giving the target event name.}
#'   \item{donor_events}{Character vector of donor event names.}
#'   \item{Y_series_list}{List of return series, with the target series first followed by the donor series.}
#'   \item{covariates_series_list}{List of donor covariate data frames.}
#'   \item{target_covariates}{Covariate data frame for the target episode.}
#'   \item{shock_time_vec}{Integer vector of shock-time indices for the target and donor series.}
#'   \item{shock_length_vec}{Integer vector of post-shock window lengths for the target and donor series.}
#' }
#'
#' @details
#' The outcome series are daily log returns for IYG. The covariates include
#' aligned macro-financial variables used for donor balancing in the
#' `SynthVolForecast()` workflow. Since the true conditional variance is latent
#' in real data, empirical evaluation uses the next-day squared return as a
#' realized variance proxy.
#'
#' The `all` donor pool contains all candidate donor episodes. The `restricted`
#' donor pool contains the 2012 Election, Brexit Poll Released, and Brexit
#' episodes. The `elections_only` donor pool contains the 2004, 2008, and 2012
#' U.S. election episodes.
#'
#' @examples
#' data("iyg_onestep_input")
#'
#' names(iyg_onestep_input)
#' names(iyg_onestep_input$inputs)
#'
#' iyg_restricted <- iyg_onestep_input$inputs$restricted
#' names(iyg_restricted)
#' iyg_restricted$donor_events
#'
#' @name iyg_onestep_input
#' @docType data
#' @keywords datasets
"iyg_onestep_input"