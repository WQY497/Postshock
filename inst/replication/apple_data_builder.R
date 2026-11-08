build_aapl_s1_input <- function(aapl_onestep_input, aapl_control_input) {
  if (!requireNamespace("xts", quietly = TRUE)) stop("Building Apple S1 requires xts.")
  apple <- aapl_onestep_input
  market_names <- c("GSPC", "IRX")
  panels <- list(apple$Z_tgt_base,
    cbind(apple$Y_series_list[[1L]], apple$X_for_dbw[[1L]]),
    cbind(aapl_control_input$Y[[1L]], aapl_control_input$X[[1L]]))
  for (panel in panels) {
    stopifnot(inherits(panel, "xts"), identical(colnames(panel), c("AAPL", market_names)),
      !anyDuplicated(zoo::index(panel)), all(is.finite(as.matrix(panel))))
  }
  for (first in seq_len(length(panels) - 1L)) {
    for (second in seq.int(first + 1L, length(panels))) {
      dates <- intersect(as.character(zoo::index(panels[[first]])), as.character(zoo::index(panels[[second]])))
      stopifnot(length(dates) > 0L,
        identical(as.matrix(panels[[first]][match(dates, as.character(zoo::index(panels[[first]]))), ]),
                  as.matrix(panels[[second]][match(dates, as.character(zoo::index(panels[[second]]))), ])))
    }
  }
  target_panel <- do.call(rbind, panels)
  target_panel <- target_panel[!duplicated(zoo::index(target_panel)), ]
  target_panel <- target_panel[order(zoo::index(target_panel)), ]
  origin <- as.Date(apple$last_obs_day)
  evaluation_date <- as.Date(apple$target_day)
  control_date <- as.Date(apple$control_shock_day)
  stopifnot(identical(origin, as.Date("2020-11-10")),
    identical(evaluation_date, as.Date("2020-11-11")),
    identical(control_date, as.Date("2020-10-09")))
  target_window <- tail(target_panel[as.Date(zoo::index(target_panel)) <= origin, ], 25L)
  stopifnot(NROW(target_window) == 25L,
    identical(as.Date(zoo::index(target_window)[1L]), as.Date("2020-10-07")),
    identical(as.matrix(target_window), as.matrix(apple$Z_tgt_base)),
    all(as.Date(zoo::index(target_window)) < evaluation_date))
  donors <- apple$Y_series_list[-1L]
  markets <- apple$X_for_dbw[-1L]
  effect_dates <- as.Date(apple$donor_effect_days)
  stopifnot(identical(unname(effect_dates), as.Date(c("2017-06-07", "2018-10-31", "2019-06-04"))))
  effect_rows <- vapply(seq_along(donors), function(donor) {
    stopifnot(NROW(donors[[donor]]) == 25L, NROW(markets[[donor]]) == 25L,
      identical(zoo::index(donors[[donor]]), zoo::index(markets[[donor]])),
      all(is.finite(as.matrix(donors[[donor]]))), all(is.finite(as.matrix(markets[[donor]]))),
      identical(colnames(markets[[donor]]), market_names))
    match(effect_dates[donor], as.Date(zoo::index(donors[[donor]])))
  }, integer(1L))
  stopifnot(!anyNA(effect_rows), all(effect_rows > 10L), all(effect_rows == 13L))
  boundaries <- effect_rows - 1L
  control_row <- match(control_date, as.Date(zoo::index(target_window)))
  evaluation_row <- match(evaluation_date, as.Date(zoo::index(target_panel)))
  stopifnot(!is.na(control_row), !is.na(evaluation_row))
  responses <- c(list(target_window[, "AAPL", drop = FALSE]), donors)
  covariates <- c(list(target_window[, market_names, drop = FALSE]), markets)
  controls <- vector("list", length(responses))
  controls[[1L]] <- list(list(time = control_row, length = 1L, shape = "point"))
  boundary_rows <- c(NROW(target_window), boundaries)
  date_records <- do.call(rbind, lapply(seq_along(responses), function(series) {
    dates <- as.Date(zoo::index(responses[[series]]))
    data.frame(series = if (series == 1L) "target" else paste0("donor", series - 1L),
      supplied_rows = length(dates), start = dates[1L], end = tail(dates, 1L),
      boundary_row = boundary_rows[series], boundary = dates[boundary_rows[series]],
      effect = if (series == 1L) evaluation_date else dates[effect_rows[series - 1L]],
      fitted_rows = if (series == 1L) length(dates) else effect_rows[series - 1L],
      matching_start = dates[boundary_rows[series] - 9L])
  }))
  list(Y_series_list = responses, covariates_series_list = covariates,
    shock_time_vec = boundary_rows, shock_length_vec = rep(1L, length(responses)),
    control_shocks = controls,
    settings = list(k = 1L, use_dbw = TRUE, dbw_indices = 1:2, covariate_indices = NULL,
      arima_order = NULL, seasonal = FALSE, user_ic_choice = "aicc",
      stepwise = TRUE, approximation = TRUE, dbw_center = TRUE, dbw_scale = TRUE,
      days_before_shocktime_vec = 10L, princ_comp_input = NULL,
      penalty_lambda = 0, penalty_normchoice = "l1"),
    dates = date_records, control_date = control_date,
    evaluation = list(date = evaluation_date,
      realized = as.numeric(target_panel[evaluation_row, "AAPL"])),
    provenance = list(schema = "postshock-apple-S1-v1",
      sources = c("aapl_onestep_input", "aapl_control_input"),
      illustration = "Retrospective one-step illustration; not a real-time-vintage backtest.",
      limitations = "Historical download timestamps and raw provider snapshots are not recorded; original object builders are incomplete. Annual CPI processing and revised provider histories are retained without changes."))
}
