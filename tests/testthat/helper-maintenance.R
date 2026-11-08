maintenance_input <- function() {
  set.seed(51026)
  rows <- 90L
  dates <- as.Date("2020-01-01") + seq_len(rows) - 1L
  responses <- lapply(seq_len(4L), function(series) {
    values <- as.numeric(stats::arima.sim(list(ar = 0.45), n = rows)) + 10
    values[71:73] <- values[71:73] + series / 2
    xts::xts(values, order.by = dates)
  })
  markets <- lapply(seq_len(4L), function(series)
    xts::xts(cbind(first = stats::rnorm(rows) + series,
                  second = stats::rnorm(rows)), order.by = dates))
  names(responses) <- names(markets) <- c("target", "donorA", "donorB", "donorC")
  list(Y_series_list = responses, covariates_series_list = markets,
       shock_time_vec = rep(70L, 4L), shock_length_vec = rep(3L, 4L),
       k = 2L, arima_order = c(1L, 0L, 0L), seasonal = FALSE,
       use_dbw = FALSE)
}

