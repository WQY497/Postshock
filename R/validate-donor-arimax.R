validate_donor_arimax <- function(response, covariates, columns, last, donor_label) {
  fail <- function(...) stop("SynthPrediction: ", donor_label, " ARIMAX: ", ..., call. = FALSE)
  if (length(last) != 1L || !is.finite(last) || last < 1 || last != floor(last)) {
    fail("invalid fitting-window endpoint.")
  }
  if (NROW(response) < last) {
    fail("response requires ", last, " original rows; available ", NROW(response), ".")
  }
  if (is.null(covariates) || NROW(covariates) < last) {
    fail("covariates require ", last, " original rows; available ", NROW(covariates), ".")
  }
  available <- seq_len(NCOL(covariates))
  names(available) <- colnames(covariates)
  selected <- tryCatch(available[columns], error = function(condition) NULL)
  if (!length(selected) || anyNA(selected) ||
      (is.numeric(columns) && any(!is.finite(columns) | columns != floor(columns)))) {
    fail("invalid or absent covariate columns: ", paste(columns, collapse = ", "), ".")
  }
  rows <- seq_len(last)
  response_window <- response[rows]
  market_window <- tryCatch(as.matrix(covariates[rows, selected, drop = FALSE]),
    error = function(condition) fail("cannot select covariate columns: ", conditionMessage(condition)))
  if (NCOL(response) != 1L || length(response_window) != last || NROW(market_window) != last) {
    fail("response must be univariate and Y/X must share the same original rows.")
  }
  actual_index <- function(object) {
    if (inherits(object, "zoo")) return(zoo::index(object)[rows])
    if (stats::is.ts(object)) return(stats::time(object)[rows])
    NULL
  }
  response_index <- actual_index(response)
  market_index <- actual_index(covariates)
  if (!is.null(response_index) && !is.null(market_index)) {
    mismatches <- which(as.character(response_index) != as.character(market_index))
    if (length(mismatches)) {
      position <- mismatches[1L]
      fail("Y/X index mismatch at original row ", position, " (response ",
           as.character(response_index[position]), "; covariate ",
           as.character(market_index[position]), "); ", length(mismatches), " mismatched rows.")
    }
  }
  check_values <- function(values, labels, dates) {
    if (!is.numeric(values) || is.complex(values)) fail(paste(labels, collapse = ", "), " must be real numeric.")
    positions <- which(!is.finite(values), arr.ind = TRUE)
    if (!NROW(positions)) return(invisible(NULL))
    display <- utils::head(seq_len(NROW(positions)), 8L)
    problems <- vapply(display, function(position) {
      row <- positions[position, 1L]
      column <- positions[position, 2L]
      value <- values[row, column]
      kind <- if (is.nan(value)) "NaN" else if (is.na(value)) "NA" else if (value > 0) "+Inf" else "-Inf"
      date <- if (is.null(dates)) "" else paste0(" (", as.character(dates[row]), ")")
      paste0(labels[column], " at original row ", row, " ", date, ": ", kind)
    }, character(1L))
    fail(paste(problems, collapse = "; "), "; ", NROW(positions), " invalid values total, ",
         NROW(positions) - length(display), " omitted.")
  }
  labels <- colnames(market_window)
  if (is.null(labels)) labels <- paste0("covariate column ", selected)
  check_values(as.matrix(response_window), "response", response_index)
  check_values(market_window, labels, market_index)
  list(response = response_window, covariates = market_window)
}
