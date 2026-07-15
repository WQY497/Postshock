#' Helper Function: Build control-shock xreg matrix
#'
#' @description
#' Internal helper to convert shock specifications into a design matrix for ARIMA xreg.
#'
#' @param last Integer; length of the training period.
#' @param control_spec List or data.frame defining shock timing, length, and shape.
#'
#' @return A numeric matrix or NULL.
#' @export
build_control_xreg <- function(last, control_spec = NULL) {
  # 1. Basic checks
  if (is.null(last) || length(last) != 1L || !is.finite(last) || last < 1) {
    stop("`last` must be a single positive integer.")
  }
  last <- as.integer(last)
  
  # 2. Return NULL if no controls
  if (is.null(control_spec)) return(NULL)
  
  # 3. Normalize input to a list of specs
  specs <- NULL
  if (is.data.frame(control_spec)) {
    if (!all(c("time") %in% names(control_spec))) {
      stop("If `control_spec` is a data.frame, it must contain a `time` column.")
    }
    # Fill defaults
    if (!("length" %in% names(control_spec))) control_spec$length <- 1L
    if (!("shape"  %in% names(control_spec))) control_spec$shape  <- "point"
    
    specs <- lapply(seq_len(nrow(control_spec)), function(i) {
      list(
        time   = control_spec$time[[i]],
        length = control_spec$length[[i]],
        shape  = control_spec$shape[[i]]
      )
    })
  } else if (is.list(control_spec)) {
    # A named list containing any recognized field is treated
    # as one control-shock specification, even if `time` is missing.
    spec_fields <- c("time", "length", "shape")
    
    is_single_spec <-
      !is.null(names(control_spec)) &&
      any(names(control_spec) %in% spec_fields)
    
    if (is_single_spec) {
      specs <- list(control_spec)
    } else {
      specs <- control_spec
    }
  } else {
    stop("`control_spec` must be NULL, a list, or a data.frame.")
  }
  
  # Filter empty specifications
  specs <- Filter(
    function(s) !is.null(s) && length(s) > 0L,
    specs
  )
  
  if (length(specs) == 0L) {
    return(NULL)
  }
  
  # Build one design-matrix column per control shock
  cols <- lapply(seq_along(specs), function(j) {
    s <- specs[[j]]
    
    if (!is.list(s)) {
      stop("Each control shock specification must be a list.")
    }
    
    t0 <- if (!is.null(s$time)) {
      as.integer(s$time)
    } else {
      NA_integer_
    }
    
    L <- if (!is.null(s$length)) {
      as.integer(s$length)
    } else {
      1L
    }
    
    shp <- if (!is.null(s$shape)) {
      tolower(as.character(s$shape))
    } else {
      "point"
    }
    
    if (length(t0) != 1L || is.na(t0) || !is.finite(t0)) {
      stop("Control shock must have a valid `time`.")
    }
    
    if (length(L) != 1L || is.na(L) || !is.finite(L) || L < 1L) {
      stop("Control shock `length` must be >= 1.")
    }
    
    if (length(shp) != 1L ||
        !shp %in% c("point", "window", "step")) {
      stop("Unknown shock shape: ", paste(shp, collapse = ", "))
    }
    
    v <- rep(0, last)
    
    if (shp == "point") {
      if (t0 >= 1L && t0 <= last) {
        v[t0] <- 1
      }
    } else if (shp == "window") {
      a <- max(1L, t0)
      b <- min(last, t0 + L - 1L)
      
      if (a <= b) {
        v[a:b] <- 1
      }
    } else if (shp == "step") {
      a <- max(1L, t0)
      
      if (a <= last) {
        v[a:last] <- 1
      }
    }
    
    v
  })
  
  C <- do.call(cbind, cols)
  if (is.null(C)) return(NULL) # Safety if cols was empty
  colnames(C) <- paste0("ctrl_shock_", seq_len(ncol(C)))
  storage.mode(C) <- "double"
  return(C)
}