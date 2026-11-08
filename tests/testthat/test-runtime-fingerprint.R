runtime_checker <- function() {
  checker <- new.env(parent = globalenv())
  sys.source(system.file("replication", "check_runtime.R", package = "postshock",
    mustWork = TRUE), checker)
  checker
}

test_that("runtime fingerprints use binary UTF-8 with fixed LF bytes", {
  checker <- runtime_checker()
  payloads <- list()
  checker$writeLines <- function(text, con, ...) {
    stop("Platform-dependent text output must not be used for fingerprints.")
  }
  checker$writeBin <- function(object, con, ...) {
    payloads[[length(payloads) + 1L]] <<- object
    base::writeBin(object, con, ...)
  }
  manifest <- checker$check_postshock_runtime()
  expect_length(payloads, nrow(manifest))
  for (position in seq_len(nrow(manifest))) {
    lines <- enc2utf8(deparse(body(get(manifest$function_name[position],
      envir = asNamespace("postshock"))), width.cutoff = 500L))
    expected <- charToRaw(paste0(paste(lines, collapse = "\n"), "\n"))
    expect_identical(payloads[[position]], expected)
    expect_false(as.raw(13L) %in% payloads[[position]])
  }
})

test_that("legacy CRLF hashes reproduce the Windows mismatch without changing bodies", {
  manifest <- runtime_checker()$check_postshock_runtime()
  temporary <- tempfile()
  on.exit(unlink(temporary))
  for (name in manifest$function_name) {
    lines <- enc2utf8(deparse(body(get(name, envir = asNamespace("postshock"))),
      width.cutoff = 500L))
    writeBin(charToRaw(paste0(paste(lines, collapse = "\r\n"), "\r\n")), temporary)
    observed <- unname(tools::md5sum(temporary))
    expect_false(identical(observed, manifest$body_md5[manifest$function_name == name]))
  }
})

test_that("runtime identity checks still reject a different function body", {
  checker <- runtime_checker()
  checker$get <- function(name, envir, ...) {
    observed <- base::get(name, envir = envir, ...)
    if (identical(name, "SynthPrediction")) {
      body(observed) <- quote(stop("Changed implementation for identity-check testing"))
    }
    observed
  }
  expect_error(checker$check_postshock_runtime(),
    "postshock implementation does not match the packaged replication manifest.",
    fixed = TRUE)
  expect_silent(runtime_checker()$check_postshock_runtime())
})

test_that("runtime identity checks still reject a mismatched version", {
  checker <- runtime_checker()
  manifest <- checker$check_postshock_runtime()
  manifest$version <- "0.0.0"
  temporary <- tempfile(fileext = ".csv")
  on.exit(unlink(temporary))
  utils::write.csv(manifest, temporary, row.names = FALSE)
  checker$system.file <- function(...) temporary
  expect_error(checker$check_postshock_runtime())
})
