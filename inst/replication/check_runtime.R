check_postshock_runtime <- function(verbose = FALSE) {
  manifest_path <- system.file("replication", "implementation.csv", package = "postshock", mustWork = TRUE)
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  stopifnot(all(manifest$version == as.character(utils::packageVersion("postshock"))))
  fingerprint <- function(name) {
    temporary <- tempfile()
    on.exit(unlink(temporary))
    lines <- enc2utf8(deparse(body(get(name, envir = asNamespace("postshock"))), width.cutoff = 500L))
    writeBin(charToRaw(paste0(paste(lines, collapse = "\n"), "\n")), temporary)
    unname(tools::md5sum(temporary))
  }
  observed <- vapply(manifest$function_name, fingerprint, character(1L))
  if (!identical(unname(observed), manifest$body_md5)) {
    stop("postshock implementation does not match the packaged replication manifest.")
  }
  if (isTRUE(verbose)) {
    cat("IMPLEMENTATION VERIFIED:", find.package("postshock"), "\n")
    print(manifest)
    print(sessionInfo())
  }
  invisible(manifest)
}
