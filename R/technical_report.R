#' Read the full Technical Report
#'
#' Locate the bundled PDF and sources for the full postshock Technical Report,
#' covering methods, simulations, and the COP, IYG, and Apple applications.
#'
#' @name technical_report
#' @details
#' The report is installed in \code{technical_report/} and can be located with
#' \code{system.file()}. It is a prebuilt document, distinct from the shorter
#' \code{postshock-workflow} vignette. Reading it requires a PDF viewer, not
#' R Markdown, LaTeX, an internet connection, or rerunning any analysis.
#'
#' The verified PDF is retained unchanged from the 0.2.1 documentation build.
#' See the package NEWS for subsequent maintenance changes, including the
#' 0.2.2 Apple diagnostic checks. The R Markdown report and two exploratory
#' simulation notebooks are included for reference; the notebooks are not
#' required to read or regenerate the main report. Regeneration requires
#' \pkg{rmarkdown}, \pkg{knitr}, \pkg{ggplot2}, Pandoc, and XeLaTeX in addition
#' to the installed package and its dependencies.
#'
#' @seealso \code{utils::vignette("postshock-workflow", package = "postshock")}
#' @examples
#' report <- system.file("technical_report", "Technical_report.pdf",
#'                       package = "postshock", mustWork = TRUE)
#' file.exists(report)
#' \dontrun{
#' utils::browseURL(report)
#' }
NULL
