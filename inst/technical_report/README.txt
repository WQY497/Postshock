Full Technical Report
=====================

Technical_report.pdf is the complete 46-page report verified on 7 October 2026
with postshock 0.2.1. Its contents are unchanged in this documentation bundle.
It is not Technical_report-apple-data-preview.pdf or a partial replacement.
See the package NEWS for later maintenance changes, including the 0.2.2 Apple
diagnostic checks; the PDF is not presented as a new 0.2.2 analysis run.

After installing the package, locate and open the PDF from R:

  report <- system.file("technical_report", "Technical_report.pdf",
                        package = "postshock", mustWork = TRUE)
  utils::browseURL(report)

The PDF can also be opened directly in any PDF viewer. No computation, network
access, R Markdown, or TeX installation is needed just to read it.
Use help("technical_report", package = "postshock") for the package help topic.

Source files included:
  Technical_report.Rmd - full report with embedded executable examples
  sim-test.Rmd - exploratory mean-simulation notebook
  sim_vol_tets.Rmd - exploratory variance-simulation notebook (original name)

The main report does not source either notebook. Old previews, chunk caches,
and intermediate figure directories are deliberately not distributed.
To regenerate the main PDF, copy the source to a writable working directory
and render it using rmarkdown, knitr, ggplot2, Pandoc, and XeLaTeX, with the
intended postshock package installed. Do not render into the package library.
The sources use frozen package data and installed replication resources;
they do not require the neighboring journal manuscript directory.
