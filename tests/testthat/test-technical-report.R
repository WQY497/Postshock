test_that("the full technical report and its sources are installed", {
  directory <- system.file("technical_report", package = "postshock",
                           mustWork = TRUE)
  required <- c("Technical_report.pdf", "Technical_report.Rmd", "sim-test.Rmd",
                "sim_vol_tets.Rmd", "README.txt")
  expect_true(all(file.exists(file.path(directory, required))))
  report <- file.path(directory, "Technical_report.pdf")
  expect_identical(readBin(report, "raw", n = 5L), charToRaw("%PDF-"))
  expect_length(utils::help("technical_report", package = "postshock"), 1L)
})
