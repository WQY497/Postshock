
# postshock

Source for `postshock` 0.2.3, available on [CRAN](https://CRAN.R-project.org/package=postshock).

<!-- badges: start -->

<!-- badges: end -->

`postshock` provides donor-adjusted methods for forecasting conditional
means and variances after structural shocks. Historical donor episodes
are weighted using information observed before each shock, and their
estimated post-shock effects are combined with forecasts from a
target-series model.

The package builds on the methodology of Lin and Eck (2021) and extends
it with practical tools for donor weighting, conditional mean
forecasting, volatility forecasting, structured donor pools, and
automatic GARCH-X model selection.

------------------------------------------------------------------------

## Main Features

### Donor balancing weights (`dbw()`)

Estimate donor weights from pre-shock covariates so that weighted donor
characteristics resemble those of the target series.

Features include:

- optional scaling and centering
- lower and upper weight bounds
- L1 and L2 regularization
- uniform-weight benchmark

------------------------------------------------------------------------

### Conditional mean forecasting (`SynthPrediction()`)

Generate donor-adjusted forecasts by

- fitting an ARIMA/SARIMA model to the target series,
- estimating donor-specific post-shock effects,
- combining donor effects using balancing weights, and
- adjusting the target forecast.

------------------------------------------------------------------------

### Structured donor pools

Support multiple donor pools representing different economic or
substantive mechanisms, allowing donor effects to be estimated and
combined within predefined groups.

------------------------------------------------------------------------

### Automatic GARCH-X model selection (`auto_garchx()`)

Automatically search across GARCH-X and GJR-GARCH specifications using
information criteria.

Supported features include

- external regressors,
- robust covariance estimation,
- backcast handling, and
- automatic order selection.

------------------------------------------------------------------------

### Conditional variance forecasting (`SynthVolForecast()`)

Forecast post-shock conditional variances by combining

- a target-series GARCH-X model,
- donor-estimated volatility adjustments, and
- balancing weights estimated from pre-shock information.

------------------------------------------------------------------------

### Covariate construction (`build_control_xreg()`)

Construct control variables and external regressors for use in
forecasting workflows.

------------------------------------------------------------------------

## Installation

Install the released version from CRAN:

``` r
install.packages("postshock")
```

------------------------------------------------------------------------

## Getting Started

Load the package

``` r
library(postshock)
```

The package vignette provides a complete workflow covering

- donor balancing weights,
- conditional mean forecasting,
- structured donor pools,
- volatility forecasting, and
- practical guidance for applied analyses.

Open the vignette with

``` r
vignette("postshock-workflow", package = "postshock")
```

### Full Technical Report

The package also includes the full Technical Report, covering the methods,
simulation examples, and COP, IYG, and Apple applications. Open the installed
PDF without rerunning any analysis:

``` r
report <- system.file("technical_report", "Technical_report.pdf",
                      package = "postshock", mustWork = TRUE)
utils::browseURL(report)
```

Use `help("technical_report", package = "postshock")` for details. The verified
PDF and its source files are bundled separately from the shorter workflow
vignette; reading the PDF does not require R Markdown or LaTeX.

------------------------------------------------------------------------

## Main Functions

| Function | Description |
|----|----|
| `dbw()` | Estimate donor balancing weights |
| `SynthPrediction()` | Donor-adjusted conditional mean forecasting |
| `SynthPredictionMultiPool()` | Forecasting with structured donor pools |
| `auto_garchx()` | Automatic GARCH-X model selection |
| `SynthVolForecast()` | Donor-adjusted conditional variance forecasting |
| `build_control_xreg()` | Construct control variables and external regressors |

------------------------------------------------------------------------

## Reference

Lin, Y. and Eck, D. J. (2021).

*Post-shock forecasting using synthetic controls.*

International Journal of Forecasting.

<https://doi.org/10.1016/j.ijforecast.2021.03.010>

------------------------------------------------------------------------

## License

MIT License.
