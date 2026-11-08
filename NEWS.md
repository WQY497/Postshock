# postshock 0.2.3

## Fixes

* Use a test-only tolerance of 1e-7 instead of 1e-10 (or exact equality) for comparisons of
  independently fitted ARIMAX coefficients and analogous refitted model or
  optimized-weight results. This addresses the MKL discrepancy reported by
  CRAN (2.04085310337 versus 2.04085310316). The tolerance follows testthat's
  numeric comparison semantics; it is not an optimizer stopping criterion.
* Add regression coverage that accepts the reported rounding difference but
  still rejects larger coefficient/weight discrepancies, non-finite results
  and mismatched convergence status. Keep all existing tests, exact input
  checks, convergence requirements and arithmetic-identity checks.
* Fitting code, optimization settings, penalty defaults, frozen data and
  Technical Report resources are unchanged. Update the replication manifest
  version while retaining all function-body fingerprints.

# postshock 0.2.2

## Documentation

* Include the previously verified full Technical Report PDF and its R Markdown
  sources under the installed `technical_report/` directory. Add a
  `technical_report` help topic and access instructions in the README and
  workflow vignette. The report is prebuilt; no extra analysis or dependencies
  are required to read it. Historical previews and rendering caches are excluded.

## Fixes

* Separate PCA integration-test contracts: explicitly use L2 regularization
  for the three-path weight comparison, retain an unpenalized example with
  a unique optimum, and verify the original one-component example's multiple
  optima without requiring a particular solver result. Failed explicit PCA
  fits remain errors. No fitting code, penalty defaults or existing numerical
  tolerances are changed.
* Make the installed Apple replication diagnostic's combined-adjustment
  reference check consistent with its unchanged donor-effect (1e-9) and
  weight (1e-7) absolute tolerances. Propagate these component error bounds,
  account for rounded reference constants, and independently verify the
  weighted sum to 1e-10. This addresses a reproduced cross-BLAS failure of
  the former 1e-9 aggregate reference check without changing forecasts,
  fitting, optimization, frozen inputs or six-decimal article checks.
* Retain version/function-body identity checks, including portable UTF-8/LF
  fingerprints. Public APIs and dependencies are unchanged.

# postshock 0.2.1

## Fixes

* Forward explicit PCA settings from `SynthPrediction()` to `dbw()`,
  including the MultiPool dispatcher. Validate dimensions against the
  available SVD directions. Explicit PCA failures, unsuccessful or unknown
  DBW status, and invalid weights now produce errors instead of silently
  substituting equal weights; the default NULL and disabled-DBW paths remain.
* Forward all shared model, matching, penalty, seasonal and control settings
  through the MultiPool dispatcher. Reject conflicting settings while
  preserving separate baseline/pool DBW switches and unnormalized pool weights.
* Require complete, finite, aligned inputs when additional donor ARIMAX
  covariates are selected. Report invalid rows and available dates rather
  than deleting covariate rows and misaligning responses. Other input-policy
  branches retain their existing behavior.
* Restrict automatic GARCH-X selection to fits with numerical convergence
  code zero, finite nonempty coefficients and finite BIC. Preserve the first
  eligible fit on exact BIC ties; report an error if none qualifies.
  Backcast refits must meet the same conditions and strictly improve BIC.
  Non-finite standard errors alone do not disqualify a point fit.
* Pass the specified or detected per-window period to automatic ARIMA's
  response frequency, with strict validation of explicit periods. Retain
  the non-seasonal ARIMA(1,1,1) donor failure fallback.
* Correct and align documentation for GARCH-X order `c(q, p, o)`:
  `q` is GARCH, `p` is ARCH, and `o` is the asymmetry order.

## Reproducibility

* Add `aapl_s1_input`, derived from the unchanged frozen Apple objects,
  with evaluation observations separate from fitting and matching inputs.
* Install the Apple input builder, shared analysis and version/function-body
  identity checks under `replication/`. The Technical Report uses these
  installed resources and executable results without an adjacent manuscript.
* Hash replication function bodies as UTF-8 bytes with fixed LF line endings,
  avoiding false implementation-mismatch errors on Windows. Retain the
  version, implementation-identity and numerical checks.

## Compatibility

Public function signatures and existing frozen data are unchanged. Previously
ignored explicit settings may now change forecasts; invalid donor ARIMAX inputs,
failed explicit PCA requests and ineligible automatic GARCH-X fits can now
produce errors. These corrections preserve the article's existing main
application specifications.
