# CKD-DD reproducibility

Analysis code for **Population-level co-burden and spatiotemporal patterns of chronic kidney disease and depressive disorders across 204 countries and territories**, submitted to *International Urology and Nephrology*.

The study uses GBD 2021 incidence estimates for 1992–2021. It describes ecological co-burden, age patterns, exploratory temporal predictability and spatial heterogeneity. A significant temporal test is not evidence of patient-level causation; the VAR diagnostic limitations are reported with the results. Improvement over a last-observation forecast does not isolate any contribution of depressive-disorder terms.

## Contents

- `scripts`: portable analysis and plotting code; no manuscript generators.
- `data/reference`: attributed city/reference and Natural Earth geometry points, without disease data.
- `reference_results`: aggregate published summaries for numerical verification.
- `docs/DATA_ACCESS.md`: official acquisition steps and required private input schemas.
- `docs/private_input_manifest.json`: hashes and schemas of the original local input snapshots.
- `docs/figure_source_map.csv`: figure-to-result mapping, using the final S2 age / S3 scenario numbering.
- `docs/analysis_map.csv`: analysis-to-supplement mapping; `docs/supplement_sheet_map.csv` maps the 78 data/provenance sheets to Tables S1–S11.
- `CITATION.cff` and `.zenodo.json`: release citation/archival metadata. No DOI is asserted until the archive is published.

## Environment

The executed analysis environment is R 4.5.0. Install the versions recorded in `docs/dependency_versions.csv` in your own R library. No package installation or online data acquisition occurs during analysis. Python helpers use Python 3.11+ with pandas and numpy; the age-reference calculation uses Python 3.12.14, numpy 2.4.4 and pandas 2.3.3 in its executed environment; input validation and manifest audits use only the standard library.

Run commands **from the repository root**. Place authorized data under `data/private`, or set `CKDDD_PRIVATE_DATA` to an external directory. `CKDDD_RESULTS` may set an external output directory; the default `results` is ignored by Git.

```sh
python scripts/validate_inputs.py
Rscript scripts/core_analysis.R
python scripts/run_age_reference.py --ckd data/private/CKD_age_2021.csv --dd data/private/DD_age_2021.csv --output results/age_reference
Rscript scripts/temporal_models.R
Rscript scripts/temporal_diagnostics.R
Rscript scripts/temporal_summary.R
Rscript scripts/arimax_scenarios.R 1
Rscript scripts/arimax_scenarios.R 2
Rscript scripts/summarize_original_arimax.R
Rscript scripts/arimax_scenarios.R 1 raw
Rscript scripts/arimax_scenarios.R 2 raw
Rscript scripts/summarize_raw_arimax.R
Rscript scripts/spatial_age_analysis.R
Rscript scripts/verify_spatial_math.R
Rscript scripts/align_moran_reporting.R
python scripts/summarize_spatial_age_evidence.py
Rscript scripts/figures_descriptive.R
Rscript scripts/figures_models.R
Rscript scripts/figure6_local_signals.R
Rscript scripts/figureS2_age_temporal.R
python scripts/check_reference_results.py
```

To reproduce the additional differenced VARX / restricted-order ARIMAX sensitivity forecast matrices and their all-model symmetric evidence summaries, also run:

```sh
Rscript scripts/forecast_sensitivities.R
python scripts/summarize_temporal_evidence.py
```

The differenced VARX sensitivity uses 1,000 residual paths per forecast and can take substantially longer than the primary analyses. It is not the primary lag-2 forecast in Figure 4. All prespecified origins are retained.

The two ARIMAX partitions are deterministic halves of the full country set; run both. The raw-scale branch is an exploratory sensitivity, not the manuscript's log1p scenario model. For an external `CKDDD_PRIVATE_DATA`/`CKDDD_RESULTS` directory, pass its explicit paths to the age helper as `--ckd`, `--dd` and `--output`; that helper deliberately uses CLI paths rather than environment variables. `figures_models.R` also reruns the descriptive plotting stage so its shared plotting functions are available.

Internal filename prefixes (S/A) identify reproducibility outputs; published supplementary table numbers are mapped explicitly in `docs/analysis_map.csv` and in Supplementary Data S1.

### Key analysis definitions

- `core_analysis.R` joins incidence and SDI, ranks countries, describes 2021 age profiles, runs the age bootstrap, produces differenced/log-level temporal sensitivities and the original-point GWR fits. Its output field `spec=primary` in **S9** is a retained machine identifier for the differenced sensitivity; it is **not** the manuscript's primary temporal model.
- `temporal_models.R` defines the primary **A1 `original_raw_endo`** model: raw CKD, DD and SDI endogenous; constant; AIC over lags 1–2. Disease-equation tests of CKD→DD and DD→CKD are distinct from system-wide tests. BH/BY adjustment pools all 408 disease-direction/country tests within a specification. HC3 uses the actual fitted VAR design and adjusts heteroskedasticity, not residual serial dependence.
- The forecast comparator is fixed lag-2 log1p VARX with contemporaneous SDI; the default automatic log1p ARIMAX and last-observation benchmarks use the same origins/horizons. No forecast clipping is applied to the reported validation. Conditional scenario intervals omit parameter and GBD-estimation uncertainty.
- Spatial models standardize variables over the 204 locations, use fixed Gaussian distance kernels selected by leave-one-out CV and great-circle distances, and compare both coordinate choices with SDI+PM2.5 and SDI-only adjustment. Local significance is approximate and conditional; the all-four-setting intersection does not itself carry a joint FDR guarantee.
- The editor-requested `run_age_reference.py` fits each country/sex/disease’s 19 rates on a disease-specific equal-country-weighted age profile from the other 203 countries plus an intercept, then correlates the two residual vectors. This is a country-average reference, not an official population-weighted GBD Global series. The fixed contract is distributed beside the script. It produces all 408 profiles, reference values/residuals and a sex-specific summary.
- Country bootstrap intervals are conditional descriptive intervals and do not propagate GBD estimation uncertainty or spatial dependence.

## Network and data policy

**No crawlers.** The entire analysis pipeline runs offline. Restricted GBD exports must be acquired manually through IHME's official tools under the user's own authorization. The optional WDI API helper is disabled without `--online`; requests are serial, cached and delayed, and respect rate limits and Retry-After. It is not invoked by any analysis command. Never commit credentials, tokens, signed URLs, private data or HTTP response caches.

## License and citation

Author-written code is MIT licensed; third-party data and dependencies retain their own terms. See `THIRD_PARTY_NOTICES.md`. Cite the versioned archived software release using `CITATION.cff`. The persistent archive identifier will be added only after the author confirms the GitHub Release and Zenodo archive.
