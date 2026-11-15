# Data access and preparation

## Restricted GBD data

Obtain estimates directly through IHME's official GBD Results tool (https://vizhub.healthdata.org/gbd-results/) using an authorized account and accepting the applicable terms. Do not use web scraping, authenticated endpoint extraction, shared credentials, signed download links or unofficial mirrors. Do not bypass download ceilings or collaborator/publication requirements. This repository has no GBD downloader.

Use **GBD 2021**, Incidence, Rate per 100,000, all 204 countries and territories, years 1992–2021, causes Chronic kidney disease and Depressive disorders. Export age-standardized estimates for Both/Female/Male, plus five-year age groups 5–9 through 90–94 and 95+ for Female/Male. The analyses use point estimates; uncertainty draws were unavailable. GBD 2023 was not analyzed. The study's access constraints do not grant permission to new users; obtain any required IHME approval for your own intended research/publication.

1. Save cause-specific standard IHME CSV exports locally. `scripts/prepare_gbd_inputs.py` filters them and rejects duplicated location/sex/year records. It does not connect to the internet.
2. Place `CKD_ASIR_1992_2021.csv`, `DD_ASIR_1992_2021.csv`, `CKD_age_2021.csv`, `DD_age_2021.csv` in `data/private` (or set `CKDDD_PRIVATE_DATA`). These preserve the IHME columns listed in `private_input_manifest.json`. Each age-standardized file has 18,360 records (204 × 3 × 30). Only incidence-rate rows belong in these inputs.
3. Prepare `CKD_DD_female_LONG_1992_2021.csv` and `CKD_DD_male_LONG_1992_2021.csv` with columns `location_name`, `age_name`, `year`, `ckd_val`, `dd_val` for the same 19 age groups. The preparation helper writes them from matched age-specific cause exports. Extra age groups in earlier archived matrices are discarded by analysis.
4. Obtain the GBD-2021-compatible SDI location-year series directly from IHME/GHDx (https://ghdx.healthdata.org/) and save `SDI.csv` with `location_name,val_1992,...,val_2021`. The helper validates this schema but does not redistribute SDI. `Population.csv` has the same wide schema and contains the archived country total population series used for the population-adjusted spatial sensitivity and illustrative selection. See the input manifest hashes to distinguish snapshots.

`scripts/validate_inputs.py` checks schemas, unique keys, 204 locations, annual coverage and matched positive rates. A fresh official export can have different file bytes; manifest hashes identify the exact study snapshot rather than asserting all re-downloads must be byte-identical.

## PM2.5 and reference points

The primary spatial models use a frozen World Bank WDI PM2.5 matrix (EN.ATM.PM25.MC.M3; μg/m³), matched to 204 study locations. The historical 1992–2020 series was retained and 2021 estimated separately for every location with non-seasonal `forecast::auto.arima`, using a mean fallback on forecast failure. Taiwan uses the China WDI series; Cook Islands, Niue and Tokelau use the Pacific island small states aggregate. The other 200 locations use their own matched WDI entities. Preserve those mapping rules and save `PM25.csv` as `location_name,val_1992,...,val_2021`. The historical snapshot is not assigned a newer WDI/GBD upstream version. A current WDI download is a new sensitivity input, not an exact replacement for the archived matrix.

An optional official WDI API helper is supplied in `scripts/fetch_wdi_pm25.py`. It is **offline by default**, requires `--online` for network access, sends serial requests at least 10 seconds apart, caches complete responses, honors Retry-After and backs off on rate limiting. It prohibits crawlers and follows the official API rules (https://datahelpdesk.worldbank.org/knowledgebase/topics/125589-developer-information). It never fetches GBD data and never overwrites the study snapshot. Users must stop if access restrictions or service terms prohibit the requested use.

City and land-geometry reference point tables are public, contain no disease estimates and include source attribution in `THIRD_PARTY_NOTICES.md`. They are supplied so the spatial sensitivity is reproducible without downloading map archives. The 10m geometry derivation can optionally be rerun from manually obtained official shapefiles; this does not alter the prespecified reference points used in the paper.

## Public and private results

`reference_results` contains aggregate statistics reported in the manuscript/supplement. Full country-level derived tables are delivered with the journal's Supplementary Data S1. This code repository does not redistribute raw GBD estimates, restricted downloads, saved authenticated responses or author-local working files. Generated results remain ignored by Git to prevent accidental publication of private inputs and row-level derivative matrices.
