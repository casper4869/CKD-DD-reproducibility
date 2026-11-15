# Run Rscript from the repository root. No network operations are performed.
root <- normalizePath('.', winslash='/', mustWork=TRUE)
if (!file.exists(file.path(root,'CITATION.cff'))) stop('Run from the repository root.')
private <- Sys.getenv('CKDDD_PRIVATE_DATA', unset=file.path(root,'data','private'))
result <- Sys.getenv('CKDDD_RESULTS', unset=file.path(root,'results'))
core_dir <- file.path(result,'core'); temporal_dir <- file.path(result,'temporal')
spatial_dir <- file.path(result,'spatial'); figure_dir <- file.path(result,'figures')
reference_dir <- file.path(root,'data','reference')
for (p in c(result, core_dir, temporal_dir, spatial_dir, figure_dir)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
