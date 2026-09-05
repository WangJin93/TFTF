# TFTF changelog

## TFTF 0.2.0 (2025-09)

### New features
- **Local cache for downloaded data** (GCAS-style, see `?get_tftf_cache_dir`):
  - `get_data()` now caches every TFTF API query on disk (default `rappdirs::user_cache_dir("TFTF")/data_temp`, one `.rds` file per query keyed by an md5 of the query parameters) plus an in-session memory layer, so identical queries are answered from the cache instead of hitting the server again.
  - `pantissue_cor_analysis()` caches per-gene expression queries against the Xena servers (`database = "toil"` for TCGA/GTEx and `"ccle"` for CCLE), so each gene is downloaded only once.
  - New control functions: `get_tftf_cache_dir()`, `set_tftf_cache_dir()`, `clear_tftf_cache()`; new arguments `use_cache`, `refresh`, `cache_dir`, `timeout`, `api_url`, `quiet` on `get_data()` and `use_cache`, `cache_dir` on `pantissue_cor_analysis()`.
  - Robustness: failed downloads return `NULL` with a warning (instead of aborting), and a forced refresh falls back to the stale cache when the server is unreachable.

### Bug fixes
- Restore the missing exported function `plot_venn()` (Venn diagram for up to 5 sets, flower/petal plot beyond), which was documented and exported but absent from the package sources; hardened for 1-2 sets and empty inputs.
- Remove duplicated function definitions (`R/intersection.R` duplicated `R/intersections.R`; `R/cor_analysis.R` duplicated `pantissue_cor_analysis()` with an outdated, buggy body).
- `TF_Target_batch()` no longer errors when no TF yields predictions (returns an empty data frame with a message).
- `intersections()` default dataset selection no longer misaligns when some datasets return no rows; empty intersections are consistently reported as `"None"`.
- `tissue_type()` failed at runtime (`tissue[,"TCGA"]` two-dimensional subscript on a list); now returns the correct TCGA (33) / GTEx (30) type vectors.
- Remove the global `options(timeout = 200)` side effect in `predict_target()`/`predict_TF()` (per-request timeouts instead).
- Feedback module of the Shiny app no longer embeds plain-text SMTP credentials (read from environment variables instead).

### Improvements
- Input validation and proper HTTP handling (timeouts, errors) in `get_data()`.
- `get_data_info()` accepts real table names (`FIMO_JASPAR`, `PWMEnrich_JASPAR`) and multiple datasets.
- `TF_Target_batch()` default `PWMEnrich.p` unified to `0.1` (same as `predict_target()`).
- `Find_TF()` no longer hard-codes the number of datasets.
- DESCRIPTION/NAMESPACE dependencies aligned (removed unused `RMySQL`, added `digest`, `rappdirs`, `httr`, `ggrepel`, `psych`, `plotrix`, `shiny`, `utils`; `stats`/`utils` symbols now explicitly imported).
- `R CMD check` is now fully clean (previously 1 ERROR / 3 WARNING / 2 NOTE): dataset documentation for all 9 bundled data objects, `flowerplot()` argument docs, ASCII-only package code, `LICENSE` declared in DESCRIPTION, duplicate case-only file removed from the Shiny app, runtime logs no longer shipped.
- Added `NEWS.md`, `.gitignore`, `.Rbuildignore`, `URL`/`BugReports` fields in DESCRIPTION.

## TFTF 0.1.0

- Initial release: transcription factor target prediction (`predict_target()`), upstream TF prediction (`predict_TF()`), batch prediction (`TF_Target_batch()`), pan-tissue correlation analysis (`pantissue_cor_analysis()` / `viz_cor_results()`), intersection and Venn/flower visualization, and the integrated Shiny application (`TFTF_app()`).
