#' @title Get TF-target data
#' @description
#' Get the TF-target targeting results by using the TFTF API (results stored in
#' a MySQL database on the server side). Downloaded data are cached locally
#' following the design of the GCAS package (see \code{get_tftf_cache_dir()}),
#' so repeated identical queries are answered from the cache instead of hitting
#' the server again.
#' @import jsonlite
#' @param table One of the dataset in "hTFtarget","KnockTF","FIMO_JASPAR","PWMEnrich_JASPAR","ENCODE","CHEA","TRRUST","GTRD","ChIP_Atlas". For correlation analysis, "cor_"+tissue type was the table name, i.e. cor_TCGA, cor_Lung.
#' @param searchType TF or Target.
#' @param gene Any TF names for searchType == TF, and any gene symbols for searchType == Target.
#' @param use_cache Logical (default \code{TRUE}). When \code{TRUE}, results of
#' previous identical queries are returned from the local cache (memory first,
#' then disk) and no request is sent to the server.
#' @param refresh Logical (default \code{FALSE}). When \code{TRUE}, the data are
#' re-downloaded from the server and the cache is updated, even if a cached
#' copy exists.
#' @param cache_dir Optional custom directory for the cache (default
#' \code{NULL}, i.e. the active cache directory, see \code{set_tftf_cache_dir()}).
#' @param timeout Time in seconds allowed for the server request (default 120).
#' @param api_url Optional custom API base URL. \code{NULL} (default) uses
#' \code{https://www.jingege.wang/TFTF/api.php}.
#' @param quiet Logical (default \code{FALSE}). Suppress progress messages
#' about cache reads/writes and downloads.
#' @examples
#' \dontrun{
#' results <- get_data(table = "ENCODE", searchType = "Target", gene = "GAPDH")
#' # the second call is served from the local cache, no request is sent:
#' results2 <- get_data(table = "ENCODE", searchType = "Target", gene = "GAPDH")
#' # force re-download:
#' results3 <- get_data(table = "ENCODE", searchType = "Target", gene = "GAPDH", refresh = TRUE)
#' }
#' @export
#'
get_data <- function(table = "ENCODE",
                     searchType = "Target",
                     gene = "GAPDH",
                     use_cache = TRUE,
                     refresh = FALSE,
                     cache_dir = NULL,
                     timeout = 120,
                     api_url = NULL,
                     quiet = FALSE){
  if (!is.character(table) || length(table) != 1 || !nzchar(table)) {
    stop("table must be a single, non-empty character string.")
  }
  if (!is.character(searchType) || length(searchType) != 1 ||
      !searchType %in% c("TF", "Target")) {
    stop("searchType must be one of 'TF' or 'Target'.")
  }
  if (!is.character(gene) || !length(gene) || all(!nzchar(gene))) {
    stop("gene must be a non-empty character vector.")
  }
  if (length(gene) > 1) gene <- paste0(gene, collapse = ",")

  base_url <- if (is.null(api_url)) {
    "https://www.jingege.wang/TFTF/api.php"
  } else {
    api_url
  }
  parts <- list(table = table, searchType = searchType, searchContent = gene)

  if (!isTRUE(use_cache)) {
    return(.tftf_api_fetch(base_url, parts, timeout = timeout, quiet = quiet))
  }
  .tftf_query_cached(
    parts,
    fetch_fun = function() .tftf_api_fetch(base_url, parts, timeout = timeout,
                                           quiet = quiet),
    use_cache = use_cache, refresh = refresh,
    cache_dir = cache_dir, quiet = quiet, mem_extra = base_url
  )
}

# send one request to the TFTF API and parse the JSON answer; returns NULL
# (with a warning) when the request fails.
.tftf_api_fetch <- function(base_url, parts, timeout = 120, quiet = FALSE) {
  url <- httr::modify_url(base_url, query = parts)
  if (!quiet) message("Fetching data from ", url)
  tryCatch({
    resp <- httr::GET(url,
                      httr::timeout(timeout),
                      httr::user_agent("TFTF (https://github.com/WangJin93/TFTF)"))
    httr::stop_for_status(resp)
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    if (is.null(txt) || !nzchar(txt)) return(NULL)
    jsonlite::fromJSON(txt)
  }, error = function(e) {
    warning("TFTF API request failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}
