# =============================================================================
# Local data cache for TFTF
#
# Downloads performed by this package (TF-target queries against the TFTF API
# via get_data(), and per-gene expression queries against the Xena servers via
# UCSCXenaShiny::query_pancan_value() in pantissue_cor_analysis()) can be
# expensive and are frequently repeated: the same TF/target is searched again
# and again, e.g. when several analysis functions are called with overlapping
# inputs or when the built-in Shiny app is used interactively.
#
# Following the design of the GCAS package (WangJin93/GCAS,
# get_expr_data()/get_OSF_data()), every network result is stored on disk in a
# per-user cache directory:
#
#   * default location: rappdirs::user_cache_dir("TFTF")/data_temp
#     (Linux:  ~/.cache/TFTF/data_temp,
#      macOS:  ~/Library/Caches/TFTF/data_temp,
#      Windows: C:\Users\<user>\AppData\Local\TFTF\TFTF\data_temp)
#   * overridable with set_tftf_cache_dir(), with the TFTF_CACHE_DIR
#     environment variable, or with a cache_dir argument on the fetching
#     functions.
#   * one .rds file per query; the file name embeds an md5 hash of the
#     normalised query parameters (mirroring GCAS's
#     paste0(dataset, "_", digest::digest(genes), ".RData") scheme).
#
# A small in-memory cache is layered on top so that repeated identical queries
# inside one R session (or one Shiny session) are answered instantly.
# =============================================================================

# version salt: bump this string whenever the cached format / semantics change
.tftf_cache_version <- "tftf-v1"

# session-level in-memory cache (keyed by the same md5 used for file names)
.tftf_mem_cache <- new.env(parent = emptyenv())

#' @title Get the TFTF data cache directory
#' @description
#' Print / return the directory used by TFTF to store locally cached downloads
#' (mirroring the cache layout of the GCAS package). The cache directory is
#' resolved in the following order:
#' \enumerate{
#'   \item the \code{cache_dir} argument passed to the fetching function;
#'   \item \code{options(TFTF.cache_dir = ...)};
#'   \item the \code{TFTF_CACHE_DIR} environment variable;
#'   \item \code{rappdirs::user_cache_dir("TFTF")} (default).
#' }
#' @return Path to the cache folder (invisibly).
#' @examples
#' \dontrun{
#' get_tftf_cache_dir()
#' }
#' @export
get_tftf_cache_dir <- function() {
  invisible(.tftf_cache_data_dir(cache_dir = NULL))
}

#' @title Set the TFTF data cache directory
#' @description
#' Override the directory used to cache downloaded data for the current R
#' session. The setting is stored in \code{options(TFTF.cache_dir)} and affects
#' all subsequent data downloads. Use \code{set_tftf_cache_dir(NULL)} to reset
#' to the default (\code{rappdirs::user_cache_dir("TFTF")}).
#' @param dir Path of the new cache directory. \code{NULL} (default) resets to
#' the default user cache directory.
#' @return The new cache directory path (invisibly).
#' @examples
#' \dontrun{
#' set_tftf_cache_dir("~/TFTF_cache")
#' set_tftf_cache_dir(NULL)  # back to default
#' }
#' @export
set_tftf_cache_dir <- function(dir = NULL) {
  if (is.null(dir)) {
    options(TFTF.cache_dir = NULL)
  } else {
    if (!is.character(dir) || length(dir) != 1 || !nzchar(dir)) {
      stop("dir must be a single, non-empty character string or NULL.")
    }
    options(TFTF.cache_dir = normalizePath(path.expand(dir), mustWork = FALSE))
  }
  invisible(.tftf_cache_data_dir(cache_dir = NULL))
}

#' @title Clear the TFTF data cache
#' @description
#' Remove cached .rds files (and the in-memory cache) produced by
#' \code{get_data()} and \code{pantissue_cor_analysis()}. This can be useful
#' when the underlying databases have been updated and you want the data to be
#' downloaded again (otherwise use \code{refresh = TRUE} on the fetching
#' function).
#' @param cache_dir Optional cache directory to clean; default \code{NULL} uses
#' the currently active cache directory.
#' @param pattern Optional regular expression; when given, only cached files
#' whose name matches are removed.
#' @param quiet Logical, whether to suppress the summary message.
#' @return Number of removed files (invisibly).
#' @examples
#' \dontrun{
#' clear_tftf_cache()
#' clear_tftf_cache(pattern = "hTFtarget")
#' }
#' @export
clear_tftf_cache <- function(cache_dir = NULL, pattern = NULL, quiet = FALSE) {
  # clear the in-memory cache first
  rm(list = ls(envir = .tftf_mem_cache, all.names = TRUE), envir = .tftf_mem_cache)

  dir <- .tftf_cache_data_dir(cache_dir)
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  if (!is.null(pattern)) {
    files <- files[grepl(pattern, basename(files))]
  }
  n <- 0L
  if (length(files) > 0) {
    n <- length(files)
    unlink(files)
  }
  if (!quiet) message(n, " cached file(s) removed from ", dir)
  invisible(n)
}

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

.tftf_cache_root <- function(cache_dir = NULL) {
  if (!is.null(cache_dir)) {
    root <- normalizePath(path.expand(cache_dir[1]), mustWork = FALSE)
  } else if (!is.null(getOption("TFTF.cache_dir"))) {
    root <- getOption("TFTF.cache_dir")
  } else if (nzchar(Sys.getenv("TFTF_CACHE_DIR"))) {
    root <- normalizePath(Sys.getenv("TFTF_CACHE_DIR"), mustWork = FALSE)
  } else {
    root <- rappdirs::user_cache_dir("TFTF")
  }
  if (!dir.exists(root)) {
    if (!dir.create(root, recursive = TRUE, showWarnings = FALSE)) {
      warning("Could not create cache directory '", root,
              "', falling back to '", file.path(tempdir(), "TFTF-cache"), "'.")
      root <- file.path(tempdir(), "TFTF-cache")
      dir.create(root, recursive = TRUE, showWarnings = FALSE)
    }
  }
  root
}

# cache files live under <root>/data_temp (same sub-folder name as GCAS)
.tftf_cache_data_dir <- function(cache_dir = NULL) {
  dir <- file.path(.tftf_cache_root(cache_dir), "data_temp")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  dir
}

# path of the cache file for a query; hash is computed over the normalised
# query parameters (+ the format version salt) so that different queries never
# collide, exactly in the spirit of GCAS's <table>_<md5>.RData naming.
.tftf_cache_file <- function(parts, cache_dir = NULL) {
  parts <- as.list(parts)
  prefix <- gsub("[^[:alnum:]_-]", "_", as.character(parts[[1]]))
  hash <- digest::digest(c(.tftf_cache_version, parts), algo = "md5")
  file.path(.tftf_cache_data_dir(cache_dir), paste0(prefix, "_", hash, ".rds"))
}

.tftf_disk_read <- function(file, quiet = FALSE) {
  if (!file.exists(file)) return(NULL)
  if (!quiet) message("Loading cached data from ", file)
  tryCatch(
    readRDS(file),
    error = function(e) {
      warning("Cached file ", file, " is corrupt and will be re-downloaded: ",
              conditionMessage(e))
      unlink(file)
      NULL
    }
  )
}

.tftf_disk_write <- function(file, value, quiet = FALSE) {
  if (!dir.exists(dirname(file))) dir.create(dirname(file), recursive = TRUE)
  saveRDS(value, file = file)
  if (!quiet) message("Data cached to ", file)
  invisible(file)
}

.tftf_mem_get <- function(key) {
  if (exists(key, envir = .tftf_mem_cache, inherits = FALSE)) {
    get(key, envir = .tftf_mem_cache, inherits = FALSE)
  } else {
    NULL
  }
}

.tftf_mem_set <- function(key, value) {
  assign(key, value, envir = .tftf_mem_cache)
  invisible(value)
}

# Generic get-or-fetch: return a cached value when available (memory first,
# then disk), otherwise run fetch_fun(), store the result in cache and return
# it. fetch_fun() must return NULL to signal a failed download.
# Disk file names are content-addressed (same query -> same file in every
# cache directory), while in-memory entries are additionally scoped by
# cache_dir and mem_extra (e.g. the API endpoint), so redirecting the cache or
# the endpoint never returns a stale in-memory hit.
.tftf_query_cached <- function(parts, fetch_fun, use_cache = TRUE, refresh = FALSE,
                               cache_dir = NULL, quiet = FALSE, mem_extra = NULL) {
  parts <- as.list(parts)
  file <- .tftf_cache_file(parts, cache_dir = cache_dir)
  mem_key <- digest::digest(list(file = file, extra = mem_extra), algo = "md5")
  tag <- basename(file)

  if (use_cache && !refresh) {
    mem <- .tftf_mem_get(mem_key)
    if (!is.null(mem)) {
      if (!quiet) message("Loading cached data from memory (", tag, ")")
      return(mem)
    }
    disk <- .tftf_disk_read(file, quiet = quiet)
    if (!is.null(disk)) {
      .tftf_mem_set(mem_key, disk)
      return(disk)
    }
  }

  value <- fetch_fun()
  if (is.null(value)) {
    # download failed -> fall back to an existing (stale) cache when refreshing
    if (refresh && file.exists(file)) {
      stale <- .tftf_disk_read(file, quiet = TRUE)
      if (!is.null(stale)) {
        warning("Download failed; returning the previously cached data instead.",
                call. = FALSE)
        return(stale)
      }
    }
    return(NULL)
  }

  if (use_cache) {
    .tftf_disk_write(file, value, quiet = quiet)
    .tftf_mem_set(mem_key, value)
  }
  value
}

# Cache wrapper around UCSCXenaShiny::query_pancan_value() (used by
# pantissue_cor_analysis()). The (possibly large) per-gene expression vector
# is downloaded only once per molecule / database / data type and then served
# from the local cache.
.tftf_xena_query_cached <- function(molecule, data_type = "mRNA", database = "toil",
                                    use_cache = TRUE, refresh = FALSE,
                                    cache_dir = NULL, quiet = TRUE) {
  parts <- list(source = "xena", database = database,
                data_type = data_type, molecule = molecule)
  if (!use_cache) {
    return(UCSCXenaShiny::query_pancan_value(molecule, data_type = data_type,
                                             database = database))
  }
  .tftf_query_cached(
    parts,
    fetch_fun = function() {
      UCSCXenaShiny::query_pancan_value(molecule, data_type = data_type,
                                        database = database)
    },
    use_cache = use_cache, refresh = refresh,
    cache_dir = cache_dir, quiet = quiet
  )
}
