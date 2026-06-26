#' @title Perform Robust Rank Aggregation (RRA) analysis on DEGs lists
#' @description
#' This function performs RRA analysis on differentially expressed genes (DEGs) lists
#' obtained from various studies. It ranks genes based on their differential expression
#' and aggregates the ranks to identify consistently regulated genes across studies.
#' @param DEGs_lists A list of DEGs data frames. Each data frame should contain at least a 'gene' column and a 'logFC' column.
#' @param top.num Numeric, the number of top genes to select based on their ranks. Default is 0, which selects all genes passing the thresholds.
#' @param rra.p Numeric, the p-value threshold for RRA. Default is 0.05.
#' @param logFC_cut Numeric, the log fold change threshold for filtering genes. Default is 1.
#' @param p_cut Numeric, the p-value threshold for filtering genes. Default is 0.05.
#' @return A list containing the number of up- and down-regulated genes and a data matrix of aggregated log fold changes.
#' @examples
#' \dontrun{
#' df1 <- get_OSF_data(table = "GSE31210", action = "geo_data")
#' results1 <- DEGs_analysis(df1)
#' df2 <- get_OSF_data(table = "GSE19188", action = "geo_data")
#' results2 <- DEGs_analysis(df2)
#' DEGs_lists <- list("GSE31210" = results1, "GSE19188" = results2)
#' RRA_results <- RRA_analysis(DEGs_lists)
#' ComplexHeatmap::pheatmap(RRA_results$RRA_results)
#' }
#' @export
RRA_analysis <- function(DEGs_lists, top.num = 0, rra.p = 0.05, logFC_cut = 1, p_cut = 0.05) {
  if (!requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    stop("The 'RobustRankAggreg' package is required for this function. Please install it with: install.packages('RobustRankAggreg')")
  }
  
  # Input validation
  if (is.null(DEGs_lists)) {
    warning("DEGs_lists is NULL")
    return(NULL)
  }
  
  if (!is.list(DEGs_lists)) {
    stop("DEGs_lists must be a list")
  }
  
  if (length(DEGs_lists) < 2) {
    stop("DEGs_lists must contain at least 2 datasets for RRA analysis")
  }
  
  # Validate each element in the list
  for (name in names(DEGs_lists)) {
    df <- DEGs_lists[[name]]
    if (!is.data.frame(df)) {
      stop(paste("Element", name, "in DEGs_lists is not a data frame"))
    }
    if (!all(c("gene", "logFC", "P.Value") %in% colnames(df))) {
      stop(paste("Element", name, "must contain 'gene', 'logFC', and 'P.Value' columns"))
    }
  }
  
  if (!is.numeric(top.num) || top.num < 0) {
    stop("top.num must be a non-negative numeric value")
  }
  
  if (!is.numeric(rra.p) || rra.p <= 0 || rra.p >= 1) {
    stop("rra.p must be a numeric value between 0 and 1")
  }
  
  if (!is.numeric(logFC_cut) || logFC_cut < 0) {
    stop("logFC_cut must be a non-negative numeric value")
  }
  
  if (!is.numeric(p_cut) || p_cut <= 0 || p_cut >= 1) {
    stop("p_cut must be a numeric value between 0 and 1")
  }
  
  # Get lists of DEGs based on logFC and p-value cutoffs
  results <- get_DEGs_list(DEGs_lists, logFC_cut = logFC_cut, p_cut = p_cut)

  # Aggregate ranks for up- and down-regulated genes
  ups <- RobustRankAggreg::aggregateRanks(results$DEG_up)
  downs <- RobustRankAggreg::aggregateRanks(results$DEG_down)

  # Calculate frequency of each gene in the DEGs lists
  calculate_freq <- function(glist, ranked_data) {
    gene_freq <- as.data.frame(table(unlist(glist)))
    ranked_data$Freq <- gene_freq[match(ranked_data$Name, gene_freq[, 1]), 2]
    return(ranked_data)
  }

  ups <- calculate_freq(results$DEG_up, ups)
  downs <- calculate_freq(results$DEG_down, downs)

  # Filter genes based on their rank and frequency
  filter_genes <- function(data, num, gse_num) {
    selected <- data[data$Score < rra.p & data$Freq == gse_num, 1]
    if (num > 0 && length(selected) >= num) {
      selected <- selected[1:num]
    }
    return(as.character(selected))
  }

  gse.num <- length(DEGs_lists)
  gs_ups <- filter_genes(ups, top.num, gse.num)
  gs_downs <- filter_genes(downs, top.num, gse.num)

  # If the number of selected genes is less than the required top number, return early
  if (length(gs_ups) < top.num || length(gs_downs) < top.num) {
    return(list(min(length(gs_ups), length(gs_downs)), NULL))
  }

  # Combine up- and down-regulated genes
  gs <- unique(c(gs_ups, gs_downs))

  # Create a data matrix and set gene names as row names
  dat <- data.frame(row.names = gs)
  for (i in seq_along(DEGs_lists)) {
    logFC_values <- DEGs_lists[[i]][match(gs, DEGs_lists[[i]]$gene), 'logFC']
    names(logFC_values) <- DEGs_lists[[i]][match(gs, DEGs_lists[[i]]$gene), 'gene']
    dat <- cbind(dat, logFC_values)
  }
  colnames(dat) <- names(DEGs_lists)

  gs_list <- c("up" = length(gs_ups), "down" = length(gs_downs))

  # Return the number of genes and the RRA results matrix
  return(list("gene number" = gs_list, "RRA_results" = dat))
}
