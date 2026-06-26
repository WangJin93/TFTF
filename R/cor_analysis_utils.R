#' @title Internal utility function for correlation analysis
#' @description This is an internal function used by other correlation analysis functions.
#' It calculates correlations between a target variable and multiple variables across datasets.
#' @import dplyr psych stats
#' @param df A dataframe containing the data to be analyzed.
#' @param target_col The column name of the target variable for correlation.
#' @param sig_cols A vector of column names to correlate with the target variable.
#' @param group_col The column name to group the data by (default: "dataset").
#' @param cor_method The correlation method to use ("pearson", "spearman", or "kendall").
#' @param min_samples The minimum number of samples required for correlation calculation (default: 4).
#' @param adjust_method The method for adjusting p-values for multiple testing. Options include "none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni". Default is "BH".
#' @param conf_level The confidence level for confidence interval calculation (default: 0.95).
#' @return A list containing correlation matrix (r), p-value matrix (p), adjusted p-value matrix (p_adj), 
#'         sample size matrix (n), t-statistic matrix (t), confidence interval lower bound matrix (ci_lower), 
#'         confidence interval upper bound matrix (ci_upper), and split data (sss).
#' @keywords internal
.calculate_correlation <- function(df, target_col, sig_cols, 
                                   group_col = "dataset", cor_method = "pearson",
                                   min_samples = 4, adjust_method = "BH",
                                   conf_level = 0.95) {
  
  # Input validation
  if (is.null(df)) {
    warning("Input dataframe is NULL")
    return(NULL)
  }
  
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  
  if (!target_col %in% colnames(df)) {
    stop(paste("target_col", target_col, "not found in dataframe"))
  }
  
  missing_cols <- setdiff(sig_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from dataframe:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  if (!group_col %in% colnames(df)) {
    stop(paste("group_col", group_col, "not found in dataframe"))
  }
  
  if (!cor_method %in% c("pearson", "spearman", "kendall")) {
    stop("cor_method must be 'pearson', 'spearman', or 'kendall'")
  }
  
  if (!adjust_method %in% c("none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni")) {
    stop("adjust_method must be one of: 'none', 'BH', 'BY', 'holm', 'hochberg', 'hommel', 'bonferroni'")
  }
  
  # Split the dataframe by the grouping column
  sss <- split(df, df[[group_col]])
  groups <- names(sss)
  nrow <- length(groups)
  ncol <- length(sig_cols)
  
  # Initialize matrices to store all statistics
  rvalue <- matrix(nrow = nrow, ncol = ncol)
  rownames(rvalue) <- groups
  colnames(rvalue) <- sig_cols
  
  pvalue <- matrix(nrow = nrow, ncol = ncol)
  rownames(pvalue) <- groups
  colnames(pvalue) <- sig_cols
  
  pvalue_adj <- matrix(nrow = nrow, ncol = ncol)
  rownames(pvalue_adj) <- groups
  colnames(pvalue_adj) <- sig_cols
  
  nvalue <- matrix(nrow = nrow, ncol = ncol)
  rownames(nvalue) <- groups
  colnames(nvalue) <- sig_cols
  
  tvalue <- matrix(nrow = nrow, ncol = ncol)
  rownames(tvalue) <- groups
  colnames(tvalue) <- sig_cols
  
  ci_lower <- matrix(nrow = nrow, ncol = ncol)
  rownames(ci_lower) <- groups
  colnames(ci_lower) <- sig_cols
  
  ci_upper <- matrix(nrow = nrow, ncol = ncol)
  rownames(ci_upper) <- groups
  colnames(ci_upper) <- sig_cols
  
  # Calculate correlation for each group
  for (i in seq_along(groups)) {
    group_data <- sss[[i]]
    if (nrow(group_data) < min_samples) {
      next
    } else {
      corr_result <- psych::corr.test(x = group_data[, target_col, drop = TRUE],
                                      y = group_data[, sig_cols, drop = FALSE],
                                      method = cor_method)
      r <- corr_result[["r"]]
      p <- corr_result[["p"]]
      rvalue[i, ] <- r
      pvalue[i, ] <- p
      
      # Apply multiple testing correction within each group
      if (adjust_method != "none") {
        pvalue_adj[i, ] <- stats::p.adjust(p, method = adjust_method)
      } else {
        pvalue_adj[i, ] <- p
      }
      
      # Calculate sample size
      n <- nrow(group_data)
      nvalue[i, ] <- n
      
      # Calculate t-statistic (only for pearson correlation)
      if (cor_method == "pearson") {
        t <- r * sqrt((n - 2) / (1 - r^2))
        tvalue[i, ] <- t
      } else {
        tvalue[i, ] <- NA
      }
      
      # Calculate confidence interval
      for (j in seq_along(r)) {
        ci <- .calculate_correlation_ci(r[j], n, conf_level)
        ci_lower[i, j] <- ci["lower"]
        ci_upper[i, j] <- ci["upper"]
      }
    }
  }
  
  # Check for datasets with insufficient samples
  insufficient_datasets <- groups[apply(rvalue, 1, function(row) all(is.na(row)))]
  if (length(insufficient_datasets) > 0) {
    warning(paste("The following datasets have insufficient samples (<", min_samples, 
                  ") and will be excluded from correlation analysis:", 
                  paste(insufficient_datasets, collapse = ", ")))
  }
  
  # Transpose all matrices, and omit any NA values
  rvalue_T <- t(na.omit(rvalue))
  pvalue_T <- t(na.omit(pvalue))
  pvalue_adj_T <- t(na.omit(pvalue_adj))
  nvalue_T <- t(na.omit(nvalue))
  tvalue_T <- t(na.omit(tvalue))
  ci_lower_T <- t(na.omit(ci_lower))
  ci_upper_T <- t(na.omit(ci_upper))
  
  # Create a list to hold the results
  plist <- list(r = rvalue_T, p = pvalue_T, p_adj = pvalue_adj_T,
                n = nvalue_T, t = tvalue_T, ci_lower = ci_lower_T, 
                ci_upper = ci_upper_T, sss = sss)
  
  return(plist)
}

#' @title Internal utility function to extract tumor subtypes
#' @description This is an internal helper function to extract and validate tumor subtypes.
#' @param subtype_data The subtype mapping data.
#' @param tumor_subtype The tumor subtypes to extract.
#' @return A vector of extracted subtypes.
#' @keywords internal
.extract_tumor_subtype <- function(subtype_data, tumor_subtype) {
  if (is.null(tumor_subtype)) {
    return(NULL)
  }
  
  if (!is.character(tumor_subtype)) {
    stop("tumor_subtype must be a character vector")
  }
  
  extracted <- extract_subset(subtype_data, tumor_subtype)
  
  if (length(extracted) == 0) {
    warning("No matching subtypes found for the given tumor_subtype")
    return(NULL)
  }
  
  return(extracted)
}

#' @title Internal utility function to determine sample type
#' @description This is an internal helper function to determine sample types (Tumor/Normal).
#' @param df The dataframe containing sample data.
#' @param tumor_subtype The tumor subtypes to consider.
#' @param subtype_col The column name for subtype information.
#' @return The dataframe with an added "type" column.
#' @keywords internal
.determine_sample_type <- function(df, tumor_subtype = NULL, subtype_col = "subtype") {
  
  if (!subtype_col %in% colnames(df)) {
    stop(paste("subtype_col", subtype_col, "not found in dataframe"))
  }
  
  if (is.null(tumor_subtype)) {
    df <- df %>%
      dplyr::mutate(type = ifelse(!!sym(subtype_col) %in% c("Normal", "Adjacent"), 
                                  "Normal", "Tumor"))
  } else {
    df <- df %>%
      dplyr::filter(!!sym(subtype_col) %in% c(tumor_subtype, "Normal", "Adjacent")) %>%
      dplyr::mutate(type = ifelse(!!sym(subtype_col) %in% tumor_subtype, "Tumor", "Normal"))
  }
  
  return(df)
}

#' @title Calculate Confidence Interval for Correlation Coefficient
#' @description Calculate confidence interval for a correlation coefficient using Fisher's z-transformation.
#' @param r The correlation coefficient.
#' @param n The sample size.
#' @param conf_level The confidence level (default: 0.95).
#' @return A vector containing the lower and upper bounds of the confidence interval.
#' @keywords internal
.calculate_correlation_ci <- function(r, n, conf_level = 0.95) {
  # Fisher's z-transformation
  z <- 0.5 * log((1 + r) / (1 - r))
  se <- 1 / sqrt(n - 3)
  
  # Calculate confidence interval in z-space
  z_alpha <- qnorm((1 + conf_level) / 2)
  z_lower <- z - z_alpha * se
  z_upper <- z + z_alpha * se
  
  # Transform back to r-space
  r_lower <- (exp(2 * z_lower) - 1) / (exp(2 * z_lower) + 1)
  r_upper <- (exp(2 * z_upper) - 1) / (exp(2 * z_upper) + 1)
  
  return(c(lower = r_lower, upper = r_upper))
}

#' @title Calculate Effect Size for Differential Expression
#' @description Calculate various effect size metrics for differential expression analysis.
#' @param logFC The log fold change.
#' @param se The standard error of the logFC.
#' @param mean_tumor The mean expression in tumor samples.
#' @param mean_normal The mean expression in normal samples.
#' @param sd_tumor The standard deviation in tumor samples.
#' @param sd_normal The standard deviation in normal samples.
#' @return A list containing effect size metrics (Cohen's d, Hedges' g, and Glass's delta).
#' @keywords internal
.calculate_effect_size <- function(logFC, se = NULL, mean_tumor = NULL, mean_normal = NULL, 
                                   sd_tumor = NULL, sd_normal = NULL) {
  results <- list()
  
  # Cohen's d approximation from logFC and SE
  if (!is.null(logFC) && !is.null(se)) {
    results$cohens_d_approx <- abs(logFC) / se
  }
  
  # Cohen's d using means and pooled SD
  if (!is.null(mean_tumor) && !is.null(mean_normal) && !is.null(sd_tumor) && !is.null(sd_normal)) {
    pooled_sd <- sqrt((sd_tumor^2 + sd_normal^2) / 2)
    results$cohens_d <- abs(mean_tumor - mean_normal) / pooled_sd
    
    # Hedges' g (corrected for small sample bias)
    results$hedges_g <- results$cohens_d * (1 - 3 / (4 * (length(logFC) * 2) - 9))
    
    # Glass's delta (using control group SD)
    results$glass_delta <- abs(mean_tumor - mean_normal) / sd_normal
  }
  
  # Interpret effect size
  if (!is.null(results$cohens_d)) {
    results$interpretation <- .interpret_effect_size(results$cohens_d)
  } else if (!is.null(results$cohens_d_approx)) {
    results$interpretation <- .interpret_effect_size(results$cohens_d_approx)
  }
  
  return(results)
}

#' @title Interpret Effect Size
#' @description Provide interpretation of effect size based on Cohen's guidelines.
#' @param effect_size The effect size (Cohen's d).
#' @return A character string describing the effect size magnitude.
#' @keywords internal
.interpret_effect_size <- function(effect_size) {
  effect_size <- abs(effect_size)
  if (effect_size < 0.2) {
    return("negligible (d < 0.2)")
  } else if (effect_size < 0.5) {
    return("small (0.2 <= d < 0.5)")
  } else if (effect_size < 0.8) {
    return("medium (0.5 <= d < 0.8)")
  } else {
    return("large (d >= 0.8)")
  }
}