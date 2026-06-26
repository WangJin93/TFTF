#' @title Meta-analysis in Multiple Datasets
#' @description
#' This function performs a meta-analysis on multiple datasets and generates a forest plot.
#' It also tests for publication bias.
#' @param results Data frame. The results data frame containing columns for dataset, n_Tumor, mean_Tumor, sd_Tumor, n_Normal, mean_Normal, and sd_Normal.
#' @param method Character. The statistical method to use (default is "wilcox").
#' @param k.min Integer. Minimum number of studies for bias test (default is 7).
#' @return A forest plot object.
#' @examples
#' \dontrun{
#' df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
#' results <- data_summary(df, tumor_subtype = "LUAD")
#' plot_meta_forest(results)
#' }
#' @export
plot_meta_forest <- function(results, method = "wilcox", k.min = 10) {
  if (!requireNamespace("meta", quietly = TRUE)) {
    stop("The 'meta' package is required for this function. Please install it with: install.packages('meta')")
  }
  # Remove rows with NA values
  results <- na.omit(results)

  # 确定研究数目
  num_studies <- nrow(results)

  # Perform meta-analysis using metacont function from the meta package
  metawsd <- meta::metacont(
    n.e = results$n_Tumor,
    mean.e = results$mean_Tumor,
    sd.e = results$sd_Tumor,
    n.c = results$n_Normal,
    mean.c = results$mean_Normal,
    sd.c = results$sd_Normal,
    data = results,
    sm = "SMD",  # Standard Mean Difference
    comb.fixed = FALSE,  # Do not combine fixed effects
    comb.random = TRUE,  # Combine random effects
    studlab = results$dataset
  )

  # Test for publication bias only if the number of studies is greater than or equal to k.min
  if (num_studies >= k.min) {
    bias_result <- meta::metabias(metawsd, k.min = k.min)

    # Print bias test results
    cat("\nMetabias Test Results:\n")
    print(bias_result)
  } else {
    cat("\nNot enough studies to perform metabias test. Minimum required is", k.min, "but only", num_studies, "studies were provided.\n")
  }

  # Generate and return forest plot
  meta::forest(metawsd, digits.se = 2, lab.e = "Tumor", lab.c = "Normal")
}

