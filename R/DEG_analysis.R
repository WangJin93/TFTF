#' @title Differential Expression Gene (DEG) Analysis
#' @description Perform differential expression gene analysis on a given dataset with comprehensive statistical metrics.
#' @import dplyr
#' @param df A dataframe containing gene expression data with sample IDs as columns.
#' @param tumor_subtype A character vector specifying the tumor subtypes to be analyzed. Default is NULL, which means all tumor subtypes will be included.
#' @param adjust_method The method for adjusting p-values for multiple testing. Options include "none", "BH" (Benjamini-Hochberg), "BY" (Benjamini-Yekutieli), "holm", "hochberg", "hommel", "bonferroni". Default is "BH".
#' @param confint Logical, whether to calculate confidence intervals for log fold changes. Default is TRUE.
#' @param ... Additional arguments passed to `lmFit`, `contrasts.fit`, and `eBayes`.
#' @return A dataframe with DEG analysis results, including log fold changes, p-values, adjusted p-values, confidence intervals, and effect size metrics.
#' @examples
#' \dontrun{
#' df <- get_OSF_data(table = "GSE74706", action = "geo_data")
#' results <- DEGs_analysis(df, tumor_subtype = c("NSCLC"))
#' }
#' @export
DEGs_analysis <- function(df, tumor_subtype = NULL, adjust_method = "BH", confint = TRUE, ...) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("The 'limma' package is required for this function. Please install it with: BiocManager::install('limma')")
  }
  # 检查 sample_subtype 是否存在且是数据框
  if (!exists("sample_subtype")) {
    stop("sample_subtype data not found. Please load it with: data('sample_subtype')")
  }
  
  if (!is.data.frame(sample_subtype)) {
    stop("sample_subtype must be a data frame")
  }
  
  if (!"subtype" %in% colnames(sample_subtype)) {
    stop("sample_subtype must contain a 'subtype' column")
  }
  
  # Input validation
  if (is.null(df)) {
    warning("Input dataframe is NULL")
    return(NULL)
  }
  
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  
  if (!"ID" %in% colnames(df)) {
    stop("df must contain an 'ID' column")
  }
  
  if (!is.null(tumor_subtype) && !is.character(tumor_subtype)) {
    stop("tumor_subtype must be a character vector")
  }
  
  if (!adjust_method %in% c("none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni")) {
    stop("adjust_method must be one of: 'none', 'BH', 'BY', 'holm', 'hochberg', 'hommel', 'bonferroni'")
  }

  # Filter sample information based on the provided sample IDs
  sample_info <- sample_subtype[sample_subtype$ID %in% colnames(df), , drop = FALSE]
  
  # 确保 sample_info 是数据框
  if (!is.data.frame(sample_info)) {
    warning("sample_info is not a data frame")
    return(NULL)
  }
  
  # 如果没有匹配的样本，返回 NULL
  if (nrow(sample_info) == 0) {
    warning("No samples found in sample_subtype matching the provided dataframe")
    return(NULL)
  }

  # Assign type based on tumor_subtype
  if (is.null(tumor_subtype)) {
    sample_info$type2 <- ifelse(sample_info$subtype %in% c("Normal", "Adjacent"), "Normal", "Tumor")
  } else {
    tumor_subtype <- extract_subset(subtype, tumor_subtype)
    
    # 确保 tumor_subtype 是字符向量
    if (!is.character(tumor_subtype)) {
      tumor_subtype <- as.character(tumor_subtype)
    }
    
    sample_info <- dplyr::filter(sample_info, subtype %in% c(tumor_subtype, "Normal", "Adjacent"))
    sample_info$type2 <- ifelse(sample_info$subtype %in% tumor_subtype, "Tumor", "Normal")
  }

  # Subset df to include only relevant sample columns
  df <- df[, c("ID", sample_info$ID)]

  # Define group and design matrix for linear modeling
  group <- sample_info$type2
  design <- model.matrix(~0 + factor(group))
  colnames(design) <- levels(factor(group))

  # Create contrast matrix for Tumor vs Normal comparison
  contrast.matrix <- limma::makeContrasts(contrasts = "Tumor - Normal", levels = design)

  # Fit linear model
  fit <- limma::lmFit(df, design, ...)

  # Apply contrasts and compute eBayes statistics
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2, ...)

  # Extract top differentially expressed genes with confidence intervals
  if (confint) {
    DEG <- limma::topTable(fit2, coef = 1, n = Inf, adjust = adjust_method, confint = 0.95, ...) %>% na.omit()
  } else {
    DEG <- limma::topTable(fit2, coef = 1, n = Inf, adjust = adjust_method, ...) %>% na.omit()
  }

  # Add effect size metrics (Cohen's d-like measure based on logFC and standard error)
  if ("se" %in% colnames(DEG)) {
    DEG <- DEG %>%
      dplyr::mutate(
        effect_size = abs(logFC) / se,  # Cohen's d approximation
        effect_size_category = cut(effect_size, 
                                  breaks = c(-Inf, 0.2, 0.5, 0.8, Inf),
                                  labels = c("negligible", "small", "medium", "large"))
      )
  }

  return(DEG)
}