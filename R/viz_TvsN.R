#' @title Visualizing genes expression
#' @description
#' Visualizing the different expression of mRNA expression data between Tumor and Normal tissues in GEO database.
#' @import dplyr ggplot2 reshape2 ggpubr
#' @param df Gene expression data obtained from get_expr_data().
#' @param df_type The type of gene expression data, one value of "single","multi_gene", and "multi_set".
#' @param tumor_subtype Tumor subtype for filtering data.
#' @param Show.P.value Whether to display the results of differential analysis, default TRUE.
#' @param Show.P.label Whether to display significance markers for differential analysis, default TRUE.
#' @param Method Methods of differential analysis, "t.test" or "limma", default "t.test".
#' @param values Color palette for normal and tumor groups. Default gcas_palette(2).
#' @param Show.n Display sample size.
#' @param Show.n.location Y-axis position displayed for sample size.
#' @examples
#' \dontrun{
#' df_single <- get_expr_data(datasets = "GSE27262",genes = c("TP53"))
#' df_multi_gene <- get_expr_data(datasets = "GSE27262",genes = c("TP53","TNS1"))
#' df_multi_set <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
#' viz_TvsN(df_single,df_type = "single")
#' viz_TvsN(df_multi_gene,df_type = "multi_gene",tumor_subtype ="LC")
#' viz_TvsN(df_multi_set,df_type = "multi_set")
#' }
#' @export
viz_TvsN <- function(df, df_type = c("single", "multi_gene", "multi_set"),
                     tumor_subtype = NULL,
                     Show.P.value = TRUE,
                     Show.P.label = TRUE,
                     Method = "t.test",
                     values = gcas_palette(2),
                     Show.n = TRUE,
                     Show.n.location = "default") {
  
  # Input validation
  if (is.null(df)) {
    warning("Input dataframe is NULL")
    return(NULL)
  }
  
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  
  df_type <- match.arg(df_type)
  
  # Determine sample type using utility function
  df <- .determine_sample_type(df, tumor_subtype)
  
  # Filter to only Normal and Tumor types
  df <- df %>% dplyr::filter(type %in% c("Normal", "Tumor"))
  
  # Define significance cutpoints and symbols
  sig_cutpoints <- c(0, 0.001, 0.01, 0.05, 1)
  sig_symbols <- c("***", "**", "*", "ns")
  
  # Helper function to calculate y position for annotations
  .get_y_position <- function(data, value_col, offset = 0.1) {
    values <- data[[value_col]]
    y_max <- max(values, na.rm = TRUE)
    y_min <- min(values, na.rm = TRUE)
    list(
      top = y_max + (y_max - y_min) * offset,
      bottom = y_min - (y_max - y_min) * offset * 2
    )
  }
  
  # Single gene visualization
  if (df_type == "single") {
    # Get gene column (third column after ID and dataset)
    gene_col <- colnames(df)[3]
    df <- df %>% 
      dplyr::rename(value = all_of(gene_col)) %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::mutate(value = as.numeric(value))
    
    # Calculate p-values if requested
    pv <- NULL
    if (Show.P.value) {
      pv <- df %>%
        ggpubr::compare_means(value ~ type, data = ., method = Method, 
                              symnum.args = list(cutpoints = sig_cutpoints, 
                                                 symbols = sig_symbols),
                              p.adjust.methods = "BH") %>%
        dplyr::select(c("p", "p.signif", "p.adj"))
    }
    
    # Calculate sample counts
    count_N <- df %>% 
      dplyr::group_by(type) %>% 
      dplyr::tally() %>%
      dplyr::mutate(n = paste("n =", n))
    
    # Create plot
    y_pos <- .get_y_position(df, "value")
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = type, y = value, fill = type)) +
      ggplot2::geom_violin(trim = FALSE, show.legend = FALSE) +
      ggplot2::geom_boxplot(width = 0.2, fill = "white", show.legend = FALSE) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Expression") +
      ggplot2::scale_fill_manual(values = values) +
      gcas_theme()
    
    # Add sample counts if requested
    if (Show.n) {
      n_location <- if (Show.n.location == "default") y_pos$bottom else Show.n.location
      p <- p +
        ggplot2::geom_text(data = count_N, 
                          ggplot2::aes(label = n, x = type, y = n_location, color = type),
                          position = ggplot2::position_dodge2(0.9),
                          size = 6, hjust = 0.5, show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = values)
    }
    
    # Add p-value annotations if requested
    if (Show.P.value && Show.P.label && !is.null(pv)) {
      p <- p + 
        ggplot2::geom_text(ggplot2::aes(x = 1.5, y = y_pos$top, label = .data$p.signif),
                          data = pv, size = 8, inherit.aes = FALSE)
    }
    
    if (Show.P.value && !Show.P.label && !is.null(pv)) {
      p <- p + 
        ggplot2::geom_text(ggplot2::aes(
          x = 1.5, y = y_pos$top,
          label = ifelse(signif(.data$p, 3) < 0.001, "P < 0.001", 
                        paste0("P = ", signif(.data$p, 3)))
        ), data = pv, size = 6, inherit.aes = FALSE)
    }
  }
  
  # Multi-gene visualization
  if (df_type == "multi_gene") {
    # Get gene columns (exclude non-expression columns)
    non_expr_cols <- c("ID", "subtype", "dataset", "tissue", "Patient.ID", "type")
    gene_cols <- setdiff(colnames(df), non_expr_cols)
    
    # Melt data for plotting
    df <- df %>%
      reshape2::melt(measure.vars = gene_cols) %>%
      dplyr::mutate(value = as.numeric(value)) %>%
      dplyr::filter(!is.na(value))
    
    # Calculate p-values if requested
    pv <- NULL
    if (Show.P.value) {
      pv <- df %>%
        ggpubr::compare_means(value ~ type, data = ., method = Method, 
                              group.by = "variable",
                              symnum.args = list(cutpoints = sig_cutpoints, 
                                                 symbols = sig_symbols),
                              p.adjust.methods = "BH") %>%
        dplyr::select(c("variable", "p", "p.signif", "p.adj"))
      message("Counting P value finished")
    }
    
    # Calculate sample counts
    count_N <- df %>% 
      dplyr::group_by(variable, type) %>% 
      dplyr::tally() %>%
      dplyr::mutate(n = paste("n =", n))
    
    # Create plot
    y_pos <- .get_y_position(df, "value")
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = variable, y = value, fill = type)) +
      ggplot2::geom_boxplot() +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Expression") +
      ggplot2::scale_fill_manual(values = values) +
      gcas_theme()
    
    # Add sample counts if requested
    if (Show.n) {
      n_location <- if (Show.n.location == "default") y_pos$bottom else Show.n.location
      p <- p +
        ggplot2::geom_text(data = count_N, 
                          ggplot2::aes(label = n, y = n_location, color = type),
                          position = ggplot2::position_dodge2(0.9),
                          size = 6, angle = 90, hjust = 0, show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = values)
    }
    
    # Add p-value annotations if requested
    if (Show.P.value && Show.P.label && !is.null(pv)) {
      p <- p + 
        ggplot2::geom_text(ggplot2::aes(x = .data$variable, y = y_pos$top, 
                                        label = .data$p.signif),
                          data = pv, size = 6, inherit.aes = FALSE)
    }
    
    if (Show.P.value && !Show.P.label && !is.null(pv)) {
      p <- p + 
        ggplot2::geom_text(ggplot2::aes(x = .data$variable, y = y_pos$top,
                                        label = as.character(signif(.data$p, 2))),
                          data = pv, inherit.aes = FALSE)
    }
  }
  
  # Multi-dataset visualization
  if (df_type == "multi_set") {
    # Get gene column (third column)
    if (ncol(df) < 3) {
      warning("Dataset does not contain expression columns")
      return(NULL)
    }
    
    gene_col <- colnames(df)[3]
    df <- df %>% 
      dplyr::rename(value = all_of(gene_col)) %>%
      dplyr::mutate(value = as.numeric(value)) %>%
      dplyr::filter(!is.na(value))
    
    # Filter out datasets that don't have both Normal and Tumor types
    dataset_types <- df %>%
      dplyr::group_by(dataset) %>%
      dplyr::summarize(n_types = dplyr::n_distinct(type))
    
    valid_datasets <- dataset_types %>%
      dplyr::filter(n_types >= 2) %>%
      dplyr::pull(dataset)
    
    if (length(valid_datasets) == 0) {
      warning("No datasets contain both Normal and Tumor types for comparison")
      return(NULL)
    }
    
    df <- df %>%
      dplyr::filter(dataset %in% valid_datasets)
    
    # Calculate p-values if requested
    pv <- NULL
    if (Show.P.value) {
      pv <- df %>%
        ggpubr::compare_means(value ~ type, data = ., method = Method, 
                              group.by = "dataset",
                              symnum.args = list(cutpoints = sig_cutpoints, 
                                                 symbols = sig_symbols),
                              p.adjust.methods = "BH") %>%
        dplyr::select(c("dataset", "p", "p.signif", "p.adj"))
    }
    
    # Calculate sample counts
    count_N <- df %>% 
      dplyr::group_by(dataset, type) %>% 
      dplyr::tally() %>%
      dplyr::mutate(n = paste("n =", n))
    
    # Create plot
    y_pos <- .get_y_position(df, "value")
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = dataset, y = value, fill = type)) +
      ggplot2::geom_boxplot() +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Expression") +
      ggplot2::scale_fill_manual(values = values) +
      gcas_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1.0),
                     plot.margin = ggplot2::unit(c(0.2, 0.2, 0.2, 2), 'cm'))
    
    # Add sample counts if requested
    if (Show.n) {
      n_location <- if (Show.n.location == "default") y_pos$bottom else Show.n.location
      p <- p +
        ggplot2::geom_text(data = count_N, 
                          ggplot2::aes(label = n, y = n_location, color = type),
                          position = ggplot2::position_dodge2(0.9),
                          size = 5, angle = 90, hjust = 0, show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = values)
    }
    
    # Add p-value annotations if requested
    if (Show.P.value && Show.P.label && !is.null(pv)) {
      p <- p + 
        ggplot2::geom_text(ggplot2::aes(x = .data$dataset, y = y_pos$top, 
                                        label = .data$p.signif),
                          data = pv, inherit.aes = FALSE)
    }
    
    if (Show.P.value && !Show.P.label && !is.null(pv)) {
      p <- p + 
        ggplot2::geom_text(ggplot2::aes(x = .data$dataset, y = y_pos$top,
                                        label = as.character(signif(.data$p, 2))),
                          data = pv, inherit.aes = FALSE)
    }
  }
  
  return(p)
}