#' @title Correlation analysis
#' @description
#'  Correlation analysis between TF and target gene in pan-tissue in "TCGA", "GTEx" or "CCLE" databases.
#' @import UCSCXenaShiny dplyr psych tibble
#' @param Gene1 Gene1 name.
#' @param Gene2 Gene2 name.
#' @param data_source Database used for analysis, including "TCGA", "GTEx" or "CCLE".
#' @param type Sample types used for analysis, including "normal","tumor".
#' @param cor_method Method for correlation analysis, "pearson" or "spearman".
#'
#' @export
#'
pantissue_cor_analysis <- function(Gene1 = "FOXM1",
                         Gene2 = "GAPDH",
                         data_source = "TCGA",
                         type = c("normal","tumor"),
                         cor_method = "pearson") {
  if (data_source != "CCLE"){
    tcga_gtex <- tcga_gtex %>%
      dplyr::group_by(.data$tissue) %>%
      dplyr::distinct(.data$sample, .keep_all = TRUE)

    t1 <- query_pancan_value(Gene1, data_type = "mRNA")
    if (is.list(t1)) t1 <- t1[[1]]
    t3 <- query_pancan_value(Gene2, data_type = "mRNA")
    if (is.list(t3)) t3 <- t3[[1]]
    if (all(is.na(t1))){
      # showModal(modalDialog(
      #   title = "Message",   easyClose = TRUE,
      cat(paste0("No results were returned, your inputed TF, ",Gene1,", does not have valid data in the TCGA database." ))
      # ))
      return(NULL)
    }
    if (all(is.na(t3))){
      # showModal(modalDialog(
      #   title = "Message",   easyClose = TRUE,
      cat(paste0("No results were returned, your inputed target, ",Gene2,", does not have valid data in the TCGA database." ))
      # ))
      return(NULL)
    }
    t2 <- t1 %>%
      as.data.frame() %>%
      dplyr::rename("tpm" = ".") %>%
      tibble::rownames_to_column(var = "sample") %>%
      dplyr::inner_join(tcga_gtex, by = "sample")

    t4 <- t3 %>%
      as.data.frame() %>%
      dplyr::rename("tpm" = ".") %>%
      tibble::rownames_to_column(var = "sample") %>%
      dplyr::inner_join(tcga_gtex, by = "sample")

    # merge
    t2 <- t2 %>% inner_join(t4[, c("sample", "tpm")], by = "sample")

    df <- data.frame(
      sample = t2$sample,
      tissue = t2$tissue,
      type = t2$type,
      type1 = t2$type1,

      type2 = t2$type2,
      gene1 = t2$tpm.x,
      gene2 = t2$tpm.y,
      stringsAsFactors = F
    )

  }else{
    t1 <- query_pancan_value(Gene1, data_type = data_type1, database = "ccle")
    if (is.list(t1)) t1 <- t1[[1]]
    t1 <- log2(t1 + 1)
    if (all(is.na(t1))) {
      message("All NAs returned, return NULL instead.")
      return(NULL)
    }

    t2 <- t1 %>%
      as.data.frame() %>%
      dplyr::rename("tpm" = ".") %>%
      tibble::rownames_to_column(var = "cell") %>%
      dplyr::inner_join(ccle_info, by = c("cell" = "CCLE_name"))

    t3 <- query_pancan_value(Gene2, data_type = data_type2, database = "ccle")
    if (is.list(t3)) t3 <- t3[[1]]
    t3 <- log2(t3 + 1)
    if (all(is.na(t3))) {
      message("All NAs returned, return NULL instead.")
      return(NULL)
    }

    t4 <- t3 %>%
      as.data.frame() %>%
      dplyr::rename("tpm" = ".") %>%
      tibble::rownames_to_column(var = "cell") %>%
      dplyr::inner_join(ccle_info, by = c("cell" = "CCLE_name"))

    t2 <- t2 %>% inner_join(t4[, c("cell", "tpm")], by = "cell")

    df <- data.frame(sample = t2$cell,
                     gene1 = t2$tpm.x,
                     gene2 = t2$tpm.y,
                     tissue = t2$Site_Primary,
                     stringsAsFactors = F)

  }
  cat("Performing correlation analysis...")

  if (data_source == "TCGA"){
    df %>%
      dplyr::filter(.data$type1 == data_source) %>%
      dplyr::filter(.data$type2 %in% type) -> df
  }
  if (data_source == "GTEx"){
    df %>%
      dplyr::filter(.data$type1 == data_source) -> df
  }


  output<-data.frame("",0,0,0)
  names(output)<-c("tissue","n","r","p")
  output<-output[-1,]
  for(type in unique(df$tissue)){
    data<- df[which(df$tissue == type),]
    if (nrow(data)==0){
      next
    }
    cor_result<- psych::corr.test(x =data$gene1, y = data$gene2,method = cor_method)
    nvalue = as.numeric(cor_result$n)
    rvalue = signif(as.numeric(cor_result$r),4)
    pvalue = signif(as.numeric(cor_result$p),4)
    output1<-data.frame(tcga=type,nvalue,rvalue,pvalue)
    names(output1)<-c("tissue","n","r","p")
    output<-rbind(output,output1)
  }
  output$logP<-log10(output$p)*(-1)
  output$sigmark <- ifelse(output$p>0.05,"No",ifelse(output$r> input$xinter,"Positive",ifelse(output$r< -input$xinter,"Negative","No")))

  return(list(cor_data = df,
              cor_result = output))
}

