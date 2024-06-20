#' @title Predict TFs targeting gene
#' @description
#'  Predict the upstream Transcription Factors regulating user inputted gene in multiple TF-target prediction databases and correlation analysis.
#' @import dplyr jsonlite tibble
#' @param datasets Including the TF-target regulatory data of 7 public TF online tools and the TF-targets prediction results analyzed by using FIMO and PWMEnrich.
#' @param target Input a gene symbol, e.g. GAPDH.
#' @param TCGA_tissue Cancer type in TCGA database, you can use tissue_type("TCGA") to abtain the tissue types.
#' @param GTEx_tissue Cancer type in GTEx database, you can use tissue_type("GTEx") to abtain the tissue types.
#' @param cor_DB The database used for the correlation analyze between TF and targets. You can use 2 databases, viz. TCGA (33 cancer types) and GTEx (31 normal tissue types).
#' @param cor_cutoff Threshold of correlation coefficient for correlation analysis.
#' @param FIMO.score Threshold of the score of the prediction TF-target results by using FIMO algorithm (bigger is better), default 10.
#' @param PWMEnrich.p Threshold of the score of the prediction TF-target results by using PWMEnrich algorithm (smaller is better), default 0.10.
#' @param cut.log2FC Threshold of log2FC for KnockTF dataset, default 1.
#' @param down.only Logic value. If true, only the downregulated genes in TF knockout/knockdown cells were returned in KnockTF dataset.
#' @param app Logic value. TRUE only used in the shiny app.
#' @examples
#' \dontrun{
#' results <- predict_TF(datasets=c("hTFtarget","KnockTF","ChIP_Atlas"), cor_DB = c("TCGA","GTEx"),target = "GAPDH")
#' }
#' @export
#'
predict_TF <- function(datasets=c("hTFtarget",
                                  "KnockTF",
                                  "FIMO_JASPAR",
                                  "PWMEnrich_JASPAR",
                                  "ENCODE",
                                  "CHEA",
                                  "TRRUST",
                                  "GTRD",
                                  "ChIP_Atlas"),
                       target = "GAPDH",
                       TCGA_tissue = "COAD",
                       GTEx_tissue = "Colon",
                       cor_DB = c("TCGA","GTEx"),
                       cor_cutoff = 0.3,
                       FIMO.score=10,
                       PWMEnrich.p =0.1,
                       cut.log2FC = 1,
                       down.only = T,
                       app = F){
  options(timeout=200)
  TF_result <- as.list(datasets)
  names(TF_result) <- datasets

  ###########HTFtarget
  if ("hTFtarget" %in% datasets){
    cat("Searching hTFtarget .... \n")
    if (isTRUE(app)){
      showNotification("Searching hTFtarget .... ",duration = 2)
      }
    # ensenbl <- idmap[idmap$gene == target ,]$ID
    # url<-paste0("https://guolab.wchscu.cn/hTFtarget/api/chipseq/targets/target?target=",ensenbl)
    # hTFtarget <- jsonlite::fromJSON(url)
    # if (is.null(hTFtarget)){
    #   hTFtarget <- data.frame()
    # }else{
    #   hTFtarget <-  hTFtarget[["regulations"]] %>%
    #     dplyr::rename("TF" = "tf_id")
    #   hTFtarget[["regulations"]][["regulation_info"]] <- hTFtarget[["regulations"]][["regulation_info"]][["regulation"]]
    # }
    hTFtarget <- get_data("hTFtarget","Target",target)

    TF_result[["hTFtarget"]] <- hTFtarget %>% na.omit()
  }
  ####KnockTF##数据
  if ("KnockTF" %in% datasets){
    cat("Searching KnockTF .... \n")
    if (isTRUE(app)){
      showNotification("Searching knocktf .... ",duration = 2)
    }
    KnockTF <- get_data("knocktf","Target",target)
    if (length(KnockTF) > 0){
      if (down.only) {
        KnockTF <- KnockTF %>%
          dplyr::filter(.,as.numeric( Log2FC) < -cut.log2FC ) %>%
          dplyr::filter(.,P_value < 0.05 )
      }else{
        KnockTF <- KnockTF %>%
          dplyr::filter(.,abs(as.numeric( Log2FC))> cut.log2FC )%>%
          dplyr::filter(.,P_value < 0.05 )
      }
      KnockTF <- merge(KnockTF,knocktf_data,by = "Sample_ID")
    }
    TF_result[["KnockTF"]] <- KnockTF %>% na.omit()
  }


  ####TRRUST###
  if ("TRRUST" %in% datasets){
    cat("Searching TRRUST .... \n")
    if (isTRUE(app)){
      showNotification("Searching TRRUST .... ",duration = 2)
    }
    TRRUST <- get_data("TRRUST","Target",target)
    TF_result[["TRRUST"]] <- TRRUST %>% na.omit()
  }
  # ###Cistrome
  # if ("Cistrome" %in% datasets){
  #   cat("Searching Cistrome .... This will take a long time... \n")
  #   # showNotification("Searching Cistrome, this will take a long time...  ",duration = 5,type = "message")
  #   GB_acc <- refgene[which(refgene$V13 == target),]$V2
  #   GB_acc <- GB_acc[which(stringr::str_detect( GB_acc,"NM_"))]
  #   Cistrome <- data.frame()
  #   for (x in GB_acc) {
  #     url <- paste0("http://dbtoolkit.Cistrome.org/?specie=hg38&keyword=",x,"&factor=factor&distance=10k")
  #     aaa <- url %>% read_html() %>% html_nodes("#resultPanel > div > div > div > div.panel-body > div:nth-child(3) > div > div") %>% html_nodes("iframe")
  #     aaa <- xml_attrs(aaa[[1]])[["src"]] %>% strsplit(.,"/") %>% .[[1]]
  #     url <-paste0("http://dbtoolkit.Cistrome.org/download/?fname=",aaa[5] %>% stringr::str_remove(.,".html"))
  #     ddd <- read.delim(url,sep = ",")
  #     ddd$GB_acc <- x
  #     Cistrome <- rbind(Cistrome,ddd)
  #   }
  #   Cistrome <- Cistrome[!duplicated(Cistrome),]
  #   TF_result[["Cistrome"]] <- Cistrome
  # }
  ####ENCODE
  if ("ENCODE" %in% datasets){
    cat("Searching ENCODE .... \n")
    showNotification("Searching ENCODE .... ",duration = 2)
    ENCODE <- get_data("ENCODE","Target",target)
    TF_result[["ENCODE"]] <- ENCODE %>% na.omit()
  }
  ###JASPAR
  if ("FIMO_JASPAR" %in% datasets){
    cat("Searching FIMO_JASPAR .... \n")
    if (isTRUE(app)){
      showNotification("Searching FIMO_JASPAR .... ",duration = 2)
    }
    Jaspar <- get_data("FIMO_JASPAR","Target",target)
    if (length(Jaspar) > 0){
      Jaspar$score <- as.numeric(Jaspar$score)
      Jaspar$p.value <- as.numeric(Jaspar$p.value)
      Jaspar$q.value <- as.numeric(Jaspar$q.value)
      Jaspar <- Jaspar[Jaspar$score > FIMO.score,]
    }
    TF_result[["FIMO_JASPAR"]] <- Jaspar %>% na.omit()
  }
  ###PWMEnrich_JASPAR
  if ("PWMEnrich_JASPAR" %in% datasets){
    cat("Searching PWMEnrich_JASPAR .... \n")
    if (isTRUE(app)){
      showNotification("Searching PWMEnrich_JASPAR .... ",duration = 2)
    }
    PWMEnrich_JASPAR <- get_data("PWMEnrich_JASPAR","Target",target)
    if (length(PWMEnrich_JASPAR) > 0){
      PWMEnrich_JASPAR$P.value <- as.numeric(PWMEnrich_JASPAR$P.value)
      PWMEnrich_JASPAR <- PWMEnrich_JASPAR[PWMEnrich_JASPAR$P.value < PWMEnrich.p,]
    }
    TF_result[["PWMEnrich_JASPAR"]] <- PWMEnrich_JASPAR %>% na.omit()
  }
  ###"CHEA"
  if ("CHEA" %in% datasets){
    cat("Searching CHEA .... \n")
    if (isTRUE(app)){
      showNotification("Searching CHEA .... ",duration = 2)
    }
    CHEA <- get_data("CHEA","Target",target)
    TF_result[["CHEA"]] <- CHEA %>% na.omit()
  }
  ###"GTRD"
  if ("GTRD" %in% datasets){
    cat("Searching GTRD .... \n")
    showNotification("Searching GTRD .... ",duration = 2)
    GTRD <- get_data("GTRD","Target",target)
    TF_result[["GTRD"]] <- GTRD %>% na.omit()
  }
  ####ChIP_Atlas##数据
  if ("ChIP_Atlas" %in% datasets){
    cat("Searching ChIP_Atlas .... \n")
    if (isTRUE(app)){
      showNotification("Searching ChIP_Atlas .... ",duration = 2)
    }
    ChIP_Atlas <- get_data("ChIP_Atlas","Target",target)
    TF_result[["ChIP_Atlas"]] <- ChIP_Atlas %>% na.omit()
  }

  ####cor_TCGA##数据
  if ("TCGA" %in% cor_DB){
    cat("Fetching correlation results of TCGA .... \n")
    if (isTRUE(app)){
      showNotification("Fetching correlation results of TCGA .... ",duration = 2)
    }
    cor_res <- get_data(paste0("cor_",TCGA_tissue),"Target",target)
    rownames(cor_res) <- NULL
    if (length(cor_res) > 0){
      cor_res <- cor_res %>% tibble::column_to_rownames("gene")%>%
      t()%>% as.data.frame()%>%
      na.omit() %>%
      rename(.,"cor"=all_of(target) ) %>%
      dplyr::mutate(.,cor = as.numeric( cor)/1000) %>%
      dplyr::filter(.,abs(cor) >= cor_cutoff) %>%
      tibble::rownames_to_column(.,"TF")
    }
    TF_result$cor_TCGA <- cor_res %>% na.omit()
  }
  ####cor_GTEx##数据
  if ("GTEx" %in% cor_DB){
    cat("Fetching correlation results of GTEx .... \n")
    if (isTRUE(app)){
      showNotification("Fetching correlation results of GTEx .... ",duration = 2)
    }
    cor_res <- get_data(paste0("cor_",GTEx_tissue),"Target",target)
    rownames(cor_res) <- NULL
    if (length(cor_res) > 0){
      cor_res <-  cor_res %>% tibble::column_to_rownames("gene")%>%
      t()%>% as.data.frame()%>%
      na.omit() %>%
      rename(.,"cor"=all_of(target) ) %>%
      dplyr::mutate(.,cor = as.numeric( cor)/1000) %>%
      dplyr::filter(.,abs(cor) >= cor_cutoff) %>%
      tibble::rownames_to_column(.,"TF")
    }
    TF_result$cor_GTEx <- cor_res %>% na.omit()
  }


  ##intersect
  results <- list()
  for (i in names(TF_result)) {
    if (length(TF_result[[i]])>0){
      results[[i]] <- TF_result[[i]] %>% as.data.frame() %>% .[,"TF"] %>% unique()
    }
  }
  TF_result[["results"]] <- results

  return(TF_result)
}
