#' predict_target(tf = "STAT3")
#' @import RMySQL
#' @import dplyr
#' @import jsonlite
#' @import tibble
#' @export
#'
predict_target <- function(datasets=c("hTFtarget",
                                      "KnockTF",
                                      "FIMO_JASPAR",
                                      "PWMEnrich_JASPAR",
                                      "ENCODE",
                                      "CHEA",
                                      "TRRUST",
                                      "GTRD",
                                      "ChIP_Atlas"),
                           cor_DB = c("TCGA","GTEx"),
                           TCGA_tissue = "COAD",
                           GTEx_tissue = "Colon",
                           cor_cutoff = 0.3,
                           tf = "STAT3",
                           FIMO.score=10,
                           PWMEnrich_JASPAR.score =10,
                           cut.log2FC = 1,
                           down.only = T){
  options(timeout=200)
  targets <- as.list(datasets)
  names(targets) <- datasets

  ###########HTFtarget
  if ("hTFtarget" %in% datasets){
    cat("Searching hTFtarget .... \n")
    # showNotification("Searching hTFtarget .... ",duration = 3,type = "message")
    url<-paste0("https://guolab.wchscu.cn/hTFtarget/api/chipseq/targets/tf?page=1&size=20000&tf=",tf)
    hTFtarget <- jsonlite::fromJSON(url)
    if(is.null(hTFtarget[["targets"]][["ensemblgene"]] )){
      hTFtarget <-  character()
    }else{
      hTFtarget <-  hTFtarget[["targets"]][["ensemblgene"]] %>%
        dplyr::rename("Target" = "name")
    }

    targets[["hTFtarget"]] <- hTFtarget %>% na.omit()
  }

  ####TRRUST###
  if ("TRRUST" %in% datasets){
    cat("Searching TRRUST .... \n")
    # showNotification("Searching TRRUST .... ",duration = 3,type = "message")
    TRRUST <- get_data("TRRUST","TF",tf)

    targets[["TRRUST"]] <- TRRUST %>% na.omit()
  }
  ####ENCODE
  if ("ENCODE" %in% datasets){
    cat("Searching ENCODE .... \n")
    # showNotification("Searching ENCODE .... ",duration = 3,type = "message")
    ENCODE <- get_data("ENCODE","TF",tf)
    targets[["ENCODE"]] <- ENCODE %>% na.omit()
  }
  ###JASPAR
  if ("FIMO_JASPAR" %in% datasets){
    cat("Searching FIMO_JASPAR .... \n")
  # showNotification("Searching FIMO_JASPAR .... ",duration = 3,type = "message")
    Jaspar <- get_data("FIMO_JASPAR","TF",tf)
    if (length(Jaspar) > 0){
      Jaspar <- Jaspar[Jaspar$score > FIMO.score,]
      Jaspar <- merge(Jaspar,refgene,by="Target")
      colnames(Jaspar)[c(1,4)] <- c("Target_refseq","Target")
    }
    targets[["FIMO_JASPAR"]] <- Jaspar %>% na.omit()
  }
  ###PWMEnrich_JASPAR
  if ("PWMEnrich_JASPAR" %in% datasets){
    cat("Searching PWMEnrich_JASPAR .... \n")
    # showNotification("Searching PWMEnrich_JASPAR .... ",duration = 3,type = "message")
    PWMEnrich_JASPAR <- get_data("PWMEnrich_JASPAR","TF",tf)
    if (length(PWMEnrich_JASPAR) > 0){
      PWMEnrich_JASPAR <- PWMEnrich_JASPAR[PWMEnrich_JASPAR$raw.score > PWMEnrich_JASPAR.score,]
      PWMEnrich_JASPAR <- merge(PWMEnrich_JASPAR,refgene,by="Target")
      colnames(PWMEnrich_JASPAR)[c(1,5)] <- c("Target_refseq","Target")
    }
    targets[["PWMEnrich_JASPAR"]] <- PWMEnrich_JASPAR %>% na.omit()
  }
  ###CHEA
  if ("CHEA" %in% datasets){
    cat("Searching CHEA .... \n")
    # showNotification("Searching CHEA .... ",duration = 3,type = "message")
    CHEA <- get_data("CHEA","TF",tf)
    targets[["CHEA"]] <- CHEA %>% na.omit()
  }
  ####KnockTF##数据
  if ("KnockTF" %in% datasets){
    cat("Searching KnockTF .... \n")
    KnockTF <- get_data("knocktf","TF",tf)
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
    targets[["KnockTF"]] <- KnockTF %>% na.omit()
  }
  ####GTRD##数据
  if ("GTRD" %in% datasets){
    cat("Searching GTRD .... \n")
    # showNotification("Searching GTRD .... ",duration = 3,type = "message")
    GTRD <- get_data("GTRD","TF",tf)
    targets[["GTRD"]] <- GTRD %>% na.omit()
  }
  ####ChIP_Atlas##数据
  if ("ChIP_Atlas" %in% datasets){
    cat("Searching ChIP_Atlas .... \n")
    # showNotification("Searching ChIP_Atlas .... ",duration = 3,type = "message")
    ChIP_Atlas <- get_data("ChIP_Atlas","TF",tf)
    targets[["ChIP_Atlas"]] <- ChIP_Atlas %>% na.omit()
  }

  ####cor_TCGA##数据
  if ("TCGA" %in% cor_DB){
    cat("Fetching correlation results of TCGA .... \n")
    # showNotification("Fetching correlation results of TCGA .... ",duration = 3,type = "message")
    cor_res <- get_data(paste0("cor_",TCGA_tissue),"TF",tf)
    if (length(cor_res) > 0){
      cor_res <- cor_res %>% na.omit() %>%
      rename(.,"cor"=tf ) %>%
      rename(.,"Target"="gene" ) %>%
      dplyr::mutate(.,cor = as.numeric( cor)/100) %>%
      dplyr::filter(.,abs(cor) >= cor_cutoff) %>%
      dplyr::mutate(.,TF = tf)
    }
    targets$cor_TCGA <- cor_res %>% na.omit()
  }
  ####cor_GTEx##数据
  if ("GTEx" %in% cor_DB){
    cat("Fetching correlation results of GTEx .... \n")
    # showNotification("Fetching correlation results of GTEx .... ",duration = 3,type = "message")
    cor_res <- get_data(paste0("cor_",GTEx_tissue),"TF",tf)
    if (length(cor_res) > 0){
      cor_res <- cor_res %>% na.omit() %>%
        rename(.,"cor"=tf ) %>%
        rename(.,"Target"="gene" ) %>%
        dplyr::mutate(.,cor = as.numeric( cor)/100) %>%
        dplyr::filter(.,abs(cor) >= cor_cutoff) %>%
        dplyr::mutate(.,TF = tf)
    }
    targets$cor_GTEx <- cor_res %>% na.omit()
  }

  ##intersect
  results <- list()
  for (i in names(targets)) {
    if (length(targets[[i]])>0){
      print(targets[[i]] %>% as.data.frame() %>% .[,"Target"])
      results[[i]] <- targets[[i]] %>% as.data.frame() %>% .[,"Target"] %>% unique()
    }
  }
  # for (i in names(results)) {
  #   length(results[[i]]) <- lapply(results, length) %>% unlist() %>% max()
  # }
  #
  # results <- as.data.frame(results)

  # ##intersect
  # inter_results <- list()
  # for (i in datasets) {
  #   if (nrow(targets[[i]]) != 0){
  #     inter_results[[i]] <- targets[[i]] %>% as.data.frame() %>% .[,"Target"] %>% unique()
  #   }else{
  #     inter_results[[i]] <- NULL
  #   }
  # }
  #
  # inter_TF <- Reduce(intersect,inter_results)
  # inter_results$intersect <- inter_TF

  targets$results <- results

  return(targets)
}


