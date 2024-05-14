#' predict_TF(target = "GAPDH")
#' @importFrom shiny validate
#' @import RMySQL
#' @import dplyr
#' @import jsonlite
#' @import tibble
#'
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
                       cor_DB = c("TCGA","GTEx"),
                       TCGA_tissue = "COAD",
                       GTEx_tissue = "Colon",
                       cor_cutoff = 0.3,
                       target = "GAPDH",
                       FIMO.score=10,
                       PWMEnrich_JASPAR.score =10,
                       cut.log2FC = 1,
                       down.only = T){
  options(timeout=200)
  TF_result <- as.list(datasets)
  names(TF_result) <- datasets

  ###########HTFtarget
  if ("hTFtarget" %in% datasets){
    cat("Searching hTFtarget .... \n")
    # # showNotification("Searching hTFtarget .... ",duration = 2)
    ensenbl <- idmap[idmap$gene == target ,]$ID
    url<-paste0("https://guolab.wchscu.cn/hTFtarget/api/chipseq/targets/target?target=",ensenbl)
    hTFtarget <- jsonlite::fromJSON(url)
    if (is.null(hTFtarget)){
      hTFtarget <- data.frame()
    }else{
      hTFtarget <-  hTFtarget[["regulations"]] %>%
        dplyr::rename("TF" = "tf_id")
      hTFtarget[["regulations"]][["regulation_info"]] <- hTFtarget[["regulations"]][["regulation_info"]][["regulation"]]
    }
    TF_result[["hTFtarget"]] <- hTFtarget[-2] %>% na.omit()
  }
  ####KnockTF##数据
  if ("KnockTF" %in% datasets){
    cat("Searching KnockTF .... \n")
    # showNotification("Searching KnockTF .... ",duration = 2)
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
    # showNotification("Searching TRRUST .... ",duration = 2,type = "message")
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
    # showNotification("Searching ENCODE .... ",duration = 2,type = "message")
    ENCODE <- get_data("ENCODE","Target",target)
    TF_result[["ENCODE"]] <- ENCODE %>% na.omit()
  }
  ###JASPAR
  if ("FIMO_JASPAR" %in% datasets){
    cat("Searching FIMO_JASPAR .... \n")
    # # showNotification("Searching FIMO_JASPAR .... ",duration = 3,type = "message")
    target2 <- refgene[refgene$symbol == target ,]$Target %>% na.omit()
    Jaspar <- get_data("FIMO_JASPAR","Target",target2)
    if (length(Jaspar) > 0){
      Jaspar <- Jaspar[Jaspar$score > FIMO.score,]
    }
    TF_result[["FIMO_JASPAR"]] <- Jaspar %>% na.omit()
  }
  ###PWMEnrich_JASPAR
  if ("PWMEnrich_JASPAR" %in% datasets){
    cat("Searching PWMEnrich_JASPAR .... \n")
    # # showNotification("Searching PWMEnrich_JASPAR .... ",duration = 3,type = "message")
    target2 <- refgene[refgene$symbol == target ,]$Target %>% na.omit()
    PWMEnrich_JASPAR <- get_data("PWMEnrich_JASPAR","Target",target2)
    if (length(PWMEnrich_JASPAR) > 0){
      PWMEnrich_JASPAR <- PWMEnrich_JASPAR[PWMEnrich_JASPAR$raw.score > PWMEnrich_JASPAR.score,]
    }
    TF_result[["PWMEnrich_JASPAR"]] <- PWMEnrich_JASPAR %>% na.omit()
  }
  ###"CHEA"
  if ("CHEA" %in% datasets){
    cat("Searching CHEA .... \n")
    # showNotification("Searching CHEA .... ",duration = 2,type = "message")
      CHEA <- get_data("CHEA","Target",target)
    TF_result[["CHEA"]] <- CHEA %>% na.omit()
  }
  ###"GTRD"
  if ("GTRD" %in% datasets){
    cat("Searching GTRD .... \n")
    # showNotification("Searching GTRD .... ",duration = 2,type = "message")
    GTRD <- get_data("GTRD","Target",target)
    TF_result[["GTRD"]] <- GTRD %>% na.omit()
  }
  ####ChIP_Atlas##数据
  if ("ChIP_Atlas" %in% datasets){
    cat("Searching ChIP_Atlas .... \n")
    # showNotification("Searching ChIP_Atlas .... ",duration = 3,type = "message")
    ChIP_Atlas <- get_data("ChIP_Atlas","Target",target)
    TF_result[["ChIP_Atlas"]] <- ChIP_Atlas %>% na.omit()
  }

  ####cor_TCGA##数据
  if ("TCGA" %in% cor_DB){
    cat("Fetching correlation results of TCGA .... \n")
    # # showNotification("Fetching correlation results of TCGA .... ",duration = 3,type = "message")
    cor_res <- get_data(paste0("cor_",TCGA_tissue),"Target",target)
    rownames(cor_res) <- NULL
    if (length(cor_res) > 0){
      cor_res <- cor_res %>% tibble::column_to_rownames("gene")%>%
      t()%>% as.data.frame()%>%
      na.omit() %>%
      rename(.,"cor"=target ) %>%
      dplyr::mutate(.,cor = as.numeric( cor)/100) %>%
      dplyr::filter(.,abs(cor) >= cor_cutoff) %>%
      tibble::rownames_to_column(.,"TF")
    }
    TF_result$cor_TCGA <- cor_res %>% na.omit()
  }
  ####cor_GTEx##数据
  if ("GTEx" %in% cor_DB){
    cat("Fetching correlation results of GTEx .... \n")
    # # showNotification("Fetching correlation results of GTEx .... ",duration = 3,type = "message")
    cor_res <- get_data(paste0("cor_",GTEx_tissue),"Target",target)
    rownames(cor_res) <- NULL
    if (length(cor_res) > 0){
      cor_res <-  cor_res %>% tibble::column_to_rownames("gene")%>%
      t()%>% as.data.frame()%>%
      na.omit() %>%
      rename(.,"cor"=target ) %>%
      dplyr::mutate(.,cor = as.numeric( cor)/100) %>%
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
