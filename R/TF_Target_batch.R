#' @title Predict the target genes of TFs
#' @description
#'  Predict the target genes of Transcription Factors in multiple TF-target prediction databases and correlation analysis.
#' @import RMySQL dplyr jsonlite tibble
#' @param datasets Including the TF-target regulatory data of 7 public TF online tools and the TF-targets prediction results analyzed by using FIMO and PWMEnrich.
#' @param tfs Transcription factor names
#' @param TCGA_tissue Cancer type in TCGA database, you can use tissue_type("TCGA") to abtain the tissue types.
#' @param GTEx_tissue Cancer type in GTEx database, you can use tissue_type("GTEx") to abtain the tissue types.
#' @param cor_DB The database used for the correlation analyze between TF and targets. You can use 2 databases, viz. TCGA (33 cancer types) and GTEx (31 normal tissue types).
#' @param cor_cutoff Threshold of correlation coefficient for correlation analysis.
#' @param FIMO.score Threshold of the score of the prediction TF-target results by using FIMO algorithm.
#' @param PWMEnrich.score Threshold of the score of the prediction TF-target results by using PWMEnrich algorithm..
#' @param cut.log2FC Threshold of log2FC for KnockTF dataset.
#' @param down.only Logic value. If true, only the downregulated genes in TF knockout/knockdown cells were returned in KnockTF dataset.
#' @param app Logic value. TRUE only used in the shiny app.
#' @examples
#' \dontrun{
#' results <- TF_Target_batch(datasets=c("hTFtarget","KnockTF"), cor_DB = c("TCGA","GTEx"),tf = c("STAT3","FOXM1"))
#' }
#' @export
#'
TF_Target_batch <- function(datasets=c("FIMO_JASPAR"),
                            tfs = c("STAT3","AHR"),
                            TCGA_tissue = "COAD",
                            GTEx_tissue = "Colon",
                            cor_DB = c("TCGA","GTEx"),
                            cor_cutoff = 0.3,
                            FIMO.score=10,
                            PWMEnrich.score =10,
                            cut.log2FC = 1,
                            down.only = T,
                            app = F){
  all_results <- do.call(rbind,
                         lapply(tfs,function(x){
                           dd <- predict_target(datasets,
                                          cor_DB = cor_DB,
                                          TCGA_tissue = TCGA_tissue,
                                          GTEx_tissue = GTEx_tissue,
                                          cor_cutoff = cor_cutoff,
                                          FIMO.score=FIMO.score,
                                          PWMEnrich.score =PWMEnrich.score,
                                          cut.log2FC = cut.log2FC,
                                          down.only = down.only,
                                          tf = x,
                                          app = app)
                           data.frame(tf = x,
                                           Target =   intersections(dd)[["intersection"]])
                         }
                         )
  )
  all_results <- all_results[which(all_results$Target != "None"),]
  return(all_results)

}
