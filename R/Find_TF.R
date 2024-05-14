#' Detect which datasets contain input TFs
#' @import dplyr
#' @param tf Transcription Factor names.
#' @examples
#' \dontrun{
#' Find_TF(c("STAT3","AHR"))
#' }
#' @export
#'
Find_TF <- function(tf = "STAT3"){
  df <- tf_list %>% dplyr::filter(TF %in% tf)
  if (nrow(df)==0){
    cat("The input TFs are not included in the TF list!")
  }else{
    df <- df %>% tibble::column_to_rownames(.,"TF")%>% t() %>% as.data.frame() %>% tibble::rownames_to_column(.,"Dataset")
    df$All_in <- lapply(1:9,function(x) ifelse("F" %in% df[x,],"F","T"))
    print(df)
    }
}
