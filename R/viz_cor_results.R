#' @title Visualization of pan-tissue correlation analysis
#' @description
#'  Visualization of pan-tissue correlation analysis using ggplot2.
#' @import ggplot2
#' @param cor_results Correlation analysis results obtained from pantissue_cor_analysis() function.
#' @param values Color of the points， default c("green","black","red").
#' @param alpha Transparency of points， default 0.6.
#' @param size Size of points， default 3
#' @param r.cut Threshold of correlation coefficient， default 0.6.
#' @param p.cut Threshold of P value， default 0.05.
#' @param max.overlaps Max overlaps of the point labels, default 10.
#' @examples
#' \dontrun{
#' cor_results <- pantissue_cor_analysis(Gene1 = "FOXM1",Gene2 = "GAPDH")
#' viz_cor_results(cor_results,
#'                 values = c("black","red"))
#' }
#' @export
#'
viz_cor_results <- function(cor_results,
                            values = c("green","black","red"),
                            alpha = 0.6,
                            size = 3,
                            r.cut = 0.3,
                            p.cut = 0.05,
                            max.overlaps = 10
                            ){
  output <- cor_results[["cor_result"]]
  output$sigmark <- ifelse(output$p>0.05,"No",ifelse(output$r> r.cut,"Positive",
                                                     ifelse(output$r< -r.cut,"Negative","No")))

  output$label=ifelse(output$sigmark != "No",as.vector(output$tissue),"")
  p <- ggplot(data = output,aes(x = r,y = logP)) +
    geom_point(alpha=alpha, size= size,aes(color=sigmark)) +
    scale_color_manual(values= values)+      ##########color
    geom_vline(xintercept=c(r.cut,-r.cut),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(as.numeric(p.cut)) ,lty=4,col="black",lwd=0.8) +
    ylab("-log10(P value)")+xlab("Correlation coefficient")+
    theme_bw()

  p<- p + ggrepel::geom_label_repel(aes(label = label),data = output,    color="black" ,max.overlaps =  max.overlaps)+
    theme(axis.title.x =element_text(size=14,    color="black" ), axis.title.y=element_text(size=14,    color="black" ),
          axis.text.x = element_text(size=14,    color="black" ), axis.text.y=element_text(size=14,    color="black"),
          legend.text=element_text(size=14),legend.title=element_text(size=14))


  p+ ggplot2::theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15)
  )
}

