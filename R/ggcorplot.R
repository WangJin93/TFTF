#' @title Plotting scatter plot
#' @description
#'  Scatter plot with sample size (n), correlation coefficient (r) and p value (p.value).
#' @import ggplot2
#' @param data A gene expression dataset with at least two genes included, rows represent samples, and columns represent gene expression in the matrix.
#' @param a Gene A
#' @param b Gene B
#' @param method Method for correlation analysis, "pearson" or "spearman".
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#'
#' @export
#'
ggcorplot <- function(data,
                      a,b,
                      method="spearman",
                      x_lab=" relative expression",
                      y_lab=" relative expression"){
  corr_eqn <- function(x,y,digits=3) {
    test <- cor.test(x,y,method=method,exact=FALSE)
    pp<-round(test$p.value,digits)
    paste(paste0("n = ",length(x)),
          paste0("r = ",round(test$estimate,digits),"(",method,")"),
          paste0("p.value ",ifelse(pp==0,"< 0.001",paste0("= ",pp))),
          sep = ", ")
  }
  data <- data[,c(a,b)]
  names(data) <- c("geneA","geneB")
  ggplot(data,aes(geneA,geneB))+
    geom_point(col="black")+
    geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="red")+
    geom_rug(col="#006fbc")+
    theme_minimal()+
    xlab(paste0(a,x_lab))+
    ylab(paste0(b,y_lab))+
    labs(title = corr_eqn(data$geneA,data$geneB))+
    theme(plot.title = element_text(size=16, hjust = 0.5),
          plot.margin = margin(1, 1, 1, 1, "cm"),axis.title= element_text(size=16, color="black"))+theme(axis.text= element_text(size=16, color="black"))
}
