#' Automatic quantification of signals for all experiments using the information located in the ROI patterns file.
#'
#' @param dataset Dataset
#' @param ppm ppm
#' @param limits Limits
#' @param method Correlation method.
#' @param ROI_separator ROI separator.

#'
#' @return STOCSY plot
#' @export STOCSY
#' @import plotly


STOCSY= function(dataset,ppm,limits,method) {
  ROI_buckets = which.min(abs(as.numeric(limits[1])-ppm)):which.min(abs(as.numeric(limits[2])-ppm))
  if (ROI_buckets[2]>ROI_buckets[1]) {
cor_values=as.vector(cor(rowSums(dataset[,ROI_buckets]),dataset,method=method))
  } else {
    cor_values=rep(0,length(ppm))
  }
Ydata = apply(dataset,2,median)
ay <- list(tickfont = list(color = "red"),overlaying = "y",side = "right",title = "correlation",range = c(0, max(Ydata)))
az = list(title = "Intensity",range = c(-1, max(Ydata)-1))
p=plot_ly()%>%
  add_lines(x=~ppm,y = ~Ydata,name='Median spectrum')%>%
  add_lines(x=~ppm,y = ~cor_values, name='Correlation',yaxis = "y2")%>%
  layout(xaxis=list(title='ppm',range=c(max(ppm),min(ppm))),yaxis=az, yaxis2 = ay)

return(p)
}
