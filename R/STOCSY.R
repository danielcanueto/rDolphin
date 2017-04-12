#' Automatic quantification of signals for all experiments using the information located in the ROI patterns file.
#'
#' @param dataset Dataset
#' @param ppm ppm
#' @param limits Left and right edges of the region to correalte with the rest of spectrum.
#' @param method Correlation method ('pearson','spearman').
#'
#' @return STOCSY interactive plot
#' @export STOCSY
#' @import plotly
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' STOCSY(imported_data$dataset,imported_data$ppm,c(1,0.995),'spearman')



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
