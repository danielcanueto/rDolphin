#' STOCSY or RANSY identification tool.
#'
#' @param dataset Dataset
#' @param ppm ppm axis
#' @param driver_peak_edges Left and right edges of the driver peak to correlate with the rest of spectrum.
#' @param method Method ('pearson','spearman','ransy').
#'
#' @return interactive plot
#' @export identification_tool
#' @import plotly
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' identification_tool(dataset,ppm,c(1,0.995),'spearman')



identification_tool= function(dataset,ppm,driver_peak_edges,method) {
  ROI_buckets = which.min(abs(as.numeric(driver_peak_edges[1])-ppm)):which.min(abs(as.numeric(driver_peak_edges[2])-ppm))

  if (method=='ransy') {
    #Loosely based on equivalent function in 'muma' R package
    driver = rowSums(dataset[, ROI_buckets]) %*% matrix(rep(1, ncol(dataset)),nrow=1)
    R = apply(dataset/driver,2,function(x)median(x,na.rm=T))/apply(dataset/driver, 2, function(x)mad(x,na.rm=T))
    R[R==Inf]=max(R[is.finite(R)])
    R [is.na(R)]=0
    cor_values=R/max(R,na.rm=T)
  } else {
    cor_values=as.vector(cor(rowSums(dataset[,ROI_buckets]),dataset,method=method))
  }


Ydata = apply(dataset,2,median)
p=plot_ly()%>%
  add_lines(x=~ppm,y = ~Ydata,name='Median spectrum')%>%
  add_lines(x=~ppm,y = ~cor_values, name='Correlation',yaxis = "y2")%>%
  layout(xaxis=list(title='ppm',range=c(max(ppm),min(ppm))),yaxis=list(title = "Intensity (arbitrary unit)"))
p2 <- plot_ly(x=~ppm)%>%add_lines(y =cor_values, name='correlation value',line = list(color = 'rgba(255, 0, 0, 1)'))%>%
  layout(xaxis=list(title='ppm',range=c(max(ppm),min(ppm))))
p <- subplot(p, p2,nrows=2,heights=c(0.75,0.25),margin=0,shareX = T)
return(p)
}
