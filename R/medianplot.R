#' Creation of plotly figure with the median spectrum for each group of spectra
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#'
#' @return Plotly figure with the median spectrum for each group of spectra
#' @export median_plot
#' @import plot_ly
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' median_plot=median_plot(imported_data)

median_plot = function(imported_data) {
  types=unique(imported_data$Metadata[,2])
  mediandataset=matrix(NA,length(types),ncol(imported_data$dataset))
  for (i in 1:length(types)) mediandataset[i,]=apply(imported_data$dataset[which(imported_data$Metadata[,2]==types[i]),,drop=F],2,median)
  p_value_bucketing=as.vector(p_values(imported_data$dataset,imported_data$Metadata))
  p_value_bucketing[p_value_bucketing>0.1]=0.1
  p_value_bucketing= matrix(p_value_bucketing,1,length(p_value_bucketing))

  p=plot_ly(x=~imported_data$ppm)
  for (i in seq(nrow(mediandataset))){
    p=p%>%add_lines(y = mediandataset[i,],name=types[i])
  }
  p=p%>%layout(xaxis=list(title='ppm',range=c(max(imported_data$ppm),min(imported_data$ppm))),yaxis=list(title = "Intensity (arbitrary unit)"))
  p2 <- plot_ly(x=~imported_data$ppm,z =p_value_bucketing, colorscale = "Greys", type = "heatmap")%>%
    layout(xaxis=list(title='ppm',range=c(max(imported_data$ppm),min(imported_data$ppm))))
  p <- subplot(p, p2,nrows=2,heights=c(0.95,0.05),margin=0,shareX = T)
  return(p)
}
