#' Creation of plotly figure with the median spectrum for each group of spectra
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#'
#' @return Plotly figure with the median spectrum for each group of spectra 
#' @export medianplot
#' @import reshape2
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' median_plot=medianplot(imported_data)

medianplot = function(imported_data) {
  types=unique(imported_data$Metadata[,2])
  mediandataset=matrix(NA,length(types),ncol(imported_data$dataset))
  for (i in 1:length(types)) mediandataset[i,]=apply(imported_data$dataset[which(imported_data$Metadata[,2]==types[i]),,drop=F],2,median)

  plotdata = data.frame(Xdata=imported_data$ppm, signals = t(mediandataset ))
  colnames(plotdata)=c('Xdata',types)
 plotdata=melt(plotdata,id = "Xdata")
  median_plot=plot_ly(data=plotdata,x=~Xdata,y=~ value,color=~variable,type='scatter',mode='lines')%>% layout(xaxis = list(range = c(imported_data$ppm[1], imported_data$ppm[length(imported_data$ppm)]),title='ppm'),yaxis = list(range = range,title='Intensity'))

  return(median_plot)
}
