#' Creation of plotly figure with the median spectrum for each group of spectra
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#'
#' @return Plotly figure with the median spectrum for each group of spectra
#' @export medianplot
#' @import reshape2
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' median_plot=medianplot(imported_data)

medianplot = function(imported_data) {
  types=unique(imported_data$Metadata[,2])
  mediandataset=matrix(NA,length(types),ncol(imported_data$dataset))
  for (i in 1:length(types)) mediandataset[i,]=apply(imported_data$dataset[which(imported_data$Metadata[,2]==types[i]),,drop=F],2,median)
  p_value_bucketing=as.vector(p_values(imported_data$dataset,imported_data$Metadata))

  ay <- list(tickfont = list(color = "red"),overlaying = "y",side = "right",title = "p value",range = c(0,max(mediandataset)))
  az = list(title = "Intensity",range = c(-1, max(mediandataset)-1))

  p=plot_ly(x=~imported_data$ppm)
  for (i in seq(nrow(mediandataset))){
    p=p%>%add_lines(y = mediandataset[i,],name=types[i])
  }
  p=p%>%add_lines(y = ~p_value_bucketing,name='p value', yaxis = "y2",line = list(color = 'rgba(255, 0, 0, 1)'))%>%
    layout(xaxis=list(title='ppm',range=c(max(imported_data$ppm),min(imported_data$ppm))),yaxis=az, yaxis2 = ay)
  return(p)
}
