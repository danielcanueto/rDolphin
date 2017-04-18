#' Type of analysis plot
#'
#' @param rows_selected Rows selected in GUi table.
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param ROI_limits Range of plot
#' @param medianplot Meidan plot
#' @param clusterplot Cluster plot
#'
#' @return Plot
#' @export type_plot
#' @import plotly
#' @import reshape2
#'
type_plot = function(imported_data,ROI_limits,rows_selected,medianplot,clusterplot) {


  if (length(rows_selected)>1&any(rows_selected<3)) {
    print('Not possible to combine cluster or median spectra with individual spectra.')
    return(NULL)
  } else if (rows_selected==1) {
    range=c(0,max(imported_data$dataset[,which.min(abs(imported_data$ppm-ROI_limits[1])):which.min(abs(imported_data$ppm-ROI_limits[2]))])+10)
    p=clusterplot %>% layout(xaxis = list(range = c(ROI_limits[1], ROI_limits[2]),title='ppm'),yaxis = list(range = range,title='Intensity'))
  }else if (rows_selected==2) {
    range=c(0,max(apply(imported_data$dataset[,which.min(abs(imported_data$ppm-ROI_limits[1])):which.min(abs(imported_data$ppm-ROI_limits[2]))],2,median)))
    p=medianplot %>% layout(xaxis = list(range = c(ROI_limits[1], ROI_limits[2]),title='ppm'),yaxis = list(range = range,title='Intensity'))
  }else if (all(rows_selected>2)){
    range=c(0,max(imported_data$dataset[rows_selected-2,which.min(abs(imported_data$ppm-ROI_limits[1])):which.min(abs(imported_data$ppm-ROI_limits[2]))]))
    plotdata = data.frame(Xdata=imported_data$ppm, t(imported_data$dataset[rows_selected-2,,drop=F]))
    plotdata <- melt(plotdata, id = "Xdata")
    p=plot_ly(data=plotdata,x=~Xdata,y=~ value,color=~variable,type='scatter',mode='lines')%>% layout(xaxis = list(range = c(ROI_limits[1], ROI_limits[2]),title='ppm'),yaxis = list(range = range,title='Intensity'))
  }

return(p)
}
