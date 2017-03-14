
type_plot = function(imported_data,ROI_limits,rows_selected,medianplot,clusterplot) {

  range=c(0,max(imported_data$dataset[,which.min(abs(imported_data$ppm-ROI_limits[1])):which.min(abs(imported_data$ppm-ROI_limits[2]))]))
  if (length(rows_selected)>1&any(rows_selected<3)) {
    print('Wrong selection.')
    return(NULL)
  } else if (rows_selected==1) {
    clusterplot %>% layout(xaxis = list(range = c(ROI_limits[1], ROI_limits[2]),title='ppm'),yaxis = list(range = range,title='Intensity'))
  }else if (rows_selected==2) {
    medianplot %>% layout(xaxis = list(range = c(ROI_limits[1], ROI_limits[2]),title='ppm'),yaxis = list(range = range,title='Intensity'))
  }else if (all(rows_selected>2)){
    plotdata = data.frame(Xdata=imported_data$ppm, t(imported_data$dataset[rows_selected-2,,drop=F]))
    plotdata <- melt(plotdata, id = "Xdata")
    plot_ly(data=plotdata,x=~Xdata,y=~ value,color=~variable,type='scatter',mode='lines')%>% layout(xaxis = list(range = c(ROI_limits[1], ROI_limits[2]),title='ppm'),yaxis = list(range = range,title='Intensity'))

  }


}
