#' Load quantification
#'
#' @param useful_data List with necessary information to load quantifications on the Shiny GUI.
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param info List with 'row' and 'column' indicating spectrum and signal to load.
#' @param ROI_data ROIs data
#'
#' @return Loaded plot, signals parameters and quality of fitting parameters of quantifications
#' @export load_quantification
#' @import reshape2
#' @import plotly
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' resulting_data=not_automatic_quant(imported_data,imported_data$final_output,c(1,4),imported_data$ROI_data[1:2,],imported_data$useful_data)
#' loaded_quantification=load_quantification(resulting_data$useful_data,imported_data,resulting_data$final_output,list(row=1,col=1),imported_data$ROI_data)



load_quantification=function(useful_data,imported_data,final_output,info,ROI_data) {
  loaded_quantification=list()
row=info$row
col=info$col
Xdata=useful_data[[row]][[col]]$Xdata
Ydata=useful_data[[row]][[col]]$Ydata
plot_data=useful_data[[row]][[col]]$plot_data
ROI_profile=useful_data[[row]][[col]]$ROI_profile

rownames(plot_data) = c("signals_sum",
                        "baseline_sum",
                        "fitted_sum",
                        as.character(paste(ROI_profile[,4],ROI_profile[,5],sep='_')),rep('additional signal',dim(plot_data)[1]-length(ROI_profile[,4])-3))

plotdata2 = data.frame(Xdata,
                       Ydata,
                       plot_data[3, ],
                       plot_data[2, ] )
plotdata3 <- reshape2::melt(plotdata2, id = "Xdata")
plotdata3$variable = c(
  rep('Original Spectrum', length(Ydata)),
  rep('Generated Spectrum', length(Ydata)),
  rep('Generated Background', length(Ydata))
)
plot_title = paste(imported_data$Experiments[row],"- ROI ",ROI_profile[1,1],"-",ROI_profile[1,2],"ppm")
colors=c(I('red'),I('blue'),I('black'),I('brown'),I('cyan'),I('green'),I('yellow'))
loaded_quantification$plot=plotly::plot_ly(plotdata3,x=~Xdata,y=~value,color=~variable,type='scatter',mode='lines',fill=NULL) %>% layout(title = plot_title,xaxis = list(range=c(Xdata[1],Xdata[length(Xdata)]),title = 'ppm'), yaxis = list(range=c(0,max(Ydata)),title = 'Intensity'))
for (i in 4:nrow(plot_data)) {
  plotdata5 =  data.frame(Xdata=Xdata, variable=rownames(plot_data)[i] ,value=plot_data[i,])

  loaded_quantification$plot=loaded_quantification$plot%>%add_trace(data=plotdata5,x=~Xdata,y=~value,name=~variable,type='scatter',mode='lines',fill='tozeroy',fillcolor=colors[i-3])
}
#Preparation of ROI parameters and of indicators of quality of quantification
loaded_quantification$ROIpar=ROI_profile
loaded_quantification$signpar=matrix(NA,2,7)
colnames(loaded_quantification$signpar)=c("intensity",	"shift",	"half bandwidth",	"gaussian",	"J_coupling",	"multiplicities",	"roof_effect")

if (!is.null(useful_data[[row]][[col]]$signals_parameters)) {
  loaded_quantification$signpar=t(useful_data[[row]][[col]]$signals_parameters)
  if (is.null(rownames(loaded_quantification$signpar))) {
  if (ROI_profile[1,3]=="Baseline Fitting") {
    rownames(loaded_quantification$signpar)=c(as.character(paste(ROI_profile[,4],ROI_profile[,5],sep='_')),rep('baseline signal',nrow(loaded_quantification$signpar)-nrow(ROI_profile)))
  } else {
    rownames(loaded_quantification$signpar)=as.character(paste(ROI_profile[,4],ROI_profile[,5],sep='_'))
}}}


dummy = which(is.na(ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(ROI_data)[1]+1
    lal=which(duplicated(ROI_data[-dummy,1:2])==F)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(ROI_data[-dummy,])[1]))

	ind=which(ROI_separator[,2]-col>=0)[1]
	loaded_quantification$ind=(ROI_separator[ind, 1]:ROI_separator[ind, 2])

	loaded_quantification$qualitypar=cbind(t(final_output$quantification[row,info$col,drop=F]),t(final_output$fitting_error[row,info$col,drop=F]),t(final_output$signal_area_ratio[row,info$col,drop=F]))
	colnames(loaded_quantification$qualitypar)=c('Quantification','fitting_error','signal/total spectrum ratio')
	rownames(loaded_quantification$qualitypar)=imported_data$signals_names[col]
return(loaded_quantification)
}
