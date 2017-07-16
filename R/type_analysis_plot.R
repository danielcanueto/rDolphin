#' Type of analysis plot
#'
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param data Data to be used to get the analysis plot
#' @param type Kind of plot wanted ('boxplot','pca','dendrogram_heatmap')
#'
#' @return Analysis plot
#' @export type_analysis_plot
#' @import plotly
#' @import missForest
#' @import heatmaply
#' @import reshape2
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' load("MTBLS242_subset_example.RData")
#' type_analysis_plot(quantification_variables$final_output$quantification,quantification_variables$final_output,imported_data,'boxplot')


type_analysis_plot = function(data,final_output,imported_data,type=c('boxplot','pca','dendrogram_heatmap')) {

  m <- list(l = 150, r = 0, b = 150, t = 0,pad = 4)
  ind=apply(data,2,function(x)!all(is.na(x)))
data=data[,ind]
if (imported_data$program_parameters$automatic_removal=='Y') {
  data[final_output$fitting_error[,ind]>imported_data$program_parameters$fitting_error_analysis_limit]=NA
  data[final_output$signal_area_ratio[,ind]<imported_data$program_parameters$signal_area_ratio_analysis_limit]=NA
  data=data[,!apply(data,2,function(x)length(which(is.na(x)))>0.5*length(x))]
}

switch(type, boxplot = {
p_value_final=p_values(data,imported_data$Metadata)
boxplotdata=as.data.frame(data)
colnames(boxplotdata)=paste(colnames(boxplotdata),'(p= ',p_value_final,')',sep='')
boxplotdata=cbind(boxplotdata,factor(imported_data$Metadata[,2]))
boxplotdata=reshape2::melt(boxplotdata)
colnames(boxplotdata)=c('Metadata','Signal','Value')
plot_ly(boxplotdata, x = ~Signal, y = ~Value, color = ~Metadata, type = "box") %>%
  layout(boxmode='group',margin=m)
}, pca = {
a=cbind(scale(data),imported_data$Metadata)
a=missForest(a)$ximp
b=prcomp(a[,-c(ncol(a)-1,ncol(a))])
carsDf2 <- data.frame(b$rotation)
carsDf <- data.frame(b$x,metadata=imported_data$Metadata)
colnames(carsDf)[length(colnames(carsDf))]='metadata'
p <- plot_ly(x=~ carsDf2$PC1,y=~ carsDf2$PC2,type='scatter',
  mode=~"markers",text = rownames(carsDf2),color='loadings',marker=list(size=8))%>% add_trace(x=~ carsDf$PC1,y=~ carsDf$PC2,
    mode=~"markers",text = rownames(carsDf),color =~ (as.factor(carsDf$metadata)),marker=list(size=11))
p <- layout(p,title="PCA scores and loadings",
  xaxis=list(title="PC1"),
  yaxis=list(title="PC2"),margin=m)
print(p)
}, dendrogram_heatmap = {
  data=as.data.frame(scale(data))
ind=apply(data,2,function(x)!all(is.na(x)))
data=data[,ind]
heatmaply::heatmaply(data[seq(nrow(data)),])
})
}
