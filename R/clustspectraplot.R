#' Creation of plotly figure of subset of representative spectra of the dataset. The plot generated helps to determine the best parameters for quantification (e.g. number of signals to fit, chemical $chemical_shift tolerance...).
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#'
#' @return Plotly figure with subset of spectra that represent different kinds of spectra that can be found on the dataset depending on chemical shift, half bandwidth, etc.
#' @export exemplars_plot
#' @import apcluster
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' plot=exemplars_plot(imported_data)




exemplars_plot = function(imported_data) {

  if (nrow(imported_data$dataset)>10) {
  scaled_roi=scale(imported_data$dataset[ , sort(colMeans(imported_data$dataset),decreasing=T,index.return=T)$ix[1:(ncol(imported_data$dataset)/3)],drop=F])
  updated_scaled_roi=scaled_roi
  rm_ind=c()
  stop=0
ind=seq(nrow(imported_data$dataset))
  while ((!is.null(rm_ind)|stop==0)&(nrow(scaled_roi)>15)) {
    stop=1
    if (length(rm_ind)>0) scaled_roi=scaled_roi[-rm_ind,]
    rm_ind=c()
    apres <- suppressWarnings(apcluster::apclusterK(apcluster::negDistMat(r=2), scaled_roi, K=min(c(dim(scaled_roi)[1]-1,10)),verbose=F))
    for (i in 1:length(apres@clusters)) {
      if (length(apres@clusters[[i]])==1) rm_ind=c(rm_ind,apres@clusters[[i]][1])
    }
    ind=apres@exemplars
  }
    updated_scaled_roi=scaled_roi[ind,,drop=F]
  spectra_lag=rep(NA,nrow(updated_scaled_roi))
  dummy=apply(updated_scaled_roi, 2, function(x)  median(x,na.rm=T))

  for (i in 1:nrow(updated_scaled_roi)) {
    d <-ccf(updated_scaled_roi[i,],dummy,type = 'covariance',plot = FALSE)
    spectra_lag[i]=d$lag[which.max(d$acf)]
  }
  visual_roi=original_roi=imported_data$dataset[ind[sort(spectra_lag,index.return=T)$ix],,drop=F]

  } else {
    visual_roi=original_roi=imported_data$dataset
    ind=seq(nrow(imported_data$dataset))
}
  # for (i in 1:nrow(original_roi))   visual_roi[i,]=original_roi[i,]+(i-1)*mean(original_roi)

  p_value_bucketing=as.vector(p_values(imported_data$dataset,imported_data$Metadata))
  p_value_bucketing[p_value_bucketing>0.1]=0.1

  p_value_bucketing= matrix(p_value_bucketing,1,length(p_value_bucketing))

  p=plot_ly(x=~imported_data$ppm)
  shade=as.character(seq(0.2,1,length.out = nrow(visual_roi)))
  for (i in seq(nrow(visual_roi))){
        p=p%>%add_lines(y = visual_roi[i,],name=imported_data$Experiments[ ind[i]],line = list(color = paste('rgba(0, 0, 255,',shade[i],')',sep='')))
    }
  p=p%>%layout(xaxis=list(title='ppm',range=c(max(imported_data$ppm),min(imported_data$ppm))),yaxis=list(title = "Intensity (arbitrary unit)"))
  p2 <- plot_ly(x=~imported_data$ppm,z =p_value_bucketing, colorscale = "Greys", type = "heatmap")%>%    layout(xaxis=list(title='ppm',range=c(max(imported_data$ppm),min(imported_data$ppm))))
  p <- subplot(p, p2,nrows=2,heights=c(0.95,0.05),margin=0,shareX = T)



  return(p)
}
