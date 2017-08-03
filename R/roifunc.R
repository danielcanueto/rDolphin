#' Helper function
#' @param ROI_data
#' @param Metadata
#' @param Experiments
#' @return dummy

roifunc <- function(ROI_data,Metadata,Experiments) {
	dummy = which(is.na(ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(ROI_data)[1]+1
    lal=which(duplicated(ROI_data[-dummy,1:2])==F)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(ROI_data[-dummy,])[1]))

    ROI_names=paste(ROI_data[ROI_separator[, 1],1],ROI_data[ROI_separator[, 1],2])
    select_options=seq_along(ROI_names)
    names(select_options)=ROI_names
    mm=matrix(NA,2,dim(Metadata)[2])
    colnames(mm)=colnames(Metadata)
    spectra=cbind(c('Exemplars','Median Spectrum per group',Experiments),rbind(mm,Metadata))
    colnames(spectra)=c('spectrum',colnames(mm))
    dummy=list(select_options=select_options,spectra=spectra)
	return(dummy)
  }
