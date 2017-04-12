#' Helper function to generate ROI data. Not to be used in console.
#'
#' @param datapath Type of valdiation to perform (1: fitting error, 2: signal area ratio, 3: chemical shift, 4: half bandwidth, 5: outliers, 6: relative intensity of signals of same metabolite)
#' @param ROI_data ROIs data
#' @param Metadata Metadata
#' @param Experiments Experiments
#' @return dummy
#' @export roifunc

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