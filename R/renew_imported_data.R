#' Adaptation of variables of imported data to changes applied to the ROI profiles
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#'
#' @return Changed imported data
#' @export renew_imported_data
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' imported_data$ROI_data=imported_data$ROI_data[-1,]
#' imported_data=renew_imported_data(imported_data)

renew_imported_data= function(imported_data) {

imported_data$signals_names=paste(imported_data$ROI_data[which(!is.na(imported_data$ROI_data[, 1])),4],imported_data$ROI_data[which(!is.na(imported_data$ROI_data[, 1])),5],sep='_')
  imported_data$signals_codes = seq_along(signals_names)

  dummy=matrix(NaN,nrow(imported_data$dataset),length(imported_data$signals_names),dimnames=list(imported_data$Experiments,imported_data$signals_names))
  imported_data$final_output = list(Area= dummy,signal_area_ratio = dummy,fitting_error = dummy, shift = dummy,intensity = dummy, half_band_width = dummy)

  #creation of list of necessary parameters to load quantifications and evaluate quality of them
  imported_data$useful_data=vector('list',length(imported_data$Experiments))
  for (i in seq_along(imported_data$useful_data)) imported_data$useful_data[[i]]=vector('list',length(imported_data$signals_codes))
  for (i in seq_along(imported_data$useful_data)) {
    for (j in seq_along(imported_data$useful_data[[i]])) {
      imported_data$useful_data[[i]][[j]]=list(Ydata=NULL,Xdata=NULL,ROI_profile=NULL,program_parameters=NULL,plot_data=NULL,FeaturesMatrix=NULL,signals_parameters=NULL,results_to_save=NULL,error1=1000000)
    }}
	
	return(imported_data)
	}