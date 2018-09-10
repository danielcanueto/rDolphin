
#' Automatic quantification of signals for all experiments using the information located in the ROI patterns file.
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param ROI_data ROIs data.
#' @param optimization By default TRUE. If TRUE, profiling quality is maximized through the analysis og signals parameters, with the tradeoff of additional computing time.
#' @param spectra_to_profile By default NA, and all spectra are profiled. If a vector of spectra is provided, profiling is limited to this vector.
#'
#' @return List with final_output (with metabolite signal relative concentrations and quality indicators) and reproducibility_data (with the necessary data to reproduce the profiling performed).
#' @export automatic_profiling
#' @import baseline
#' @import caret
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' # Not run:
#' # profiling_data=automatic_profiling(imported_data,imported_data$ROI_data)


automatic_profiling = function(imported_data, ROI_data,optimization=TRUE,spectra_to_profile=NULL) {

  signals_names=make.names(paste(ROI_data[,4],ROI_data[,5],sep='_'))
	dummy=matrix(NaN,nrow(imported_data$dataset),length(signals_names),dimnames=list(imported_data$Experiments,signals_names))
  final_output = list(quantification= dummy,signal_area_ratio = dummy,fitting_error = dummy, chemical_shift = dummy,intensity = dummy, half_bandwidth = dummy)

  #creation of list of necessary parameters to load quantifications and evaluate quality of them
  reproducibility_data=vector('list',length(imported_data$Experiments))
  for (i in seq_along(reproducibility_data)) reproducibility_data[[i]]=vector('list',length(signals_names))
  for (i in seq_along(reproducibility_data)) {
    for (j in seq_along(reproducibility_data[[i]])) {
      reproducibility_data[[i]][[j]]=list(Ydata=NULL,Xdata=NULL,ROI_profile=imported_data$ROI_data[j,],program_parameters=NULL,plot_data=NULL,FeaturesMatrix=NULL,signals_parameters=NULL,results_to_save=NULL,error1=1000000)
    }}
  print('Be patient. Gonna take a while. You should be writing, meanwhile.')

  #Splitting of ROI data into individual ROIs to be quantified
	dummy = which(is.na(ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(ROI_data)[1]+1
    lal=which(duplicated(ROI_data[-dummy,1:2])==F)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(ROI_data[-dummy,])[1]))

  baselinedataset=baseline::baseline.rollingBall(imported_data$dataset,5,5)$baseline

  #For every ROI
  for (ROI_index in seq_along(ROI_separator[, 1])) {


    #Preparation of ROI parameters
    ROI_profile = ROI_data[ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2],]
    ROI_buckets = which.min(abs(as.numeric(ROI_profile[1, 1])-imported_data$ppm)):which.min(abs(as.numeric(ROI_profile[1, 2])-imported_data$ppm))
 if (length(ROI_buckets)<20) {
	print ("Ignoring ROI as width is too small")
	next
	}
	if (ROI_buckets[1]>ROI_buckets[2]) ROI_buckets=rev(ROI_buckets)


    #Preparation of program parameters to be sued during fitting, with some variables added to ease interpretability of code
    program_parameters=imported_data$program_parameters
    program_parameters$freq = imported_data$freq
    program_parameters$ROI_buckets = ROI_buckets
    program_parameters$buck_step = imported_data$buck_step

    Xdata = imported_data$ppm[ROI_buckets]
    fitting_type = as.character(ROI_profile[1, 3])
    if (length(grep("Clean",fitting_type))==1) {
      program_parameters$clean_fit="Y"
    } else {
      program_parameters$clean_fit="N"
    }
	signals_to_quantify = which(ROI_profile[, 5] >= 1)
	signals_codes = (ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2])


    print(paste(ROI_profile[1,1], ROI_profile[1,2], sep = '-'))
    print(paste('ROI',ROI_index,'of',nrow(ROI_separator)))

    #Quantification for every spectrum
    pb   <- txtProgressBar(1, nrow(imported_data$dataset), style=3)
    if (is.null(spectra_to_profile)) spectra_to_profile=1:nrow(imported_data$dataset)

    for (spectrum_index in spectra_to_profile) {
      output=profiling_func(spectrum_index,signals_codes,
                            imported_data,
                            ROI_buckets,fitting_type,
                            program_parameters,Xdata,Ydata,
                            final_output,
                            reproducibility_data,
                            ROI_profile,baselinedataset,
                            signals_to_quantify,pb)
      final_output=output$final_output
      reproducibility_data=output$reproducibility_data
    }
}

  tryCatch({
  if (optimization==TRUE&length(spectra_to_profile)>20&nrow(ROI_data)>20) {
  optimized_profiling_data=automatic_profiling_improv(imported_data,final_output,reproducibility_data,ROI_data)
  nn=optimized_profiling_data$final_output$fitting_error-final_output$fitting_error
  no=boxplot.stats(nn)$stats[5]
  ind=which(nn>no)
  for (i in 1:length(optimized_profiling_data$final_output)) {
    optimized_profiling_data$final_output[[i]][ind]=final_output[[i]][ind]
  }
  for (i in 1:nrow(optimized_profiling_data$final_output$fitting_error)) {
    ind=which(nn[i,]>no)
    if (length(ind)>0) optimized_profiling_data$reproducibility_data[[i]][ind]=reproducibility_data[[i]][ind]
  }
  final_output=optimized_profiling_data$final_output
  reproducibility_data=optimized_profiling_data$reproducibility_data
  }
},error=function(e)NA)
  profiling_data=list(final_output=final_output,reproducibility_data=reproducibility_data)

  print("Done!")

  return(profiling_data)
}

