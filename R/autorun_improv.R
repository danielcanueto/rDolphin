
#' Improvement of automatic profiling using the inforamtion collected on a previous implementation.
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param reproducibility_data List with necessary information to load quantifications on the Shiny GUI.
#' @param ROI_data ROIs data.
#' @param improvement_option If "correction", quantifications are updated taking into account the predicted signal parameters. If "reimplemetation", profiling is repeated using the prediction information.
#' @param level How extensive should be the improvement? If "all", all quantifications are changed. If "outliers", quantifications whoss signal parameters behave as outliers are changed (please take into account that only the quantifications will be updated). If a number is introduced, the quantifications with a higher fitting error than the number specified are repeated.
#'
#' @return List with updated final_output and reproducibility_data variables.
#' @export automatic_profiling_improv
#' @import baseline
#' @import robustbase
#' @import caret
#' @import missRanger
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' # Not run:
#' # load(file.path(system.file(package = "rDolphin"),"extdata","MTBLS242_subset_profiling_data.RData"))
#' # profiling_data_2=automatic_profiling_improv(imported_data,profiling_data$final_output,profiling_data$reproducibility_data,imported_data$ROI_data)


#TODO: Choose criteria to repeat only individual quantification and all signals of all spectra.

automatic_profiling_improv = function(imported_data, final_output,reproducibility_data,ROI_data,improvement_option='reimplementation',level='outliers') {
print("Starting maximization of profiling data quality using information of original profiling...")
  print("Now estimating the predicted signal parameters with confidence intervals...")
  sink("aux");
  predicted_info=rf_pred(final_output$half_bandwidth,final_output$fitting_error)
predicted_width=as.matrix(predicted_info$predicted_matrix)

  max_width=as.matrix(predicted_info$upper_bound_matrix)
  min_width=as.matrix(predicted_info$lower_bound_matrix)
ind=which(is.na(predicted_width[1,]))
predicted_width[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,8]))
min_width[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,8]*0.75))
max_width[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,8]*1.25))

predicted_info=rf_pred(final_output$chemical_shift,final_output$fitting_error)
predicted_shift=as.matrix(predicted_info$predicted_matrix)
max_shift=as.matrix(predicted_info$upper_bound_matrix)
min_shift=as.matrix(predicted_info$lower_bound_matrix)

  ind=which(is.na(predicted_shift[1,]))
  if (length(ind)>0) {
    predicted_shift[,ind]=as.matrix(t(replicate(nrow(predicted_width),ROI_data[ind,6])))
  max_shift[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,6]+ROI_data[ind,7]))
  min_shift[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,6]-ROI_data[ind,7]))
  }

  predicted_info=rf_pred_intensity(final_output$intensity,ROI_data[,4],final_output$fitting_error)
  predicted_intensity=as.matrix(predicted_info$predicted_matrix)
  max_intensity=as.matrix(predicted_info$upper_bound_matrix)
  min_intensity=as.matrix(predicted_info$lower_bound_matrix)

  ind=which(is.na(predicted_intensity[1,]))
  max_intensity[!is.finite(max_intensity)]=NA
  min_intensity[!is.finite(min_intensity)]=NA
  min_intensity[min_intensity<0]=0
  sink(NULL);
  print("Beginning the maximization of profiling data quality according to the option selected...")


  quantifications_to_repeat=matrix(0,nrow(predicted_width),ncol(predicted_width))
  if (level=="all") quantifications_to_repeat[,]=1
	if (is.numeric(level)) quantifications_to_repeat[which(final_output$fitting_error>level)]=1
	if (level=="outliers") {
	  outlier_indicator=sapply(which(!is.na(predicted_shift)),
	             function(x)findInterval(final_output$chemical_shift[x],
	                                     c(min_shift[x],max_shift[x])))
	  if (length(outlier_indicator)>0) quantifications_to_repeat[which(!is.na(predicted_shift))][sapply(outlier_indicator,function(x)x==0|x==2)]=1
	  outlier_indicator=sapply(which(!is.na(predicted_width)),
	             function(x)findInterval(final_output$half_bandwidth[x],
	                                     c(min_width[x],max_width[x])))
	  if (length(outlier_indicator)>0) quantifications_to_repeat[which(!is.na(predicted_width))][sapply(outlier_indicator,function(x)x==0|x==2)]=1
	  outlier_indicator=sapply(which(!is.na(predicted_intensity)),
	             function(x)findInterval(final_output$intensity[x],
	                                     c(min_intensity[x],max_intensity[x])))
	  if (length(outlier_indicator)>0) quantifications_to_repeat[which(!is.na(predicted_intensity))][sapply(outlier_indicator,function(x)x==0|x==2)]=1
}

if (improvement_option=='reimplementation') {  #Splitting of ROI data into individual ROIs to be quantified
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
    if (length(ROI_buckets)<5) next
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


	index_to_use_3=which(rowSums(quantifications_to_repeat[,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2],drop=F])>0)
	pb   <- txtProgressBar(1, nrow(imported_data$dataset), style=3)

    #Quantification for every spectrum
    for (spectrum_index in index_to_use_3) {

         #Preparation of necessary variables to store figures and information of the fitting
      Ydata = as.numeric(imported_data$dataset[spectrum_index, ROI_buckets])

      #If the quantification is through integration with or without baseline
      if (fitting_type == "Clean Sum" ||
          fitting_type == "Baseline Sum") {
        dummy = integration(program_parameters$clean_fit, Xdata,Ydata,program_parameters$buck_step)

        results_to_save=dummy$results_to_save
        #Generation of useful variables specific of every quantification
        reproducibility_data[[spectrum_index]][[signals_codes]]$ROI_profile=ROI_profile
        reproducibility_data[[spectrum_index]][[signals_codes]]$plot_data=dummy$plot_data
        reproducibility_data[[spectrum_index]][[signals_codes]]$Xdata=Xdata
        reproducibility_data[[spectrum_index]][[signals_codes]]$Ydata=Ydata
        reproducibility_data[[spectrum_index]][[signals_codes]]$results_to_save=results_to_save
        reproducibility_data[[spectrum_index]][[signals_codes]]$error1=results_to_save$fitting_error

        #If the quantification is through fitting with or without baseline
      } else if (fitting_type == "Clean Fitting" || fitting_type ==
          "Baseline Fitting") {

         #Adaptation of the info of the parameters into a single matrix and preparation (if necessary) of the background signals that will conform the baseline
        FeaturesMatrix = fitting_prep_2(Xdata,
          Ydata,
          ROI_profile[, 5:11,drop=F],
          program_parameters,baselinedataset[spectrum_index,ROI_buckets],max_shift,min_shift,max_intensity,
          min_intensity,max_width,min_width,spectrum_index,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2])

        #Calculation of the parameters that will achieve the best fitting
        dummy = fittingloop(FeaturesMatrix,
          Xdata,
          Ydata,
          program_parameters)
        signals_parameters=dummy$signals_parameters
		Xdata_2=imported_data$ppm
		Ydata_2 = as.numeric(imported_data$dataset[spectrum_index, ])
        #Fitting of the signals
        multiplicities=c(FeaturesMatrix[,11],rep(1,(length(signals_parameters)/5)-dim(FeaturesMatrix)[1]))
        roof_effect=c(FeaturesMatrix[,12],rep(0,(length(signals_parameters)/5)-dim(FeaturesMatrix)[1]))
        fitted_signals = signal_fitting(signals_parameters,
          Xdata_2,multiplicities,roof_effect,program_parameters$freq)
               dim(signals_parameters) = c(5, length(signals_parameters)/5)
        rownames(signals_parameters) = c('intensity','$chemical_shift','half_bandwidth','gaussian','J_coupling')
        signals_parameters=rbind(signals_parameters,multiplicities,roof_effect)

        #Generation of output data about the fitting and of the necessary variables for the generation ofa figure
        dummy = output_generator(signals_to_quantify,fitted_signals,Ydata_2,Xdata_2,signals_parameters,multiplicities,program_parameters$buck_step)
        output_data=dummy$output_data
        error1=ifelse(is.nan(dummy$error1),3000,dummy$error1)

        #Generation of the dataframe with the final output variables
        results_to_save = data.frame(
          chemical_shift = output_data$chemical_shift,
          quantification = output_data$quantification,
          signal_area_ratio = output_data$signal_area_ratio,
          fitting_error = output_data$fitting_error,
          intensity = output_data$intensity,
          half_bandwidth = output_data$half_bandwidth
        )

        #Generation of the figure data
        plot_data = rbind(output_data$signals_sum,output_data$baseline_sum,output_data$fitted_sum,output_data$signals)
        plot_data = plot_data[,ROI_buckets]

         rownames(plot_data) = c("signals_sum","baseline_sum","fitted_sum",make.names(paste(ROI_profile[,4],ROI_profile[,5],sep='_')),rep('additional signal',dim(plot_data)[1]-length(ROI_profile[,4])-3))

        #Generation of useful variables specific of every quantification
        for (i in seq_along(signals_codes)) {
          reproducibility_data[[spectrum_index]][[signals_codes[i]]]$ROI_profile=ROI_profile
          reproducibility_data[[spectrum_index]][[signals_codes[i]]]$program_parameters=program_parameters
          reproducibility_data[[spectrum_index]][[signals_codes[i]]]$plot_data=plot_data
          reproducibility_data[[spectrum_index]][[signals_codes[i]]]$error1=error1
          reproducibility_data[[spectrum_index]][[signals_codes[i]]]$FeaturesMatrix=FeaturesMatrix
          reproducibility_data[[spectrum_index]][[signals_codes[i]]]$signals_parameters=signals_parameters
          reproducibility_data[[spectrum_index]][[signals_codes[i]]]$Xdata=Xdata
          reproducibility_data[[spectrum_index]][[signals_codes[i]]]$Ydata=Ydata
          reproducibility_data[[spectrum_index]][[signals_codes[i]]]$results_to_save=results_to_save
          }
     }

      #Generation of output variables specific of every quantification
      final_output = save_output(spectrum_index,signals_codes,results_to_save,imported_data$buck_step,final_output)
      setTxtProgressBar(pb, spectrum_index)

      }

  }


 } else if (improvement_option=='correction') {

prova_intensity=predicted_intensity
prova_intensity[,apply(predicted_intensity,2,function(x)all(is.na(x)))]=final_output$intensity[,apply(predicted_intensity,2,function(x)all(is.na(x)))]
prova_shift=predicted_shift
prova_shift[,apply(predicted_shift,2,function(x)all(is.na(x)))]=final_output$chemical_shift[,apply(predicted_shift,2,function(x)all(is.na(x)))]
prova_width=predicted_width
prova_width[,apply(predicted_width,2,function(x)all(is.na(x)))]=final_output$half_bandwidth[,apply(predicted_width,2,function(x)all(is.na(x)))]

tec=sapply(seq(length(prova_intensity)),function(x)sum(peakpvoigt(c(prova_intensity[x],prova_shift[x],prova_width[x]*0.5/600.2,0),imported_data$ppm))*imported_data$buck_step)
dim(tec)=dim(prova_intensity)
tec[,apply(tec,2,function(x)all(is.na(x)))]=final_output$quantification[,apply(tec,2,function(x)all(is.na(x)))]
	for (i in 1:ncol(final_output$quantification)) {
	index_to_use_3=which(quantifications_to_repeat[,i]>0)
	final_output$quantification[index_to_use_3,i]=tec[index_to_use_3,i]
}
 }
  print("Done!")
  profiling_data=list(final_output=final_output,reproducibility_data=reproducibility_data,
                                predicted_shift=predicted_shift,predicted_width=predicted_width,
                                predicted_intensity=predicted_intensity,max_width=max_width,
                                min_width=min_width,max_shift=max_shift,min_shift=min_shift,
                                max_intensity=max_intensity,min_intensity=min_intensity)
  return(profiling_data)
}

