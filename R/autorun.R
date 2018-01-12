
#' Automatic quantification of signals for all experiments using the information located in the ROI patterns file.
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param useful_data List with necessary information to load quantifications on the Shiny GUI.
#' @param ROI_data ROIs data.
#' @param optimization By default TRUE. If TRUE, profiling quality is maximized through the analysis og signals parameters, with the tradeoff of additional computing time.
#' @param spectra_to_fit By default NA, and all spectra are profiled. If a vector of spectra is provided, profiling is limited to this vector.
#'
#' @return List with updated final_output and useful_data variables.
#' @export autorun
#' @import baseline
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' # Not run:
#' # profiling_data=autorun(imported_data,imported_data$final_output,imported_data$useful_data,imported_data$ROI_data)


autorun = function(imported_data, final_output,useful_data,ROI_data,optimization=TRUE,spectra_to_fit=NA) {

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



    #Quantification for every spectrum
    pb   <- txtProgressBar(1, nrow(imported_data$dataset), style=3)
    if (is.na(spectra_to_fit)) spectra_to_fit=1:nrow(imported_data$dataset)
    for (spectrum_index in spectra_to_fit) {

      #Preparation of necessary variables to store figures and information of the fitting
      Ydata = as.numeric(imported_data$dataset[spectrum_index, ROI_buckets])

      #If the quantification is through integration with or without baseline
      if (fitting_type == "Clean Sum" ||
          fitting_type == "Baseline Sum") {
        #Fitting error is calculated through the comparison with the median spectrum, so singals interfering with the integration can be controlled
        # program_parameters$clean_fit = ifelse(fitting_type == "Clean Sum", "Y",
        #                                       "N")
        # program_parameters$freq=imported_data$freq
        # baseline_int = fitting_prep_integration(Xdata,Ydata,program_parameters,baselinedataset[spectrum_index, ROI_buckets])
        # Ydatamedian=as.numeric(apply(imported_data$dataset[, ROI_buckets,drop=F],2,median))
        dummy = integration(program_parameters$clean_fit, Xdata,Ydata,program_parameters$buck_step)

        results_to_save=dummy$results_to_save
        #Generation of useful variables specific of every quantification
        useful_data[[spectrum_index]][[signals_codes]]$ROI_profile=ROI_profile
        useful_data[[spectrum_index]][[signals_codes]]$plot_data=dummy$plot_data
        useful_data[[spectrum_index]][[signals_codes]]$Xdata=Xdata
        useful_data[[spectrum_index]][[signals_codes]]$Ydata=Ydata
        useful_data[[spectrum_index]][[signals_codes]]$results_to_save=results_to_save
        useful_data[[spectrum_index]][[signals_codes]]$error1=results_to_save$fitting_error

        #If the quantification is through fitting with or without baseline
      } else if (fitting_type == "Clean Fitting" || fitting_type ==
          "Baseline Fitting") {


        #Adaptation of the info of the parameters into a single matrix and preparation (if necessary) of the background signals that will conform the baseline
        FeaturesMatrix = fitting_prep(Xdata,
          Ydata,
          ROI_profile[, 5:11,drop=F],
          program_parameters,baselinedataset[spectrum_index,ROI_buckets])

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
          useful_data[[spectrum_index]][[signals_codes[i]]]$ROI_profile=ROI_profile
          useful_data[[spectrum_index]][[signals_codes[i]]]$program_parameters=program_parameters
          useful_data[[spectrum_index]][[signals_codes[i]]]$plot_data=plot_data
          useful_data[[spectrum_index]][[signals_codes[i]]]$error1=error1
          useful_data[[spectrum_index]][[signals_codes[i]]]$FeaturesMatrix=FeaturesMatrix
          useful_data[[spectrum_index]][[signals_codes[i]]]$signals_parameters=signals_parameters
          useful_data[[spectrum_index]][[signals_codes[i]]]$Xdata=Xdata
          useful_data[[spectrum_index]][[signals_codes[i]]]$Ydata=Ydata
          useful_data[[spectrum_index]][[signals_codes[i]]]$results_to_save=results_to_save
          }
     }

      #Generation of output variables specific of every quantification
      final_output = save_output(spectrum_index,signals_codes,results_to_save,imported_data$buck_step,final_output)

      setTxtProgressBar(pb, spectrum_index)
      }

  }


  if (optimization==TRUE) {
  optimized_profiling_data=autorun_improv(imported_data,final_output,useful_data,ROI_data)
  nn=optimized_profiling_data$final_output$fitting_error-final_output$fitting_error
  no=boxplot.stats(nn)$stats[5]
  ind=which(nn>no)
  for (i in 1:length(optimized_profiling_data$final_output)) {
    optimized_profiling_data$final_output[[i]][ind]=final_output[[i]][ind]
  }
  for (i in 1:nrow(optimized_profiling_data$final_output$fitting_error)) {
    ind=which(nn[i,]>no)
    if (length(ind)>0) optimized_profiling_data$useful_data[[i]][ind]=useful_data[[i]][ind]
  }
  final_output=optimized_profiling_data$final_output
  useful_data=optimized_profiling_data$useful_data
  }

  profiling_data=list(final_output=final_output,useful_data=useful_data)

  print("Done!")

  return(profiling_data)
}
