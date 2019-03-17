
  profiling_func=function(spectrum_index,signals_codes,imported_data,ROI_buckets,
                          fitting_type,
                          program_parameters,Xdata,Ydata,final_output,
                          reproducibility_data,
                          ROI_profile,baselinedataset,signals_to_quantify,pb=NULL,reimplementation=F,
                         max_shift=NULL,min_shift=NULL,
                         max_intensity=NULL,
                         min_intensity=NULL,
                         max_width=NULL,
                         min_width=NULL,
                         signal_index=NULL) {
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
    if (reimplementation==F) FeaturesMatrix = fitting_prep(Xdata,
                                  Ydata,
                                  ROI_profile[, 5:11,drop=F],
                                  program_parameters,baselinedataset[spectrum_index,ROI_buckets])
								  
	if (reimplementation==T)  FeaturesMatrix = fitting_prep_2(Xdata,
          Ydata,
          ROI_profile[, 5:11,drop=F],
          program_parameters,baselinedataset[spectrum_index,ROI_buckets],max_shift,min_shift,max_intensity,
          min_intensity,max_width,min_width,spectrum_index,signal_index)

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

  if (!is.null(pb)) setTxtProgressBar(pb, spectrum_index)

  profiling_data=list(final_output=final_output,reproducibility_data=reproducibility_data)
  return(profiling_data)
  }

