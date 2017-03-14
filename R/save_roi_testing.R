
save_roi_testing=function(provisional_data,imported_data,finaloutput,useful_data) {


if (provisional_data$fitting_type == "Clean Sum" ||
    provisional_data$fitting_type == "Baseline Sum") {

  provisional_data$useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes]]$ROI_profile=ROI_profile
  provisional_data$useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes]]$plot_data=dummy$plot_data
  provisional_data$useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes]]$Xdata=Xdata
  provisional_data$useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes]]$Ydata=Ydata
  provisional_data$useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes]]$results_to_save=results_to_save
  provisional_data$useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes]]$error1=results_to_save$fitting_error



} else if (provisional_data$fitting_type == "Clean Fitting" || provisional_data$fitting_type ==
    "Baseline Fitting") {



  provisional_data$program_parameters$signals_to_quantify=NULL

  for (i in seq_along(provisional_data$signals_codes)) {
    useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$ROI_profile=provisional_data$ROI_profile
    useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$program_parameters=provisional_data$program_parameters
    # useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$fitted_signals=provisional_data$fitted_signals
    useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$plot_data=provisional_data$plot_data
    useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$Xdata=provisional_data$Xdata
    useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$Ydata=provisional_data$Ydata
    useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$results_to_save=provisional_data$results_to_save
    useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$signals_parameters=provisional_data$signals_parameters
    useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$FeaturesMatrix=provisional_data$FeaturesMatrix
    useful_data[[provisional_data$spectrum_index]][[provisional_data$signals_codes[i]]]$error1=provisional_data$error1


  }
  finaloutput = save_output(
    provisional_data$spectrum_index,
    provisional_data$signals_codes,
    provisional_data$results_to_save,
    imported_data$buck_step,
    finaloutput)

}



  # tryCatch({write_info(imported_data$export_path, finaloutput)}, error = function(err) {
  #   print('Not possible to overwrite a csv file open with Microsoft Excel')
  # })

    dummy=list(finaloutput=finaloutput,useful_data=useful_data)
    return(dummy)

  }
