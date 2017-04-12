#' Helper function to import data in GUI. Not to be used in console.
#'
#' @param datapath
#' @param dummy 
#'
#' @return dummy
#' @export helperimport



	helperimport = function(datapath,dummy) {


	dummy$imported_data = tryCatch({suppressWarnings(import_data(datapath))}, error = function(e) {
	stop('Import of data not possible. Please revise the introduced parameters.')
	# return(dummy)
	})
    # reactiveprogramdata$originaldataset=imported_data$dataset
	    #plot of quantification in model spectrum with current roi profiles
	dummy2=tryCatch({autorun_model_spectrum(dummy$imported_data)}, error = function(e) {
	print('Automatic quantification of model spectrum not possible.')
	return(dummy)
	})
	dummy$autorun_plot=dummy2$p
	dummy$total_signals_parameters=dummy2$total_signals_parameters
	# dummy$indicators=dummy2$indicators


	#plots of representative spectra and median spectra per group to help setting the right ROI parameters

    dummy$clusterplot=tryCatch({clustspectraplot(dummy$imported_data)  }, error = function(e) {
	print('Generation of subsets or representative spectra not possible.')
	return(dummy)
	})
    dummy$medianplot=tryCatch({medianplot(dummy$imported_data)
	  }, error = function(e) {
	print('Generation of median spectra not possible.')
	return(dummy)
	})

    #Subsetting of ROIs is prepared
    #Names of ROIS and cluster and median spectra are prepared
	roi_variables=tryCatch({roifunc(dummy$imported_data$ROI_data,dummy$imported_data$Metadata,dummy$imported_data$Experiments)
      }, error = function(e) {
	print('Generation of ROI data not possible.')
	return(dummy)
	})
	dummy$select_options=roi_variables$select_options
	dummy$spectra=roi_variables$spectra
	dummy$beginning=T
	dummy$jres_plot=tryCatch(twod_data(dummy$imported_data$jres_path), error = function(e) NA)
	return(dummy)
	}
