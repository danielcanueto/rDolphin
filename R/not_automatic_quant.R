
#' Quantification of individual ROIs with selected spectra.
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param ind Experiment to quantify.
#' @param ROI_profile Information of signals to fit in ROI
#' @param useful_data List with necessary information to load quantifications on the Shiny GUI.
#' @param interface Is the function being used with the Shiny GUI or not? By default F.
#'
#' @return Output depends on if the function is being used or not in the interface. If not in the interface, list with updated final_output and useful_data variables. If in the interface, necessary parameters to evaluate quality of the quantification before confiramtion by the user.
#' @export not_automatic_quant
#' @import baseline
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' resulting_data=not_automatic_quant(imported_data,imported_data$final_output,c(1,4),imported_data$ROI_data[1:2,],imported_data$useful_data)


not_automatic_quant = function(imported_data, final_output,ind,ROI_profile,useful_data,interface=F) {

  resulting_data=list(final_output=final_output,useful_data=useful_data)

  if (identical(ind,seq(nrow(imported_data$dataset)))) pb <- txtProgressBar(1, length(ind), style=3)

  ROI_buckets = which.min(abs(as.numeric(ROI_profile[1, 1])-imported_data$ppm)):which.min(abs(as.numeric(ROI_profile[1, 2])-imported_data$ppm))
  Xdata= as.numeric(imported_data$ppm[ROI_buckets])
  program_parameters=imported_data$program_parameters
  program_parameters$freq = imported_data$freq
  program_parameters$ROI_buckets = ROI_buckets
  program_parameters$buck_step = imported_data$buck_step
  fitting_type = as.character(ROI_profile[1, 3])
  if (length(grep("Clean",fitting_type))==1) {
    program_parameters$clean_fit="Y"
  } else {
    program_parameters$clean_fit="N"
  }
  signals_to_quantify = which(ROI_profile[, 5] >0)
  signals_codes = signals_names = rep(NA,length(signals_to_quantify))
  j = 1
  for (i in signals_to_quantify) {
    k = which(imported_data$signals_names == paste(ROI_profile[i,
      4],ROI_profile[i,5],sep='_'))

    signals_codes[j] = imported_data$signals_codes[k]
    signals_names[j] = as.character(imported_data$signals_names[k])
    j = j + 1
  }

  for (spectrum_index in ind) {
    dummy=imported_data$dataset[spectrum_index,(ROI_buckets[1]-5):(ROI_buckets[length(ROI_buckets)]+5)]
    baseline=baseline.rollingBall(rbind(dummy,dummy),5,5)$baseline[1,]
    baseline=baseline[6:(5+length(ROI_buckets))]
    # print(paste("Spectrum ",spectrum_index))

    Ydata = as.numeric(imported_data$dataset[spectrum_index, ROI_buckets])

    experiment_name = imported_data$Experiments[[spectrum_index]]

    # If the quantification is through integration with or without baseline
    if (fitting_type == "Clean Sum" ||
        fitting_type == "Baseline Sum") {
      # program_parameters$clean_fit = ifelse(fitting_type == "Clean Sum", "Y",
      #                                       "N")
      # program_parameters$freq=imported_data$freq
      # baseline_int = fitting_prep_integration(Xdata,Ydata,program_parameters,baseline)
      # Ydatamedian=as.numeric(apply(imported_data$dataset[, ROI_buckets,drop=F],2,median))

     dummy = integration(program_parameters$clean_fit, Xdata,Ydata,interface='T')

      results_to_save=dummy$results_to_save
      p=dummy$p
      plot_data=dummy$plot_data

      # resulting_data$integration_parameters=integration_parameters
      #Generation of output variables specific of every quantification
	if (identical(ind,seq(nrow(imported_data$dataset)))) {
        resulting_data$useful_data[[spectrum_index]][[signals_codes]]$ROI_profile=ROI_profile
        # resulting_data$useful_data[[spectrum_index]][[signals_codes]]$integration_parameters=integration_parameters
        resulting_data$useful_data[[spectrum_index]][[signals_codes]]$plot_data=dummy$plot_data
        resulting_data$useful_data[[spectrum_index]][[signals_codes]]$Xdata=Xdata
        resulting_data$useful_data[[spectrum_index]][[signals_codes]]$Ydata=Ydata
        resulting_data$useful_data[[spectrum_index]][[signals_codes]]$results_to_save=results_to_save
        resulting_data$useful_data[[spectrum_index]][[signals_codes]]$error1=results_to_save$fitting_error

        final_output = save_output(
          spectrum_index,
          signals_codes,
          results_to_save,
          imported_data$buck_step,
          final_output
        )

      }

      #If the quantification is through fitting with or without baseline
    } else if (fitting_type == "Clean Fitting" || fitting_type ==
               "Baseline Fitting") {

      program_parameters$freq=imported_data$freq

      FeaturesMatrix = fitting_prep(Xdata,
                                    Ydata,
                                    ROI_profile[, 5:11,drop=F],
                                    program_parameters,baseline)
      #Calculation of the parameters that will achieve the best fitting
      dummy = fittingloop(FeaturesMatrix,
                                       Xdata,
                                       Ydata,
                                       program_parameters)
      signals_parameters=dummy$signals_parameters
      multiplicities=c(FeaturesMatrix[,11],rep(1,(length(signals_parameters)/5)-dim(FeaturesMatrix)[1]))
      roof_effect=c(FeaturesMatrix[,12],rep(0,(length(signals_parameters)/5)-dim(FeaturesMatrix)[1]))

      signals_parameters_2=signals_parameters
      multiplicities_2=multiplicities
      roof_effect_2=roof_effect
      #Fitting of the signals
      dim(signals_parameters) = c(5, length(signals_parameters)/5)
      rownames(signals_parameters) = c(
        'intensity',
        'shift',
        'half_band_width',
        'gaussian',
        'J_coupling'
      )

      Xdata_2=imported_data$ppm
	   program_parameters$signals_to_quantify=signals_to_quantify
      Ydata_2 = as.numeric(imported_data$dataset[spectrum_index, ])
      # signals_parameters_2=unlist(signals_parameters_2)
      # multiplicities_2=unlist(multiplicities_2)
      # roof_effect_2=unlist(roof_effect_2)
      # signals_parameters_2=unlist(signals_parameters)
      multiplicities_2=unlist(multiplicities)
      roof_effect_2=unlist(roof_effect)

      fitted_signals = signal_fitting(signals_parameters_2,
                                         Xdata_2,multiplicities_2,roof_effect_2,program_parameters$freq)
      # signals_parameters=as.matrix(signals_parameters)

      # dim(signals_parameters_2) = c(5, length(signals_parameters_2)/5)
      # rownames(signals_parameters_2) = c(
      #   'intensity',
      #   'shift',
      #   'half_band_width',
      #   'gaussian',
      #   'J_coupling'
      # )


      #Generation of output data about the fitting and of the necessary variables for the generation ofa figure
      dummy = output_generator(
        signals_to_quantify,
        fitted_signals,
        Ydata_2,
        Xdata_2,
        signals_parameters,multiplicities,ROI_buckets
      )
      output_data=dummy$output_data
      error1=dummy$error1

      if (any(output_data$fitting_error>0.05)==T) {
        dummy = fittingloop(FeaturesMatrix,
          Xdata,
          Ydata,
          program_parameters)
        signals_parameters=dummy$signals_parameters
        multiplicities=c(FeaturesMatrix[,11],rep(1,(length(signals_parameters)/5)-dim(FeaturesMatrix)[1]))
        roof_effect=c(FeaturesMatrix[,12],rep(0,(length(signals_parameters)/5)-dim(FeaturesMatrix)[1]))

        signals_parameters_2=signals_parameters
        # multiplicities_2=multiplicities
        # roof_effect_2=roof_effect
        #Fitting of the signals
        dim(signals_parameters) = c(5, length(signals_parameters)/5)
        rownames(signals_parameters) = c(
          'intensity',
          'shift',
          'half_band_width',
          'gaussian',
          'J_coupling'
        )

        # signals_parameters_2=unlist(signals_parameters_2)
        # multiplicities_2=unlist(multiplicities_2)
        # roof_effect_2=unlist(roof_effect_2)
        # signals_parameters_2=unlist(signals_parameters)
        multiplicities_2=unlist(multiplicities)
        roof_effect_2=unlist(roof_effect)

        fitted_signals = signal_fitting(signals_parameters_2,
          Xdata_2,multiplicities_2,roof_effect_2,program_parameters$freq)
        # signals_parameters=as.matrix(signals_parameters)

        # dim(signals_parameters_2) = c(5, length(signals_parameters_2)/5)
        # rownames(signals_parameters_2) = c(
        #   'intensity',
        #   'shift',
        #   'half_band_width',
        #   'gaussian',
        #   'J_coupling'
        # )

        #Generation of output data about the fitting and of the necessary variables for the generation ofa figure
        dummy = output_generator(
          signals_to_quantify,
          fitted_signals,
          Ydata_2,
          Xdata_2,
          signals_parameters,multiplicities,ROI_buckets
        )

        #If new deconvolution has improved previous one
        if (dummy$error1<error1) {
          output_data=dummy$output_data
          error1=dummy$error1
        }}

      #Generation of the dataframe with the final output variables
      results_to_save = data.frame(
        shift = output_data$shift,
        quantification = output_data$quantification,
        signal_area_ratio = output_data$signal_area_ratio,
        fitting_error = output_data$fitting_error,
        intensity = output_data$intensity,
        half_band_width = output_data$half_band_width
      )

      plot_data = rbind(
        output_data$signals_sum,
        output_data$baseline_sum,
        output_data$fitted_sum,
        output_data$signals
      )
      rownames(plot_data) = c("signals_sum",
        "baseline_sum",
        "fitted_sum",
        as.character(paste(ROI_profile[,4],ROI_profile[,5],sep='_')),rep('additional signal',dim(plot_data)[1]-length(ROI_profile[,4])-3))

      plotdata2 = data.frame(Xdata=Xdata_2,
        Ydata=Ydata_2,
        plot_data[3, ],
        plot_data[2, ] )
      plotdata3 <- melt(plotdata2, id = "Xdata")
      plotdata3$variable = c(
        rep('Original Spectrum', length(Ydata_2)),
        rep('Generated Spectrum', length(Ydata_2)),
        rep('Generated Background', length(Ydata_2))
      )
      plot_title = paste(imported_data$Experiments[spectrum_index],"- ROI ",ROI_profile[1,1],"-",ROI_profile[1,2],"ppm")
colors=c(I('red'),I('blue'),I('black'),I('brown'),I('cyan'),I('green'),I('yellow'))
      p=plot_ly(plotdata3,x=~Xdata,y=~value,color=~variable,type='scatter',mode='lines',fill=NULL) %>% layout(title = plot_title,xaxis = list(range=c(Xdata[1],Xdata[length(Xdata)]),title = 'ppm'), yaxis = list(range=c(0,max(Ydata)),title = 'Intensity'))
        for (i in 4:nrow(plot_data)) {
          plotdata5 =  data.frame(Xdata=Xdata_2, variable=rownames(plot_data)[i] ,value=plot_data[i,])

        p=p %>%add_trace(data=plotdata5,x=~Xdata,y=~value,name=~variable,type='scatter',mode='lines',fill='tozeroy',fillcolor=colors[i-3])
}


    signals_parameters=rbind(signals_parameters,multiplicities,roof_effect)
    if (fitting_type == "Clean Fitting") {
      colnames(signals_parameters)=paste(ROI_profile[,4],ROI_profile[,5],sep='_')
    } else {
      colnames(signals_parameters)=c(paste(ROI_profile[,4],ROI_profile[,5],sep='_'),paste('baseline_signal',seq(ncol(signals_parameters)-nrow(ROI_profile)),sep='_'))
    }

    # if (resulting_data$useful_data[[spectrum_index]][[signals_codes[1]]]$error1>0.8*error1) {

    # }
if (identical(ind,seq(nrow(imported_data$dataset)))| interface ==F)  {
	# if (resulting_data$useful_data[[spectrum_index]][[signals_codes[1]]]$error1>error1) {
      for (i in seq_along(signals_codes)) {
        resulting_data$useful_data[[spectrum_index]][[signals_codes[i]]]$ROI_profile=ROI_profile
        resulting_data$useful_data[[spectrum_index]][[signals_codes[i]]]$program_parameters=program_parameters
        resulting_data$useful_data[[spectrum_index]][[signals_codes[i]]]$plot_data=plot_data[,ROI_buckets]
        resulting_data$useful_data[[spectrum_index]][[signals_codes[i]]]$FeaturesMatrix=FeaturesMatrix
        resulting_data$useful_data[[spectrum_index]][[signals_codes[i]]]$signals_parameters=signals_parameters
        resulting_data$useful_data[[spectrum_index]][[signals_codes[i]]]$error1=error1
        resulting_data$useful_data[[spectrum_index]][[signals_codes[i]]]$Xdata=Xdata
        resulting_data$useful_data[[spectrum_index]][[signals_codes[i]]]$Ydata=Ydata
        resulting_data$useful_data[[spectrum_index]][[signals_codes[i]]]$results_to_save=results_to_save


      }
        resulting_data$final_output = save_output(
          spectrum_index,
          signals_codes,
          results_to_save,
          imported_data$buck_step,
          resulting_data$final_output)
      } else {
	  resulting_data$program_parameters=program_parameters
    resulting_data$results_to_save=results_to_save
    resulting_data$ROI_profile=ROI_profile
    resulting_data$Ydata=Ydata
    resulting_data$plot_data=plot_data[,ROI_buckets]
    resulting_data$FeaturesMatrix=FeaturesMatrix
    resulting_data$error1=error1
    resulting_data$signals_parameters=signals_parameters
    resulting_data$Xdata=Xdata
	}

    }

	    if (identical(ind,seq(nrow(imported_data$dataset))))  setTxtProgressBar(pb, spectrum_index)
    # if (fitting_type == "Clean Sum" || fitting_type == "Baseline Sum"&&identical(ind,seq(nrow(imported_data$dataset)))) {
    #
    #   dummy=integration_error(ROI_data,useful_data,final_output,signals_codes)
    #   resulting_data$useful_data=dummy$useful_data
    #   resulting_data$final_output=dummy$final_output

    }
	if (interface == T) {
		resulting_data$p=p
		resulting_data$results_to_save=results_to_save
		resulting_data$spectrum_index=spectrum_index
		resulting_data$signals_codes=signals_codes
		resulting_data$fitting_type=fitting_type
		# resulting_data$signals_names=signals_names
	}
  # }
  return(resulting_data)
}
