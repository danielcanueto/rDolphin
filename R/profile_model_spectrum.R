#' Automatic quantification in the model spectrum of the dataset using the information located in the ROI patterns file, with associated p values for every bin. The reulting information can be useful to analyze possible new ROIs to profile or modification of existing ones.
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param ROI_data ROIs data.

#'
#' @return plotly figure with performed fitting for every ROI and p value for every bin
#' @export profile_model_spectrum
#' @import baseline
#' @import minpack.lm
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' model_spectrum_plot=profile_model_spectrum(imported_data,imported_data$ROI_data)



profile_model_spectrum = function(imported_data, ROI_data,spectrum_index=NA) {

  print('Preparing quantifications in a model spectrum with the current ROI Profiles. A figure with the performed quantifications will be shown, as well as a chemometric model with the metadata given. Then you can go to change the parameters of the ROI Profiles')

#Splitting of ROI data into individual ROIs to be quantified
	dummy = which(is.na(ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(ROI_data)[1]+1
    lal=which(duplicated(ROI_data[-dummy,1:2])==F)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(ROI_data[-dummy,])[1]))

  # indicators=matrix(NA,nrow(ROI_data),2,dimnames=list(imported_data$signals_names)))
  total_signals_parameters=matrix(NA,nrow(ROI_data),9,dimnames=list(imported_data$signals_names))
  colnames(total_signals_parameters)=c("intensity",	" chemical shift",	"half_band_width",	"gaussian %",	"J coupling",	"multiplicities",	"roof_effect","fitting error","signal / total area ratio")
if (is.na(spectrum_index)) {
  quartile_spectrum = as.numeric(apply(imported_data$dataset, 2, function(x)
    quantile(x, 0.75,na.rm=T)))
  spectrum_index = which.min(apply(imported_data$dataset, 1, function(x)
    sqrt(mean((x - quartile_spectrum) ^ 2
      ,na.rm=T))))
}
  baseline=baseline::baseline.rollingBall(rbind(imported_data$dataset[spectrum_index,],imported_data$dataset[spectrum_index,]),5,5)$baseline[1,]

  plotdata = data.frame(Xdata=as.numeric(imported_data$ppm),Ydata = as.numeric(imported_data$dataset[spectrum_index,]))
  fitted_data=rep(0,length(imported_data$ppm))
  pb   <- txtProgressBar(1, nrow(ROI_separator), style=3)

  for (ROI_index in seq_along(ROI_separator[, 1])) {
    #Loading of every ROI parameters
    ROI_profile = ROI_data[ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2],]
    ROI_buckets = which.min(abs(as.numeric(ROI_profile[1, 1])-imported_data$ppm)):which.min(abs(as.numeric(ROI_profile[1, 2])-imported_data$ppm))
    signals_to_quantify = which(ROI_profile[, 5] >= 1)
    signals_codes = (ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2])

    #Preparation of necessary parameters
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
    if (length(ROI_buckets)<5) next
    if (ROI_buckets[1]>ROI_buckets[2]) ROI_buckets=rev(ROI_buckets)



    Xdata= as.numeric(imported_data$ppm[ROI_buckets])
    Ydata = as.numeric(imported_data$dataset[spectrum_index, ROI_buckets])



    # If the quantification is through integration with or without baseline
    if (fitting_type == "Clean Sum" ||
        fitting_type == "Baseline Sum") {

      integration_variables = integration(program_parameters$clean_fit, Xdata,Ydata,program_parameters$buck_step)

      total_signals_parameters[signals_codes,]=c(integration_variables$results_to_save$intensity,integration_variables$results_to_save$shift,rep(NA,5),integration_variables$results_to_save$fitting_error,integration_variables$results_to_save$signal_area_ratio)


      #preparation of output
      fitted_data[ROI_buckets]= integration_variables$plot_data[3,]


    } else if (fitting_type == "Clean Fitting" || fitting_type ==
        "Baseline Fitting") {


      #Adaptation of the info of the parameters into a single matrix and preparation (if necessary) of the background signals that will conform the baseline
      FeaturesMatrix = fitting_prep(Xdata,
        Ydata,
        ROI_profile[, 5:11,drop=F],
        program_parameters,baseline[ROI_buckets])

      #Calculation of the parameters that will achieve the best fitting
      dummy = fittingloop(FeaturesMatrix,
        Xdata,
        Ydata,
        program_parameters)
      signals_parameters=dummy$signals_parameters


      #Fitting of the signals
      multiplicities=FeaturesMatrix[,11]
      roof_effect=FeaturesMatrix[,12]


      fitted_signals = signal_fitting(signals_parameters,
        Xdata,multiplicities,roof_effect,program_parameters$freq)
      dim(signals_parameters) = c(5, length(signals_parameters)/5)
      rownames(signals_parameters) = c(
        'intensity',
        'shift',
        'half_band_width',
        'gaussian',
        'J_coupling'
      )

      signals_parameters=rbind(signals_parameters,multiplicities,roof_effect)


      dummy = output_generator(
        signals_to_quantify,
        fitted_signals,
        Ydata,
        Xdata,
        signals_parameters,multiplicities,program_parameters$buck_step
      )
      output_data=dummy$output_data

      fitted_data[ROI_buckets]=output_data$fitted_sum
    total_signals_parameters[signals_codes,]=cbind(t(signals_parameters[,seq(nrow(ROI_profile))]),output_data$fitting_error,output_data$signal_area_ratio)


    }
    setTxtProgressBar(pb, ROI_index)

  }
  p_value_bucketing=as.vector(p_values(imported_data$dataset,imported_data$Metadata))

  ay <- list(tickfont = list(color = "red"),overlaying = "y",side = "right",title = "p value",range = c(0, max(as.numeric(imported_data$dataset[spectrum_index, ]))))
  az = list(title = "Intensity (arbitrary unit)",range = c(-1, max(as.numeric(imported_data$dataset[spectrum_index, ]))-1))
  p=plot_ly()%>%
    add_lines(x=~imported_data$ppm,y = ~as.numeric(imported_data$dataset[spectrum_index, ]),name='Model spectrum')%>%
    add_lines(x=~imported_data$ppm,y = ~fitted_data,name='Fitted spectrum',fill='tozeroy')%>%
  add_lines(x=~imported_data$ppm,y = ~p_value_bucketing,name='p value', yaxis = "y2",line = list(color = 'rgba(255, 0, 0, 1)'))%>%
    layout(xaxis=list(title='ppm',range=c(max(imported_data$ppm),min(imported_data$ppm))),yaxis=az, yaxis2 = ay)
  dummy=list(p=p,total_signals_parameters=round(total_signals_parameters,4),ROI_data=ROI_data)
  return(dummy)
  print("DOne!")
}
