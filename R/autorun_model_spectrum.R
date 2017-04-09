#' Automatic quantification in the model spectrum of the dataset using the information located in the ROI patterns file, with associated p values for every bin. The reulting information can be useful to analyze possible new ROIs to profile or modification of existing ones.
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#'
#' @return plotly figure with performed fitting for every ROI and p value for every bin
#' @export autorun_model_spectrum
#' @import baseline
#' @import minpack.lm
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' model_spectrum_plot=autorun_model_spectrum(imported_data)



autorun_model_spectrum = function(imported_data) {

  print('Preparing quantifications in a model spectrum with the current ROI Profiles. A figure with the performed quantifications will be shown, as well as a chemometric model with the metadata given. Then you can go to change the parameters of the ROI Profiles')


  baselinedataset=baseline.rollingBall(imported_data$dataset,5,5)$baseline
  # indicators=matrix(NA,nrow(imported_data$ROI_data),2,dimnames=list(imported_data$signals_names)))
  total_signals_parameters=matrix(NA,nrow(imported_data$ROI_data),9,dimnames=list(imported_data$signals_names))
  colnames(total_signals_parameters)=c("intensity",	" chemical shift",	"half_band_width",	"gaussian %",	"J coupling",	"multiplicities",	"roof_effect","fitting error","signal / total area ratio")

  quartile_spectrum = as.numeric(apply(imported_data$dataset, 2, function(x)
    quantile(x, 0.75,na.rm=T)))
  spectrum_index = which.min(apply(imported_data$dataset, 1, function(x)
    sqrt(mean((x - quartile_spectrum) ^ 2
      ,na.rm=T))))
  plotdata = data.frame(Xdata=as.numeric(imported_data$ppm),Ydata = as.numeric(imported_data$dataset[spectrum_index,]))
  fitted_data=rep(0,length(imported_data$ppm))
  pb   <- txtProgressBar(1, nrow(imported_data$ROI_separator), style=3)

  for (ROI_index in seq_along(imported_data$ROI_separator[, 1])) {
    #Loading of every ROI parameters
    ROI_profile = imported_data$ROI_data[imported_data$ROI_separator[ROI_index, 1]:imported_data$ROI_separator[ROI_index, 2],]
    ROI_buckets = which.min(abs(as.numeric(ROI_profile[1, 1])-imported_data$ppm)):which.min(abs(as.numeric(ROI_profile[1, 2])-imported_data$ppm))
    signals_to_quantify = which(ROI_profile[, 5] >= 1)
    signals_codes = (imported_data$ROI_separator[ROI_index, 1]:imported_data$ROI_separator[ROI_index, 2])[signals_to_quantify]

    #Preparation of necessary parameters
    program_parameters=imported_data$program_parameters
    program_parameters$freq = imported_data$freq
    program_parameters$ROI_buckets = ROI_buckets
    program_parameters$buck_step = imported_data$buck_step


    fitting_type = as.character(ROI_profile[1, 3])
    signals_to_quantify = which(ROI_profile[, 5] >= 1)
    ROI_buckets = which.min(abs(as.numeric(ROI_profile[1, 1])-imported_data$ppm)):which.min(abs(as.numeric(ROI_profile[1, 2])-imported_data$ppm))
    if (length(ROI_buckets)<5) next
    if (ROI_buckets[1]>ROI_buckets[2]) ROI_buckets=rev(ROI_buckets)



    Xdata= as.numeric(imported_data$ppm[ROI_buckets])
    Ydata = as.numeric(imported_data$dataset[spectrum_index, ROI_buckets])



    # If the quantification is through integration with or without baseline
    if (fitting_type == "Clean Sum" ||
        fitting_type == "Baseline Sum") {
      baseline = ifelse(fitting_type == "Clean Sum", rep(0,length(Xdata)),seq(min(Ydata[1:5]), min(Ydata[(length(Xdata) - 4):length(Xdata)]), len=length(Xdata)))

      Ydatamedian=as.numeric(apply(imported_data$dataset[, ROI_buckets,drop=F],2,median))

      clean_fit = ifelse(fitting_type == "Clean Sum", "Y", "N")
      integration_variables = integration(clean_fit, Xdata,Ydata,Ydatamedian)

      total_signals_parameters[signals_codes,]=c(integration_variables$results_to_save$intensity,integration_variables$results_to_save$shift,rep(NA,5),integration_variables$results_to_save$fitting_error,integration_variables$results_to_save$signal_area_ratio)

      #integration ad chechk that there are no negative values
      # integrated_signal = Ydata - baseline
      # integrated_signal[integrated_signal<0]=0
      #preparation of output
      fitted_data[ROI_buckets]= integration_variables$plot_data[3,]


    } else if (fitting_type == "Clean Fitting" || fitting_type ==
        "Baseline Fitting") {
      is_roi_testing = "N"
      program_parameters$clean_fit='N'

      initial_fit_parameters = ROI_profile[, 5:11,drop=F]

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

      #Fitting of the signals
      multiplicities=FeaturesMatrix[,11]
      roof_effect=FeaturesMatrix[,12]
      # dummy = ROI_data[which(is.na(ROI_data[, 1])),]
      # dummy2= list()
      # for (i in 1:dim(ROI_profile)[1]) {
      #
      #   if (length(which(dummy[,4] == ROI_profile[i,4]))>0 & ROI_profile[i,5] == 1) {
      #     dummy2[[length(dummy2)+1]]=which(dummy[,4] == ROI_profile[i,4])
      #     for (j in 1:length(dummy2[[i]])) {
      #
      #       cc= signals_parameters[(5*i-4):(5*i)]
      #       cc[5]=dummy[dummy2[[i]][j],][9]
      #       cc[1]=dummy[dummy2[[i]][j],][12]*cc[1]
      #       cc[2]=as.numeric(dummy[dummy2[[i]][j],][5])+(as.numeric(cc[2])-as.numeric(ROI_profile[i,6]))
      #       # cc[5]=dummy[dummy2[[i]][j],][9]
      #       # cc[1]=dummy[dummy2[[i]][j],][12]*cc[1]
      #       # cc[2]=dummy[dummy2[[i]][j],][5]-(cc[2]-ROI_profile[i,5])
      #       signals_parameters=c(signals_parameters,cc)
      #       multiplicities=c(multiplicities,dummy[dummy2[[i]][j],][8])
      #       roof_effect=c(roof_effect,dummy[dummy2[[i]][j],][10])
      #     }
      #   }
      # }
      # Xdata=imported_data$ppm
      # multiplicities=unlist(multiplicities)
      # roof_effect=unlist(roof_effect)

      fitted_signals = signal_fitting(signals_parameters,
        Xdata,multiplicities,roof_effect,Ydata,program_parameters$freq)
      dim(signals_parameters) = c(5, length(signals_parameters)/5)
      rownames(signals_parameters) = c(
        'intensity',
        'shift',
        'half_band_width',
        'gaussian',
        'J_coupling'
      )

      # Ydata = as.numeric(imported_data$dataset[spectrum_index, ])
      signals_parameters=rbind(signals_parameters,multiplicities,roof_effect)



      dummy = output_generator(
        signals_to_quantify,
        fitted_signals,
        Ydata,
        Xdata,
        signals_parameters,multiplicities
      )
      output_data=dummy$output_data
      # results_to_save = data.frame(
      #   shift = output_data$shift,
      #   Area = output_data$Area,
      #   signal_area_ratio = output_data$signal_area_ratio,
      #   fitting_error = output_data$fitting_error,
      #   intensity = output_data$intensity,
      #   half_band_width = output_data$half_band_width
      # )

      # plotdata2$Ydata[ROI_buckets]=output_data$signals_sum[ROI_buckets]
      # plotdata3$Ydata[ROI_buckets]=output_data$fitted_sum[ROI_buckets]
      fitted_data[ROI_buckets]=output_data$fitted_sum


    }
    setTxtProgressBar(pb, ROI_index)
    # indicators[signals_codes,]=unlist(results_to_save)
    total_signals_parameters[signals_codes,]=cbind(t(signals_parameters[,signals_to_quantify]),output_data$fitting_error[signals_to_quantify],output_data$signal_area_ratio[signals_to_quantify])

  }
  p_value_bucketing=as.vector(p_values(imported_data$dataset,imported_data$Metadata))

  ay <- list(tickfont = list(color = "red"),overlaying = "y",side = "right",title = "p value",range = c(0, max(as.numeric(imported_data$dataset[spectrum_index, ]))))
  az = list(title = "Intensity",range = c(-1, max(as.numeric(imported_data$dataset[spectrum_index, ]))-1))
  p=plot_ly()%>%
    add_lines(x=~imported_data$ppm,y = ~as.numeric(imported_data$dataset[spectrum_index, ]),,name='Model spectrum')%>%
    add_lines(x=~imported_data$ppm,y = ~fitted_data,name='Fitted spectrum',fill='tozeroy')%>%
  add_lines(x=~imported_data$ppm,y = ~p_value_bucketing,name='p value', yaxis = "y2")%>%
    layout(xaxis=list(title='ppm',range=c(max(imported_data$ppm),min(imported_data$ppm))),yaxis=az, yaxis2 = ay)
  dummy=list(p=p,total_signals_parameters=round(total_signals_parameters,4),ROI_data=imported_data$ROI_data)
  return(dummy)
}
