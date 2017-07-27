
#' Automatic quantification of signals for all experiments using the information located in the ROI patterns file.
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param useful_data List with necessary information to load quantifications on the Shiny GUI.
#' @param ROI_data ROIs data.
#'
#' @return List with updated final_output and useful_data variables.
#' @export autorun
#' @import baseline
#' @import robustbase

#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' # Not run:
#' # quantification_variables=autorun(imported_data,imported_data$final_output,imported_data$useful_data,imported_data$ROI_data)


autorun_shift = function(imported_data, final_output,useful_data,ROI_data) {

  print('Be patient. Gonna take a while. You should be writing, meanwhile.')

    final_output$shift[final_output$shift==Inf]=NA
  original_shift=final_output$shift
  adjusted_shift_indicator=apply(original_shift,2,function(x)sort(table(x),decreasing=TRUE)[1])
  original_shift[,adjusted_shift_indicator>0.5*nrow(original_shift)]=NA
  analyzed_signals=apply(original_shift,2,function(x)! all(is.na(x)))
  original_shift[,analyzed_signals]=missForest::missForest(original_shift[,analyzed_signals])$ximp
  subset_shift=original_shift
  predicted_shift=matrix(NA,nrow(original_shift),ncol(original_shift))
  possible_predictors=seq(ncol(original_shift))
  for (i in 1:ncol(subset_shift)) {
    if (all(is.na(subset_shift[,i]))) next
    r2_all_rlr=rep(NA,ncol(original_shift))
    for (j in possible_predictors) r2_all_rlr[j]=tryCatch({summary(robustbase::lmrob(subset_shift[,i]~original_shift[,j,drop=F], k.max = 2000))$adj.r.squared},error=function(e)NA)
    r2_all_rlr[i]=NA
    if (max(r2_all_rlr,na.rm=T)<0.2) next
    model=robustbase::lmrob(subset_shift[,i]~original_shift[,which(r2_all_rlr>0.2),drop=F], k.max = 2000)
    predicted_shift[,i]=predict(model)
  }

  outliers=matrix(NA,nrow(original_shift),ncol(original_shift))
  for (i in 1:ncol(subset_shift)) {
    lol=summary(robustbase::lmrob(subset_shift[,i]~predicted_shift[,i],max.it = 2000))$residuals
    outliers[which(lol %in% boxplot.stats(lol)$out),i]=1
  }


  for (i in 1:ncol(subset_shift)) {
    if (all(is.na(subset_shift[,i]))) next

    r2_all_rlr=rep(NA,ncol(original_shift))
    for (j in possible_predictors) r2_all_rlr[j]=tryCatch({summary(robustbase::lmrob(subset_shift[,i]~original_shift[,j,drop=F], k.max = 2000))$adj.r.squared},error=function(e)NA)
    r2_all_rlr[i]=NA

    if (max(r2_all_rlr,na.rm=T)<0.2) next
    las=which(is.na(outliers[,i])==F)
    las2=which(apply(outliers[las,,drop=F],2,function(x)all(is.na(x)))==T)
    int=intersect(which(r2_all_rlr>0.2),las2)
    if (length(int)==0) next
    predicted_shift[las,i]=predict(robustbase::lmrob(subset_shift[,i]~original_shift[,int,drop=F], k.max = 2000))[las]

  }

  matrnew=matrix(NA,nrow(original_shift),ncol(original_shift))
  for (i in 1:ncol(subset_shift)) {
    lol=summary(robustbase::lmrob(subset_shift[,i]~predicted_shift[,i],max.it = 2000))$residuals
    matrnew[which(lol %in% boxplot.stats(lol)$out),i]=1
  }

  for (i in 1:ncol(subset_shift)) {

    r2_all_rlr=rep(NA,ncol(original_shift))
    for (j in possible_predictors) r2_all_rlr[j]=tryCatch({summary(robustbase::lmrob(subset_shift[,i]~original_shift[,j,drop=F], k.max = 2000))$adj.r.squared},error=function(e)NA)
    r2_all_rlr[i]=NA

    if (max(r2_all_rlr,na.rm=T)<0.2) next
    las=which(is.na(matrnew[,i])==F)
    las2=which(apply(matrnew[las,,drop=F],2,function(x)all(is.na(x)))==T)
    int=intersect(which(r2_all_rlr>0.2),las2)
    if (length(int)==0) next
    predicted_shift[las,i]=predict(robustbase::lmrob(subset_shift[,i]~original_shift[,int,drop=F], k.max = 2000))[las]

    # provisional_predictors=provisional_predictors[duplicated(provisional_predictors)==F]
    # predictors_train[[i]]=provisional_predictors
  }

  ll=apply(predicted_shift-original_shift,2,function(x)boxplot.stats(x)$stats[5])
  ind=which(is.na(ll))
  ll[ind]=ROI_data[ind,7]
  if (length(ind)>0) predicted_shift[,ind]=t(replicate(nrow(original_shift),ROI_data[ind,6]))
  # ind=which(la>0.5*nrow(shift2))
  # if (length(ind)>0) {
  #   predicted_shift[,ind]=t(replicate(nrow(original_shift),ind2[ind]))
  #   ll[ind]=ROI_data[ind,7]
  # }


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
    for (spectrum_index in 1:nrow(imported_data$dataset)) {

      #Preparation of necessary variables to store figures and information of the fitting
      Ydata = as.numeric(imported_data$dataset[spectrum_index, ROI_buckets])

      #If the quantification is through integration with or without baseline
      if (fitting_type == "Clean Sum" ||
          fitting_type == "Baseline Sum") {
#
#         far=(ROI_profile[,1]- ROI_profile[,2])-ROI_profile[,7]*2
#         far=far/2+  ll[ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]]
#         ROI_profile[,1:2]=mean(as.numeric(ROI_profile[,1:2]))+c(far,-far)
#
#         ROI_profile[,6]=predicted_shift[spectrum_index,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]]
#         ROI_profile[,7]=ll[ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]]


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

        far=(ROI_profile[,1]- ROI_profile[,2])/2-(ROI_profile[,7]-max(ll[ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]],na.rm=T))
        ROI_profile[,1]=(ROI_profile[,1]- ROI_profile[,2])/2+far
        ROI_profile[,2]=(ROI_profile[,1]- ROI_profile[,2])/2-far
        ROI_profile[,6]=predicted_shift[spectrum_index,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]]
        ROI_profile[,7]=ll[ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]]
        limits_intensity=replicate(nrow(ROI_profile),c(0,max(Ydata)))
        widthtol= rep(program_parameters$widthtolerance,nrow(ROI_profile))

         #Adaptation of the info of the parameters into a single matrix and preparation (if necessary) of the background signals that will conform the baseline
        FeaturesMatrix = fitting_prep_4(Xdata,
          Ydata,
          ROI_profile[, 5:11,drop=F],
          program_parameters,baselinedataset[spectrum_index,ROI_buckets],limits_intensity,widthtol)

        #Calculation of the parameters that will achieve the best fitting
        dummy = fittingloop2(FeaturesMatrix,
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
        rownames(signals_parameters) = c('intensity','shift','half_band_width','gaussian','J_coupling')
        signals_parameters=rbind(signals_parameters,multiplicities,roof_effect)

        #Generation of output data about the fitting and of the necessary variables for the generation ofa figure
        dummy = output_generator(signals_to_quantify,fitted_signals,Ydata_2,Xdata_2,signals_parameters,multiplicities,program_parameters$buck_step)
        output_data=dummy$output_data
        error1=ifelse(is.nan(dummy$error1),3000,dummy$error1)

        # #If any of the qunatificaitons has more than 5% fitting error, try again the deconvolution
        # if (any(output_data$fitting_error>0.05)==T) {
        # dummy = fittingloop(FeaturesMatrix,Xdata,Ydata,program_parameters)
        # signals_parameters=dummy$signals_parameters

		# #Fitting of the signals
        # multiplicities=c(FeaturesMatrix[,11],rep(1,(length(signals_parameters)/5)-dim(FeaturesMatrix)[1]))
        # roof_effect=c(FeaturesMatrix[,12],rep(0,(length(signals_parameters)/5)-dim(FeaturesMatrix)[1]))
        # fitted_signals = signal_fitting(signals_parameters,
          # Xdata_2,multiplicities,roof_effect,program_parameters$freq)
        # dim(signals_parameters) = c(5, length(signals_parameters)/5)
		# signals_parameters=rbind(signals_parameters,multiplicities,roof_effect)
        # rownames(signals_parameters) = c('intensity','shift','half_band_width','gaussian','J_coupling','multiplicities','roof effect')
        # if (fitting_type == "Clean Fitting") {
          # colnames(signals_parameters)=paste(ROI_profile[,4],ROI_profile[,5],sep='_')
        # } else {
          # colnames(signals_parameters)=c(paste(ROI_profile[,4],ROI_profile[,5],sep='_'),paste('baseline_signal',seq(ncol(signals_parameters)-nrow(ROI_profile)),sep='_'))
        # }
        # #Generation of output data about the fitting and of the necessary variables for the generation ofa figure
        # dummy = output_generator(signals_to_quantify,fitted_signals,Ydata_2,Xdata_2,signals_parameters,multiplicities,program_parameters$buck_step)

        # #If new deconvolution has improved previous one
        # if (dummy$error1<error1) {
          # output_data=dummy$output_data
          # error1=dummy$error1
        # }}

        #Generation of the dataframe with the final output variables
        results_to_save = data.frame(
          shift = output_data$shift,
          quantification = output_data$quantification,
          signal_area_ratio = output_data$signal_area_ratio,
          fitting_error = output_data$fitting_error,
          intensity = output_data$intensity,
          half_band_width = output_data$half_band_width
        )

        #Generation of the figure data
        plot_data = rbind(output_data$signals_sum,output_data$baseline_sum,output_data$fitted_sum,output_data$signals)
        plot_data = plot_data[,ROI_buckets]

         rownames(plot_data) = c("signals_sum","baseline_sum","fitted_sum",as.character(paste(ROI_profile[,4],ROI_profile[,5],sep='_')),rep('additional signal',dim(plot_data)[1]-length(ROI_profile[,4])-3))

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
  print("Now optimizing quantifications...")


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
    if (fitting_type == "Clean Sum" ||
        fitting_type == "Baseline Sum") next


    if (length(grep("Clean",fitting_type))==1) {
      program_parameters$clean_fit="Y"
    } else {
      program_parameters$clean_fit="N"
    }
    signals_to_quantify = which(ROI_profile[, 5] >= 1)
    signals_codes = (ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2])



    print(paste(ROI_profile[1,1], ROI_profile[1,2], sep = '-'))
    hal=final_output$fitting_error[,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2],drop=F]
    hal2=apply(final_output$fitting_error,1,function(x)median(x,na.rm = T))
    hal3=hal/replicate(ncol(hal),hal2)

    index_to_use_3=unique(unlist(apply(hal3,2,function(x) which(x %in% boxplot.stats(x)$out==T))))
    print(length(index_to_use_3))
    # index_to_use_3=seq(nrow(imported_data$dataset))
    #Quantification for every spectrum
    # pb   <- txtProgressBar(1, nrow(imported_data$dataset), style=3)
    for (spectrum_index in index_to_use_3) {

      #Preparation of necessary variables to store figures and information of the fitting
      Ydata = as.numeric(imported_data$dataset[spectrum_index, ROI_buckets])
      far=(ROI_profile[,1]- ROI_profile[,2])/2-(ROI_profile[,7]-max(ll[ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]],na.rm=T))
      ROI_profile[,1]=(ROI_profile[,1]- ROI_profile[,2])/2+far
      ROI_profile[,2]=(ROI_profile[,1]- ROI_profile[,2])/2-far
      ROI_profile[,6]=predicted_shift[spectrum_index,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]]
      ROI_profile[,7]=ll[ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]]
      limits_intensity=replicate(nrow(ROI_profile),c(0,max(Ydata)))
      widthtol= rep(program_parameters$widthtolerance,nrow(ROI_profile))

      #Adaptation of the info of the parameters into a single matrix and preparation (if necessary) of the background signals that will conform the baseline
      FeaturesMatrix = fitting_prep_4(Xdata,
                                      Ydata,
                                      ROI_profile[, 5:11,drop=F],
                                      program_parameters,baselinedataset[spectrum_index,ROI_buckets],limits_intensity,widthtol)
      #Calculation of the parameters that will achieve the best fitting
      dummy = fittingloop2(FeaturesMatrix,
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
      rownames(signals_parameters) = c('intensity','shift','half_band_width','gaussian','J_coupling')
      signals_parameters=rbind(signals_parameters,multiplicities,roof_effect)

      #Generation of output data about the fitting and of the necessary variables for the generation ofa figure
      dummy = output_generator(signals_to_quantify,fitted_signals,Ydata_2,Xdata_2,signals_parameters,multiplicities,program_parameters$buck_step)
      if(mean(dummy$output_data$fitting_error[signals_to_quantify])>mean(final_output$fitting_error[spectrum_index,signals_codes[signals_to_quantify]]))next


      output_data=dummy$output_data
      error1=ifelse(is.nan(dummy$error1),3000,dummy$error1)



      #Generation of the dataframe with the final output variables
      results_to_save = data.frame(
        shift = output_data$shift,
        quantification = output_data$quantification,
        signal_area_ratio = output_data$signal_area_ratio,
        fitting_error = output_data$fitting_error,
        intensity = output_data$intensity,
        half_band_width = output_data$half_band_width
      )
      #Generation of the figure data
      plot_data = rbind(output_data$signals_sum,output_data$baseline_sum,output_data$fitted_sum,output_data$signals)
      plot_data = plot_data[,ROI_buckets]

      rownames(plot_data) = c("signals_sum","baseline_sum","fitted_sum",as.character(paste(ROI_profile[,4],ROI_profile[,5],sep='_')),rep('additional signal',dim(plot_data)[1]-length(ROI_profile[,4])-3))

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



      #Generation of output variables specific of every quantification
      final_output = save_output(spectrum_index,(ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]),results_to_save,imported_data$buck_step,final_output)

      # setTxtProgressBar(pb, spectrum_index)
    }

  }
  print("Done!")
  quantification_variables=list(final_output=final_output,useful_data=useful_data)
  return(quantification_variables)
}
