
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
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' # Not run:
#' # profiling_data=autorun(imported_data,imported_data$final_output,imported_data$useful_data,imported_data$ROI_data)


autorun2 = function(imported_data, final_output,useful_data,ROI_data) {

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
    for (spectrum_index in 1:nrow(imported_data$dataset)) {

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



    hal=final_output$fitting_error[,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2],drop=F]
    hal2=apply(final_output$fitting_error,1,function(x)median(x,na.rm = T))
    hal3=hal/replicate(ncol(hal),hal2)

    index_to_use_3=unique(unlist(apply(hal3,2,function(x) which(x %in% boxplot.stats(x)$out==T))))

    for (spectrum_index in index_to_use_3) {

      ROI_profile_2=ROI_profile
      #Preparation of necessary variables to store figures and information of the fitting
      Ydata = as.numeric(imported_data$dataset[spectrum_index, ROI_buckets])

      #If the quantification is through integration with or without baseline

      #Adaptation of the info of the parameters into a single matrix and preparation (if necessary) of the background signals that will conform the baseline
      FeaturesMatrix = fitting_prep(Xdata,
                                    Ydata,
                                    ROI_profile_2[, 5:11,drop=F],
                                    program_parameters,baselinedataset[spectrum_index,ROI_buckets])

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
      if (fitting_type == "Clean Fitting") {
        colnames(signals_parameters)=paste(ROI_profile[,4],ROI_profile[,5],sep='_')
      } else {
        colnames(signals_parameters)=c(paste(ROI_profile[,4],ROI_profile[,5],sep='_'),paste('baseline_signal',seq(ncol(signals_parameters)-nrow(ROI_profile)),sep='_'))
      }
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
        useful_data[[spectrum_index]][[signals_codes[i]]]$ROI_profile=ROI_profile_2
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
  profiling_data=list(final_output=final_output,useful_data=useful_data)
  return(profiling_data)
}



fittingloop2 = function(FeaturesMatrix,Xdata,Ydata,program_parameters) {

  #Preallocation of output and setting of necessary variables for loop
  signals_parameters=rep(0,length(as.vector(t(FeaturesMatrix[, seq(1, 9, 2), drop = F]))))
  iterrep = 0
  fitting_maxiterrep = program_parameters$fitting_maxiterrep
  signals_to_fit = which(FeaturesMatrix[, 11] != 0)
  paramprov=rep(0,nrow(FeaturesMatrix)*5)





  #Necessary information to incorporate additional signals if necessary
  range_ind = round(program_parameters$additional_signal_ppm_distance / program_parameters$buck_step)



  #Function where to find a minimum
  residFun <-
    function(par, observed, xx,multiplicities,roof_effect,freq,bins)
      observed[bins] - colSums(signal_fitting(par, xx,multiplicities,roof_effect,freq))[bins]


  # Loop to control if additional signals are incorporated, until a maximum of iterations specified bt fitting_maxiterrep.
  # If at the last fitting the improvement was lesser than 25% respective to the previous fitting,
  # iterrep becomes equal to fitting_maxiterrep and the loop is stooped
  while (iterrep <= fitting_maxiterrep) {

    iter = 0
    errorprov = error1 = 3000
    worsterror = 0
    dummy = error2 = 3000
    multiplicities=FeaturesMatrix[,11]
    roof_effect=FeaturesMatrix[,12]
    fitted_signals = signal_fitting(as.vector(t(FeaturesMatrix[,c(2,3,6,7,9)])),
                                    Xdata,multiplicities,roof_effect,program_parameters$freq)
    bins=c()
    for (i in signals_to_fit) {
      sorted_bins=sort(fitted_signals[i,]/sum(fitted_signals[i,]),decreasing=T,index.return=T)
      if(length(sorted_bins$x)>0) {
        bins2= sorted_bins$ix[1:which.min(abs(cumsum(sorted_bins$x)-0.9))]
        distance=diff(FeaturesMatrix[i,3:4])/program_parameters$buck_step
        distance2=round(diff(FeaturesMatrix[i,9:10])/program_parameters$buck_step/program_parameters$freq)

        bins2=unique(as.vector(sapply(seq(distance),function(x)bins2-x)))
        bins2=unique(bins2,min(bins2)-distance2,max(bins2)+distance2)
        bins2=bins2[bins2>0&bins2<length(Xdata)]

        # aa=peakdet(fitted_signals[i,],0.00001)$maxtab$pos
        # aa=as.vector(sapply(seq(distance),function(x)aa-x))
        # #
        #  lol=sapply(seq(distance),function(x)min(Ydata[bins2]-fitted_signals[i,(bins2-x)]))
        # FeaturesMatrix[i,2]=FeaturesMatrix[i,2]+lol


        bins=unique(c(bins,bins2))
      }}

    #Depending on the complexity of the ROI, more or less iterations are performed
    if (is.numeric(program_parameters$fitting_maxiter)) {
      fitting_maxiter = program_parameters$fitting_maxiter
    } else {
      if (nrow(FeaturesMatrix)> 8 |
          any(FeaturesMatrix[, 4] - FeaturesMatrix[, 3] > 0.01)) {
        fitting_maxiter = 20
      } else if ((nrow(FeaturesMatrix)> 5 &&
                  nrow(FeaturesMatrix)< 9)) {
        fitting_maxiter = 14
      } else {
        fitting_maxiter = 8
      }
    }


    #Conditions to keep the loop:
    # -The error is bigger than the specified in program_parameters
    # -There is no fitting with more than 66.7% improvement from the worst solution
    # -The loop has not arrived the specified maximum of iterations
    while (error1 > program_parameters$errorprov &error1 > (1 / 3 * worsterror) & iter < fitting_maxiter) {
      #Initialization of parameters to optimize. In every iteration the initialization will be different
      lb = as.vector(t(FeaturesMatrix[, seq(1, 9, 2), drop = F]))
      ub = as.vector(t(FeaturesMatrix[, seq(2, 10, 2), drop = F]))
      set.seed(iter);s0 = lb + (ub - lb) * runif(length(ub))
      order1=order(rowMeans(FeaturesMatrix[signals_to_fit ,3:4,drop=F])[signals_to_fit])

      # aaa=iter%%3/3
      # bbb=ifelse((iter+1)%%3/3==0,1,(iter+1)%%3/3)
      # s0[which(seq_along(s0)%%5==2)]=lb[which(seq_along(s0)%%5==2)] + (ub[which(seq_along(s0)%%5==2)] - lb[which(seq_along(s0)%%5==2)]) * runif(1,min=aaa,max=bbb)

      #During the first two iterations, find the peaks on the region of the spectrum. If the number of peaks is the same that the expected on the ROI and the location is similar, the signals are located where there are the peaks with minimum shift tolerance.
      peaks_xdata = peakdet(c(Ydata[1],diff(Ydata)), program_parameters$peakdet_minimum*0.1*max(1e-10,max(Ydata)),Xdata)

      if (iter<4&length(peaks_xdata$maxtab$val)>0) {
        peaks_bindata = peakdet(c(Ydata[1],diff(Ydata)), program_parameters$peakdet_minimum*0.1*max(1e-10,max(Ydata)))
        peaks=peaks_xdata$maxtab$pos[sort(peaks_xdata$maxtab$val,decreasing=T,index.return=T)$ix[1:sum(multiplicities[signals_to_fit])]]
        peaks_compare=rowMeans(FeaturesMatrix[signals_to_fit ,3:4,drop=F])
        for (i in 1:length(peaks_compare)) {
          ind=sort(abs(peaks-peaks_compare[i]),index.return=T)$ix[1:multiplicities[i]]
          if (!is.na(mean(peaks[ind]))&&mean(peaks[ind])>FeaturesMatrix[i,3]&&mean(peaks[ind])<FeaturesMatrix[i,4]) {
            s0[which(seq_along(s0)%%5==2)[i]]=mean(peaks[ind])
            lb[which(seq_along(s0)%%5==2)[i]]=mean(peaks[ind])-0.001
            ub[which(seq_along(s0)%%5==2)[i]]=mean(peaks[ind])+0.001
          }
        }

        #Main optimization
        set.seed(iter);nls.out <-
          minpack.lm::nls.lm(
            par = s0,
            fn = residFun,
            observed = Ydata,
            xx = Xdata,
            multiplicities=multiplicities,
            roof_effect=roof_effect,
            freq=program_parameters$freq,
            lower = lb,
            upper = ub,
            control = minpack.lm::nls.lm.control(
              factor = program_parameters$factor,
              maxiter = program_parameters$nls_lm_maxiter,
              ftol = program_parameters$ftol,
              ptol = program_parameters$ptol
            )
          )

        # #Procedure to calculate the fititng error in all the ROI
        #An adapted MSE error is calculated, and the parameters of the optimization with less MSE are stored
        iter = iter + 1
        order2=order(coef(nls.out)[which(seq_along(coef(nls.out))%%5==2)][signals_to_fit])

        errorprov = (sqrt(nls.out$deviance / length(Ydata))) * 100 / (max(Ydata) -min(Ydata))
        if (is.nan(errorprov) || is.na(errorprov)) errorprov = error1
        if (errorprov < error1 && identical(order1,order2)) {
          error1 = errorprov
          paramprov=coef(nls.out)
        } else if (errorprov > worsterror) {
          worsterror = errorprov
        }
      } else {
        #If in the first two iterations the procedure of finding peaks is not effective enough, the irignal chemical shift and chemical shift tolerance of every signal is maintained
        set.seed(iter);nls.out <-
          minpack.lm::nls.lm(
            par = s0,
            fn = residFun,
            observed = Ydata,
            xx = Xdata,
            multiplicities=multiplicities,
            roof_effect=roof_effect,
            freq=program_parameters$freq,
            lower = lb,
            upper = ub,
            control = minpack.lm::nls.lm.control(
              factor = program_parameters$factor,
              maxiter = program_parameters$nls_lm_maxiter,
              ftol = program_parameters$ftol,
              ptol = program_parameters$ptol
            )

          )

        iter = iter + 1

        order2=order(coef(nls.out)[which(seq_along(coef(nls.out))%%5==2)][signals_to_fit])
        # #Procedure to calculate the fititng error in all the ROI
        #An adapted MSE error is calculated, and the parameters of the optimization with less MSE are stored
        errorprov = (sqrt(nls.out$deviance / length(Ydata))) * 100 / (max(Ydata) -min(Ydata))
        if (is.nan(errorprov) || is.na(errorprov))errorprov = error1
        if (errorprov < error1 && identical(order1,order2)) {
          error1 = errorprov
          paramprov=coef(nls.out)
        } else if (errorprov > worsterror) {
          worsterror = errorprov
        }
      }}
    signals_parameters = paramprov

    fitted_signals = signal_fitting(signals_parameters,
                                    Xdata,multiplicities,roof_effect,program_parameters$freq)

    bins=c()
    for (ind in signals_to_fit) {
      sorted_bins=sort(fitted_signals[ind,]/sum(fitted_signals[ind, ]),decreasing=T,index.return=T)
      if(length(sorted_bins$x)>0) bins= sorted_bins$ix[1:which.min(abs(cumsum(sorted_bins$x)-0.75))]

    }
    if (length(bins)==0) bins=seq_along(Ydata)

    #Correction of half_band_width and j-coupling
    iter = 0
    error22=error2=error1
    errorprov = error1=3000
    #Only half_band_width and j-coupling will have different lower und upper bounds.
    change_indexes=which(seq_along(lb)%%5!=3 & seq_along(lb)%%5!=4 & seq_along(lb)%%5!=0)
    lb[change_indexes]=ub[change_indexes]=signals_parameters[change_indexes]
    while (iter < 3) {
      set.seed(iter);s0 = lb + (ub - lb) * runif(length(ub))
      nls.out <-
        minpack.lm::nls.lm(
          par = s0,
          fn = residFun,
          observed = Ydata,
          xx = Xdata,
          multiplicities=multiplicities,
          roof_effect=roof_effect,
          freq=program_parameters$freq,
          lower = lb,
          upper = ub,
          control = minpack.lm::nls.lm.control(
            factor = program_parameters$factor,
            maxiter = program_parameters$nls_lm_maxiter,
            ftol = program_parameters$ftol,
            ptol = program_parameters$ptol
          )

        )
      iter = iter + 1
      # #Procedure to calculate the fititng error in all the ROI
      #An adapted MSE error is calculated, and the parameters of the optimization with less MSE are stored
      errorprov = (sqrt(nls.out$deviance / length(Ydata))) * 100 / (max(Ydata) -
                                                                      min(Ydata))
      if (is.nan(errorprov) || is.na(errorprov)) errorprov = error1
      if (errorprov < error1) {
        error1 = errorprov
        paramprov=coef(nls.out)
      } else if (errorprov > worsterror) {
        worsterror = errorprov
      }
      if (error1 < error2) {
        error2 = error1
        signals_parameters = paramprov
      }
    }

    #If half_band_width and j-coup change improves fitting


    iterrep = iterrep + 1



    #If the fitting seems to be still clearly improvable through the addition of signals
    if (iterrep <= fitting_maxiterrep& error22 < (program_parameters$additional_signal_improvement * dummy) &
        (error22 > program_parameters$additional_signal_percentage_limit)&length(peaks_xdata$maxtab$pos)>sum(multiplicities[signals_to_fit])) {
      # print('Trying to improve initial fit adding peaks')

      #I find peaks on the residuals
      residual_peaks = tryCatch(peakdet(c(nls.out$fvec[1],diff(nls.out$fvec)), program_parameters$peakdet_minimum*max(1e-10,max(Ydata))),error= function(e){
        dummy=list(signals_parameters=signals_parameters,error1=error1)
        return(dummy)
      })

      if (is.null(residual_peaks$maxtab) == F) {
        #Preparation of information of where signals of interest are located
        dummy=multiplicities[signals_to_fit]%%2
        dummy[dummy==0]=2
        additional_signal_matrix = matrix(paramprov,nrow(FeaturesMatrix),5,byrow = TRUE)
        points_to_avoid = abs(rbind(matrix(Xdata,length(signals_to_fit),length(Xdata),byrow = TRUE) - matrix(
          additional_signal_matrix[signals_to_fit, 2] - (additional_signal_matrix[signals_to_fit, 5]/dummy)/program_parameters$freq,length(signals_to_fit),length(Xdata)),
          matrix(Xdata,length(signals_to_fit),length(Xdata),byrow = TRUE) - matrix(additional_signal_matrix[signals_to_fit, 2] + (additional_signal_matrix[signals_to_fit, 5]/dummy)/program_parameters$freq,length(signals_to_fit),length(Xdata))))
        points_to_avoid = apply(points_to_avoid, 1, which.min)
        seq_range = c()
        for (i in (-range_ind):range_ind) seq_range = append(seq_range, points_to_avoid - i)

        #Finding of posible additional signals to incorporate if there are not in zones where the signals o interest are located
        residual_peaks=cbind(residual_peaks$maxtab$pos, residual_peaks$maxtab$val)[residual_peaks$maxtab$pos %in% peaks_bindata$maxtab$pos,,drop=F]
        valid_residual_peaks = matrix(NA, 0, 2)
        if (nrow(residual_peaks)>0) {
          for (i in seq(nrow(residual_peaks))) {
            if (any(abs(points_to_avoid - residual_peaks[i, 1]) < range_ind) == F) valid_residual_peaks = rbind(valid_residual_peaks, residual_peaks[i, ])
          }
          #Selection of more intense additional signals
          if (nrow(valid_residual_peaks) > program_parameters$signals_to_add) {
            ad = sort(valid_residual_peaks[, 2],decreasing = T,index.return = T)$ix
            valid_residual_peaks = valid_residual_peaks[ad[1:min(program_parameters$signals_to_add, length(ad))], , drop = F]
          }
          #Creation of rows to incorporate to FeaturesMatrix
          if (nrow(valid_residual_peaks)>0) {
            dummy = t(replicate(nrow(valid_residual_peaks),FeaturesMatrix[1,]))
            dummy[, 2] = Ydata[valid_residual_peaks[, 1]]
            dummy[, 3] = Xdata[valid_residual_peaks[, 1]] - 0.001
            dummy[, 4] = Xdata[valid_residual_peaks[, 1]] + 0.001
            dummy[, 5]=min(FeaturesMatrix[,5])
            dummy[, 6]=min(FeaturesMatrix[,6])
            dummy[, c(9,10,12)] = rep(0, nrow(valid_residual_peaks))
            dummy[, 11] = rep(1, nrow(valid_residual_peaks))
            FeaturesMatrix = rbind(FeaturesMatrix, dummy)
          } else {
            iterrep = fitting_maxiterrep +1
          }
        } else  {
          iterrep = fitting_maxiterrep +1
        }}
    } else {
      iterrep = fitting_maxiterrep +1
    }

  }
  optim_parameters=list(signals_parameters=signals_parameters,error1=error1)
  return(optim_parameters)
}

