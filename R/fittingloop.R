#Workflow of optimization of metabolite and baseline signals parameters. There are several optimization iterations ("iter") and the one with less error is chosen. The maximum number of iterations depends on the complexity of the ROI. Then the half bandwidth and j coupling are optimized separately as they have diferent behavior with the least squares algorithm. After the first fitting there can be further ones, depending on the need to add additional signals to adapt the signals and ROI information provided by the user to the concrete characteristics of the spectrum.


fittingloop = function(FeaturesMatrix,Xdata,Ydata,program_parameters) {

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
      s0 = lb + (ub - lb) * runif(length(ub))
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
# for (i in signals_to_fit)    {
#   aa=peakdet(fitted_signals[i,],0.00001)$maxtab$pos
#   if (length(aa)==0) next
#   bb=min(Ydata[aa]-colSums(fitted_signals[,aa,drop=F]))
#   if (bb<0) signals_parameters[which(seq_along(lb)%%5==1)[i]]=max(signals_parameters[which(seq_along(lb)%%5==1)[i]]+bb,0)
# }

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
      s0 = lb + (ub - lb) * runif(length(ub))
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
