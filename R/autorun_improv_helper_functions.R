

fitting_prep_2 = function(Xdata,Ydata,initial_fit_parameters,program_parameters,created_baseline,max_shift,
                          min_shift,max_intensity,
                          min_intensity,max_width,min_width,spectrum_index,signal_index) {
  Ydata[Ydata<0]=0
  min_intensity[spectrum_index,signal_index][is.na(min_intensity[spectrum_index,signal_index])]=0
  max_intensity[spectrum_index,signal_index][is.na(max_intensity[spectrum_index,signal_index])]=max(Ydata)

  colnames(initial_fit_parameters) = c(
    "quantification_or_not",
    "positions",
    "shift_tolerance",
    "widths",
    "multiplicities",
    "Jcoupling",
    "roof_effect"
  )
  signals_to_fit = length(initial_fit_parameters$positions)
  ROIlength = length(Xdata)


  #Calculation of number of background signals, if baseline fitting is performed
  BGSigNum = ifelse(program_parameters$clean_fit == 'N', max(round(abs(Xdata[1] -
                                                                         Xdata[ROIlength]) * program_parameters$BGdensity), 3), 0)

  #Preallocation of parameters to optimize into a matrix of features
  FeaturesMatrix = matrix(NA, (signals_to_fit + BGSigNum), 12)
  colnames(FeaturesMatrix) = c(
    'minimum_intensity',
    'maximum_intensity',
    'shift_left_limit',
    'shift_right_limit',
    'minimum_width',
    'maximum_width',
    'minimum_gaussian',
    'maximum_gaussian',
    'minimum_J_coupling',
    'maximum_J_coupling',
    'multiplicities',
    'roof_effect'
  )

  #Parameters of signals to fit
  FeaturesMatrix[1:signals_to_fit, 1] = min_intensity[spectrum_index,signal_index]
  FeaturesMatrix[1:signals_to_fit, 2] = max_intensity[spectrum_index,signal_index]
  FeaturesMatrix[1:signals_to_fit, 3] = min_shift[spectrum_index,signal_index]
  FeaturesMatrix[1:signals_to_fit, 4] = max_shift[spectrum_index,signal_index]
  FeaturesMatrix[1:signals_to_fit, 5] = min_width[spectrum_index,signal_index]
  FeaturesMatrix[1:signals_to_fit, 6] = max_width[spectrum_index,signal_index]
  FeaturesMatrix[1:signals_to_fit, 7] = 0
  FeaturesMatrix[1:signals_to_fit, 8] = program_parameters$gaussian
  FeaturesMatrix[1:signals_to_fit, 9] = initial_fit_parameters$Jcoupling -
    program_parameters$j_coupling_variation
  FeaturesMatrix[1:signals_to_fit, 10] = initial_fit_parameters$Jcoupling +
    program_parameters$j_coupling_variation
  FeaturesMatrix[1:signals_to_fit, 11] = initial_fit_parameters$multiplicities
  FeaturesMatrix[1:signals_to_fit, 12] = initial_fit_parameters$roof_effect


  FeaturesMatrix[initial_fit_parameters$multiplicities==1, 9:10] = 0

  #Finding of maximum intensity and chemical shift tolerance of every background signal
  if (BGSigNum>0) {
    BGSigrightlimits = seq(Xdata[1]-0.005, Xdata[ROIlength]+0.005, length = BGSigNum) -
      0.005
    BGSigleftlimits = BGSigrightlimits + 0.01

    peaks = peakdet(Ydata, program_parameters$peakdet_minimum*max(1e-10,max(Ydata)))
    left = which(peaks$mintab$pos < ROIlength / 5)
    right = which(peaks$mintab$pos > 4 * ROIlength / 5)
    dummy = round(seq(1, ROIlength, length = 2 * BGSigNum - 1))
    BGleftlimits = dummy[c(1, seq(2, length(dummy) - 1, 2))]
    BGrightlimits = dummy[c(seq(2, length(dummy) - 1, 2), length(dummy))]
    BGSig_maximums = replicate(BGSigNum, NA)
    for (ss in 1:BGSigNum)
      BGSig_maximums[ss] = min(Ydata[BGleftlimits[ss]:BGrightlimits[ss]])

    BG_width=max(min(initial_fit_parameters$widths,na.rm=T)*program_parameters$BG_width_factor,program_parameters$BG_width)
    #Parameters of background signals
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 1] = 0
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 2] = BGSig_maximums
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 3] = BGSigrightlimits
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 4] = BGSigleftlimits
    # FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 5] = (1.5 /
    #                                                                     program_parameters$freq) * 10
    # FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 6] = (1.5 /
    #                                                                      program_parameters$freq) * 15
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 5] = BG_width*(1-program_parameters$BG_width_tolerance)
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 6] = BG_width*(1+program_parameters$BG_width_tolerance)
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 7] = 0
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 8] = program_parameters$BG_gaussian_percentage
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 9] = 0
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 10] = 0 #j coupling makes no sense with backgorund signals
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 11] = 0 #arbitrary number used to signal later background signals
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix), 12] = 0



    # optimization of baseline parameters , to be sure that the algorithm doesn ot try ti fot spurious signals as basleine
    FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix),2] = fittingloop_bg(FeaturesMatrix[(signals_to_fit + 1):nrow(FeaturesMatrix),],
                                                                                 Xdata,
                                                                                 created_baseline,
                                                                                 program_parameters)$BG_intensities


  }


  return(FeaturesMatrix)
}
rf_pred_intensity=function(initial_matrix,met_names,fitting_error) {
  fitting_error=fitting_error
  modified_matrix=initial_matrix
  colnames(modified_matrix)=make.names(colnames(modified_matrix))
  dummy=which(apply(modified_matrix,2,function(x)! all(is.na(x)))==T)
  for (i in dummy) modified_matrix[fitting_error[,i] %in%
                                     boxplot.stats(fitting_error[,i])$out,i]=NA
  dummy2=missRanger::missRanger(as.data.frame(modified_matrix[,dummy]))
  modified_matrix[,dummy]=as.matrix(dummy2)
  analyzed_signals=which(apply(modified_matrix,2,function(x)length(which(is.na(x))))<0.1*nrow(modified_matrix))
  predicted_matrix=lower_bound_matrix=upper_bound_matrix=as.data.frame(matrix(NA,nrow(modified_matrix),ncol(modified_matrix)))
  ctrl <- caret::trainControl(method = "boot632",number=18,savePredictions="all")
  for (i in analyzed_signals) {
    sed=intersect(which(met_names==met_names[i]),analyzed_signals)
    if (length(sed)==1) next
    tel=cbind(modified_matrix[,setdiff(sed,i)],prcomp(scale(modified_matrix[,sed]))$x[,1])
    ind3=which(abs(cor(modified_matrix[,i],tel,method='spearman'))>0.7)
    if (length(ind3)==0) next
    tel=tel[,ind3]
    training_data=data.frame(y=modified_matrix[,i],scale(tel))
    plsFit <- caret::train(y ~ .,data = training_data,method = "rf",trControl = ctrl)
    iii=plsFit$pred
    ff=sapply(seq(nrow(training_data)),function(x)quantile(rnorm(1000,mean=mean(iii$pred[iii$rowIndex==x],na.rm=T),
                                                                 sd=sd(iii$pred[iii$rowIndex==x],na.rm=T)),c(0.025,0.5,0.975),na.rm=T))

    predicted_matrix[,i]=ff[2,]
    lower_bound_matrix[,i]=ff[1,]
    upper_bound_matrix[,i]=ff[3,]

  }
  if (all(is.na(predicted_matrix))==F) {
  predicted_matrix=missRanger::missRanger(predicted_matrix)
  lower_bound_matrix=missRanger::missRanger(lower_bound_matrix)
  upper_bound_matrix=missRanger::missRanger(upper_bound_matrix)
  }
  output=list(predicted_matrix=predicted_matrix,lower_bound_matrix=lower_bound_matrix,upper_bound_matrix=upper_bound_matrix)
  return(output)
}

rf_pred=function(initial_matrix,fitting_error) {
  fitting_error=fitting_error
  modified_matrix=initial_matrix
  colnames(modified_matrix)=make.names(colnames(modified_matrix))
  modified_matrix=jitter(modified_matrix,0.00001)
  dummy=which(apply(modified_matrix,2,function(x) all(is.finite(x)))==T)

  for (i in dummy) {
    modified_matrix[fitting_error[,i] %in%
                      boxplot.stats(fitting_error[,i])$out,i]=NA
  }
  dummy2=missRanger::missRanger(as.data.frame(modified_matrix[,dummy]))
  modified_matrix[,dummy]=as.matrix(dummy2)
  analyzed_signals=which(apply(modified_matrix,2,function(x) all(is.finite(x)))==T)
  predicted_matrix=lower_bound_matrix=upper_bound_matrix=as.data.frame(matrix(NA,nrow(modified_matrix),ncol(modified_matrix)))
  ctrl <- caret::trainControl(method = "boot632",number=18,savePredictions="all")
  for (i in analyzed_signals) {
    training_data=data.frame(y=modified_matrix[,i],scale(modified_matrix[,setdiff(analyzed_signals,i)]))
    plsFit <- caret::train(y ~ .,data = training_data,method = "lasso",trControl = ctrl)
    ind=which(varImp(plsFit)$importance$Overall>30)
    if (length(ind==1)) ind=order(varImp(plsFit)$importance$Overall,decreasing=T)[1:2]
    training_data=data.frame(y=modified_matrix[,i],scale(modified_matrix[,setdiff(analyzed_signals,i)][,ind]))
    plsFit <- caret::train(y ~ .,data = training_data,method = "ranger",trControl = ctrl)
    iii=plsFit$pred
    ff=sapply(seq(nrow(training_data)),function(x)quantile(rnorm(1000,mean=mean(iii$pred[iii$rowIndex==x],na.rm=T),
                                                                 sd=sd(iii$pred[iii$rowIndex==x],na.rm=T)),c(0.025,0.5,0.975),na.rm=T))
    predicted_matrix[,i]=ff[2,]
    lower_bound_matrix[,i]=ff[1,]
    upper_bound_matrix[,i]=ff[3,]
  }


  predicted_matrix=missRanger::missRanger(predicted_matrix)
  lower_bound_matrix=missRanger::missRanger(lower_bound_matrix)
  upper_bound_matrix=missRanger::missRanger(upper_bound_matrix)

  output=list(predicted_matrix=predicted_matrix,lower_bound_matrix=lower_bound_matrix,upper_bound_matrix=upper_bound_matrix)
  return(output)
}
