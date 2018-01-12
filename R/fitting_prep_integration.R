
fitting_prep_integration = function(Xdata,Ydata,program_parameters,created_baseline) {
  #Created by Daniel Canueto 30/08/2016
  #Preparation of parameters to optimize to achieve the best fitting


  Ydata[Ydata<0]=0

  ROIlength = length(Xdata)


  #Calculation of number of background signals, if baseline fitting is performed
  BGSigNum = ifelse(program_parameters$clean_fit == 'N', max(round(abs(Xdata[1] -
                                                                        Xdata[ROIlength]) * program_parameters$BGdensity), 3), 0)

  #Preallocation of parameters to optimize into a matrix of features
  FeaturesMatrix = matrix(NA, (BGSigNum), 12)

  #Finding of maximum intensity and $chemical_shift tolerance of every background signal
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


    #Parameters of background signals
    FeaturesMatrix[, 1] = 0
    FeaturesMatrix[, 2] = BGSig_maximums
    FeaturesMatrix[, 3] = BGSigrightlimits
    FeaturesMatrix[, 4] = BGSigleftlimits
    FeaturesMatrix[, 5] = program_parameters$BG_width-program_parameters$BG_width*program_parameters$BG_width_tolerance
    FeaturesMatrix[, 6] = program_parameters$BG_width+program_parameters$BG_width*program_parameters$BG_width_tolerance
    FeaturesMatrix[, 7] = 0
    FeaturesMatrix[, 8] = program_parameters$BG_gaussian_percentage
    FeaturesMatrix[, 9] = 0
    FeaturesMatrix[, 10] = 0 #j coupling makes no sense with backgorund signals
    FeaturesMatrix[, 11] = 0 #arbitrary number used to signal later background signals
    FeaturesMatrix[, 12] = 0



    # optimization of baseline parameters , to be sure that the algorithm doesn ot try ti fot spurious signals as basleine
    baseline = fittingloop_bg(FeaturesMatrix,
                                Xdata,
                                created_baseline,
                                program_parameters)$baseline


  }


  return(baseline)
}
