
fittingloop_bg = function(FeaturesMatrix, Xdata, Ydata, program_parameters) {

  residFun <-
    function(par, observed, xx,multiplicities,roof_effect,freq)
      observed - colSums(signal_fitting(par, xx,multiplicities,roof_effect,observed,freq))

  # Loop to control if additional signals are incorporated, until a maximum of iterations specified bt fitting_maxiterrep.
  # If at the last fitting the improvement was lesser than 25% respective to the previous fitting,
  # iterrep becomes equal to fitting_maxiterrep and the loop is stooped
  lb = as.vector(t(FeaturesMatrix[, seq(1, 9, 2), drop = F]))
  ub = as.vector(t(FeaturesMatrix[, seq(2, 10, 2), drop = F]))
  multiplicities=FeaturesMatrix[,11]
  roof_effect=FeaturesMatrix[,12]

  s0 = lb + (ub - lb) * runif(length(ub))

    nls.out <-
      nls.lm(
        par = s0,
        fn = residFun,
        observed = Ydata,
        xx = Xdata,
        multiplicities=multiplicities,
        roof_effect=roof_effect,
        lower = lb,
        upper = ub,
        freq=program_parameters$freq,
        control = nls.lm.control(
          factor = program_parameters$factor,
          maxiter = program_parameters$nls_lm_maxiter,
          ftol = program_parameters$ftol,
          ptol = program_parameters$ptol
        )
      )
      BG_intensities = coef(nls.out)[which(seq(length(coef(nls.out)))%%5==1)]
      
  return(BG_intensities)
}
