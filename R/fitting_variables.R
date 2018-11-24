#' Preparation of parameters involved with the program. A list of manual changes to parameters can be provided.
#'
#' @return Preparation of parameters involved with the program
#' @export fitting_variables
#' @param parameter_change List of parameters ot change, if desired.
#'
#' @examples
#' program_parameters=fitting_variables()
#' program_parameters=fitting_variables(list(fitting_maxiterrep=1,fitting_error_analysis_limit=0.05))



fitting_variables = function(parameter_change =list()) {


  program_parameters=list()

  program_parameters$spectrum_borders=c(12,-0.5)


  program_parameters$BGdensity=70 #Density of signals to prepare abaseline below the signals to fit
program_parameters$widthtolerance=0.25 #Allowed Variability of halfwidth
program_parameters$gaussian=0.1 #Allowed Variability of gaussian percentage
program_parameters$j_coupling_variation=0.2 #Allowed Variability of j-coupling
program_parameters$BG_gaussian_percentage=0.2 #Allowed gaussian percentage of baseline signals
program_parameters$BG_width=8
program_parameters$BG_width_factor=5
program_parameters$BG_width_tolerance=0.1 #Allowed Variability of half bandwidth of baseline signals

#Parameters related to the fitting loop
program_parameters$errorprov=3 #Percentage limit of fitting error on the ROI. If a lower percentage is reached, the solution is considered optimized and the optimization ends
program_parameters$fitting_maxiter=NA #The number of maximum iterations of optimization can be specified

#Parameters related to the nls algorithm. Read documentation of nls.lm package for details.
program_parameters$nls_lm_maxiter=500
program_parameters$ftol=1e-6
program_parameters$ptol=1e-6
program_parameters$factor=0.01

#Parameters related to the addition of other signals post fitting
program_parameters$additional_signal_ppm_distance=0.002 #Allowed distance for added peaks
program_parameters$signals_to_add = 2 #Allowed distance for added peaks to the ROI
program_parameters$fitting_maxiterrep = 0 #Allowed tries to add peaks
program_parameters$additional_signal_improvement=0.75 #Improvement of ROI by adding peaks. If previous addition did not achieve less than 75% of fitting error, process of addition of peaks is stopped
program_parameters$additional_signal_percentage_limit=3 #If fititng erorr is less tha nthis percentage, addition of peaks is not performed
program_parameters$peakdet_minimum=0.01 #Limit to find a peak as relevant enough to add it.


#Parameters related to criteria to automatically discard quantifications during analysis
program_parameters$automatic_removal='Y'
program_parameters$fitting_error_analysis_limit=0.15
program_parameters$signal_area_ratio_analysis_limit=10

if (length(parameter_change)>0) {
for (i in seq(length(parameter_change))) program_parameters[[which(names(program_parameters)==names(parameter_change)[i])]]=parameter_change[[i]]
}

return(program_parameters)
}
