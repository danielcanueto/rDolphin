
#' Improvement of automatic profiling using the inforamtion collected on a previous implementation.
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param reproducibility_data List with necessary information to load quantifications on the Shiny GUI.
#' @param ROI_data ROIs data.
#' @param improvement_option If "correction", quantifications are updated taking into account the predicted signal parameters. If "reimplemetation", profiling is repeated using the prediction information.
#' @param level How extensive should be the improvement? If "all", all quantifications are changed. If "outliers", quantifications whoss signal parameters behave as outliers are changed (please take into account that only the quantifications will be updated). If a number is introduced, the quantifications with a higher fitting error than the number specified are repeated.
#'
#' @return List with final_output (with metabolite signal relative concentrations and quality indicators) and reproducibility_data (with the necessary data to reproduce the profiling performed).
#' @export automatic_profiling_improv
#' @import baseline
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' # Not run:
#' # load(file.path(system.file(package = "rDolphin"),"extdata","MTBLS242_subset_profiling_data.RData"))
#' # profiling_data_2=automatic_profiling_improv(imported_data,profiling_data$final_output,profiling_data$reproducibility_data,imported_data$ROI_data)


#TODO: Choose criteria to repeat only individual quantification and all signals of all spectra.

automatic_profiling_improv = function(imported_data, final_output,reproducibility_data,ROI_data,improvement_option='reimplementation',level='outliers') {
print("Starting maximization of profiling data quality using information of original profiling...")
  print("Now estimating the predicted signal parameters with prediction intervals...")
  capture.output(predicted_info<-signparpred(final_output$half_bandwidth,fitting_error=final_output$fitting_error))

  predicted_width=as.matrix(predicted_info$predicted_matrix)

  max_width=as.matrix(predicted_info$upper_bound_matrix)
  min_width=as.matrix(predicted_info$lower_bound_matrix)
ind=which(is.na(predicted_width[1,]))
if (length(ind)>0) {
predicted_width[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,8]))
min_width[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,8]*0.75))
max_width[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,8]*1.25))
}
capture.output(predicted_info<-signparpred(final_output$chemical_shift,fitting_error=final_output$fitting_error))

predicted_shift=as.matrix(predicted_info$predicted_matrix)
max_shift=as.matrix(predicted_info$upper_bound_matrix)
min_shift=as.matrix(predicted_info$lower_bound_matrix)

  ind=which(is.na(predicted_shift[1,]))
  if (length(ind)>0) {
    predicted_shift[,ind]=as.matrix(t(replicate(nrow(predicted_width),ROI_data[ind,6])))
  max_shift[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,6]+ROI_data[ind,7]))
  min_shift[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,6]-ROI_data[ind,7]))
  }

  capture.output(predicted_info<-signparpred(final_output$intensity,fitting_error=final_output$fitting_error,met_names=ROI_data[,4]))


  predicted_intensity=as.matrix(predicted_info$predicted_matrix)
  max_intensity=as.matrix(predicted_info$upper_bound_matrix)
  min_intensity=as.matrix(predicted_info$lower_bound_matrix)

  ind=which(is.na(predicted_intensity[1,]))
  max_intensity[!is.finite(max_intensity)]=NA
  min_intensity[!is.finite(min_intensity)]=NA
  min_intensity[min_intensity<0]=0
  print("Beginning the maximization of profiling data quality according to the option selected...")



  quantifications_to_repeat=matrix(0,nrow(predicted_width),ncol(predicted_width))
  if (level=="all") quantifications_to_repeat[,]=1
	if (is.numeric(level)) quantifications_to_repeat[which(final_output$fitting_error>level)]=1
	if (level=="outliers") {
	  tryCatch({

	  outlier_indicator=sapply(which(!is.na(predicted_shift)),
	             function(x)findInterval(final_output$chemical_shift[x],
	                                     c(min_shift[x],max_shift[x])))
	  if (length(outlier_indicator)>0) quantifications_to_repeat[which(!is.na(predicted_shift))][sapply(outlier_indicator,function(x)x==0|x==2)]=1
	  outlier_indicator=sapply(which(!is.na(predicted_width)),
	             function(x)findInterval(final_output$half_bandwidth[x],
	                                     c(min_width[x],max_width[x])))
	  if (length(outlier_indicator)>0) quantifications_to_repeat[which(!is.na(predicted_width))][sapply(outlier_indicator,function(x)x==0|x==2)]=1
	  outlier_indicator=sapply(which(!is.na(predicted_intensity)),
	             function(x)findInterval(final_output$intensity[x],
	                                     c(min_intensity[x],max_intensity[x])))
	  if (length(outlier_indicator)>0) quantifications_to_repeat[which(!is.na(predicted_intensity))][sapply(outlier_indicator,function(x)x==0|x==2)]=1

	  }, error=function(e)quantifications_to_repeat[,]=1)
	  }


if (improvement_option=='reimplementation') {  #Splitting of ROI data into individual ROIs to be quantified
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
    if (length(ROI_buckets)<20) {
	print ("Ignoring ROI as width is too small")
	next
	}
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


	index_to_use_3=which(rowSums(quantifications_to_repeat[,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2],drop=F])>0)
	pb   <- txtProgressBar(1, nrow(imported_data$dataset), style=3)

    #Quantification for every spectrum
    for (spectrum_index in index_to_use_3) {

         #Preparation of necessary variables to store figures and information of the fitting
      output=profiling_func(spectrum_index,signals_codes,
                            imported_data,
                            ROI_buckets,fitting_type,
                            program_parameters,Xdata,Ydata,
                            final_output,
                            reproducibility_data,
                            ROI_profile,baselinedataset,
                            signals_to_quantify,pb)
      final_output=output$final_output
      reproducibility_data=output$reproducibility_data
      }

  }


 } else if (improvement_option=='correction') {

prova_intensity=predicted_intensity
prova_intensity[,apply(predicted_intensity,2,function(x)all(is.na(x)))]=final_output$intensity[,apply(predicted_intensity,2,function(x)all(is.na(x)))]
prova_shift=predicted_shift
prova_shift[,apply(predicted_shift,2,function(x)all(is.na(x)))]=final_output$chemical_shift[,apply(predicted_shift,2,function(x)all(is.na(x)))]
prova_width=predicted_width
prova_width[,apply(predicted_width,2,function(x)all(is.na(x)))]=final_output$half_bandwidth[,apply(predicted_width,2,function(x)all(is.na(x)))]

tec=sapply(seq(length(prova_intensity)),function(x)sum(peakpvoigt(c(prova_intensity[x],prova_shift[x],prova_width[x]*0.5/600.2,0),imported_data$ppm))*imported_data$buck_step)
dim(tec)=dim(prova_intensity)
tec[,apply(tec,2,function(x)all(is.na(x)))]=final_output$quantification[,apply(tec,2,function(x)all(is.na(x)))]
	for (i in 1:ncol(final_output$quantification)) {
	index_to_use_3=which(quantifications_to_repeat[,i]>0)
	final_output$quantification[index_to_use_3,i]=tec[index_to_use_3,i]
}
 }
  print("Done!")
  profiling_data=list(final_output=final_output,reproducibility_data=reproducibility_data,
                                predicted_shift=predicted_shift,predicted_width=predicted_width,
                                predicted_intensity=predicted_intensity,max_width=max_width,
                                min_width=min_width,max_shift=max_shift,min_shift=min_shift,
                                max_intensity=max_intensity,min_intensity=min_intensity)
  return(profiling_data)
}



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
