
#' Improvement of automatic profiling using the inforamtion collected on a previous implementation.
#'
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param useful_data List with necessary information to load quantifications on the Shiny GUI.
#' @param ROI_data ROIs data.
#' @param improvement_option If "correction", quantifications are updated taking into account the predicted signal parameters. If "reimplemetation", profiling is repeated using the prediction information.
#' @param level How extensive should be the improvement? If "all", all quantifications are changed. If "outliers", quantifications whoss signal parameters behave as outliers are changed (please take into account that only the quantifications will be updated). If a number is introduced, the quantifications with a higher fitting error than the number specified are repeated.
#'
#' @return List with updated final_output and useful_data variables.
#' @export autorun_improv
#' @import baseline
#' @import robustbase
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' # Not run:
#' # quantification_variables=autorun(imported_data,imported_data$final_output,imported_data$useful_data,imported_data$ROI_data)
#' # quantification_variables_2=autorun_improv(imported_data,quantification_variables$final_output,quantification_variables$useful_data,imported_data$ROI_data,"correction","outliers")


#TODO: Choose criteria to repeat only individual quantification and all signals of all spectra.

autorun_improv = function(imported_data, final_output,useful_data,ROI_data,improvement_option,level) {
print("Starting maximization of profiling quality using information of original profiling...")
predicted_info=rf_pred(final_output$half_band_width)
predicted_width=as.matrix(predicted_info$predicted_matrix)

  max_width=as.matrix(predicted_info$upper_bound_matrix)
  min_width=as.matrix(predicted_info$lower_bound_matrix)
ind=which(is.na(predicted_width[1,]))
predicted_width[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,8]))
min_width[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,8]*0.75))
max_width[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,8]*1.25))

predicted_info=rf_pred(final_output$shift)
predicted_shift=as.matrix(predicted_info$predicted_matrix)
max_shift=as.matrix(predicted_info$upper_bound_matrix)
min_shift=as.matrix(predicted_info$lower_bound_matrix)

  ind=which(is.na(predicted_shift[1,]))
  if (length(ind)>0) {
    predicted_shift[,ind]=as.matrix(t(replicate(nrow(predicted_width),ROI_data[ind,6])))
  max_shift[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,6]+ROI_data[ind,7]))
  min_shift[,ind]=t(replicate(nrow(predicted_width),ROI_data[ind,6]-ROI_data[ind,7]))
  }

  predicted_info=rf_pred_intensity(final_output$intensity,ROI_data[,4])
  predicted_intensity=as.matrix(predicted_info$predicted_matrix)
  max_intensity=as.matrix(predicted_info$upper_bound_matrix)
  min_intensity=as.matrix(predicted_info$lower_bound_matrix)

  ind=which(is.na(predicted_intensity[1,]))
  max_intensity[!is.finite(max_intensity)]=NA
  min_intensity[!is.finite(min_intensity)]=NA
  min_intensity[min_intensity<0]=0

  quantifications_to_repeat=matrix(0,nrow(predicted_width),ncol(predicted_width))
  if (level=="all") quantifications_to_repeat[,]=1
	if (is.numeric(level)) quantifications_to_repeat[which(final_output$fitting_error>level)]=1
	if (level=="outliers") {
pal1=pal2=pal3=matrix(1,nrow(predicted_width),ncol(predicted_width))
for (i in 1:ncol(final_output$shift)) {
indexes=intersect(which(!is.na(predicted_shift[,i])),which(!is.na(final_output$shift[,i])))
if (length(indexes)>0) pal1[indexes,i]=robustbase::lmrob(predicted_shift[indexes,i]~final_output$shift[indexes,i])$rweights
indexes=intersect(which(!is.na(predicted_width[,i])),which(!is.na(final_output$width[,i])))
if (length(indexes)>0) pal2[indexes,i]=robustbase::lmrob(predicted_width[indexes,i]~final_output$half_band_width[indexes,i])$rweights
indexes=intersect(which(!is.na(predicted_intensity[,i])),which(!is.na(final_output$intensity[,i])))
if (length(indexes)>0) pal3[indexes,i]=robustbase::lmrob(predicted_intensity[indexes,i]~final_output$intensity[indexes,i])$rweights
}
pal=pal1+pal2+pal3
for (i in 1:ncol(quantifications_to_repeat)) quantifications_to_repeat[pal[,i] %in% boxplot.stats(pal[,i])$out,i]=1

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


	index_to_use_3=which(rowSums(quantifications_to_repeat[,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2]])>0)

    #Quantification for every spectrum
    for (spectrum_index in index_to_use_3) {

         #Preparation of necessary variables to store figures and information of the fitting
      Ydata = as.numeric(imported_data$dataset[spectrum_index, ROI_buckets])

      #If the quantification is through integration with or without baseline
      if (fitting_type == "Clean Sum" ||
          fitting_type == "Baseline Sum") {
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
        FeaturesMatrix = fitting_prep_2(Xdata,
          Ydata,
          ROI_profile[, 5:11,drop=F],
          program_parameters,baselinedataset[spectrum_index,ROI_buckets],max_shift,min_shift,max_intensity,
          min_intensity,max_width,min_width,spectrum_index,ROI_separator[ROI_index, 1]:ROI_separator[ROI_index, 2])

        #Calculation of the parameters that will achieve the best fitting
        dummy = fittingloop(FeaturesMatrix,
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

      }

  }


 } else if (improvement_option=='correction') {

prova_intensity=predicted_intensity
prova_intensity[,apply(predicted_intensity,2,function(x)all(is.na(x)))]=final_output$intensity[,apply(predicted_intensity,2,function(x)all(is.na(x)))]
prova_shift=predicted_shift
prova_shift[,apply(predicted_shift,2,function(x)all(is.na(x)))]=final_output$shift[,apply(predicted_shift,2,function(x)all(is.na(x)))]
prova_width=predicted_width
prova_width[,apply(predicted_width,2,function(x)all(is.na(x)))]=final_output$half_band_width[,apply(predicted_width,2,function(x)all(is.na(x)))]

tec=sapply(seq(length(prova_intensity)),function(x)sum(peakpvoigt(c(prova_intensity[x],prova_shift[x],prova_width[x]*0.5/600.2,0),imported_data$ppm))*imported_data$buck_step)
dim(tec)=dim(prova_intensity)
tec[,apply(tec,2,function(x)all(is.na(x)))]=final_output$quantification[,apply(tec,2,function(x)all(is.na(x)))]
	for (i in 1:ncol(final_output$quantification)) {
	index_to_use_3=which(quantifications_to_repeat[,i]>0)
	final_output$quantification[index_to_use_3,i]=tec[index_to_use_3,i]
}
 }
  print("Done!")
  quantification_variables=list(final_output=final_output,useful_data=useful_data,
                                predicted_shift=predicted_shift,predicted_width=predicted_width,
                                predicted_intensity=predicted_intensity,max_width=max_width,
                                min_width=min_width,max_shift=max_shift,min_shift=min_shift,
                                max_intensity=max_intensity,min_intensity=min_intensity)
  return(quantification_variables)
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

  #Finding of maximum intensity and shift tolerance of every background signal
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
rf_pred_intensity=function(initial_matrix,lolo) {
  samples2=sample(nrow(initial_matrix),nrow(initial_matrix)*0.5)
  lol2=initial_matrix
  analyzed_signals=which(apply(lol2,2,function(x)! all(is.na(x)))==T)
  lol2=missForest::missForest(lol2[,analyzed_signals])$ximp

  hent=hent2=hent3=as.data.frame(matrix(NA,nrow(lol2),ncol(lol2)))
  stop=0
  while(stop==0) {
    samples2=replicate(30,sample(nrow(lol2),0.4*nrow(lol2)))
    if (all(table(samples2)>=5)) stop=1
  }
  for (i in seq(ncol(lol2))) {
    sed=which(lolo[analyzed_signals]==lolo[analyzed_signals][i])
    if (length(sed)==1) next
    tel=prcomp(scale(lol2[,sed]))$x
    st=matrix(NA,nrow(lol2),30)

    for (j in 1:30) {
      lol3=data.frame(y=lol2[,i],tel)
      model=ranger::ranger(y ~.,lol3[-samples2[,j],])
      st[samples2[,j],j]=predict(model,lol3[samples2[,j],])$predictions
    }
    st=t(sapply(seq(nrow(hent)),function(x)rnorm(1000,mean(st[x,],na.rm=T),sd(st[x,],na.rm=T))))
    hent[,i]=apply(st,1,function(x)quantile(x,0.5,na.rm=T))
    hent2[,i]=apply(st,1,function(x)quantile(x,0.025,na.rm=T))
    hent3[,i]=apply(st,1,function(x)quantile(x,0.975,na.rm=T))
  }




  if (length(analyzed_signals)<ncol(initial_matrix)) {
    predicted_matrix=lower_bound_matrix=upper_bound_matrix=matrix(NA,nrow(initial_matrix),ncol(initial_matrix))
    predicted_matrix[,analyzed_signals]=as.matrix(hent)
    lower_bound_matrix[,analyzed_signals]=as.matrix(hent2)
    upper_bound_matrix[,analyzed_signals]=as.matrix(hent3)

  } else {
    predicted_matrix=hent
    lower_bound_matrix=hent2
    upper_bound_matrix=hent3

  }

  output=list(predicted_matrix=predicted_matrix,lower_bound_matrix=lower_bound_matrix,upper_bound_matrix=upper_bound_matrix)
  return(output)
}
rf_pred=function(initial_matrix) {

  lol2=initial_matrix
  analyzed_signals=which(apply(lol2,2,function(x)! all(is.na(x)))==T)
  lol2=missForest::missForest(lol2[,analyzed_signals])$ximp

  hent=hent2=hent3=as.data.frame(matrix(NA,nrow(lol2),ncol(lol2)))
  tel=prcomp(scale(lol2))$x[,1:20]
  stop=0
  while(stop==0) {
    samples2=replicate(30,sample(nrow(lol2),0.4*nrow(lol2)))
    if (all(table(samples2)>=5)) stop=1
  }
  for (i in seq(ncol(lol2))) {
    st=matrix(NA,nrow(lol2),30)

    for (j in 1:30) {
      lol3=data.frame(y=lol2[,i],tel)
    model=ranger::ranger(y ~.,lol3[-samples2[,j],])
            st[samples2[,j],j]=predict(model,lol3[samples2[,j],])$predictions
    }
    st=t(sapply(seq(nrow(hent)),function(x)rnorm(1000,mean(st[x,],na.rm=T),sd(st[x,],na.rm=T))))
    hent[,i]=apply(st,1,function(x)quantile(x,0.5,na.rm=T))
    hent2[,i]=apply(st,1,function(x)quantile(x,0.025,na.rm=T))
    hent3[,i]=apply(st,1,function(x)quantile(x,0.975,na.rm=T))
  }




  if (length(analyzed_signals)<ncol(initial_matrix)) {
    predicted_matrix=lower_bound_matrix=upper_bound_matrix=matrix(NA,nrow(initial_matrix),ncol(initial_matrix))
    predicted_matrix[,analyzed_signals]=as.matrix(hent)
    lower_bound_matrix[,analyzed_signals]=as.matrix(hent2)
    upper_bound_matrix[,analyzed_signals]=as.matrix(hent3)

  } else {
    predicted_matrix=hent
    lower_bound_matrix=hent2
    upper_bound_matrix=hent3

  }
  output=list(predicted_matrix=predicted_matrix,lower_bound_matrix=lower_bound_matrix,upper_bound_matrix=upper_bound_matrix)
  return(output)
}
