

#' Import of variables stored in the parameters file and of the dataset to quantify
#'
#' @param parameters_path CSV file woth input paths and desired pre-processing steps
#'
#' @return Imported data of experiment
#' @export import_data
#' @import data.table
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")

import_data = function(parameters_path) {
  #Created by Daniel Canueto 30/08/2016
  print('Importing necessary data to begin the profiling. The more spectra and narrower bucketing, the more time I need.')


  #List of parameters to use to create the dataset
  params = list()

  #Import of parameters from the csv file
  import_profile = read.delim(
    parameters_path,
    sep = ',',
    header = T,
    stringsAsFactors = F
  )
  import_profile = as.data.frame(sapply(import_profile, function(x)
    gsub("\\\\", "/", x)))

  #Getting the names of experiments, signals and ROIs to quantify and use
  metadata_path = as.character(import_profile[3, 2])

  dummy = read.delim(
    metadata_path,
    sep = ',',
    header = T,
    stringsAsFactors = F
  )
  Experiments=as.character(dummy[,1])
  Experiments = as.vector(Experiments[Experiments != ''])
  Metadata=dummy[,-1,drop=F]

  #Import of ROI profiles and generation of names and codes of signals
  ROI_data=try(read.csv(as.character(import_profile[6, 2]), stringsAsFactors = F),silent=T)
  signals_names=paste(ROI_data[which(!is.na(ROI_data[, 1])),4],ROI_data[which(!is.na(ROI_data[, 1])),5],sep='_')
  signals_codes = 1:length(signals_names)


#Other necessary variables
  freq = as.numeric(as.character(import_profile[10, 2]))
  biofluid=import_profile[12, 2]
  jres_path=as.character(import_profile[13, 2])
  try(source(as.character(import_profile[14, 2])),silent=T)
  if (!exists("program_parameters")) program_parameters=fitting_variables()

#Creation of repository adapted to biofluid
  repository=data.frame(data.table::fread(file.path(system.file(package = "rDolphin"),"extdata","HMDB_Repository.csv")))
  biofluid_column=which(gsub('.conc.','',colnames(repository))==biofluid)
  repository=repository[repository[,biofluid_column]!=0,]
  repository=repository[order(repository[,biofluid_column],decreasing = T),c(1:3,5:7,biofluid_column+c(0,12))]



  #Kind of normalization
  normalization = import_profile[7, 2]
  pqn='N'

  params$norm_AREA = 'N'
  params$norm_PEAK = 'N'
  params$norm_left_ppm = program_parameters$spectrum_borders[1]
  params$norm_right_ppm = program_parameters$spectrum_borders[2]
  if (normalization == 1) {
    #Eretic
    params$norm_AREA = 'Y'
    params$norm_left_ppm = 11.53
    params$norm_right_ppm = 10.47
  } else if (normalization == 2) {
    #TSP
    params$norm_AREA = 'Y'
    params$norm_left_ppm = 0.1
    params$norm_right_ppm = -0.1
  } else if (normalization == 3) {
    #Creatinine (intensity, not area, maybe dangerous for rats because of oxalacetate)
    params$norm_PEAK = 'Y'
    params$norm_left_ppm = 3.10
    params$norm_right_ppm = 3
  } else if (normalization == 4) {
    #Spectrum AreA
    params$norm_AREA = 'Y'
  } else if (normalization == 0) {
    #No normailzation

  } else if (normalization == 5) {
    #PQN normailzation
    params$norm_AREA = 'Y'
    pqn='Y'
  }

  #Alignment
  alignment = import_profile[8, 2]
  params$glucose_alignment = 'N'
  params$tsp_alignment = 'N'
  params$peak_alignment = 'N'
  params$ref_peak_pos = 8.452
  if (alignment == 1) {
    #Glucose
    params$glucose_alignment = 'Y'
  } else if (alignment == 2) {
    #TSP
    params$tsp_alignment = 'Y'
  } else if (alignment == 3) {
    #Formate
    params$peak_alignment = 'Y'
  }

  #Suppresion regions
  suppression = as.character(import_profile[9, 2])
  if (suppression == '') {
    params$disol_suppression = 'N'
  } else {
    params$disol_suppression = 'Y'
    params$disol_suppression_ppm = as.numeric(strsplit(suppression, '-|;')[[1]])
    params$disol_suppression_ppm = matrix(params$disol_suppression_ppm,length(params$disol_suppression_ppm) /2,2,byrow = T)
  }

  #Variables only necessary for reading Bruker files
  bruker_path = import_profile[1, 2]
  expno = as.character(import_profile[4, 2])
  processingno = as.character(import_profile[5, 2])

  #Variables only necessary for reading dataset in csv format
  dataset_path = as.character(import_profile[2, 2])

  #If data comes from csv dataset
  if (bruker_path == '' || expno == '' || processingno == '') {
    if (dataset_path != '') {
      #Reading of dataset file (ideally with data.table::fread of data.table package, but seems that the package is not compatible with R 3.3.1). Maybe now it is possible.
      imported_data = list()
      dummy = data.matrix(data.table::fread(dataset_path, sep = ',',header=F))
      notnormalizeddataset=imported_data$dataset=dummy[-1,]
      imported_data$dataset[is.na(imported_data$dataset)]=0
	  imported_data$ppm = round(dummy[1,],4)
      if (imported_data$ppm[1]<imported_data$ppm[2]) {
        imported_data$dataset=t(apply(imported_data$dataset,1,rev))
        imported_data$ppm=rev(imported_data$ppm)
      }
	  colnames(imported_data$dataset) = imported_data$ppm
	  rownames(imported_data$dataset) = Experiments
	  params$buck_step = ifelse(
	    as.character(import_profile[11, 2]) == '',
	    abs(imported_data$ppm[1] - imported_data$ppm[length(imported_data$ppm)]) /
	      length(imported_data$ppm),
	    as.numeric(as.character(import_profile[11, 2]))
	  )

	#TODO: revise alignment and normalization when coming data from csv.
      if (alignment == 1) {
        #Glucose
        lmn=apply(imported_data$dataset,1,function(x)JTPcalibrateToGlucose(x,imported_data$ppm)$deltaPPM)
        } else if (alignment == 2) {
        #TSP
          lmn=apply(imported_data$dataset,1,function(x)JTPcalibrateToTSP(x,imported_data$ppm)$deltaPPM)
          } else if (alignment == 3) {
        #Formate
            lmn=apply(imported_data$dataset,1,function(x)JTPcalibrateToPeak(x,imported_data$ppm,params$ref_peak_pos)$deltaPPM)
          }
            if (alignment!=0&&nrow(imported_data$dataset)>1) {

              imported_data$ppm=imported_data$ppm-median(lmn)

              spectra_lag=round((lmn-median(lmn))/params$buck_step)

      so=(1+max(abs(spectra_lag))):(length(imported_data$ppm)-max(abs(spectra_lag)))
      for (i in 1:dim(imported_data$dataset)[1])   imported_data$dataset[i,so-spectra_lag[i]]=imported_data$dataset[i,so]
            }
	  norm_factor=rep(1,nrow(imported_data$dataset))
      if (params$norm_AREA == 'Y') {
        # for (i in 1:dim(imported_data$dataset)[1])
        #   norm_factor[i]=mean(rowSums(imported_data$dataset[,which.min(abs(imported_data$ppm-params$norm_left_ppm)):which.min(abs(imported_data$ppm-params$norm_right_ppm))]))/sum(imported_data$dataset[i,which.min(abs(imported_data$ppm-params$norm_left_ppm)):which.min(abs(imported_data$ppm-params$norm_right_ppm))])
        norm_factor=rowSums(imported_data$dataset[,which.min(abs(imported_data$ppm-params$norm_left_ppm)):which.min(abs(imported_data$ppm-params$norm_right_ppm))])
      } else if (params$norm_PEAK == 'Y') {
        norm_factor=apply(imported_data$dataset[,which.min(abs(imported_data$ppm-params$norm_left_ppm)):which.min(abs(imported_data$ppm-params$norm_right_ppm))],1,max)
      }
	  imported_data$dataset=imported_data$dataset/norm_factor
    } else {
      print('Problem when creating the dataset. Please revise the parameters.')
      return()
    }

  } else {

    #Reading of Bruker files
    params$dir = bruker_path
    params$expno = expno
    params$processingno = processingno
    params$buck_step = as.numeric(as.character(import_profile[11,2]))
    imported_data = Metadata2Buckets(Experiments, params,program_parameters$spectrum_borders)
    Metadata=Metadata[!Experiments %in% imported_data$not_loaded_experiments,]
    Experiments=Experiments[!Experiments %in% imported_data$not_loaded_experiments]
    norm_factor=imported_data$norm_factor
    notnormalizeddataset=imported_data$dataset*norm_factor


    # if (dim(imported_data$dataset)==2) dummy=NA
  }

  imported_data$dataset[is.na(imported_data$dataset)]=min(abs(imported_data$dataset)[abs(imported_data$dataset)>0])

  #Region suppression
    if (params$disol_suppression == 'Y') {
      ind=c()
      for (i in seq(nrow(params$disol_suppression_ppm))) ind=c(ind,which.min(abs(imported_data$ppm-params$disol_suppression_ppm[i,1])):which.min(abs(imported_data$ppm-params$disol_suppression_ppm[i,2])))
      imported_data$dataset=imported_data$dataset[,-ind,drop=F]
      imported_data$ppm=imported_data$ppm[-ind]
    }

# Possibility of removing zones without interesting information. Added in the future.
 #snr=apply(imported_data$dataset,1,function(x)stats::mad(x,na.rm=T))

  # ind=seq(1,ncol(imported_data$dataset),round(ncol(imported_data$dataset)/(0.1/params$buck_step)))
  # flan=rep(NA,length(ind))
  # for (i in 1:length(ind)) flan[i]=tryCatch(colMeans(imported_data$dataset[,ind[i]:(ind[i]+1)]), error=function(e) NA)
  #
  # snr=apply(imported_data$dataset[,ind[which.min(flan)]:ind[which.min(flan)+1]],1,function(x)stats::mad(x,na.rm=T))
  #
  # dfg=matrix(0,nrow(imported_data$dataset),ncol(imported_data$dataset))
  # for (i in 1:nrow(imported_data$dataset)) {
  #   dfg[i,which(imported_data$dataset[i,]>snr[i])]=1
  # }
  #
  # dfi=which(apply(dfg,2,sum)>0.5*nrow(imported_data$dataset))
  # dfj=c()
  # for (i in 1:length(dfi)) dfj=unique(c(dfj,round((dfi[i]-0.02/params$buck_step):(dfi[i]+0.02/params$buck_step))))
  # dfj = dfj[dfj > 0 & dfj <= ncol(imported_data$dataset)]
  #imported_data$dataset=imported_data$dataset[,dfj,drop=F]
  #imported_data$ppm=imported_data$ppm[dfj]

  #If pqn is desired
  #TODO: specify which are the control samples or if there are no control samples.
  if (pqn=='Y'&&nrow(imported_data$dataset)>1) {
    treated=t(imported_data$dataset[,which(apply(imported_data$dataset,2,median)>median(apply(imported_data$dataset,2,median))),drop=F])
    reference <- apply(treated,1,function(x)median(x,na.rm=T))
    quotient <- treated/reference
    quotient.median <- apply(quotient,2,function(x)median(x,na.rm=T))
    imported_data$dataset <- imported_data$dataset/quotient.median
    norm_factor=norm_factor*quotient.median

  }

  #Adaptation of data to magnitudes similar to 1. To be removed in the future.
  secondfactor=quantile(imported_data$dataset,0.9,na.rm=T)
  imported_data$dataset=imported_data$dataset/secondfactor
  norm_factor=norm_factor*secondfactor


  imported_data$dataset=imported_data$dataset[,which(apply(imported_data$dataset,2,function(x) all(is.na(x)))==F),drop=F]
  imported_data$dataset[is.na(imported_data$dataset)]=0

  imported_data$ppm=imported_data$ppm[which(!is.na(imported_data$ppm))]
  # if (imported_data$ppm[1]<imported_data$ppm[2]) {
    # imported_data$ppm=rev(imported_data$ppm)
    # imported_data$dataset=t(apply(imported_data$dataset,1,rev))
  # }
  #Storage of parameters needed to perform the fit in a single variable to return.

  imported_data$buck_step = params$buck_step
  imported_data$metadata_path = metadata_path
  imported_data$parameters_path = parameters_path
  imported_data$signals_names = signals_names
  imported_data$signals_codes = signals_codes
  imported_data$Experiments = Experiments
  imported_data$ROI_data = ROI_data
  imported_data$freq = freq
  imported_data$Metadata=Metadata
  imported_data$repository=repository
  imported_data$jres_path=jres_path
  imported_data$program_parameters=program_parameters
  imported_data$norm_factor=norm_factor


  #creation of list with the different final outputs
  dummy=matrix(NaN,nrow(imported_data$dataset),length(imported_data$signals_names),dimnames=list(imported_data$Experiments,imported_data$signals_names))
  imported_data$final_output = list(quantification= dummy,signal_area_ratio = dummy,fitting_error = dummy, shift = dummy,intensity = dummy, half_band_width = dummy)

  #creation of list of necessary parameters to load quantifications and evaluate quality of them
  imported_data$useful_data=vector('list',length(imported_data$Experiments))
  for (i in seq_along(imported_data$useful_data)) imported_data$useful_data[[i]]=vector('list',length(imported_data$signals_codes))
  for (i in seq_along(imported_data$useful_data)) {
    for (j in seq_along(imported_data$useful_data[[i]])) {
      imported_data$useful_data[[i]][[j]]=list(Ydata=NULL,Xdata=NULL,ROI_profile=NULL,program_parameters=NULL,plot_data=NULL,FeaturesMatrix=NULL,signals_parameters=NULL,results_to_save=NULL,error1=1000000)
    }}

	# dummy = which(is.na(imported_data$ROI_data[, 1]))
    # if (length(dummy)==0) dummy=dim(imported_data$ROI_data)[1]+1
    # lal=which(duplicated(imported_data$ROI_data[-dummy,1:2])==F)
    # imported_data$ROI_separator = cbind(lal, c(lal[-1] - 1, dim(imported_data$ROI_data[-dummy,])[1]))
  # export_path=file.path(dirname(as.character(import_profile[6, 2])),"input_data")
  #Useful data about conditions of import of data. TO BE REARRANGED
  # dir.create(export_path,showWarnings = FALSE)
  # fwrite(as.data.frame(imported_data$params),file=file.path(imported_data$export_path, 'initial_params.csv'))
  # fwrite(as.data.frame(imported_data$dataset),file=file.path(export_path, 'initial_dataset.csv'))
  # fwrite(as.data.frame(notnormalizeddataset),file=file.path(export_path, 'notnormalizeddataset.csv'))
  # fwrite(as.data.frame(norm_factor),file=file.path(export_path, 'norm_factor.csv'))
  # if ("not_loaded_experiments" %in% names(imported_data)&&length(imported_data$not_loaded_experiments)>0)
  #   fwrite(as.data.frame(imported_data$not_loaded_experiments),file=file.path(export_path, 'not_loaded_experiments.csv'),col.names = F)
  print('Done!')
  return(imported_data)

}
