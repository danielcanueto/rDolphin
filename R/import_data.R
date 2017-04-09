

#' Import of variables stored in the parameters file and of the dataset to quantify
#'
#' @param parameters_path CSV file woth input paths and desired pre-processing steps
#'
#' @return Imported data of experiment
#' @export import_data
#' @import rio
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")

import_data = function(parameters_path) {
  #Created by Daniel Canueto 30/08/2016
  print('Importing necessary data to begin the profiling. The more spectra and narrower bucketing, the more time I need.')


  #List of parameters to use to create the dataset
  params = list()

  #Import fo parameters from the csv file
  # TO DO: stringsasfactors=F
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
  Experiments=dummy[,1]
  Experiments = as.vector(Experiments[Experiments != ''])
  Metadata=dummy[,-1,drop=F]
  biofluid=import_profile[12, 2]
  profile_folder_path = as.character(import_profile[6, 2])

  ROI_data=try(read.csv(profile_folder_path, stringsAsFactors = F),silent=T)
  # if (class(ROI_data)=='try-error') {
# if  (biofluid=='Urine'){ ROI_data=read.csv('C:/Users/user/Documents/Dolphin/R/Urine_library.csv', stringsAsFactors = F) } else if (biofluid=='Blood') { ROI_data=read.csv('C:/Users/user/Documents/Dolphin/R/Blood_library.csv', stringsAsFactors = F) } else {
  # ROI_data=read.csv('C:/Users/user/Documents/Dolphin/R/Urine_library.csv', stringsAsFactors = F)
  # print('Loading urine library. Please prepare a library adapted to your dataset.')
# }}
  signals_names=paste(ROI_data[which(!is.na(ROI_data[, 1])),4],ROI_data[which(!is.na(ROI_data[, 1])),5],sep='_')

  signals_codes = 1:length(signals_names)



  #Preparing the structure of experiments and signals where to store the output
  export_path = dirname(parameters_path)
#Other necessary variables
  freq = as.numeric(as.character(import_profile[10, 2]))

  jres_path=as.character(import_profile[13, 2])
  try(source(as.character(import_profile[14, 2])),silent=T)

  if (!exists("program_parameters")) program_parameters=fitting_variables()


  repository=as.data.frame(rio::import(file.path(system.file(package = "Dolphin"),"extdata","HMDB_Repository.csv")))
  der=which(colnames(repository)==biofluid)
  fa=repository[,der]
  # repository=repository[repository[,der]>0&!is.na(repository[,der])&!is.nan(repository[,der]),]
  repository=repository[!is.na(fa),]
  fa=repository[,der]

  repository=repository[fa>0,]

  fa=repository[,der]

  repository=repository[sort(fa,decreasing = T,index.return=T)$ix,c(1:3,5:7,der)]



  #Kind of normalization
  #TO DO: add PQN (but before standardize a way to find the regions to have into account)
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
  } else if (normalization == 5) {
    #No normailzation

  } else if (normalization == 6) {
    #No normailzation
    params$norm_AREA = 'Y'
    pqn='Y'
  }

  #Alignment
  alignment = import_profile[8, 2]
  params$glucose_alignment = 'N'
  params$tsp_alignment = 'N'
  params$peak_alignment = 'N'
  params$ref_pos = 8.452
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

  if (bruker_path == '' || expno == '' || processingno == '') {
    if (dataset_path != '') {
      #Reading of dataset file (ideally with fread of data.table package, bu seems that the package is not compatible with R 3.3.1)
      imported_data = list()
      dummy = rio::import(dataset_path, sep = ',',header=F,colClasses='numeric')
      pa=dim(dummy[-1,])

      imported_data$dataset=as.numeric(as.matrix(dummy[-1,]))
      imported_data$dataset[is.na(imported_data$dataset)]=0
      dim(imported_data$dataset)=pa
      colnames(imported_data$dataset) = dummy[1,]
      imported_data$ppm = round(as.numeric(dummy[1,]),4)
      rownames(imported_data$dataset) = Experiments
      if (imported_data$ppm[1]<imported_data$ppm[2]) {
        imported_data$dataset=t(apply(imported_data$dataset,1,rev))
        imported_data$ppm=rev(imported_data$ppm)
      }
      if (alignment == 1) {
        #Glucose
        limi=c(5.5,5.1)
      } else if (alignment == 2) {
        #TSP
        limi=c(0.1,-0.1)
      } else if (alignment == 3) {
        #Formate
        limi=c(8.48,8.42)
      }
      if (alignment!=4&&nrow(imported_data$dataset)>1) {
        # alignment_ind=which.min(abs(imported_data$ppm-limi[1])):which.min(abs(imported_data$ppm-limi[2]))
      spectra_lag=rep(NA,nrow(imported_data$dataset))
      for (i in 1:dim(imported_data$dataset)[1]) {
        d <-
          ccf(imported_data$dataset[i, ],
            apply(imported_data$dataset[, which.min(abs(imported_data$ppm-limi[1])):which.min(abs(imported_data$ppm-limi[2]))], 2, median),
            type = 'covariance',
            plot = FALSE)
        spectra_lag[i]=d$lag[which.max(d$acf)]
      }
      so=(1+max(abs(spectra_lag))):(length(imported_data$ppm)-max(abs(spectra_lag)))
      for (i in 1:dim(imported_data$dataset)[1])   imported_data$dataset[i,so-spectra_lag[i]]=imported_data$dataset[i,so]
      }
      if (params$norm_AREA == 'Y') {
        for (i in 1:dim(imported_data$dataset)[1])
          imported_data$dataset[i,]=imported_data$dataset[i,]*mean(rowSums(imported_data$dataset[,which.min(abs(imported_data$ppm-params$norm_left_ppm)):which.min(abs(imported_data$ppm-params$norm_right_ppm))]))/sum(imported_data$dataset[i,which.min(abs(imported_data$ppm-params$norm_left_ppm)):which.min(abs(imported_data$ppm-params$norm_right_ppm))])
      } else if (params$norm_PEAK == 'Y') {
        for (i in 1:dim(imported_data$dataset)[1])
          imported_data$dataset[i,]=imported_data$dataset[i,]*mean(apply(imported_data$dataset[,which.min(abs(imported_data$ppm-params$norm_left_ppm)):which.min(abs(imported_data$ppm-params$norm_right_ppm))],1,max))/sum(imported_data$dataset[i,which.min(abs(imported_data$ppm-params$norm_left_ppm)):which.min(abs(imported_data$ppm-params$norm_right_ppm))])
      }


      params$buck_step = ifelse(
        as.character(import_profile[11, 2]) == '',
        abs(imported_data$ppm[1] - imported_data$ppm[length(imported_data$ppm)]) /
          length(imported_data$ppm),
        as.numeric(as.character(import_profile[11, 2]))
      )
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
    if (dim(imported_data$dataset)==2) dummy=NA
  }

  imported_data$dataset[is.na(imported_data$dataset)]=min(abs(imported_data$dataset)[abs(imported_data$dataset)>0])


  #snr=apply(imported_data$dataset,1,function(x)stats::mad(x,na.rm=T))

    if (params$disol_suppression == 'Y') {
      ind=c()
      for (i in seq(nrow(params$disol_suppression_ppm))) ind=c(ind,which.min(abs(imported_data$ppm-params$disol_suppression_ppm[i,1])):which.min(abs(imported_data$ppm-params$disol_suppression_ppm[i,2])))
      imported_data$dataset=imported_data$dataset[,-ind,drop=F]
      imported_data$ppm=imported_data$ppm[-ind]
    }

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
  if (pqn=='Y'&&nrow(imported_data$dataset)>1) {
    treated=t(imported_data$dataset[,which(apply(imported_data$dataset,2,median)>median(apply(imported_data$dataset,2,median))),drop=F])
    reference <- apply(treated,1,function(x)median(x,na.rm=T))
    quotient <- treated/reference
    quotient.median <- apply(quotient,2,function(x)median(x,na.rm=T))
    imported_data$dataset <- imported_data$dataset/quotient.median
  }

  imported_data$dataset=imported_data$dataset/quantile(imported_data$dataset,0.9,na.rm=T)

  imported_data$dataset=  imported_data$dataset[,which(apply(imported_data$dataset,2,function(x) all(is.na(x)))==F),drop=F]

  imported_data$ppm=imported_data$ppm[which(!is.na(imported_data$ppm))]
  if (imported_data$ppm[1]<imported_data$ppm[2]) {
    imported_data$ppm=rev(imported_data$ppm)
    imported_data$dataset=t(apply(imported_data$dataset,1,rev))
  }
  imported_data$dataset[is.na(imported_data$dataset)]=0
  #Storage of parameters needed to perform the fit in a single variable to return.

  imported_data$buck_step = params$buck_step
  imported_data$profile_folder_path = profile_folder_path
  imported_data$metadata_path = metadata_path
  imported_data$parameters_path = parameters_path
  imported_data$signals_names = signals_names
  imported_data$signals_codes = signals_codes
  imported_data$Experiments = setdiff(Experiments, imported_data$not_loaded_experiments)
  imported_data$ROI_data = ROI_data
  imported_data$freq = freq
  imported_data$Metadata=Metadata
  imported_data$repository=repository
  imported_data$jres_path=jres_path
  imported_data$program_parameters=program_parameters
  imported_data$export_path=export_path
  colnames(imported_data$dataset) = imported_data$ppm
  rownames(imported_data$dataset) = imported_data$Experiments

  #creation of list with the different final outputs
  dummy=matrix(NaN,nrow(imported_data$dataset),length(imported_data$signals_names),dimnames=list(imported_data$Experiments,imported_data$signals_names))
  imported_data$finaloutput = list(Area= dummy,signal_area_ratio = dummy,fitting_error = dummy, shift = dummy,intensity = dummy, half_band_width = dummy)

  #creation of list of necessary parameters to load quantifications and evaluate quality of them
  imported_data$useful_data=vector('list',length(imported_data$Experiments))
  for (i in seq_along(imported_data$useful_data)) imported_data$useful_data[[i]]=vector('list',length(imported_data$signals_codes))
  for (i in seq_along(imported_data$useful_data)) {
    for (j in seq_along(imported_data$useful_data[[i]])) {
      imported_data$useful_data[[i]][[j]]=list(Ydata=NULL,Xdata=NULL,ROI_profile=NULL,program_parameters=NULL,plot_data=NULL,FeaturesMatrix=NULL,signals_parameters=NULL,results_to_save=NULL,error1=1000000)
    }}

	dummy = which(is.na(imported_data$ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(imported_data$ROI_data)[1]+1
    lal=which(duplicated(imported_data$ROI_data[-dummy,1:2])==F)
    imported_data$ROI_separator = cbind(lal, c(lal[-1] - 1, dim(imported_data$ROI_data[-dummy,])[1]))

  #Useful data about conditions of import of data. TO BE REARRANGED
  dir.create(imported_data$export_path, showWarnings = FALSE)
  write.csv(
    as.data.frame(imported_data$params),
    file.path(imported_data$export_path, 'initial_params.csv'),
    row.names = F
  )
  write.csv(imported_data$dataset,
    file.path(imported_data$export_path, 'initial_dataset.csv'),row.names=F)
  if ("not_loaded_experiments" %in% names(imported_data))
  write.table(
      imported_data$not_loaded_experiments,
      file.path(imported_data$export_path, 'not_loaded_experiments.csv'),
      row.names = F,
      col.names = F
    )
  print('Done!')
  return(imported_data)

}
