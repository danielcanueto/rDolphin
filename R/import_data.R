

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
  tryCatch({import_profile = read.delim(
    parameters_path,
    sep = ',',
    header = T,
    stringsAsFactors = F
  )},error=function(cond) {
    message("The Parameters file could not be read.")
    return(NA)})
  import_profile = as.data.frame(sapply(import_profile, function(x)
    gsub("\\\\", "/", x)))

  #Getting the names of experiments, signals and ROIs to quantify and use
  metadata_path = as.character(import_profile[5, 2])

  tryCatch({dummy = read.delim(
    metadata_path,
    sep = ',',
    header = T,
    stringsAsFactors = F
  )},error=function(cond) {
    message("The Metadata file could not be read")
    return(NA)})

  Experiments=as.character(dummy[,1])
  Experiments = as.vector(Experiments[Experiments != ''])
  Metadata=dummy[,-1,drop=F]


  #Import of ROI profiles and generation of names and codes of signals
  tryCatch({ROI_data=read.csv(as.character(import_profile[6, 2]), stringsAsFactors = F)
  ROI_data=ROI_data[order(ROI_data[,6]),1:12]
  },error=function(cond) {
  message("The ROI Profiles file could not be read or rightly processed.")
  return(NA)})


  signals_names=make.names(paste(ROI_data[which(!is.na(ROI_data[, 1])),4],ROI_data[which(!is.na(ROI_data[, 1])),5],sep='_'))
if (any(duplicated(signals_names))==T) {
  print("Revise duplicated signal IDs.")
  return()
}
#Other necessary variables
  freq = as.numeric(as.character(import_profile[10, 2]))
  biofluid=import_profile[12, 2]
  jres_path=as.character(import_profile[13, 2])
  program_parameters=fitting_variables()

  if (as.character(import_profile[14, 2])!="") {
    source(as.character(import_profile[14, 2]))
    program_parameters=fitting_variables()

}
#Creation of repository adapted to biofluid
  repository=data.frame(data.table::fread(file.path(system.file(package = "rDolphin"),"extdata","HMDB_Repository.csv")))
  biofluid_column=which(gsub('.times','',colnames(repository))==biofluid)
  if (length(biofluid_column)==0) {
        times=rowSums(repository[,25:36])
    repository=cbind(repository[order(times,decreasing = T),c(1:3,5:7)],rep(NA,length(times)),times[order(times,decreasing = T)])
    colnames(repository)[c(7,8)]=c('Conc','Times')
  } else {
  repository=repository[repository[,biofluid_column]!=0,]
  repository=repository[order(repository[,biofluid_column],decreasing = T),c(1:3,5:7,biofluid_column+c(-12,0))]
}


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
  expno = as.character(import_profile[2, 2])
  processingno = as.character(import_profile[3, 2])

  #Variables only necessary for reading dataset in csv format
  dataset_path = as.character(import_profile[4, 2])

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
      for (i in 1:dim(imported_data$dataset)[1])   imported_data$dataset[i,so+spectra_lag[i]]=imported_data$dataset[i,so]
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
    tryCatch({imported_data = Metadata2Buckets(Experiments, params,program_parameters$spectrum_borders)
    },error=function(cond) {
      message("The Bruker files could not be read.")
      return(NA)})

    Metadata=Metadata[!Experiments %in% imported_data$not_loaded_experiments,]
    Experiments=Experiments[!Experiments %in% imported_data$not_loaded_experiments]
    norm_factor=imported_data$norm_factor
    notnormalizeddataset=imported_data$dataset*norm_factor


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

  # #Adaptation of data to magnitudes similar to 1. To be removed in the future.
  # secondfactor=quantile(imported_data$dataset,0.9,na.rm=T)
  # imported_data$dataset=imported_data$dataset/secondfactor
  # norm_factor=norm_factor*secondfactor


  imported_data$dataset=imported_data$dataset[,which(apply(imported_data$dataset,2,function(x) all(is.na(x)))==F),drop=F]
  imported_data$dataset[is.na(imported_data$dataset)]=0

  imported_data$ppm=imported_data$ppm[which(!is.na(imported_data$ppm))]

  #Storage of parameters needed to perform the fit in a single variable to return.

  imported_data$buck_step = params$buck_step
  imported_data$Experiments = Experiments
  imported_data$ROI_data = ROI_data
  imported_data$freq = freq
  imported_data$Metadata=Metadata
  imported_data$repository=repository
  imported_data$jres_path=jres_path
  imported_data$program_parameters=program_parameters
  imported_data$norm_factor=norm_factor


  #creation of list with the different final outputs
  dummy=matrix(NaN,nrow(imported_data$dataset),length(signals_names),dimnames=list(imported_data$Experiments,signals_names))
  imported_data$final_output = list(quantification= dummy,signal_area_ratio = dummy,fitting_error = dummy, chemical_shift = dummy,intensity = dummy, half_bandwidth = dummy)

  #creation of list of necessary parameters to load quantifications and evaluate quality of them
  imported_data$reproducibility_data=vector('list',length(imported_data$Experiments))
  for (i in seq_along(imported_data$reproducibility_data)) imported_data$reproducibility_data[[i]]=vector('list',length(signals_names))
  for (i in seq_along(imported_data$reproducibility_data)) {
    for (j in seq_along(imported_data$reproducibility_data[[i]])) {
      imported_data$reproducibility_data[[i]][[j]]=list(Ydata=NULL,Xdata=NULL,ROI_profile=imported_data$ROI_data[j,],program_parameters=NULL,plot_data=NULL,FeaturesMatrix=NULL,signals_parameters=NULL,results_to_save=NULL,error1=1000000)
    }}

print('Done!')
  return(imported_data)

}


Metadata2Buckets <- function(Experiments, params, spectrum_borders) {

  CURRENT = list()
  RAW = list()
  not_loaded_experiments = read_spectra =  norm_factor=c()


  left_spectral_border = ifelse(exists("left_spectral_border", where = params),
                                params$left_spectral_border,
    spectrum_borders[1])
  right_spectral_border = ifelse(exists("right_spectral_border", where = params),
                                 params$right_spectral_border,
    spectrum_borders[2])

  RAW$norm_PEAK_left_ppm = ifelse(params$norm_PEAK == "Y", params$norm_left_ppm,
    0)
  RAW$norm_PEAK_right_ppm = ifelse(params$norm_PEAK == "Y", params$norm_left_ppm,
    0)
  RAW$norm_AREA_left_ppm = ifelse(params$norm_AREA == "Y", params$norm_left_ppm,
                                  0)
  RAW$norm_AREA_right_ppm = ifelse(params$norm_AREA == "Y", params$norm_right_ppm,
                                   0)

  tsp_alignment = params$tsp_alignment
  ref_peak_pos = ifelse(params$peak_alignment == "Y", params$ref_peak_pos,
                        0)

  k2 = 1
  maxspec = length(Experiments)
  pb   <- txtProgressBar(1, maxspec, style=3)

  for (k in 1:maxspec) {

    filename = paste(params$dir, Experiments[k], params$expno, "pdata", params$processingno, sep = "/")
    partname = paste(params$dir, Experiments[k], params$expno, sep = "/")
    storedpars = topspin_read_spectrum2(partname, filename,spectrum_borders[2]-0.1, spectrum_borders[1]+0.1)
    if (all(is.nan(storedpars$real)) == 0) {
      CURRENT$minppm = storedpars$OFFSET - storedpars$SW
      CURRENT$maxppm = storedpars$OFFSET
      CURRENT$step = storedpars$SW / (length(storedpars$real) - 1)
      CURRENT$ppm = seq(CURRENT$maxppm, CURRENT$minppm,-CURRENT$step)

     tmp = (storedpars$real * ((2 ^ storedpars$NC_proc) / storedpars$RG))
	     if (left_spectral_border>CURRENT$maxppm) left_spectral_border=RAW$norm_AREA_left_ppm=floor(CURRENT$maxppm*10)/10
     if (right_spectral_border<CURRENT$minppm) right_spectral_border=RAW$norm_AREA_right=ceiling(CURRENT$minppm*10)/10
           CURRENT$left_spectral_border = 1 + round(-(left_spectral_border -
                                                   CURRENT$maxppm) / CURRENT$step)
      CURRENT$norm_PEAK_left = 1 + round(-(RAW$norm_PEAK_left_ppm -
                                             CURRENT$maxppm) / CURRENT$step)
      CURRENT$norm_PEAK_right = 1 + round(-(RAW$norm_PEAK_right_ppm -
                                              CURRENT$maxppm) / CURRENT$step)
      CURRENT$norm_AREA_left = 1 + round(-(RAW$norm_AREA_left_ppm -
                                             CURRENT$maxppm) / CURRENT$step)
      CURRENT$norm_AREA_right = 1 + round(-(RAW$norm_AREA_right_ppm -
                                              CURRENT$maxppm) / CURRENT$step)
      norm_PEAK_max = max(tmp[CURRENT$norm_PEAK_left:CURRENT$norm_PEAK_right])
      norm_PEAK_max_position = which.max(tmp[CURRENT$norm_PEAK_left:CURRENT$norm_PEAK_right])
      norm_AREA = sum(tmp[CURRENT$norm_AREA_left:CURRENT$norm_AREA_right])
      total_AREA_mean = mean(tmp[1:length(tmp)])

      RAW$len = length(storedpars$real)
      RAW$ofs = storedpars$OFFSET
      RAW$sw = storedpars$SW
      RAW$ncproc = storedpars$NC_proc
      RAW$minppm = storedpars$OFFSET - storedpars$SW
      RAW$maxppm = storedpars$OFFSET
      RAW$size = length(storedpars$real)
      RAW$step = RAW$sw / (RAW$size - 1)
      RAW$ppm = seq(RAW$maxppm, RAW$minppm,-RAW$step)
      RAW$buck_step = ifelse(params$buck_step < RAW$step, RAW$step,
                             params$buck_step)
      RAW$ppm_bucks = seq(left_spectral_border,
                          right_spectral_border,-RAW$buck_step)
      RAW$len_bucks = length(RAW$ppm_bucks)
      RAW$norm_PEAK_max = norm_PEAK_max
      RAW$total_AREA_mean = total_AREA_mean
      RAW$norm_AREA = norm_AREA
      RAW$differential_norm_AREA = 0

      if (params$norm_PEAK == "Y" &
          norm_PEAK_max > 0)
        tmp = tmp / norm_PEAK_max
      if (params$norm_AREA == "Y" &
          norm_AREA > 0)
        tmp = tmp / norm_AREA

      # calculem valors despres de normalitzar
      norm_PEAK_max2 = max(tmp[CURRENT$norm_PEAK_left:CURRENT$norm_PEAK_right])
      norm_AREA2 = sum(tmp[CURRENT$norm_AREA_left:CURRENT$norm_AREA_right])
      # referenciem
      LeftEreticBefore = ifelse(CURRENT$ppm[1] > 11.8, min(CURRENT$ppm[which(CURRENT$ppm > 11.7)]), NaN)

      if (params$glucose_alignment == "Y") {
        JTP = JTPcalibrateToGlucose(tmp, CURRENT$ppm)
      } else if (params$tsp_alignment == "Y") {
        JTP = JTPcalibrateToTSP(tmp, CURRENT$ppm)
      } else if (params$peak_alignment == "Y") {
       JTP = JTPcalibrateToPeak(tmp, CURRENT$ppm, ref_peak_pos)
      } else {
	JTP=list(ppm=CURRENT$ppm)
	}
      CURRENT$ppm = JTP$ppm
      # si hi llegim zona d'eretic llavors tornem a posar Eretic a 11 ppm
      if (!is.nan(LeftEreticBefore)) {
        LeftEreticAfter = min(CURRENT$ppm[which(CURRENT$ppm > 11.7)])
        RightEreticAfter = min(CURRENT$ppm[which(CURRENT$ppm > 10.8)])
        RightEreticBefore = LeftEreticBefore + RightEreticAfter - LeftEreticAfter
        tmp[LeftEreticAfter:RightEreticAfter] = tmp[LeftEreticBefore:RightEreticBefore]
      }



      fill2end = 0
      tmp_count = CURRENT$left_spectral_border
      tmp_buck = rep(0, RAW$len_bucks)
      jumped_bucket = 0
      for (tmp_buck_count in 2:RAW$len_bucks) {
        items = suma = 0
        while (CURRENT$ppm[tmp_count] > RAW$ppm_bucks[tmp_buck_count]) {
          items = items + 1
          suma = suma + tmp[tmp_count]
          tmp_count = tmp_count + 1
          if (tmp_count > length(CURRENT$ppm)) {
            fill2end = 1
            break
          }
        }

	# faig una modificacio per a corregir buckets sense dades.
        if (items > 0) {
          tmp_buck[tmp_buck_count - 1] = suma / items
          while (jumped_bucket > 0) {
            tmp_buck[tmp_buck_count - 1 - jumped_bucket] = suma / items
            jumped_bucket = jumped_bucket - 1
          }
        } else {
          jumped_bucket = jumped_bucket + 1
        }
        if (fill2end == 1) {
          while (tmp_buck_count <= RAW$len_bucks) {
            tmp_buck[tmp_buck_count] = tmp_buck[tmp_buck_count -
                                                  1]
            tmp_buck_count = tmp_buck_count + 1
          }
          break
        }
      }
      while (length(tmp_buck) < RAW$len_bucks)
        tmp_buck = c(tmp_buck, tmp_buck[length(tmp_buck)])


      if (params$norm_PEAK == "Y" &
          norm_PEAK_max > 0) {
        normvalue=norm_PEAK_max
      } else if (params$norm_AREA == "Y" &
          norm_AREA > 0) {
        normvalue=norm_AREA
    } else {
      normvalue=1
    }
      if (k2 == 1) dataset = matrix(NA, 0, RAW$len_bucks)

      dataset=rbind(dataset,tmp_buck[1:RAW$len_bucks])
      norm_factor = append(norm_factor,normvalue)
         read_spectra = append(read_spectra, as.character(Experiments[k]))
         k2 = k2 + 1

         setTxtProgressBar(pb, k)


    } else {
      # llista d'objectes no inclosos
      not_loaded_experiments = append(not_loaded_experiments, Experiments[k])
      print(paste(Experiments[k],"not loaded",sep=" "))
    }
  }
  rownames(dataset) = read_spectra
  colnames(dataset) = as.character(RAW$ppm_bucks)
  finaldata = list(dataset=dataset,ppm=RAW$ppm_bucks,not_loaded_experiments = not_loaded_experiments,norm_factor=norm_factor)


  return(finaldata)
}



JTPcalibrateToGlucose <- function(realSpectra, ppm){
  JTP=list()

  # Locate the target region (n ppm either side of 5.233)
  regionMask = (ppm < 5.733) & (ppm > 4.98)
  # nosaltres tenim la part positiva de l'espectre a l'esquerra
  maskOffset = which(ppm > 5.733)
  maskOffset = maskOffset[length(maskOffset)]

  # Take the approximate second derivative
  realSpectra2 = diff(realSpectra[regionMask], differences = 2)

  # Find the two lowest points, corresponding to the top of the two sharpest
  # peaks.
  min1Index = which.min(realSpectra2) ##ok
  min1Index = maskOffset + min1Index
  while (realSpectra[min1Index+1]>realSpectra[min1Index]) {
    min1Index = min1Index + 1
  }
  while (realSpectra[min1Index-1]>realSpectra[min1Index]) {
    min1Index = min1Index - 1
  }

  peak1PPM = ppm[min1Index];

  # Having found the first peak, flatten it so we can locate the second.
  peak1Mask = (ppm[regionMask] > (peak1PPM - 0.004)) & (ppm[regionMask] < (peak1PPM + 0.004));
  realSpectra2[peak1Mask] = 0

  min2Index = which.min(realSpectra2); ##ok
  min2Index = maskOffset + min2Index
  while (realSpectra[min2Index+1]>realSpectra[min2Index]) {
    min2Index = min2Index + 1
  }
  peak2PPM = ppm[min2Index]

  # Reference to the midpoint of peak 1 and 2.

  JTP$deltaPPM = mean(c(peak1PPM, peak2PPM)) - 5.233
  JTP$ppm = ppm - JTP$deltaPPM
  return(JTP)
}




JTPcalibrateToTSP <- function(realSpectra, ppm){
  JTP=list()
  # Locate the target region
  regionMask = (ppm < 0.2) & (ppm > -0.2)
  # nosaltres tenim la part positiva de l'espectre a l'esquerra
  maskOffset = which(ppm > 0.2)
  maskOffset = maskOffset[length(maskOffset)]

  # Take the approximate second derivative
  TSP_max_position = which.max(realSpectra[regionMask])
  deltaPPM = round(TSP_max_position + maskOffset)
  JTP$deltaPPM = ppm[deltaPPM]

  JTP$ppm = ppm - JTP$deltaPPM;

  return(JTP)
}



JTPcalibrateToPeak2 <- function(realSpectra, ppm, PeakPos, winsize){
  JTP=list()
  # Locate the target region (n ppm either side of PeakPos)
  regionMask = (ppm < PeakPos+winsize) & (ppm > PeakPos-winsize)
  # nosaltres tenim la part positiva de l'espectre a l'esquerra
  maskOffset = which(ppm > PeakPos+winsize)
  maskOffset = maskOffset[length(maskOffset)]

  Peak_max_position = which.max(realSpectra[regionMask])
  deltaPPM = round(Peak_max_position + maskOffset)
  JTP$deltaPPM = ppm[deltaPPM] - PeakPos

  JTP$ppm = ppm - JTP$deltaPPM

  return(JTP)
}




JTPcalibrateToPeak <- function(realSpectra, ppm, PeakPos){
  JTP=list()
  # Locate the target region (n ppm either side of PeakPos)
  regionMask = (ppm < PeakPos+0.1) & (ppm > PeakPos-0.1)
  #regionMask = (ppm < PeakPos+0.05) & (ppm > PeakPos-0.05);
  # nosaltres tenim la part positiva de l'espectre a l'esquerra
  maskOffset = ifelse(length(which(ppm > PeakPos+0.1))>0,length(which(ppm > PeakPos+0.1)),ppm[1])


  # Take the approximate second derivative
  realSpectra2 = diff(realSpectra[regionMask], differences = 2)

  # Find the lowest point, corresponding to the top of the peak.
  min1Index = which.min(realSpectra2) ##ok
  min1Index = maskOffset + min1Index
  while (realSpectra[min1Index+1]>realSpectra[min1Index]) {
    min1Index = min1Index + 1
  }
  while (realSpectra[min1Index+1]>realSpectra[min1Index]) {
    min1Index = min1Index - 1
  }

  JTP$deltaPPM = ppm[min1Index] - PeakPos

  JTP$ppm = ppm - JTP$deltaPPM

  return(JTP)
}


## Internal function for parsing Bruker acquisition files
## inFile - string; directory containing the necessary acquisition files
## params - string; desired parameters to return, if missing will return
##           relevant paramaters for 1D or 2D file
## note: all values are returned as string arguments
## returns values for the designated acquisition parameters
topspin_read_spectrum2 <- function(partname, filename, minppm, maxppm){
  acquPar=parseAcqus(partname)
  pars=parseProcs(filename)
  storedpars <-append(pars, acquPar)
  if (is.null(storedpars$SI)) {
    storedpars$SI     = 0
    storedpars$OFFSET    = 0
    storedpars$SW     = 0
    storedpars$real   = NaN
    storedpars$NC_proc = 0
    storedpars$XDIM = 0
    return(storedpars)
  }

  storedpars$RG = as.numeric(storedpars$RG)
  storedpars$OFFSET=as.numeric(storedpars$OFFSET)
  storedpars$SI=as.numeric(storedpars$SI)
  storedpars$NC_proc=as.numeric(storedpars$NC_proc)
  storedpars$SW=as.numeric(storedpars$SW)
  storedpars$XDIM=as.numeric(storedpars$XDIM)


  minppmindex = floor((storedpars$OFFSET-minppm)/storedpars$SW*(storedpars$SI-1))
  maxppmindex = ceiling((storedpars$OFFSET-maxppm)/storedpars$SW*(storedpars$SI-1))

  if (maxppmindex<0) {
    maxppmindex = 0
    minppmindex = storedpars$SI-1
  }
  realfile = paste(filename, '1r', sep='/')

  if (file.exists(realfile)==0) {
    sprintf('Error: file %s not found', realfile)
    storedpars$SI     = 0
    storedpars$OFFSET    = 0
    storedpars$SW     = 0
    storedpars$real   = NaN
    storedpars$NC_proc = 0
    storedpars$XDIM = 0
    return(storedpars)

  }


  readCon <- file(realfile, 'rb')
  data <- try(readBin(readCon, size=4, what='integer', n=storedpars$SI, endian=storedpars$BYTORDP),
              silent=TRUE)
  storedpars$real=data[(maxppmindex+1):(minppmindex+1)]
  close(readCon)

  storedpars$OFFSET = (storedpars$OFFSET - storedpars$SW/(storedpars$SI-1)*maxppmindex);
  storedpars$SW = (storedpars$OFFSET - storedpars$SW/(storedpars$SI-1)*maxppmindex) - (storedpars$OFFSET - storedpars$SW/(storedpars$SI-1)*minppmindex)
  storedpars$SI=length(storedpars$real)
  return(storedpars)
}





parseAcqus <- function(inDir, params){

  ## Designate parameters if not provided
  if (missing(params))
    params <- c('NC', 'RG', 'OVERFLW', 'NS', 'DATE')
  paramVar <- paste('##$', params, sep='')

  ## Search inDir for necessary acquisition parameter files
  acqus <- list.files(inDir, full.names=TRUE, pattern='^acqus$')[1]
  if (is.na(acqus)) acqus <- list.files(inDir, full.names=TRUE, pattern='^acqu$')[1]
  if (is.na(acqus)) {
    paste('Could not find acquisition parameter files ("acqu" or "acqus")',
          ' in:\n"', inDir, '".', sep='')
    return()
  }
  acqu2s <- list.files(inDir, full.names=TRUE, pattern='^acqu2s')[1]
  if (is.na(acqu2s))
    acqu2s <- list.files(inDir, full.names=TRUE, pattern='^acqu2$')[1]
  if (is.na(acqu2s)) files <- acqus else files <- c(acqus, acqu2s)
  files <- acqus

  ## Search acquisition files for designated parameters
  acquPar <- NULL
  for (i in seq_along(files)){

    ## Determine paramater/value separator
    for (paramSep in c('= ', '=', ' =', ' = ')){
      splitText <- strsplit(readLines(files[i]), paramSep)
      parNames <- sapply(splitText, function(x) x[1])
      parVals <- sapply(splitText, function(x) x[2])
      matches <- match(paramVar, parNames)
      if (any(is.na(matches)))
        next
      else
        break
    }

    ## Return an error if any parameters can not be found
    if (any(is.na(matches)))
      stop(paste('One or more of the following parameters could not be found: ',
                 paste("'", params[which(is.na(matches))], "'", sep='',
                       collapse=', '), ' in:\n"', files[i], sep=''))
    acquPar <- rbind(acquPar, parVals[matches])
  }

  ## Format the data
  colnames(acquPar) <- params
  acquPar <- data.frame(acquPar, stringsAsFactors=FALSE)
  if (!is.null(acquPar$NUC1)){
    for (i in seq_along(acquPar$NUC1)){
      acquPar$NUC1[i] <- unlist(strsplit(unlist(strsplit(acquPar$NUC1[i],
                                                         '<'))[2], '>'))
    }
  }
  if (!is.na(acqu2s))
    rownames(acquPar) <- c('w2', 'w1')
  return(acquPar)
}


## Internal function for parsing Bruker processing files
parseProcs <- function(inDir, params){

  ## Designate parameters if not provided
  if (missing(params))
    params <- c('SW_p','SF','SI','OFFSET','NC_proc','BYTORDP','XDIM')
  paramVar <- paste('##$', params, sep='')

  ## Search inDir for necessary processing parameter files
  procs <- list.files(inDir, full.names=TRUE, pattern='^procs$')[1]
  if (is.na(procs))
    procs <- list.files(inDir, full.names=TRUE, pattern='^proc$')[1]
  if (is.na(procs)) {
    paste('Could not find processing parameter files ("proc" or "procs")',
          ' in:\n"', inDir, '".', sep='')
    return()
  }
  proc2s <- list.files(inDir, full.names=TRUE, pattern='^proc2s$')[1]
  if (is.na(proc2s))
    proc2s <- list.files(inDir, full.names=TRUE, pattern='^proc2$')[1]
  if (is.na(proc2s))
    files <- procs
  else
    files <- c(procs, proc2s)

  ## Search processing files for designated parameters
  pars <- NULL
  for (i in seq_along(files)){

    ## Determine paramater/value separator
    for (paramSep in c('= ', '=', ' =', ' = ')){
      splitText <- strsplit(readLines(files[i]), paramSep)
      parNames <- sapply(splitText, function(x) x[1])
      parVals <- sapply(splitText, function(x) x[2])
      matches <- match(paramVar, parNames)
      if (any(is.na(matches)))
        next
      else
        break
    }

    ## Return an error if any parameters can not be found
    if (any(is.na(matches)))
      stop(paste('One or more of the following parameters could not be found: ',
                 paste("'", params[which(is.na(matches))], "'", sep='',
                       collapse=', '), ' in:\n"', files[i], sep=''))
    pars <- rbind(pars, parVals[matches])
  }

  ## Format the data
  colnames(pars) <- params
  pars <- data.frame(pars, stringsAsFactors=FALSE)
  if (!is.null(pars$BYTORDP))
    pars$BYTORDP <- ifelse(as.numeric(pars$BYTORDP), 'big', 'little')
  if (!is.na(proc2s))
    rownames(pars) <- c('w2', 'w1')
  pars$SW=as.numeric(pars$SW_p)/as.numeric(pars$SF)
  pars=pars[3:8]
  return(pars)
}

