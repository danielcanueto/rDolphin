#' Creation of matrix to validate quantifications, with information of difference with predicted shift, signal to total area ratio, fitting error, and difference with expected relative intensity
#'
#' @param finaloutput List with quantifications and indicators of quality of quantification.
#' @param validation_type Type of valdiation to perform (1: fitting error, 2: signal area ratio, 3: chemical shift, 4: half bandwidth, 5: outliers, 6: relative intensity of signals of same metabolite)
#' @param ROI_data ROIs data
#' @param metadata Dataset metadata
#'
#' @return Matrix with data required
#' @export validation
#' @import robust
#' @import robustbase
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' quantification_variables=autorun(imported_data,imported_data$finaloutput,imported_data$useful_data)
#' validation_data=validation(quantification_variables$finaloutput,5,imported_data$ROI_data,imported_data$Metadata)
#' DT::datatable(round(validation_data$alarmmatrix,4),selection = list(mode = 'single', target = 'cell')) %>% formatStyle(colnames(validation_data$alarmmatrix), backgroundColor = styleInterval(validation_data$brks, validation_data$clrs))



validation = function(finaloutput,validation_type,ROI_data,metadata) {


    alarmmatrix=matrix(NA,dim(finaloutput$shift)[1],dim(finaloutput$shift)[2],dimnames=list(rownames(finaloutput$shift),colnames(finaloutput$shift)))

   #Analysis of fitting error of quantifications

  if (validation_type==1) {
  alarmmatrix=finaloutput$fitting_error
  brks <- seq(0.01,0.19,0.01)
  clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

	#Analysis of which quantifications have too low signal to toal area ratio

} else if (validation_type==2) {
  alarmmatrix=finaloutput$signal_area_ratio
  brks <- seq(0,20,length.out=19)
  clrs <- round(seq(40, 255, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

	#Analysis of which quantifications deviate too much from expected shift, according to prediction with linear model of signals with similar behavior

} else if (validation_type==3) {

  ind=which(apply(finaloutput$shift,2, function(x) all(is.na(x)))==F) #find quantified signals
  shift_corrmatrix=cor(finaloutput$shift,use='pairwise.complete.obs',method='spearman')

  for (i in ind) {
    similar_signals=finaloutput$shift[,unique(c(i,ind[sort(abs(shift_corrmatrix[,i]),decreasing=T,index.return=T)$ix][1:3]))]
   j=is.na(rowMeans(similar_signals)) #find signals with similar behavior
   #Create linear models with two most simila signals and predict shift
     lm_similar_signals=tryCatch({lmrob(similar_signals[!j,1] ~ similar_signals[!j,2],control = lmrob.control(maxit.scale=5000))},error= function(e) {lm(similar_signals[!j,1] ~ similar_signals[!j,2])},warning= function(e) {lm(similar_signals[!j,1] ~ similar_signals[!j,2])})
    prediction_similar_signal_1=suppressWarnings(predict(lm_similar_signals, interval='prediction'))
    lm_similar_signals=tryCatch({lmrob(similar_signals[!j,1] ~ similar_signals[!j,3],control = lmrob.control(maxit.scale=5000))},error= function(e) {lm(similar_signals[!j,1] ~ similar_signals[!j,3])},warning= function(e) {lm(similar_signals[!j,1] ~ similar_signals[!j,3])})
    prediction_similar_signal_2=suppressWarnings(predict(lm_similar_signals, interval='prediction'))
  alarmmatrix[!j,i]=apply(rbind(finaloutput$shift[!j,i]-prediction_similar_signal_1[,1],finaloutput$shift[!j,i]-prediction_similar_signal_2[,1]),2,min)
  }
  brks <-c(-seq(max(abs(alarmmatrix),na.rm=T), 0, length.out=10),seq(0, max(abs(alarmmatrix),na.rm=T), length.out=10)[-1])
  clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

  	#Analysis of which quantifications deviate too much from expected half_band_width, according to prediction with linear model of spectra with similar behavior

  } else if (validation_type==4) {
    ind=which(apply(finaloutput$half_band_width,2, function(x) all(is.na(x)))==F)#find signals with quantified half_band_width
  medianwidth=apply(finaloutput$half_band_width,2,function(x)median(x,na.rm=T))
  for (i in 1:dim(finaloutput$half_band_width)[1]) {
     #Create linear model with most similar half_band_width and predict them
    lm_similar_spectrum=tryCatch({lmrob(as.numeric(finaloutput$half_band_width[i,]) ~ medianwidth,control = lmrob.control(maxit.scale=5000))},error= function(e) {lm(as.numeric(finaloutput$half_band_width[i,]) ~ medianwidth)},warning= function(e) {lm(as.numeric(finaloutput$half_band_width[i,]) ~ medianwidth)})
    prediction_similar_spectrum=suppressWarnings(predict(lm_similar_spectrum, interval='prediction'))
    alarmmatrix[i,ind][!is.na(finaloutput$half_band_width[i,ind])]=finaloutput$half_band_width[i,ind][!is.na(finaloutput$half_band_width[i,ind])]-prediction_similar_spectrum[,1]
  }
  brks <-c(-seq(max(abs(alarmmatrix),na.rm=T), 0, length.out=10),seq(0, max(abs(alarmmatrix),na.rm=T), length.out=10)[-1])
  clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

  #Analysis of outliers for every class and of their magnitude

  } else if (validation_type==5) {

 metadata_types=unique(metadata[,2])
    for (k in 1:length(metadata_types)) {
    iqr_data=apply(finaloutput$Area[metadata[,2]==metadata_types[k],],2,function(x) IQR(x,na.rm=T))
    quartile_data=rbind(apply(finaloutput$Area[metadata[,2]==metadata_types[k],],2,function(x) quantile(x,0.25,na.rm=T)),apply(finaloutput$Area[metadata[,2]==metadata_types[k],],2,function(x) quantile(x,0.75,na.rm=T)))
    for (i in which(metadata[,2]==metadata_types[k])) {
      for (j in which(!is.na(iqr_data))) {
        if (!is.na(finaloutput$Area[i,j]) && finaloutput$Area[i,j]>quartile_data[1,j]&&finaloutput$Area[i,j]<quartile_data[2,j]) {
          alarmmatrix[i,j]=0
        } else if (!is.na(finaloutput$Area[i,j]) &&finaloutput$Area[i,j]<quartile_data[1,j]) {
          alarmmatrix[i,j]=abs(finaloutput$Area[i,j]-quartile_data[1,j])/iqr_data[j]
        } else if (!is.na(finaloutput$Area[i,j]) &&finaloutput$Area[i,j]>quartile_data[2,j]) {
          alarmmatrix[i,j]=abs(finaloutput$Area[i,j]-quartile_data[2,j])/iqr_data[j]
        }
      }
    }}
  brks <- quantile(alarmmatrix, probs = seq(.05, .95, .05), na.rm = TRUE)
  clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

    #Analysis of difference with expected intensity comparing with another signal from the same metabolite

  } else if (validation_type==6) {

    relative_intensity = ROI_data[,12]

    alarmmatrix=finaloutput$intensity
    alarmmatrix[,]=NA
    ind=unique(ROI_data[,4][duplicated(ROI_data[,4])])
    for (i in seq_along(ind)) {
      ab=which(ROI_data[,4]==ind[i])
      ab2=ab[which.min(colMeans(finaloutput$fitting_error[,ab],na.rm=T))]
      if (length(ab2)>0) alarmmatrix[,ab]=(finaloutput$intensity[,ab]/finaloutput$intensity[,ab2])*relative_intensity[ab]

  }

    brks <-c(seq(min(alarmmatrix,na.rm=T), 1, length.out=10),seq(1, max(alarmmatrix,na.rm=T), length.out=10)[-1])
    clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
    {paste0("rgb(255,", ., ",", ., ")")}
  }



validationdata=list(alarmmatrix=alarmmatrix,brks=brks,clrs=clrs)

return(validationdata)
}
