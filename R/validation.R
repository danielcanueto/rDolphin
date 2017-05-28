#' Creation of matrix to validate quantifications, with information of difference with predicted shift, signal to total area ratio, fitting error, and difference with expected relative intensity
#'
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param validation_type Type of valdiation to perform (1: fitting error, 2: signal area ratio, 3: chemical shift, 4: half bandwidth, 5: outliers, 6: relative intensity of signals of same metabolite)
#' @param ROI_data ROIs data
#' @param metadata Dataset metadata
#'
#' @return Matrix with data required
#' @export validation
#' @import randomForest
#' @import robustbase
#'
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' load("MTBLS242_subset_example.RData")
#' validation_data=validation(quantification_variables$final_output,5,imported_data$ROI_data,imported_data$Metadata)



validation = function(final_output,validation_type,ROI_data,metadata) {
print("Updating the chosen validation method...")
if (is.null(final_output)) return(NULL)
    alarmmatrix=matrix(NA,dim(final_output$shift)[1],dim(final_output$shift)[2],dimnames=list(rownames(final_output$shift),colnames(final_output$shift)))

   #Analysis of fitting error of quantifications

  if (validation_type==1) {
  alarmmatrix=final_output$fitting_error
  brks <- seq(0.01,0.19,0.01)
  clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

	#Analysis of which quantifications have too low signal to toal area ratio

} else if (validation_type==2) {
  alarmmatrix=final_output$signal_area_ratio
  brks <- seq(0,20,length.out=19)
  clrs <- round(seq(40, 255, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

	#Analysis of which quantifications deviate too much from expected shift, according to prediction with linear model of signals with similar behavior

} else if (validation_type==3) {
  ind=which(apply(final_output$shift,2,function(x)length(which(is.na(x))))<0.5*nrow(final_output$shift))
  for (i in ind) {
    shift3=data.frame(y=final_output$shift[,i],d=final_output$shift[,setdiff(ind,i)])
    alarmmatrix[,i]=predict(randomForest::randomForest(y~.,data=shift3, importance =TRUE),shift3)
  }
 alarmmatrix=alarmmatrix-final_output$shift

  brks <-c(-seq(max(abs(alarmmatrix),na.rm=T), 0, length.out=10),seq(0, max(abs(alarmmatrix),na.rm=T), length.out=10)[-1])
  clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

  	#Analysis of which quantifications deviate too much from expected half_band_width, according to prediction with linear model of spectra with similar behavior

  } else if (validation_type==4) {
  #   ind=which(apply(final_output$half_band_width,2, function(x) all(is.na(x)))==F)#find signals with quantified half_band_width
  # medianwidth=apply(final_output$half_band_width,2,function(x)median(x,na.rm=T))
  # for (i in 1:dim(final_output$half_band_width)[1]) {
  #    #Create linear model with most similar half_band_width and predict them
  #   lm_similar_spectrum=tryCatch({lmrob(as.numeric(final_output$half_band_width[i,]) ~ medianwidth,control = lmrob.control(maxit.scale=5000))},error= function(e) {lm(as.numeric(final_output$half_band_width[i,]) ~ medianwidth)},warning= function(e) {lm(as.numeric(final_output$half_band_width[i,]) ~ medianwidth)})
  #   prediction_similar_spectrum=suppressWarnings(predict(lm_similar_spectrum, interval='prediction'))
  #   alarmmatrix[i,ind][!is.na(final_output$half_band_width[i,ind])]=final_output$half_band_width[i,ind][!is.na(final_output$half_band_width[i,ind])]-prediction_similar_spectrum[,1]
  # }
    ind=which(apply(final_output$half_band_width,2,function(x)length(which(is.na(x))))<0.5*nrow(final_output$half_band_width))
    for (i in ind) {
      shift3=data.frame(y=final_output$half_band_width[,i],d=final_output$half_band_width[,setdiff(ind,i)])
      alarmmatrix[,i]=predict(randomForest::randomForest(y ~.,data=shift3, importance =TRUE),shift3)
    }
    alarmmatrix=alarmmatrix-final_output$half_band_width

  brks <-c(-seq(max(abs(alarmmatrix),na.rm=T), 0, length.out=10),seq(0, max(abs(alarmmatrix),na.rm=T), length.out=10)[-1])
  clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

  #Analysis of outliers for every class and of their magnitude

  } else if (validation_type==5) {
    ind=which(apply(final_output$intensity,2,function(x)length(which(is.na(x))))<0.5*nrow(final_output$intensity))
    for (i in ind) {
      intensity3=data.frame(y=final_output$intensity[,i],d=final_output$intensity[,setdiff(ind,i)])
      alarmmatrix[,i]=predict(randomForest::randomForest(y~.,data=intensity3, importance =TRUE),intensity3)
    }
    alarmmatrix=alarmmatrix/final_output$intensity

    brks <-c(rev(1/1.1^seq(9)),1,1.1^seq(9))
    clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

    #Analysis of difference with expected intensity comparing with another signal from the same metabolite

  } else if (validation_type==6) {

  #   relative_intensity = ROI_data[,12]
  #
  #   alarmmatrix=final_output$intensity
  #   alarmmatrix[,]=NA
  #   ind=unique(ROI_data[,4][duplicated(ROI_data[,4])])
  #   for (i in seq_along(ind)) {
  #     ab=which(ROI_data[,4]==ind[i])
  #     ab2=ab[which.min(colMeans(final_output$fitting_error[,ab],na.rm=T))]
  #     if (length(ab2)>0) alarmmatrix[,ab]=(final_output$intensity[,ab]/final_output$intensity[,ab2])*relative_intensity[ab]
  #
  # }
    shift=final_output$quantification
    # shift=shift[,-which(is.na(shift[1,]))]
    # shift=data.matrix(shift[,-1])

    for (i in seq(ncol(shift))) {
      mm=rep(NA,ncol(shift))
      for (j in 1:ncol(shift)) mm[j]=tryCatch(suppressWarnings(summary(lmrob(shift[,i]~shift[,j]))$sigma),error=function(e)NaN)
      if (all(is.na(mm))) next
      def=predict(lmrob(shift[,i]~shift[,order(mm)[1:20]],max.it = 1000))
      alarmmatrix[,i]=shift[,i]-def
    }
    brks <-c(-seq(max(abs(alarmmatrix),na.rm=T), 0, length.out=10),seq(0, max(abs(alarmmatrix),na.rm=T), length.out=10)[-1])
    clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
    {paste0("rgb(255,", ., ",", ., ")")}

    # brks <- c(seq(0.25, 1,length.out = 10), seq(1, 4, length.out = 10)[-1])
    # clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
    # {paste0("rgb(255,", ., ",", ., ")")}
  }



validationdata=list(alarmmatrix=alarmmatrix,brks=brks,clrs=clrs)
print("Done!")

return(validationdata)
}
