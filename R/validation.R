#' Creation of matrix to validate quantifications, with information of difference with predicted shift, signal to total area ratio, fitting error, and difference with expected relative intensity
#'
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param alarmmatrix List with quantifications and indicators of quality of quantification.
#' @param validation_type Type of valdiation to perform (1: fitting error, 2: signal area ratio, 3: chemical shift, 4: half bandwidth, 5: outliers, 6: relative intensity of signals of same metabolite)
#'
#' @return Matrix with data required
#' @export validation
#' @import ranger
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' load("MTBLS242_subset_example.RData")
#' validation_data=validation(quantification_variables$final_output,quantification_variables$alarmmatrix,5)



validation = function(final_output,alarmmatrix,validation_type) {
print("Updating the chosen validation method...")
if (is.null(alarmmatrix)) {
  alarmmatrix=final_output
for (i in seq(length(alarmmatrix))) alarmmatrix[[i]][,]=NA
}
  tec=apply(final_output$half_band_width,2,function(x)!all(is.na(x)))
  final_output$shift[final_output$shift==Inf]=NA
  final_output$shift=suppressWarnings(tryCatch(missForest::missForest(lol)$ximp,error=function(e) final_output$shift))
  final_output$half_band_width[,tec]=suppressWarnings(tryCatch(missForest::missForest(lol)$ximp,error=function(e) final_output$half_band_width[,tec]))
  final_output$intensity=suppressWarnings(tryCatch(missForest::missForest(lol)$ximp,error=function(e) final_output$intensity))


  indexes=sapply(seq(ncol(final_output$signal_area_ratio)),function(x)identical(final_output$signal_area_ratio[,x],alarmmatrix$signal_area_ratio[,x]))
indexes=which(indexes==F)
  ind=which(apply(final_output$shift,2,function(x)length(which(is.na(x))))<0.5*nrow(final_output$shift))
  for (i in intersect(ind,indexes)) {

    shift3=data.frame(y=final_output$shift[,i],d=final_output$shift[,setdiff(ind,i)])
    alarmmatrix$shift[,i]=predict(ranger::ranger(y~.,data=shift3,mtry=3),shift3)$predictions-final_output$shift[,i]
  }
  ind=which(apply(final_output$half_band_width,2,function(x)length(which(is.na(x))))<0.5*nrow(final_output$half_band_width))
  for (i in intersect(ind,indexes)) {
    shift3=data.frame(y=final_output$half_band_width[,i],d=final_output$half_band_width[,setdiff(ind,i)])
    alarmmatrix$half_band_width[,i]=predict(ranger::ranger(y ~.,data=shift3,mtry=3),shift3)$predictions/final_output$half_band_width[,i]
  }
  ind=which(apply(final_output$intensity,2,function(x)length(which(is.na(x))))<0.5*nrow(final_output$intensity))
  for (i in intersect(ind,indexes)) {
    shift3=data.frame(y=final_output$intensity[,i],d=final_output$intensity[,setdiff(ind,i)])
    alarmmatrix$intensity[,i]=predict(ranger::ranger(y~.,data=shift3,mtry=3),shift3)$predictions/final_output$intensity[,i]
  }

  alarmmatrix$fitting_error=final_output$fitting_error
  alarmmatrix$signal_area_ratio=final_output$signal_area_ratio

  if (validation_type==1) {
    shownmatrix=alarmmatrix$fitting_error
  brks <- seq(0.01,0.19,0.01)
  clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

	#Analysis of which quantifications have too low signal to toal area ratio

} else if (validation_type==2) {
  shownmatrix=alarmmatrix$signal_area_ratio
  brks <- seq(0,20,length.out=19)
  clrs <- round(seq(40, 255, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

	#Analysis of which quantifications deviate too much from expected shift, according to prediction with linear model of signals with similar behavior

} else if (validation_type==3) {
  shownmatrix=alarmmatrix$shift
  brks <-c(-seq(max(abs(shownmatrix),na.rm=T), 0, length.out=10),seq(0, max(abs(shownmatrix),na.rm=T), length.out=10)[-1])
  clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}
  	#Analysis of which quantifications deviate too much from expected half_band_width, according to prediction with linear model of spectra with similar behavior

  } else if (validation_type==4) {
    shownmatrix=alarmmatrix$half_band_width
    brks <- c(seq(0.25, 1,length.out = 10), seq(1, 4, length.out = 10)[-1])
    clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

  #Analysis of outliers for every class and of their magnitude

  } else if (validation_type==5) {
    shownmatrix=alarmmatrix$intensity
     brks <- c(seq(0.25, 1,length.out = 10), seq(1, 4, length.out = 10)[-1])
    clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
    {paste0("rgb(255,", ., ",", ., ")")}

}
validationdata=list(alarmmatrix=alarmmatrix,shownmatrix=shownmatrix,brks=brks,clrs=clrs)
print("Done!")

return(validationdata)
}
