#' Estimation of quality indicators, with signal to total area ratio, fitting error, and difference between expected and found signal parameter values.
#'
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param alarmmatrix List of previous indicators.
#' @param validation_type Type of validation to perform (1: fitting error, 2: signal total area ratio, 3: expected chemical $chemical_shift, 4: expected half bandwidth, 5: expected relative intensity)
#'
#' @return List with estimated indicators, matrix to show in GUI, and breaks and colors to prepare a different color for each cell of the matrix.
#' @export validation
#' @import ranger
#' @import missRanger
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' load("MTBLS242_subset_profiling_data.RData")
#' validation_data=validation(profiling_data$final_output,profiling_data$alarmmatrix,5)



validation = function(final_output,alarmmatrix,validation_type) {
print("Updating the chosen validation method...")
  sink(type='message');
  #If some required information does not exist yet, it is created
  if (is.null(alarmmatrix)) {
  alarmmatrix=final_output
for (i in seq(length(alarmmatrix))) alarmmatrix[[i]][,]=NA
}
 if (validation_type=="0") validation_type=1
  for (i in 1:length(final_output)) colnames(final_output[[i]])=make.names(colnames(final_output[[i]]))


#Check of which signals have different values
  modified_signals_indexes=sapply(seq(ncol(final_output$signal_area_ratio)),function(x)identical(final_output$signal_area_ratio[,x],alarmmatrix$signal_area_ratio[,x]))
  modified_signals_indexes=which(modified_signals_indexes==F)

  if (length(modified_signals_indexes)>0) {

  #Imputation of missing values
  final_output$half_bandwidth=suppressWarnings(tryCatch(as.matrix(missRanger::missRanger(as.data.frame(final_output$half_bandwidth))),error=function(e) final_output$half_bandwidth))
  final_output$chemical_shift[final_output$chemical_shift==Inf]=NA
  final_output$chemical_shift=suppressWarnings(tryCatch(as.matrix(missRanger::missRanger(as.data.frame(final_output$chemical_shift))),error=function(e) final_output$chemical_shift))
  final_output$intensity=suppressWarnings(tryCatch(as.matrix(missRanger::missRanger(as.data.frame(final_output$intensity))),error=function(e) final_output$intensity))


  #Chemical $chemical_shift value prediction
  valid_indexes=which(apply(final_output$chemical_shift,2,function(x)length(which(is.na(x))))<0.5*nrow(final_output$chemical_shift))
  for (i in intersect(valid_indexes,modified_signals_indexes)) {
    shift3=data.frame(y=final_output$chemical_shift[,i],d=final_output$chemical_shift[,setdiff(valid_indexes,i)])
    alarmmatrix$chemical_shift[,i]=predict(ranger::ranger(y~.,data=shift3,mtry=3),shift3)$predictions-final_output$chemical_shift[,i]
  }

  #Half bandwidth value prediction
  valid_indexes=which(apply(final_output$half_bandwidth,2,function(x)length(which(is.na(x))))<0.5*nrow(final_output$half_bandwidth))
  for (i in intersect(valid_indexes,modified_signals_indexes)) {
    shift3=data.frame(y=final_output$half_bandwidth[,i],d=final_output$half_bandwidth[,setdiff(valid_indexes,i)])
    alarmmatrix$half_bandwidth[,i]=predict(ranger::ranger(y ~.,data=shift3,mtry=3),shift3)$predictions/final_output$half_bandwidth[,i]
  }

  #Intensity value prediction
  valid_indexes=which(apply(final_output$intensity,2,function(x)length(which(is.na(x))))<0.5*nrow(final_output$intensity))
  for (i in intersect(valid_indexes,modified_signals_indexes)) {
    shift3=data.frame(y=final_output$intensity[,i],d=final_output$intensity[,setdiff(valid_indexes,i)])
    alarmmatrix$intensity[,i]=predict(ranger::ranger(y~.,data=shift3,mtry=3),shift3)$predictions/final_output$intensity[,i]
  }

  alarmmatrix$fitting_error=final_output$fitting_error
  alarmmatrix$signal_area_ratio=final_output$signal_area_ratio
}

  #Preparation of matrix to show in GUI and red shade for each cell
  if (validation_type==1) {
    shownmatrix=alarmmatrix$fitting_error
  brks <- seq(0.01,0.19,0.01)
  clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}


} else if (validation_type==2) {
  shownmatrix=alarmmatrix$signal_area_ratio
  brks <- seq(0,20,length.out=19)
  clrs <- round(seq(40, 255, length.out = length(brks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}


} else if (validation_type==3) {
  shownmatrix=alarmmatrix$chemical_shift
  brks <-c(-seq(max(abs(shownmatrix),na.rm=T), 0, length.out=10),seq(0, max(abs(shownmatrix),na.rm=T), length.out=10)[-1])
  clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}

  } else if (validation_type==4) {
    shownmatrix=alarmmatrix$half_bandwidth
    brks <- c(seq(0.25, 1,length.out = 10), seq(1, 4, length.out = 10)[-1])
    clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")}


  } else if (validation_type==5) {
    shownmatrix=alarmmatrix$intensity
     brks <- c(seq(0.25, 1,length.out = 10), seq(1, 4, length.out = 10)[-1])
    clrs <- round(c(seq(40, 255, length.out = (length(brks) + 1)/2),seq(255, 40, length.out = (length(brks) + 1)/2)), 0) %>%
    {paste0("rgb(255,", ., ",", ., ")")}

  }
  alarmmatrix=lapply(alarmmatrix,as.matrix)
  shownmatrix=as.matrix(shownmatrix)
validationdata=list(alarmmatrix=alarmmatrix,shownmatrix=shownmatrix,brks=brks,clrs=clrs)
sink(NULL);
print("Done!")

return(validationdata)
}
