
#' Prediction of signal parameter information with confidence intervals.
#'
#' @param initial_matrix Matrix of signal paramter values.
#' @param fitting_error By default NULL. Fitting error associated to the signal parameter values.
#' @param met_names By default NA. In case of predicting intensity, it enables prediction based  on only signals from the same metabolite. Vector of metabolite associated to each signal whose itensity is predicted.
#' @return List with predicted values and confidence intervals.
#' @export signparpred
#' @import caret
#' @import missRanger
#' @import Boruta
#' @examples
#' # Not run:
#' # imported_data=import_data(file.path(system.file(package = "rDolphin"),"extdata","Parameters_MTBLS242_15spectra_5groups.csv"))
#' # load(file.path(system.file(package = "rDolphin"),"extdata","MTBLS242_subset_profiling_data.RData"))
#' # chemical_shift_pred=signparpred(profiling_data$final_output$chemical_shift,profiling_data$final_output$fitting_error)
#' # intensity_pred=signparpred(profiling_data$final_output$intensity,profiling_data$final_output$fitting_error,imported_data$ROI_data[,4])


signparpred=function(initial_matrix,fitting_error=NULL,met_names=NA) {
  modified_matrix=initial_matrix
  colnames(modified_matrix)=make.names(colnames(modified_matrix))
  modified_matrix=jitter(modified_matrix,0.000001*mean(modified_matrix,na.rm=T))
  dummy=which(apply(modified_matrix,2,function(x) all(is.finite(x)))==T)
  if (!is.null(fitting_error)) {
  for (i in dummy) {
    modified_matrix[fitting_error[,i] %in%
                      boxplot.stats(fitting_error[,i])$out,i]=NA
  }}
  dummy2=missRanger::missRanger(as.data.frame(modified_matrix[,dummy]))
  modified_matrix[,dummy]=as.matrix(dummy2)
  analyzed_signals=which(apply(modified_matrix,2,function(x) all(is.finite(x)))==T)

  processed_matrix=data.frame(modified_matrix[,analyzed_signals],
  prcomp(scale(modified_matrix[,analyzed_signals]))$x[,1:5])
  processed_matrix=predict(caret::preProcess(processed_matrix,c("center","scale","nzv")),processed_matrix)
  predicted_matrix=lower_bound_matrix=upper_bound_matrix=as.data.frame(matrix(NA,nrow(modified_matrix),ncol(modified_matrix)))
  ctrl <- caret::trainControl(method = "boot632",number=18,savePredictions="all")
  for (i in analyzed_signals) {
    if (is.na(met_names)) {
      nam=which(colnames(processed_matrix)==colnames(modified_matrix)[i])
      training_data=data.frame(y=modified_matrix[,i],processed_matrix[,-nam])
      boruta_output <- Boruta::Boruta(y ~ ., data=training_data, doTrace=0)
      ind=colnames(training_data) %in%
        c("y",Boruta::getSelectedAttributes(boruta_output, withTentative = TRUE))
      if (length(which(ind==T))==1) next

      plsFit <- tryCatch(caret::train(y ~ .,data = training_data[,ind],method = "ranger",trControl = ctrl),
                         error=function(e)
                           caret::train(y ~ .,data = training_data[,ind],method = "rf", trControl = ctrl))

    } else {
      sed=intersect(which(met_names==met_names[i]),analyzed_signals)
      if (length(sed)==1) next
      tel=cbind(modified_matrix[,setdiff(sed,i)],prcomp(scale(modified_matrix[,sed]))$x[,1])
      training_data=data.frame(y=modified_matrix[,i],scale(tel))
      boruta_output <- Boruta::Boruta(y ~ ., data=training_data, doTrace=0)
      ind=colnames(training_data) %in%
        c("y",Boruta::getSelectedAttributes(boruta_output, withTentative = TRUE))
      if (!is.null(fitting_error)) {
        if (all(is.na(fitting_error[,i]))) {
        weights=rep(1,nrow(training_data))
      } else {
        weights=1/fitting_error[,i]
      }}
      plsFit <- tryCatch(caret::train(y ~ .,data = training_data[,ind],method = "ranger", weights = weights,trControl = ctrl),
                         error=function(e)
                         caret::train(y ~ .,data = training_data[,ind],method = "rf", weights = weights,trControl = ctrl))
    }


    iii=plsFit$pred
    ff=sapply(seq(nrow(training_data)),function(x)quantile(rnorm(1000,mean=mean(iii$pred[iii$rowIndex==x],na.rm=T),
                                                                 sd=sd(iii$pred[iii$rowIndex==x],na.rm=T)),c(0.025,0.5,0.975),na.rm=T))
    predicted_matrix[,i]=ff[2,]
    lower_bound_matrix[,i]=ff[1,]
    upper_bound_matrix[,i]=ff[3,]
  }


  predicted_matrix=missRanger::missRanger(predicted_matrix)
  lower_bound_matrix=missRanger::missRanger(lower_bound_matrix)
  upper_bound_matrix=missRanger::missRanger(upper_bound_matrix)

  output=list(predicted_matrix=predicted_matrix,lower_bound_matrix=lower_bound_matrix,upper_bound_matrix=upper_bound_matrix)
  return(output)
}
