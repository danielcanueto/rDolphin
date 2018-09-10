
#' Prediction of signal parameter information with prediction intervals.
#'
#' @param initial_matrix Matrix of signal paramter values.
#' @param fitting_error By default NULL. Fitting error associated to the signal parameter values.
#' @param met_names By default NA. In case of predicting intensity, it enables prediction based  on only signals from the same metabolite. Vector of metabolite associated to each signal whose itensity is predicted.
#' @return List with predicted values and prediction intervals.
#' @export signparpred
#' @import caret
#' @import missRanger
#' @import Boruta
#' @examples
#' # Not run:
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' load("MTBLS242_subset_profiling_data.RData")
#' # chemical_shift_pred=signparpred(profiling_data$final_output$chemical_shift,profiling_data$final_output$fitting_error)
#' # intensity_pred=signparpred(profiling_data$final_output$intensity,profiling_data$final_output$fitting_error,imported_data$ROI_data[,4])


signparpred=function(initial_matrix,fitting_error=NULL,met_names=NA) {

  ctrl <- caret::trainControl(method = "boot632",savePredictions="final")

  colnames(initial_matrix)=make.names(colnames(initial_matrix))

  initial_matrix=jitter(initial_matrix,0.000001*mean(initial_matrix,na.rm=T))

  dummy=which(apply(initial_matrix,2,function(x) all(is.finite(x)))==T)
  if (!is.null(fitting_error)) {
  for (i in dummy) {
    initial_matrix[fitting_error[,i] %in%
                      boxplot.stats(fitting_error[,i])$out,i]=NA
  }}
  set.seed(1);initial_matrix[,dummy]=as.matrix(missRanger::missRanger(as.data.frame(initial_matrix[,dummy])))

  analyzed_signals=which(apply(initial_matrix,2,function(x) all(is.finite(x)))==T)

  features=data.frame(initial_matrix[,analyzed_signals],
  prcomp(scale(initial_matrix[,analyzed_signals]))$x[,1:5])
  features=predict(caret::preProcess(features,c("center","scale","nzv")),features)

predicted_matrix=lower_bound_matrix=upper_bound_matrix=initial_matrix
predicted_matrix[,]=lower_bound_matrix[,]=upper_bound_matrix[,]=NA


  for (i in analyzed_signals) {
    if (is.na(met_names)) {
      idx=which(colnames(features)==colnames(initial_matrix)[i])
      training_data=data.frame(y=initial_matrix[,i],features[,-idx])
      set.seed(1);boruta_output <- Boruta::Boruta(y ~ ., data=training_data, doTrace=0)
      selected_features=colnames(training_data) %in%
        c("y",Boruta::getSelectedAttributes(boruta_output, withTentative = TRUE))
      if (length(which(selected_features==T))==1) next
      tgrid=expand.grid(.mtry = 0.5*length(which(selected_features==T)),
                        .splitrule="extratrees",
                        .min.node.size=5)
      set.seed(1);model <- tryCatch(caret::train(y ~ .,
                                                  data = training_data[,selected_features],
                                                  method = "ranger",
                                                  tuneGrid=tgrid,
                                                  trControl = ctrl),
                         error=function(e)
                           caret::train(y ~ .,data = training_data[,selected_features],method = "rf", trControl = ctrl))

    } else {
      features=intersect(which(met_names==met_names[i]),analyzed_signals)
      if (length(features)==1) next
      predictors=cbind(initial_matrix[,setdiff(features,i)],prcomp(scale(initial_matrix[,features]))$x[,1])
      training_data=data.frame(y=initial_matrix[,i],scale(predictors))
      set.seed(1);boruta_output <- Boruta::Boruta(y ~ ., data=training_data, doTrace=0)
      selected_features=colnames(training_data) %in%
        c("y",Boruta::getSelectedAttributes(boruta_output, withTentative = TRUE))
      # if (!is.null(fitting_error)) {
      #   if (all(is.na(fitting_error[,i]))) {
      #   weights=rep(1,nrow(training_data))
      # } else {
      #   weights=1/fitting_error[,i]
      # }}
      tgrid=expand.grid(.mtry = length(which(selected_features==T))-1,
                        .splitrule="extratrees",
                        .min.node.size=5)
      set.seed(1);model <- tryCatch(caret::train(y ~ .,
                                                  data = training_data[,selected_features],
                                                  method = "ranger",
                                                  # weights = weights,
                                                  tuneGrid = tgrid,
                                                  trControl = ctrl),
                         error=function(e)
                         caret::train(y ~ .,data = training_data[,selected_features],method = "rf", weights = weights,trControl = ctrl))
    }

    set.seed(1);predictions=sapply(seq(nrow(training_data)),
      function(x)
        quantile(
          rnorm(1000,mean=mean(model$pred$pred[model$pred$rowIndex==x],na.rm=T),
                sd=sd(model$pred$pred[model$pred$rowIndex==x],na.rm=T)),
          c(0.025,0.5,0.975),na.rm=T))
    predicted_matrix[,i]=predictions[2,]
    lower_bound_matrix[,i]=predictions[1,]
    upper_bound_matrix[,i]=predictions[3,]
  }


  set.seed(1);predicted_matrix=as.matrix(missRanger::missRanger(as.data.frame(predicted_matrix)))
  set.seed(1);lower_bound_matrix=as.matrix(missRanger::missRanger(as.data.frame(lower_bound_matrix)))
  set.seed(1);upper_bound_matrix=as.matrix(missRanger::missRanger(as.data.frame(upper_bound_matrix)))

  output=list(predicted_matrix=predicted_matrix,lower_bound_matrix=lower_bound_matrix,upper_bound_matrix=upper_bound_matrix)
  return(output)
}
