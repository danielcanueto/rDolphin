prediction_rf = function(initial_matrix) {

original_matrix=initial_matrix
predicted_matrix=matrix(NA,nrow(original_matrix),ncol(original_matrix))
adjusted_matrix_indicator=apply(original_matrix,2,function(x)sort(table(x),decreasing=TRUE)[1])
original_matrix[,adjusted_matrix_indicator>0.5*nrow(original_matrix)]=NA
analyzed_signals=apply(original_matrix,2,function(x)! all(is.na(x)))
original_matrix[,analyzed_signals]=missForest::missForest(original_matrix[,analyzed_signals])$ximp
subset_matrix=original_matrix
possible_predictors=seq(ncol(original_matrix))
# samples <- sample(NROW(subset_matrix), NROW(subset_matrix) * .5)
samples=seq(NROW(subset_matrix))
ctrl <- trainControl("cv", number = 5)
for (i in 1:ncol(subset_matrix)) {
  iris <- as.data.frame(cbind(subset_matrix[,i],subset_matrix[,-i]))
  colnames(iris)=paste("V",seq(ncol(iris)),sep="_")
rffit <- train(V_1~.,data=iris, method = "rf", trControl = ctrl, tuneLength = 5)
  predicted_matrix[,i] <- predict(rffit, iris)
}

return(predicted_matrix)
}
