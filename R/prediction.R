prediction = function(initial_matrix) {

original_matrix=initial_matrix
predicted_matrix=matrix(NA,nrow(original_matrix),ncol(original_matrix))
adjusted_matrix_indicator=apply(original_matrix,2,function(x)sort(table(x),decreasing=TRUE)[1])
original_matrix[,adjusted_matrix_indicator>0.5*nrow(original_matrix)]=NA
analyzed_signals=apply(original_matrix,2,function(x)! all(is.na(x)))
original_matrix[,analyzed_signals]=missForest::missForest(original_matrix[,analyzed_signals])$ximp
subset_matrix=original_matrix
possible_predictors=seq(ncol(original_matrix))
for (i in 1:ncol(subset_matrix)) {
  if (all(is.na(subset_matrix[,i]))) next
  r2_all_rlr=rep(NA,ncol(original_matrix))
  for (j in possible_predictors) r2_all_rlr[j]=tryCatch({summary(robustbase::lmrob(subset_matrix[,i]~original_matrix[,j,drop=F], k.max = 2000))$adj.r.squared},error=function(e)NA)
  r2_all_rlr[i]=NA
  if (max(r2_all_rlr,na.rm=T)<0.2) next
  model=robustbase::lmrob(subset_matrix[,i]~original_matrix[,which(r2_all_rlr>0.2),drop=F], k.max = 2000)
  predicted_matrix[,i]=predict(model)
}

outliers=matrix(NA,nrow(original_matrix),ncol(original_matrix))
for (i in 1:ncol(subset_matrix)) {
  lol=summary(robustbase::lmrob(subset_matrix[,i]~predicted_matrix[,i],max.it = 2000))$residuals
  outliers[which(lol %in% boxplot.stats(lol)$out),i]=1
}


for (i in 1:ncol(subset_matrix)) {
  if (all(is.na(subset_matrix[,i]))) next

  r2_all_rlr=rep(NA,ncol(original_matrix))
  for (j in possible_predictors) r2_all_rlr[j]=tryCatch({summary(robustbase::lmrob(subset_matrix[,i]~original_matrix[,j,drop=F], k.max = 2000))$adj.r.squared},error=function(e)NA)
  r2_all_rlr[i]=NA

  if (max(r2_all_rlr,na.rm=T)<0.2) next
  las=which(is.na(outliers[,i])==F)
  las2=which(apply(outliers[las,,drop=F],2,function(x)all(is.na(x)))==T)
  int=intersect(which(r2_all_rlr>0.2),las2)
  if (length(int)==0) next
  predicted_matrix[las,i]=predict(robustbase::lmrob(subset_matrix[,i]~original_matrix[,int,drop=F], k.max = 2000))[las]

}

matrnew=matrix(NA,nrow(original_matrix),ncol(original_matrix))
for (i in 1:ncol(subset_matrix)) {
  lol=summary(robustbase::lmrob(subset_matrix[,i]~predicted_matrix[,i],max.it = 2000))$residuals
  matrnew[which(lol %in% boxplot.stats(lol)$out),i]=1
}

for (i in 1:ncol(subset_matrix)) {
  if (all(is.na(subset_matrix[,i]))) next

  r2_all_rlr=rep(NA,ncol(original_matrix))
  for (j in possible_predictors) r2_all_rlr[j]=tryCatch({summary(robustbase::lmrob(subset_matrix[,i]~original_matrix[,j,drop=F], k.max = 2000))$adj.r.squared},error=function(e)NA)
  r2_all_rlr[i]=NA

  if (max(r2_all_rlr,na.rm=T)<0.2) next
  las=which(is.na(matrnew[,i])==F)
  las2=which(apply(matrnew[las,,drop=F],2,function(x)all(is.na(x)))==T)
  int=intersect(which(r2_all_rlr>0.2),las2)
  if (length(int)==0) next
  predicted_matrix[las,i]=predict(robustbase::lmrob(subset_matrix[,i]~original_matrix[,int,drop=F], k.max = 2000))[las]

  # provisional_predictors=provisional_predictors[duplicated(provisional_predictors)==F]
  # predictors_train[[i]]=provisional_predictors
}
return(predicted_matrix)
}
