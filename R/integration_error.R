#' Helper function to generate ROI data. Not to be used in console.
#'
#' @param ROI_data ROI_data
#' @param useful_data useful_data
#' @param final_output final_output
#' @param ind ind
#' @return dummy
#' @export integration_error

integration_error <- function(ROI_data,useful_data,final_output,ind) {
for (ii in ind) {
  all4=t(as.data.frame(lapply(useful_data,function(x)x[[ii]]$plot_data[1,])))
  spectra_lag=rep(NA,dim(all4)[1])
  for (i in 1:dim(all4)[1]) {
    d <-
      ccf(all4[i, ],
          apply(all4, 2, median),
          type = 'covariance',
          plot = FALSE)
    spectra_lag[i]=d$lag[which.max(d$acf)]
  }
  so=(1+max(abs(spectra_lag))):(ncol(all4)-max(abs(spectra_lag)))
  for (i in 1:dim(all4)[1])   all4[i,so-spectra_lag[i]]=all4[i,so]

fitted_median=apply(all4,2,median)
sorted_bins=sort(fitted_median/sum(fitted_median),decreasing=T,index.return=T)
if(length(sorted_bins$x)>0) {
  bins= sorted_bins$ix[1:which.min(abs(cumsum(sorted_bins$x)-0.9))]
} else {
  bins=seq_along(fitted_median)
}

ple=apply(all4[,bins],1,function(x)summary(lm(fitted_median[bins]~x))$sigma/max(fitted_median[bins]))
for (j in seq(nrow(final_output$fitting_error))) useful_data[[j]][[ii]]$results_to_save$fitting_error=ple[j]
final_output$fitting_error[,ii]=ple
}
dummy=list(final_output=final_output,useful_data=useful_data)
}
