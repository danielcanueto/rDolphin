#' Import of variables stored in the parameters file and of the dataset to quantify
#'
#' @param info Experiment and signal where to remove information
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param finaloutput List with quantifications and indicators of quality of quantification.
#'
#' @return Imported data of experiment
#' @export remove_quant



remove_quant=function(info,imported_data,finaloutput) {
  ind1=info$row
  ind2=info$col


finaloutput$Area[ind1,ind2]=finaloutput$shift[ind1,ind2]=finaloutput$Area[ind1,ind2]=finaloutput$half_band_width[ind1,ind2]=finaloutput$signal_area_ratio[ind1,ind2]=finaloutput$fitting_error[ind1,ind2]=finaloutput$intensity[ind1,ind2]=NA

return(finaloutput)
}
