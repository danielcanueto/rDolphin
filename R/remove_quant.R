#' Import of variables stored in the parameters file and of the dataset to quantify
#'
#' @param info Experiment and signal where to remove information
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param final_output List with quantifications and indicators of quality of quantification.
#'
#' @return Imported data of experiment
#' @export remove_quant



remove_quant=function(info,imported_data,final_output) {
  ind1=info$row
  ind2=info$col


final_output$quantification[ind1,ind2]=final_output$shift[ind1,ind2]=final_output$quantification[ind1,ind2]=final_output$half_band_width[ind1,ind2]=final_output$signal_area_ratio[ind1,ind2]=final_output$fitting_error[ind1,ind2]=final_output$intensity[ind1,ind2]=NA

return(final_output)
}
