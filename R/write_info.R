

#' Writing of quantification and quality information of fitting of signals
#'
#' @param export_path Export path where the RData and the associated data is stored.
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param ROI_data ROIs data.
#'
#' @return RData with session and 'associated_data' folder with CSVs with quantification and quality information of fitting of signals.
#' @export write_info
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' quantification_variables=autorun(imported_data,imported_data$final_output,imported_data$useful_data,imported_data$ROI_data)
#' write_info('output_info',quantification_variables$final_output,imported_data$ROI_data)

write_info = function(export_path, final_output,ROI_data) {
  dir.create(export_path)
  write.csv(final_output$Area,
    file.path(export_path,
      "Area.csv"))
  write.csv(final_output$shift,
    file.path(export_path,
      "shift.csv"))
  write.csv(final_output$half_band_width,
    file.path(export_path,
      "half_band_width.csv"))
  write.csv(
    final_output$signal_area_ratio,
    file.path(export_path,
      "signal_area_ratio.csv")
  )
  write.csv(
    final_output$fitting_error,
    file.path(export_path,
      "fitting_error.csv")
  )
  write.csv(
    final_output$intensity,
    file.path(export_path,
      "intensity.csv")
  )
  write.csv(ROI_data,file.path(export_path,"ROI_profiles_used.csv"),row.names=F)

}

