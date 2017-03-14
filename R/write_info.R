

#' Writing of quantification and quality information of fitting of signals
#'
#' @param export_path Export path where the RData and the associated data is stored.
#' @param finaloutput List with quantifications and indicators of quality of quantification.
#' @param ROI_data ROI data.
#'
#' @return RData with session and 'associated_data' folder with CSVs with quantification and quality information of fitting of signals.
#' @export write_info
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' quantification_variables=autorun(imported_data,imported_data$finaloutput,imported_data$useful_data)
#' write_info('output_info',quantification_variables$finaloutput,imported_data$ROI_data)

write_info = function(export_path, finaloutput,ROI_data) {
  dir.create(export_path)
  write.csv(finaloutput$Area,
    file.path(export_path,
      "Area.csv"))
  write.csv(finaloutput$shift,
    file.path(export_path,
      "shift.csv"))
  write.csv(finaloutput$half_band_width,
    file.path(export_path,
      "half_band_width.csv"))
  write.csv(
    finaloutput$signal_area_ratio,
    file.path(export_path,
      "signal_area_ratio.csv")
  )
  write.csv(
    finaloutput$fitting_error,
    file.path(export_path,
      "fitting_error.csv")
  )
  write.csv(
    finaloutput$intensity,
    file.path(export_path,
      "intensity.csv")
  )
  write.csv(ROI_data,file.path(export_path,"ROI_profiles_used.csv"),row.names=F)

}

