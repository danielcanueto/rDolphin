
remove_quant=function(info,imported_data,finaloutput) {
  ind1=info$row
  ind2=info$col


finaloutput$Area[ind1,ind2]=finaloutput$shift[ind1,ind2]=finaloutput$Area[ind1,ind2]=finaloutput$half_band_width[ind1,ind2]=finaloutput$signal_area_ratio[ind1,ind2]=finaloutput$fitting_error[ind1,ind2]=finaloutput$intensity[ind1,ind2]=NA

# tryCatch({write_info(imported_data$export_path, finaloutput)}, error = function(err) {
#   print('Not possible to overwrite a csv file open with Microsoft Excel')
# })

return(finaloutput)
}
