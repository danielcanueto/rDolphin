save_output=function(spectrum_index,signals_codes,results_to_save,buck_step,finaloutput) {
  #Created by Daniel Canueto 30/08/2016
  #Save quantification, shift, fitting error and signal area ratio.
  
             finaloutput$Area[spectrum_index,signals_codes]=results_to_save$Area*buck_step
             finaloutput$fitting_error[spectrum_index,signals_codes]=results_to_save$fitting_error 
             finaloutput$signal_area_ratio[spectrum_index,signals_codes]=results_to_save$signal_area_ratio
             finaloutput$shift[spectrum_index,signals_codes]=results_to_save$shift
             finaloutput$intensity[spectrum_index,signals_codes]=results_to_save$intensity
             finaloutput$half_band_width[spectrum_index,signals_codes]=results_to_save$half_band_width
             
      return(finaloutput)         
}