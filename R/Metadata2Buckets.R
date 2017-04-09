
Metadata2Buckets <- function(Experiments, params, spectrum_borders) {

  CURRENT = list()
  RAW = list()
  not_loaded_experiments = read_spectra = c()


  left_spectral_border = ifelse(exists("left_spectral_border", where = params),
                                params$left_spectral_border,
    spectrum_borders[1])
  right_spectral_border = ifelse(exists("right_spectral_border", where = params),
                                 params$right_spectral_border,
    spectrum_borders[2])

  RAW$norm_PEAK_left_ppm = ifelse(params$norm_PEAK == "Y", params$norm_left_ppm,
    0)
  RAW$norm_PEAK_right_ppm = ifelse(params$norm_PEAK == "Y", params$norm_left_ppm,
    0)
  RAW$norm_AREA_left_ppm = ifelse(params$norm_AREA == "Y", params$norm_left_ppm,
                                  0)
  RAW$norm_AREA_right_ppm = ifelse(params$norm_AREA == "Y", params$norm_right_ppm,
                                   0)

  tsp_alignment = params$tsp_alignment
  ref_peak_pos = ifelse(params$peak_alignment == "Y", params$ref_peak_pos,
                        0)

  k2 = 1
  maxspec = length(Experiments)

  for (k in 1:maxspec) {
    filename = paste(params$dir, "nmr", Experiments[k], params$expno, "pdata", params$processingno, sep = "/")
    partname = paste(params$dir, "nmr", Experiments[k], params$expno, sep = "/")
    storedpars = topspin_read_spectrum2(partname, filename,spectrum_borders[2]-0.1, spectrum_borders[1]+0.1)
    if (all(is.nan(storedpars$real)) == 0) {
      CURRENT$minppm = storedpars$OFFSET - storedpars$SW
      CURRENT$maxppm = storedpars$OFFSET
      CURRENT$step = storedpars$SW / (length(storedpars$real) - 1)
      CURRENT$ppm = seq(CURRENT$maxppm, CURRENT$minppm,-CURRENT$step)

     tmp = (storedpars$real * ((2 ^ storedpars$NC_proc) / storedpars$RG))
	# left_spectral_border =  floor(min(CURRENT$maxppm,left_spectral_border)*10)/10
	#     	right_spectral_border =  ceiling(max(CURRENT$minppm,right_spectral_border)*10)/10

     if (left_spectral_border>CURRENT$maxppm | right_spectral_border<CURRENT$minppm) {
       print('Current incorrect spectrum borders. Please prepare a modified version of fitting_variables with appropiate spectrum borders and introduce its path on the Parameters CSV file.')
       finaldata=NA
       return(finaldata)
     }
      CURRENT$left_spectral_border = 1 + round(-(left_spectral_border -
                                                   CURRENT$maxppm) / CURRENT$step)
      CURRENT$norm_PEAK_left = 1 + round(-(RAW$norm_PEAK_left_ppm -
                                             CURRENT$maxppm) / CURRENT$step)
      CURRENT$norm_PEAK_right = 1 + round(-(RAW$norm_PEAK_right_ppm -
                                              CURRENT$maxppm) / CURRENT$step)
      CURRENT$norm_AREA_left = 1 + round(-(RAW$norm_AREA_left_ppm -
                                             CURRENT$maxppm) / CURRENT$step)
      CURRENT$norm_AREA_right = 1 + round(-(RAW$norm_AREA_right_ppm -
                                              CURRENT$maxppm) / CURRENT$step)
      norm_PEAK_max = max(tmp[CURRENT$norm_PEAK_left:CURRENT$norm_PEAK_right])
      norm_PEAK_max_position = which.max(tmp[CURRENT$norm_PEAK_left:CURRENT$norm_PEAK_right])
      norm_AREA = sum(tmp[CURRENT$norm_AREA_left:CURRENT$norm_AREA_right])
      total_AREA_mean = mean(tmp[1:length(tmp)])

      RAW$len = length(storedpars$real)
      RAW$ofs = storedpars$OFFSET
      RAW$sw = storedpars$SW
      RAW$ncproc = storedpars$NC_proc
      RAW$minppm = storedpars$OFFSET - storedpars$SW
      RAW$maxppm = storedpars$OFFSET
      RAW$size = length(storedpars$real)
      RAW$step = RAW$sw / (RAW$size - 1)
      RAW$ppm = seq(RAW$maxppm, RAW$minppm,-RAW$step)
      RAW$buck_step = ifelse(params$buck_step < RAW$step, RAW$step,
                             params$buck_step)
      # redieixo la part dreta de l'espectre de Matlab RAW$ppm_bucks =
      RAW$ppm_bucks = seq(left_spectral_border,
                          right_spectral_border,-RAW$buck_step)
      RAW$len_bucks = length(RAW$ppm_bucks)
      RAW$norm_PEAK_max = norm_PEAK_max
      RAW$total_AREA_mean = total_AREA_mean
      RAW$norm_AREA = norm_AREA
      RAW$differential_norm_AREA = 0
      # }

      if (params$norm_PEAK == "Y" &
          norm_PEAK_max > 0)
        tmp = tmp / norm_PEAK_max
      if (params$norm_AREA == "Y" &
          norm_AREA > 0)
        tmp = tmp / norm_AREA

      # calculem valors despres de normalitzar
      norm_PEAK_max2 = max(tmp[CURRENT$norm_PEAK_left:CURRENT$norm_PEAK_right])
      norm_AREA2 = sum(tmp[CURRENT$norm_AREA_left:CURRENT$norm_AREA_right])
      # referenciem
      LeftEreticBefore = ifelse(CURRENT$ppm[1] > 11.8, min(CURRENT$ppm[which(CURRENT$ppm > 11.7)]), NaN)

      if (params$glucose_alignment == "Y") {
        JTP = JTPcalibrateToGlucose(tmp, CURRENT$ppm)
      } else if (params$tsp_alignment == "Y") {
        JTP = JTPcalibrateToTSP(tmp, CURRENT$ppm)
      } else if (params$peak_alignment == "Y") {
        if (nrow(as.matrix(ref_peak_pos)) == maxspec) {
          JTP = JTPcalibrateToPeak2(tmp, CURRENT$ppm, 5.3, 0.4)
          JTP = JTPcalibrateToPeak(tmp, CURRENT$ppm, ref_peak_pos[k])
        } else {
          JTP = JTPcalibrateToPeak(tmp, CURRENT$ppm, ref_peak_pos)
	}
      } else {
	JTP=list(ppm=CURRENT$ppm)
	}
      CURRENT$ppm = JTP$ppm
      # si hi llegim zona d'eretic llavors tornem a posar Eretic a 11 ppm
      if (!is.nan(LeftEreticBefore)) {
        LeftEreticAfter = min(CURRENT$ppm[which(CURRENT$ppm > 11.7)])
        RightEreticAfter = min(CURRENT$ppm[which(CURRENT$ppm > 10.8)])
        RightEreticBefore = LeftEreticBefore + RightEreticAfter - LeftEreticAfter
        tmp[LeftEreticAfter:RightEreticAfter] = tmp[LeftEreticBefore:RightEreticBefore]
      }



      fill2end = 0
      tmp_count = CURRENT$left_spectral_border
      tmp_buck = rep(0, RAW$len_bucks)
      jumped_bucket = 0
      for (tmp_buck_count in 2:RAW$len_bucks) {
        items = suma = 0
        while (CURRENT$ppm[tmp_count] > RAW$ppm_bucks[tmp_buck_count]) {
          items = items + 1
          suma = suma + tmp[tmp_count]
          tmp_count = tmp_count + 1
          if (tmp_count > length(CURRENT$ppm)) {
            fill2end = 1
            break
          }
        }

	# faig una modificacio per a corregir buckets sense dades.
        if (items > 0) {
          tmp_buck[tmp_buck_count - 1] = suma / items
          while (jumped_bucket > 0) {
            tmp_buck[tmp_buck_count - 1 - jumped_bucket] = suma / items
            jumped_bucket = jumped_bucket - 1
          }
        } else {
          jumped_bucket = jumped_bucket + 1
        }
        if (fill2end == 1) {
          while (tmp_buck_count <= RAW$len_bucks) {
            tmp_buck[tmp_buck_count] = tmp_buck[tmp_buck_count -
                                                  1]
            tmp_buck_count = tmp_buck_count + 1
          }
          break
        }
      }
      while (length(tmp_buck) < RAW$len_bucks)
        tmp_buck = c(tmp_buck, tmp_buck[length(tmp_buck)])




      if (k2 == 1) dataset = matrix(NA, maxspec, RAW$len_bucks)
      dataset[k2,] = tmp_buck[1:RAW$len_bucks]
      k2 = k2 + 1
      read_spectra = append(read_spectra, as.character(Experiments[k]))


    } else {
      # llista d'objectes no inclosos
      not_loaded_experiments = append(not_loaded_experiments, Experiments[k])
    }
  }
  # finaldata$params = params
  # finaldata$RAW = RAW
  # dataset = dataset[complete.cases(dataset),,drop=F]
  rownames(dataset) = read_spectra
  colnames(dataset) = as.character(RAW$ppm_bucks)
  finaldata = list(dataset=dataset,ppm=RAW$ppm_bucks,not_loaded_experiments = not_loaded_experiments)


  return(finaldata)
}
