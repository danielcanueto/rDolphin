
integration = function(clean_fit, Xdata, Ydata,interface='F') {

  #preallocation of results_to_save
  results_to_save = list(
    shift = NA,
    quantification = NA,
    signal_area_ratio = NA,
    fitting_error = NA,
    intensity = NA,
    half_band_width = NA
  )

  baseline = rep(0,length(Xdata))
  #preparation of baseline, if specified by the user
  if (clean_fit == 'N') baseline = seq(min(Ydata[1:5]), min(Ydata[(length(Xdata) - 4):length(Xdata)]), len =length(Xdata))
baseline[which((baseline-Ydata)>0)]=Ydata[which((baseline-Ydata)>0)]
  #integration ad chechk that there are no negative values
  integrated_signal = Ydata - baseline
    integrated_signal[integrated_signal<0]=0
  #preparation of results_to_save
  results_to_save$quantification = sum(integrated_signal)
  results_to_save$intensity = max(integrated_signal)

  cumulative_area = cumsum(integrated_signal) / sum(integrated_signal)

  if (all(is.na(cumulative_area))) {
    p1=1
    p2=length(cumulative_area)
  } else {
	  p1 = which(cumulative_area< 0.05)[length(which(cumulative_area< 0.05))]
	  p2 = which(cumulative_area > 0.95)[1]
  }

  results_to_save$signal_area_ratio = tryCatch((sum(integrated_signal[p1:p2]) / sum(Ydata[p1:p2])) *
    100,error = function(e) NaN, silent=T)

  # d=ccf(Ydata,Ydatamedian,type = 'covariance',plot = FALSE)
  # so=(1+max(abs(d$lag[which.max(d$acf)]))):(length(Ydata)-max(abs(d$lag[which.max(d$acf)])))
  # results_to_save$fitting_error=summary(lm(Ydatamedian[so]~Ydata [so]))$sigma/max(Ydatamedian[so])
  results_to_save$half_band_width = NaN
#
#   peaks = peakdet(integrated_signal, 0.2*max(0.000001,max(integrated_signal)))
#   results_to_save$shift = mean(Xdata[peaks$maxtab$pos])
  wer=which.min(abs(cumulative_area-0.5))
  if (cumulative_area[wer]>0.5) {
    wer=c(wer-1,wer)
  } else {
    wer=c(wer,wer+1)
  }
  wer2=(0.5-cumulative_area[wer[1]])/diff(cumulative_area[wer])
  results_to_save$shift = Xdata[wer[1]]-wer2*diff(Xdata[wer])

p=''
if (interface=='T') {
plotdata = data.frame(Xdata, signal = integrated_signal)
  plotdata2 = data.frame(Xdata, Ydata)
  plotdata3 <- melt(plotdata2, id = "Xdata")
  plotdata3$variable = rep('Original Spectrum', length(Ydata))
  plotdata4 = data.frame(Xdata, integrated_signal)
  plotdata5 = melt(plotdata4, id = "Xdata")
  p=plot_ly(plotdata,x = ~Xdata, y = ~signal, type = 'scatter', color= 'Signal',mode = 'lines', fill = 'tozeroy') %>% add_trace(data=plotdata3,x=~Xdata,y=~value,color=~variable,type='scatter',mode='lines',fill=NULL) %>%
    layout(xaxis = list(range=c(Xdata[1],Xdata[length(Xdata)]),title = 'ppm'),
      yaxis = list(range=c(0,max(Ydata)),title = 'Intensity'))
  }
  plot_data=rbind(integrated_signal,baseline,integrated_signal+baseline,integrated_signal)

  integration_variables=list(results_to_save=results_to_save,p=p,plot_data=plot_data)
  return(integration_variables)
}
