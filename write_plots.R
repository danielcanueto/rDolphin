#' Import of variables stored in the parameters file and of the dataset to quantify
#'
#' @param path Path where plots are stored inside a created 'plots' folder
#' @param finaloutput List with quantifications and indicators of quality of quantification.
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param useful_data List with necessary information to load quantifications on the Shiny GUI.#'
#' @return Plots in pdf files
#' @export write_plots
#' @import gridExtra
#' @import grid
#' @import ggplot2
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' quantification_variables=autorun(imported_data,imported_data$finaloutput,imported_data$useful_data)
#' write_plots('quantification_plots',quantification_variables$finaloutput,imported_data,quantification_variables$useful_data)


write_plots = function(path,finaloutput,imported_data,useful_data) {
  path=paste(path,'plots/',sep='/')
  dir.create(path)
  ind3=which(apply(finaloutput$shift,2, function(x) all(is.na(x)))==F)
  print('Be patient. This could take a while. Take another cup of coffee, meanwhile')
  p <- vector(mode = "list", length = nrow(imported_data$dataset))
  pb   <- txtProgressBar(1, length(ind3), style=3)

  for (ind2 in ind3) {
    for (ind in 1:nrow(imported_data$dataset)) {
      Xdata=try(useful_data[[ind]][[ind2]]$Xdata,silent=T)
      Ydata=useful_data[[ind]][[ind2]]$Ydata
      plot_data=useful_data[[ind]][[ind2]]$plot_data
      ROI_profile=useful_data[[ind]][[ind2]]$ROI_profile
      plotdata2 = data.frame(Xdata=Xdata,
        Ydata=Ydata,
        plot_data[3, ] ,
        plot_data[2, ])
      colnames(plotdata2)=c('Xdata','Ydata',"fitted_sum","baseline_sum")

      plotdata3 <- melt(plotdata2, id = "Xdata")
      plotdata3$variable = c(
        rep('Original Spectrum', length(Ydata)),
        rep('Generated Spectrum', length(Ydata)),
        rep('Generated Background', length(Ydata))
      )
      plotdata4 = data.frame(Xdata, (t(plot_data[-c(1, 2, 3), , drop = F])))

      colnames(plotdata4)=c('Xdata',rownames(plot_data)[-c(1, 2, 3)])
      plotdata5 = melt(plotdata4, id = "Xdata")
      r=which(paste(ROI_profile[,4],ROI_profile[,5],sep='_')==imported_data$signals_names[ind2])
      plotdata = data.frame(Xdata, signals = plot_data[3 + r, ] )



      p[[ind]]=ggplot() +
        geom_line(data = plotdata3,
          aes(
            x = Xdata,
            y = value,
            colour = variable,
            group = variable
          )) +
        geom_line(data = plotdata5,
          aes(
            x = Xdata,
            y = value,
            colour = 'Surrounding signals',
            group = variable
          )) +
        geom_area(
          data = plotdata,
          aes(
            x = Xdata,
            y = signals,
            fill = 'Quantified Signal'
          )
        ) + theme(legend.position = "none", text = element_text(size=5)) + ggtitle(paste(imported_data$Experiments[ind]," - fitting error ",round(finaloutput$fitting_error[ind,ind2],3)," - signal/area ratio ",round(finaloutput$signal_area_ratio[ind,ind2],3),sep=''))+
        scale_x_reverse() + labs(x='ppm',y='Intensity')
    }
    grid.arrange(rectGrob(), rectGrob())
    ml <- marrangeGrob(p, top = imported_data$signals_names[ind2],nrow=3, ncol=1)
    ggsave(paste(path,imported_data$signals_names[ind2],".pdf",sep=''),  ml)
    setTxtProgressBar(pb, ind2)

  }
  print("Done!")

}
