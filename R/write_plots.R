#' Import of variables stored in the parameters file and of the dataset to quantify
#'
#' @param path Path where plots are stored inside a created 'plots' folder
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param useful_data List with necessary information to load quantifications on the Shiny GUI.
#' @param signals_to_plot Vector of indexes of signals in ROI data to plot. By default, NA and all figures are outputted.
#' @return Plots in pdf files
#' @export write_plots
#' @import gridExtra
#' @import grid
#' @import ggplot2
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' load("MTBLS242_subset_profiling_data.RData")
#' #Not run:
#' #write_plots('quantification_plots',profiling_data$final_output,imported_data,profiling_data$useful_data)


write_plots = function(path,final_output,imported_data,useful_data,signals_to_plot=NA) {
  dir.create(file.path(path,'plots'))
  if (is.na(signals_to_plot)) signals_to_plot=which(apply(final_output$quantification,2, function(x) all(is.na(x)))==F)
  print('Be patient. This could take a while. Take another cup of coffee, meanwhile')
  p <- vector(mode = "list", length = nrow(final_output$quantification))
  if (max(signals_to_plot)>1) pb   <- txtProgressBar(1, max(signals_to_plot), style=3)

  for (ind2 in signals_to_plot) {
    for (ind in 1:nrow(final_output$quantification)) {
      Xdata=useful_data[[ind]][[ind2]]$Xdata
      if(is.null(Xdata)) next
      Ydata=useful_data[[ind]][[ind2]]$Ydata
      plot_data=useful_data[[ind]][[ind2]]$plot_data
      ROI_profile=useful_data[[ind]][[ind2]]$ROI_profile
      plotdata2 = data.frame(Xdata=Xdata,
        Ydata=Ydata,
        plot_data[3, ] ,
        plot_data[2, ])
      colnames(plotdata2)=c('Xdata','Ydata',"fitted_sum","baseline_sum")

      plotdata3 <- reshape2::melt(plotdata2, id = "Xdata")
      plotdata3$variable = c(
        rep('Original Spectrum', length(Ydata)),
        rep('Generated Spectrum', length(Ydata)),
        rep('Generated Background', length(Ydata))
      )
      plotdata4 = data.frame(Xdata, (t(plot_data[-c(1, 2, 3), , drop = F])))

      colnames(plotdata4)=c('Xdata',rownames(plot_data)[-c(1, 2, 3)])
      plotdata5 = reshape2::melt(plotdata4, id = "Xdata")
      r=which(make.names(paste(ROI_profile[,4],ROI_profile[,5],sep='_'))==imported_data$signals_names[ind2])
      if (length(r)==0) next
      plotdata = data.frame(Xdata, signals = plot_data[3 + r,] )



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
        ) + theme(legend.position = "none", text = element_text(size=5)) + ggtitle(paste(imported_data$Experiments[ind]," - fitting error ",round(final_output$fitting_error[ind,ind2],3)," - signal/area ratio ",round(final_output$signal_area_ratio[ind,ind2],3),sep=''))+
        scale_x_reverse() + labs(x='ppm',y="Intensity (arbitrary unit)")
    }
      gridExtra::grid.arrange(grid::rectGrob(), grid::rectGrob())
      ml <- gridExtra::marrangeGrob(p, top = imported_data$signals_names[ind2],nrow=3, ncol=1)
      ggplot2::ggsave(file.path(path,paste(imported_data$signals_names[ind2],".pdf",sep='')),  ml)

    setTxtProgressBar(pb, ind2)

  }
  print("Done!")

}
