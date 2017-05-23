#' Import of variables stored in the parameters file and of the dataset to quantify
#'
#' @param path Path where plots are stored inside a created 'plots' folder
#' @param final_output List with quantifications and indicators of quality of quantification.
#' @param imported_data List with typical elements necessary to perform quantification of ROIs.
#' @param useful_data List with necessary information to load quantifications on the Shiny GUI.
#' @return Plots in pdf files
#' @export write_plots
#' @import gridExtra
#' @import grid
#' @import ggplot2
#'
#' @examples
#' setwd(paste(system.file(package = "rDolphin"),"extdata",sep='/'))
#' load("MTBLS242_subset_example.RData")
#' #Not run:
#' #write_plots('quantification_plots',quantification_variables$final_output,imported_data,quantification_variables$useful_data)


write_plots = function(path,final_output,imported_data,useful_data) {
  path=file.path(path,'plots')
  dir.create(path)
  ind3=which(apply(final_output$shift,2, function(x) all(is.na(x)))==F)
  print('Be patient. This could take a while. Take another cup of coffee, meanwhile')
  p <- vector(mode = "list", length = nrow(imported_data$dataset))
  pb   <- txtProgressBar(1, max(ind3), style=3)

  for (ind2 in ind3) {
    for (ind in 1:nrow(imported_data$dataset)) {
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
      if (length(r)==0) next
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
        ) + theme(legend.position = "none", text = element_text(size=5)) + ggtitle(paste(imported_data$Experiments[ind]," - fitting error ",round(final_output$fitting_error[ind,ind2],3)," - signal/area ratio ",round(final_output$signal_area_ratio[ind,ind2],3),sep=''))+
        scale_x_reverse() + labs(x='ppm',y='Intensity')
    }
    gridExtra::grid.arrange(grid::rectGrob(), grid::rectGrob())
    ml <- gridExtra::marrangeGrob(p, top = imported_data$signals_names[ind2],nrow=3, ncol=1)
    ggplot2::ggsave(file.path(path,paste(imported_data$signals_names[ind2],".pdf",sep='')),  ml)
    setTxtProgressBar(pb, ind2)

  }
  print("Done!")

}
