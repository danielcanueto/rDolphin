

#' Alignment of signals of the dataset through the CluPA algorithm provided by the speaq package.
#'
#' @param dataset The 1D NMR dataset where to align signals
#' @param buck_step The bucketing of the dataset (e.g. 0.001)
#'
#' @return aligneddataset The dataset with the signals aligned
#' @export alignment
#' @import speaq
#' @import MassSpecWavelet
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' aligned_data=alignment(imported_data$dataset,imported_data$buck_step)


alignment=function(dataset,buck_step) {
    print('Be patient. Gonna take a while. You should be writing, meanwhile.')

    peakList <- detectSpecPeaks(dataset,
      nDivRange = c(128),
      scales = seq(1, 16, 2),
      baselineThresh = quantile(dataset,0.60,na.rm=T),
      SNR.Th = -1,
      verbose=FALSE
    );
    resFindRef<- findRef(peakList);
    refInd <- resFindRef$refInd;

    maxShift = 0.025/buck_step;
    aligneddataset <- dohCluster(dataset,
      peakList = peakList,
      refInd = refInd,
      maxShift = maxShift,
      acceptLostPeak = TRUE, verbose=FALSE);

    print('Done!')
    aligneddataset[is.na(aligneddataset)]=0
    return(aligneddataset)

}
