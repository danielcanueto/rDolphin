server = function(input, output,session) {

  #Increase of maximum memory size that can be uploaded
  options(shiny.maxRequestSize=1000*1024^2)
  options(warn =-1)
   session$onSessionEnded(stopApp)

  #Setting of reactive parameters in the Shiny GUI
  reactiveROItestingdata <- reactiveValues(signpar = matrix(NA,2,7,dimnames = list(seq(2),c("intensity",	"chemical shift",	"half bandwidth",	"gaussian %",	"J coupling",	"multiplicities",	"roof effect"))),qualitypar = matrix(NA,2,3,dimnames=list(seq(2),c('Quantification','fitting_error','signal/total area ratio'))))
  reactivequantdata <- reactiveValues(method1=NA)
 reactiveprogramdata <- reactiveValues(ROIdata_subset=NA,ind=NA,beginning=FALSE,dataset=NA,final_output=list(),reproducibility_data=list(),imported_data=NA,p_value_final=NA,ROI_data=NA,ROI_data_check=NA,info=c(),select_options=NA,new_roi_profile=NA,p=NA,bgColScales=NA,automatic_profiling_plot=NA,ROI_names=NA,clusterplot=NA,median_plot=NA,jres_plot=NA)
reac=reactiveValues(cho=NA)
  ## FIRST TAB REACTIVE OUTPUTS
 observe({
   dummy=ifelse(reactiveprogramdata$beginning ==F,0,1)

      toggle(condition = dummy, selector = "#mynavlist li a[data-value=tab2]")
      toggle(condition = dummy, selector = "#mynavlist li a[data-value=tab3]")
      toggle(condition = dummy, selector = "#mynavlist li a[data-value=tab4]")
      # toggle(condition = dummy, selector = "#mynavlist li a[data-value=tab5]")
      toggle(condition = dummy, selector = "#mynavlist li a[data-value=tab6]")
 })


  #Read of input provided by user
  observeEvent(input$file1, {
    reactiveprogramdata$inFile <- input$file1
    if (is.null(reactiveprogramdata$inFile)) {
      return(NULL)
    }


	tryCatch({
	reactiveprogramdata$imported_data = tryCatch({suppressWarnings(import_data(reactiveprogramdata$inFile$datapath))}, error = function(e) {
	  print('Import of data not posible with current input')
	  return(NULL)
	})
	reset("file1")

  reactiveprogramdata$final_output=reactiveprogramdata$imported_data$final_output
	reactiveprogramdata$reproducibility_data=reactiveprogramdata$imported_data$reproducibility_data
	reactiveprogramdata$ROI_data=reactiveprogramdata$ROI_data_check=reactiveprogramdata$imported_data$ROI_data
	reactiveprogramdata$list=seq(nrow(reactiveprogramdata$ROI_data))
	reactiveprogramdata$imported_data$final_output=reactiveprogramdata$imported_data$reproducibility_data=reactiveprogramdata$imported_data$ROI_data=NULL
	colnames(reactiveprogramdata$ROI_data)=c("ROI left edge (ppm)","ROI right edge (ppm)","Quantification Mode","Metabolite","Quantification Signal","Chemical shift (ppm)","Chemical shift tolerance (ppm)","Half bandwidth (Hz)","Multiplicity","J coupling (Hz)","Roof effect")
	reactiveprogramdata$validation_data=list(alarm_matrix=reactiveprogramdata$final_output)
	reactiveprogramdata$validation_data=validation(reactiveprogramdata$final_output,1,reactiveprogramdata$validation_data$alarm_matrix)
	dummy=tryCatch({profile_model_spectrum(reactiveprogramdata$imported_data,reactiveprogramdata$ROI_data)}, error = function(e) {
	print('Automatic quantification of model spectrum not possible.')
	})
	reactiveprogramdata$automatic_profiling_plot=dummy$p
	reactiveprogramdata$total_signals_parameters=dummy$total_signals_parameters
	# dummy$indicators=dummy2$indicators

	print('Generating additional information...')

	#plots of representative spectra and median spectra per group to help setting the right ROI parameters

    reactiveprogramdata$clusterplot=tryCatch({exemplars_plot(reactiveprogramdata$imported_data)  }, error = function(e) {
	print('Generation of subsets or representative spectra not possible.')
	})
    reactiveprogramdata$median_plot=tryCatch({median_plot(reactiveprogramdata$imported_data)
	  }, error = function(e) {
	print('Generation of median spectra not possible.')
	})

    #Subsetting of ROIs is prepared
    #Names of ROIS and cluster and median spectra are prepared
	dummy=NULL
	dummy=tryCatch({roifunc(reactiveprogramdata$ROI_data,reactiveprogramdata$imported_data$Metadata,reactiveprogramdata$imported_data$Experiments)
  }, error = function(e) {
	print('Generation of Regions of Interest not possible. Please explain the issue in the Github page.')
	return(NULL)
	})
	if (!is.null(dummy)) {
	reactiveprogramdata$select_options=dummy$select_options
	reactiveprogramdata$spectra=dummy$spectra
	reactiveprogramdata$beginning =TRUE

	reactiveprogramdata$jres_plot=tryCatch(twod_data(reactiveprogramdata$imported_data$jres_path), error = function(e) NA)


	 #Variables that can change during the use of the GUI are separated from 'imported_data'.
	output$moreControls <- renderUI({
	  if (reactiveprogramdata$beginning==T) selectInput("select",label=NULL,choices = reactiveprogramdata$select_options,selected = 1)
	})
	shinyjs::show('automatic_profiling')
	shinyjs::show('alignment')
	shinyjs::show('model_spectrum')
	#When the session is prepared, the tabs and some inputs become active
	print('Done!')


	}}, error = function(e) {
	  print('Error. Please explain the issue in the Github page.')
	  reactiveprogramdata$imported_data=NA
	  return(reactiveprogramdata$imported_data)
	})

	})



  #Loading of previous session

    #Read of input provided by user
  observeEvent(input$file2, {
    reactiveprogramdata$inFile2 <- input$file2
    if (is.null(reactiveprogramdata$inFile2))
      return(NULL)
    print("Uploading saved session.")
    #Session is loaded in 'savedreactiveddata' variable and passed to 'reactiveprogramdata'.
    tryCatch({

      load(reactiveprogramdata$inFile2$datapath)
      plo=names(sapply(savedreactivedata, names))
      for (i in 1:length(plo)) reactiveprogramdata[[plo[i]]]=savedreactivedata[[plo[i]]]
      p=plot_ly(x=reactiveprogramdata$imported_data$ppm,y=reactiveprogramdata$dataset)
      reactiveprogramdata$plot$dependencies=reactiveprogramdata$automatic_profiling_plot$dependencies=reactiveprogramdata$clusterplot$dependencies=reactiveprogramdata$median_plot$dependencies=p$dependencies
      reactiveprogramdata$validation_data=validation(reactiveprogramdata$final_output,1,reactiveprogramdata$validation_data$alarm_matrix)

      rm(savedreactivedata)
reset("file2")
}, error = function(e) {
      print('Not possible to load the session. Please revise your choice.')
      return(NULL)
    })
	colnames(reactiveprogramdata$ROI_data)=c("ROI left edge (ppm)","ROI right edge (ppm)","Quantification Mode","Metabolite","Quantification Signal","Chemical shift (ppm)","Chemical shift tolerance (ppm)","Half bandwidth (Hz)","Multiplicity","J coupling (Hz)","Roof effect")
if (is.null(reactiveprogramdata$validation_data)) reactiveprogramdata$validation_data=validation(reactiveprogramdata$final_output,input$select_validation,reactiveprogramdata$validation_data$alarm_matrix)

    reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data
    reactiveprogramdata$list=seq(nrow(reactiveprogramdata$ROI_data))

    #Names of ROIS are prepared
# 	dummy=tryCatch({roifunc(reactiveprogramdata$ROI_data,reactiveprogramdata$imported_data$Metadata,reactiveprogramdata$imported_data$Experiments)
#   }, error = function(e) {
# 	print('Generation of Regions of Interest not possible. Please explain the issue in the Github page.')
# 	return(dummy=NULL)
# 	})
    dummy=NULL
	dummy=roifunc(reactiveprogramdata$ROI_data,reactiveprogramdata$imported_data$Metadata,reactiveprogramdata$imported_data$Experiments)

	if (!is.null(dummy)) {
	reactiveprogramdata$select_options=dummy$select_options
	reactiveprogramdata$spectra=dummy$spectra
	shinyjs::show('automatic_profiling')
	shinyjs::show('alignment')
	shinyjs::show('model_spectrum')

	reactiveprogramdata$beginning =TRUE
	output$moreControls <- renderUI({
	  if (reactiveprogramdata$beginning==T)  {
	    print("Done!")
	    selectInput("select",label=NULL,choices = reactiveprogramdata$select_options,selected = 1)
	  }
	})

	}
  })




  #Choice and storage of data associated to session
  observeEvent(input$save, {
    tryCatch({
      print('Saving information...')
    savedreactivedata=isolate(reactiveValuesToList(reactiveprogramdata))
    save(savedreactivedata, file=paste(input$caption,".RData",sep=''))
      write_info(input$caption, reactiveprogramdata$final_output, reactiveprogramdata$ROI_data)
      print('Done!')
    },error=function(e) {print('Not possible to generate the output the session files. Please revise the path given.')})
  })

observeEvent(input$folder, {
  tryCatch({
    write_plots(input$caption,reactiveprogramdata$final_output,reactiveprogramdata$reproducibility_data)},
    error= function(e) {       print('Not possible to generate the plot folder. Please check that you have permissions for the path specified.')

})})



    output$sp = DT::renderDataTable(
    reactiveprogramdata$total_signals_parameters , selection = list(selected = NULL),server = TRUE)


   #Automatic quantification of all ROIs in all spectra
  observeEvent(input$automatic_profiling, {
    tryCatch({
    profiling_data = automatic_profiling(reactiveprogramdata$imported_data, reactiveprogramdata$ROI_data)
    reactiveprogramdata$final_output=profiling_data$final_output
    reactiveprogramdata$reproducibility_data=profiling_data$reproducibility_data

    reactiveprogramdata$validation_data=validation(reactiveprogramdata$final_output,input$select_validation,reactiveprogramdata$validation_data$alarm_matrix)

    },
    error = function(e) {
      print('Error. Please explain the issue in the Github page if necessary.')
      profiling_data=NA
      return(profiling_data)
    })
  })
  #Alignment of signals
  tryCatch({observeEvent(input$alignment, {
    reactiveprogramdata$imported_data$dataset= alignment(reactiveprogramdata$imported_data$dataset,reactiveprogramdata$imported_data$buck_step)
      reactiveprogramdata$clusterplot=exemplars_plot(reactiveprogramdata$imported_data)
      reactiveprogramdata$median_plot=median_plot(reactiveprogramdata$imported_data)
      }	)}, error = function(e) {
	print('Error during alignment. Please explain the issue in the Github page if necessary.')
	return(NULL)

  })
   #Automatic quantification of all ROIs in all spectra
  observeEvent(input$model_spectrum, {
    tryCatch({
    dummy=profile_model_spectrum(reactiveprogramdata$imported_data,reactiveprogramdata$ROI_data)
    reactiveprogramdata$automatic_profiling_plot=dummy$p
    reactiveprogramdata$total_signals_parameters=dummy$total_signals_parameters
    }, error = function(e) {
      print('Automatic quantification of model spectrum not possible.')
    })
  })

  #Peak analysis. UNSTABLE!!!
  # observeEvent(input$peak_analysis, {
  #   if (is.null(reactiveprogramdata$alignment_check)) {
  #     print('Before analysing peaks, I have to align them. Then I\'ll analyze them')
  #     dummy=alignment(reactiveprogramdata$imported_data$dataset,reactiveprogramdata$imported_data$buck_step)
  #     peak_analysis(dummy,reactiveprogramdata$imported_data$ppm,reactiveprogramdata$imported_data$freq,reactiveprogramdata$imported_data$export_path,reactiveprogramdata$imported_data$Metadata,reactiveprogramdata$imported_data$repository,reactiveprogramdata$originaldataset)
  #   } else {
  #     peak_analysis(reactiveprogramdata$imported_data$dataset,reactiveprogramdata$imported_data$ppm,reactiveprogramdata$imported_data$freq,reactiveprogramdata$imported_data$export_path,reactiveprogramdata$imported_data$Metadata,reactiveprogramdata$imported_data$repository,reactiveprogramdata$originaldataset)
  #   }
  # })

  #Plot where quantification of model spectrum amnd p values for every bucket are stored
  output$automatic_profiling_plot <- renderPlotly({
    if (reactiveprogramdata$beginning==FALSE) return()
    reactiveprogramdata$automatic_profiling_plot
  })



  ## SECOND TAB REACTIVE OUTPUTS

  #Selection of ROI


  # ROI parameters are loaded when a ROI is selected
  observeEvent(input$select, {
    if (reactiveprogramdata$beginning==FALSE) return()
	#Splitting of ROI data into individual ROIs to be quantified
    tryCatch({
      reac$cho=NA
    dummy = which(is.na(reactiveprogramdata$ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(reactiveprogramdata$ROI_data)[1]+1
    lal=which(duplicated(reactiveprogramdata$ROI_data[-dummy,1:2])==FALSE)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(reactiveprogramdata$ROI_data[-dummy,])[1]))
      reactiveprogramdata$ROIdata_subset=reactiveprogramdata$ROI_data[ROI_separator[as.numeric(input$select), 1]:ROI_separator[as.numeric(input$select), 2],]
    reactiveROItestingdata$signpar <- rbind(rep(NA,7),rep(NA,7))
    colnames(reactiveROItestingdata$signpar)=c("intensity",	"chemical shift",	"half bandwidth",	"gaussian %",	"J coupling",	"multiplicities",	"roof effect")
    reactiveROItestingdata$qualitypar <- rbind(rep(NA,3),rep(NA,3))
    rownames(reactiveROItestingdata$qualitypar)=NULL
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification (arbitrary unit)','fitting_error','signal/total area ratio')

    # Plot is prepared
    ROI_limits=c(reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,1]))],reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,2]))])
    ind=ifelse(is.null(input$x1_rows_selected),1,input$x1_rows_selected)
    dummy=type_plot(reactiveprogramdata$imported_data,ROI_limits,ind,reactiveprogramdata$median_plot,reactiveprogramdata$clusterplot)
    if (!is.null(dummy)) reactiveprogramdata$plot=dummy
    }, error = function(e) {
      print('Error. Please explain the issue on Github website.')
      return()
    })
	#Setting of reactivity to edition of parameters prepared. Probably improvable, but it is a delicate matter. Now that I can debug it easily, it can be 'cleaned' in the future.
    reactiveprogramdata$info=c()


	#Analysis of ROI edition. If edition is not correct (for example, there are characters in a numeric input), the edition is rejected and shown with red colour. If correct, the change is accepted with green colour.
	#TODO: it seems sometiems the edition fails if the change was too quick. Revise possible ways to control it.
    output$ROIdata <-   DT::renderDataTable(  reactiveprogramdata$ROIdata_subset, selection = 'none', rownames = FALSE,editable=T)

    proxy_ROIdata = dataTableProxy('ROIdata')

    observeEvent(input$ROIdata_cell_edit, {
      info2 = input$ROIdata_cell_edit
      i2 = info2$row
      j2 = info2$col + 1
      v2 = info2$value
      # if (!is.na(reac$cho)) {
        reactiveprogramdata$ROIdata_subset[i2, j2] <<- DT:::coerceValue(v2, reactiveprogramdata$ROIdata_subset[i2, j2])
      replaceData(proxy_ROIdata, reactiveprogramdata$ROIdata_subset, resetPaging = FALSE, rownames = FALSE)
      # }
      reac$cho=1
      })
  })

  #Selection of spectra, or of cluster or median plots
  observeEvent(input$x1_rows_selected, {
    if (reactiveprogramdata$beginning==FALSE) return()
    tryCatch({
      reset("fit_selection_cell_clicked")

    if (reactiveprogramdata$beginning ==TRUE) {
	dummy = which(is.na(reactiveprogramdata$ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(reactiveprogramdata$ROI_data)[1]+1
    lal=which(duplicated(reactiveprogramdata$ROI_data[-dummy,1:2])==FALSE)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(reactiveprogramdata$ROI_data[-dummy,])[1]))

      if (input$select=="") {
       print("Please load the session again and wait until all the loading process is finished")
        return()
        }

  reactiveprogramdata$ROIdata_subset=reactiveprogramdata$ROI_data[ROI_separator[as.numeric(input$select), 1]:ROI_separator[as.numeric(input$select), 2],]
    }
	#Reset of parameters
    reactiveprogramdata$stop=reactiveprogramdata$stop2=0

    reactiveprogramdata$info=c()
    reactiveROItestingdata$signpar <- rbind(rep(NA,7),rep(NA,7))
    colnames(reactiveROItestingdata$signpar)=c("intensity",	"chemical_shift",	"half bandwidth",	"gaussian",	"J_coupling",	"multiplicities",	"roof effect")
    reactiveROItestingdata$qualitypar <- rbind(rep(NA,3),rep(NA,3))
    rownames(reactiveROItestingdata$qualitypar)=NULL

    colnames(reactiveROItestingdata$qualitypar)=c('Quantification (arbitrary unit)','fitting_error','signal/total area ratio')
    ROI_limits=c(reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,1]))],reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,2]))])

	#Generation of plot
    dummy=type_plot(reactiveprogramdata$imported_data,ROI_limits,input$x1_rows_selected,reactiveprogramdata$median_plot,reactiveprogramdata$clusterplot)
    if (!is.null(dummy)) reactiveprogramdata$plot=dummy
    },error=function(e) {
      print("Problem. Please explain the issue in the Github page")
    })

     })

  #Individual quantification
  observeEvent(input$action, {
    tryCatch({
	#Only when there is a spectrum selected and with number >2 (to avoid cluster and median plot) the quantifications is performed, with prior adaptation to row of dataset
    if(length(reactiveprogramdata$info)==0) reactiveprogramdata$ind=input$x1_rows_selected-2
    if (length(reactiveprogramdata$ind)!=1|reactiveprogramdata$ind<1) {
      print('Select one valid spectrum')
      return(NULL)
    }

	#The automatic quantification
    reactivequantdata$method1 <- tryCatch({individual_profiling(reactiveprogramdata$imported_data, reactiveprogramdata$final_output, reactiveprogramdata$ind,reactiveprogramdata$ROIdata_subset,reactiveprogramdata$reproducibility_data,interface=TRUE)}, warning = function(w) {},error=function(e) {
      print("There was a problem. Check the compatibility between the ROI to fit and the current ROI profiles, or try setting a wider ROI.")
      return(NULL)
      })

	#Update of tables of tab
    if (!is.null(reactivequantdata$method1$results_to_save)) {
      reactiveprogramdata$plot=reactivequantdata$method1$p
      #reactivequantdata$stop3=1
      reactiveROItestingdata$qualitypar=cbind(reactivequantdata$method1$results_to_save$quantification,round(reactivequantdata$method1$results_to_save$fitting_error,reactivequantdata$method1$results_to_save$signal_area_ratio),4)
      colnames(reactiveROItestingdata$qualitypar)=c('Quantification (arbitrary unit)','Fitting Error','Signal/total area ratio')
      ind=which(reactiveprogramdata$ROIdata_subset[,5]>0)+3

      rownames(reactiveROItestingdata$qualitypar)=rownames(reactivequantdata$method1$plot_data)[ind]
      if (!is.null(reactivequantdata$method1$signals_parameters)) reactiveROItestingdata$signpar <- t(reactivequantdata$method1$signals_parameters)
      reactiveprogramdata$stop=0
         # reactiveprogramdata$roi=1
    }
    },error=function(e) {
      print("Problem. Please explain the issue in the Github page")
    })
  })



  #Quantification of all spectra in the ROI:
  observeEvent(input$automatic_profiling_signal, {
    tryCatch({
    dummy <- individual_profiling(reactiveprogramdata$imported_data, reactiveprogramdata$final_output, seq(nrow(reactiveprogramdata$imported_data$dataset)),reactiveprogramdata$ROIdata_subset,reactiveprogramdata$reproducibility_data,interface=TRUE)
    reactiveprogramdata$final_output=dummy$final_output
    reactiveprogramdata$reproducibility_data=dummy$reproducibility_data
    reactiveprogramdata$validation_data=validation(reactiveprogramdata$final_output,input$select_validation,reactiveprogramdata$validation_data$alarm_matrix)

    },error=function(e) {print('Error. Please explain the issue on the Github website')})
    })



  #Remove quantification or save quantification or ROI profile edited
  observeEvent(input$remove_q, {
    tryCatch({
    if (!is.null(reactiveprogramdata$imported_data$signals_names[reactiveprogramdata$info$col])) {
      ind=which(reactiveprogramdata$ROI_data[,4]==reactiveprogramdata$imported_data$signals_names[reactiveprogramdata$info$col])
    } else {
      ind=as.numeric(input$select)
    }
    reactiveprogramdata$final_output <- remove_quant(reactiveprogramdata$info,reactiveprogramdata$imported_data, reactiveprogramdata$final_output)
    },error=function(e) {print('Error. Please explain the issue on the Github website')})

    })

  observeEvent(input$save_results, {
    tryCatch({
    if (is.null(reactivequantdata$method1$Ydata)) return(NULL)

    dummy=save_roi_testing(reactivequantdata$method1,reactiveprogramdata$imported_data, reactiveprogramdata$final_output,reactiveprogramdata$reproducibility_data)
    reactiveprogramdata$final_output=dummy$final_output
    reactiveprogramdata$reproducibility_data=dummy$reproducibility_data
    reactiveprogramdata$validation_data=validation(reactiveprogramdata$final_output,input$select_validation,reactiveprogramdata$validation_data$alarm_matrix)

    },error=function(e) {print('Error. Please explain the issue on the Github website')})
      })



    #Save edition of ROI profile
  observeEvent(input$save_profile, {
    tryCatch({
  dummy = which(is.na(reactiveprogramdata$ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(reactiveprogramdata$ROI_data)[1]+1
    lal=which(duplicated(reactiveprogramdata$ROI_data[-dummy,1:2])==FALSE)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(reactiveprogramdata$ROI_data[-dummy,])[1]))
    if (length(reactiveprogramdata$info$col)>0) {

      ind=which(ROI_separator[,2]-reactiveprogramdata$info$col>=0)[1]
    } else {
      ind=as.numeric(input$select)
    }

    reactiveprogramdata$ROI_data[ROI_separator[ind, 1]:ROI_separator[ind, 2],]=reactiveprogramdata$ROI_data_check[ROI_separator[ind, 1]:ROI_separator[ind, 2],]=reactiveprogramdata$ROIdata_subset
    ROI_names=paste(reactiveprogramdata$ROI_data[ROI_separator[, 1],1],reactiveprogramdata$ROI_data[ROI_separator[, 1],2])
    names(reactiveprogramdata$select_options)=ROI_names
    print("ROI Profile saved.")
    },error=function(e) {print('Error. Please explain the issue on the Github website')})
    })



  #Spectra table.
  #Spectra table.
  # output$x1 = tryCatch({DT::renderDataTable(reactiveprogramdata$spectra , selection = list(mode = 'multiple', selected = 1),server = TRUE)},error=function(e){})
  tryCatch({
  output$x1 = DT::renderDataTable(datatable(reactiveprogramdata$spectra , selection = list(mode = 'multiple', selected = 1))%>% formatStyle(0,target="row",color=styleEqual(1:2,c("red","blue")))
  )},error=function(e){})



	#Plotly figure where to analyze peak shape fitting
  output$plot <- renderPlotly({
    if (reactiveprogramdata$beginning==FALSE) return()
    if (reactiveprogramdata$beginning==T) reactiveprogramdata$plot
  })



  #Table where to analyze quantifications
  output$qualitypar = DT::renderDataTable(reactiveROItestingdata$qualitypar)


  #Repository table
  observe({
    suppressWarnings(
    if (!is.na(reactiveprogramdata$imported_data))  {
  output$repository = DT::renderDataTable(
    reactiveprogramdata$imported_data$repository[which(reactiveprogramdata$imported_data$repository[,3]>reactiveprogramdata$ROIdata_subset[1,2]&reactiveprogramdata$imported_data$repository[,3]<reactiveprogramdata$ROIdata_subset[1,1]),] , server = TRUE)
    })
  })
  observe({
    suppressWarnings(
      if (!is.na(reactiveprogramdata$imported_data))  {
        output$repository2 = DT::renderDataTable(
      reactiveprogramdata$imported_data$repository , selection = list(mode = 'single', selected = NULL),filter='top', server = TRUE)
      })
  })


  #Direct edition of parameters before quantification
  output$directedition <-   DT::renderDataTable(reactiveROItestingdata$signpar, selection = 'none', rownames = T,editable=T)

  proxy_directedition = dataTableProxy('directedition')

  observeEvent(input$directedition_cell_edit, {
    info = input$directedition_cell_edit
    i = info$row
    j = info$col
    v = info$value
    reactiveROItestingdata$signpar[i, j] <<- DT:::coerceValue(v, reactiveROItestingdata$signpar[i, j])
    replaceData(proxy_directedition, reactiveROItestingdata$signpar, resetPaging = FALSE, rownames = FALSE)
  })
  #Quantification after direct edition of paramters
  observeEvent(input$direct_edition, {
    if (all(is.na(reactiveROItestingdata$signpar))) {
      print("You can only perform direct edition of line shape fitted ROIs")
      return(NULL)
    }
    tryCatch({
    reactivequantdata$method1=signals_int(reactiveprogramdata$imported_data, reactiveprogramdata$final_output,reactiveprogramdata$ind,reactiveROItestingdata$signpar,reactiveprogramdata$ROIdata_subset)
    reactiveprogramdata$plot=reactivequantdata$method1$p
    #reactivequantdata$stop3=1
    reactiveROItestingdata$qualitypar=cbind(reactivequantdata$method1$results_to_save$quantification,reactivequantdata$method1$results_to_save$fitting_error,reactivequantdata$method1$results_to_save$signal_area_ratio)
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification (arbitrary unit)','Fitting Error','Signal/total area ratio')
    rownames(reactiveROItestingdata$qualitypar)=rownames(reactivequantdata$method1$plot_data)[4:(3+nrow(reactiveROItestingdata$qualitypar))]
    # ind=which(rownames(reactiveROItestingdata$qualitypar)=='additional signal')
    # reactiveprogramdata$bgColScales = rep(c("","info"),times=c(length(rownames(reactiveROItestingdata$qualitypar))-length(ind),length(ind)))
    if (!is.null(reactivequantdata$method1$signals_parameters)) reactiveROItestingdata$signpar <- t(reactivequantdata$method1$signals_parameters)
  },error=function(e) {
    print("Prepare a valid input")
  })
 })



  #2D plot
  output$jres_plot <- renderUI({
    if(is.na(reactiveprogramdata$jres_plot)) return()
    plotlyOutput("jres_plot2",height='250px')
  })
  output$jres_plot2 <- renderPlotly({
    print(reactiveprogramdata$jres_plot)
  })




  ## THIRD TAB REACTIVE OUTPUTS

  #Creation of table to check quantifications with the parameter chosen by the user
  observeEvent(input$select_validation, {
    if (reactiveprogramdata$beginning==FALSE) return()
	if (as.numeric(input$select_validation)>0) {
    tryCatch({
      reactiveprogramdata$validation_data=validation(reactiveprogramdata$final_output,input$select_validation,reactiveprogramdata$validation_data$alarm_matrix)
    output$fit_selection = DT::renderDataTable({ datatable(round(reactiveprogramdata$validation_data$shown_matrix,4),selection = list(mode = 'single', target = 'cell')) %>% formatStyle(colnames(reactiveprogramdata$validation_data$shown_matrix), backgroundColor = styleInterval(reactiveprogramdata$validation_data$brks, reactiveprogramdata$validation_data$clrs))
    })},error=function(e) {
      print("Not enough data to model it.")
    })
	}
  })

  #Loading of quantification
  observeEvent(input$fit_selection_cell_clicked, {
if (length(input$fit_selection_cell_clicked)<1) return()
    #Checks and setting of parameters
    tryCatch({reactiveprogramdata$info=input$fit_selection_cell_clicked
    reactiveprogramdata$ind=reactiveprogramdata$info$row
    reac$cho=NA

    updateSelectInput(session, "select",selected = NULL)

    # if (length(reactiveprogramdata$info$row)!=1) return(NULL)

    dummy=load_quantification(reactiveprogramdata$reproducibility_data,reactiveprogramdata$imported_data,reactiveprogramdata$final_output,reactiveprogramdata$info,reactiveprogramdata$ROI_data)

      reactiveprogramdata$plot=dummy$plot
      reactiveROItestingdata$signpar=dummy$signpar
      reactiveprogramdata$ROIdata_subset=dummy$ROIpar
      reactiveROItestingdata$qualitypar=dummy$qualitypar

    #Redirect to quantification tab
    updateTabsetPanel(session, "mynavlist",selected = "tab3")
    },error=function(e) {
      print("Select a valid quantification.")
    })
  })


  ## FOURTH TAB REACTIVE OUTPUTS


  observeEvent(input$roi_profile_option, {
  tryCatch({
      output$roi_profiles_plot=renderPlotly({
        if (input$roi_profile_option==1) {
          reactiveprogramdata$clusterplot
          } else if (input$roi_profile_option==2) {
            reactiveprogramdata$median_plot
            }
        })
  },error=function(e) {print('Error. Please explain the issue on the Github website')})


    })

  #Add and remove signals and save changes
  observeEvent(input$add_hmdb_signal, {
    tryCatch({
    dummy=reactiveprogramdata$imported_data$repository[input$repository2_rows_selected,]
    dummy=c(dummy[,3]+0.02,dummy[,3]-0.02,'Baseline Fitting',dummy[,1],1,dummy[,3],median(reactiveprogramdata$ROI_data_check[,7]),median(reactiveprogramdata$ROI_data_check[,8]),dummy[,4],dummy[,5],0,dummy[,2])
    if (dummy[9]=='d') {
      dummy[9]=2
    } else if (dummy[9]=='t') {
      dummy[9]=3
    } else {
      dummy[9]=1
    }

    if (is.na(as.numeric(dummy[10])))  dummy[10]=0
    dummy=as.list(dummy)
    dummy[-c(3,4,12)]=as.numeric(dummy[-c(3,4,12)])

    reactiveprogramdata$ROI_data_check=rbind(reactiveprogramdata$ROI_data_check,dummy)
    }, error = function(e) {
      print('Error. Please revise that you have chosen a signal.')
    })
  })
  observeEvent(input$open_hmdb_url, {
    tryCatch({
      browseURL(paste("http://www.hmdb.ca/metabolites/",reactiveprogramdata$imported_data$repository[input$repository2_rows_selected,2],"#spectra",sep=''))    }, error = function(e) {
      print('Not possible to load HMDB url.')
    })
  })
  observeEvent(input$add_signal, {
    tryCatch({
      dummy=c(rep(NA,4),apply(reactiveprogramdata$ROI_data_check[,5:11],2,function(x)median(x,na.rm=T)),NA)
      reactiveprogramdata$ROI_data_check=rbind(reactiveprogramdata$ROI_data_check,dummy)
    }, error = function(e) {
      print('Error. Please explain the issue in the Github page.')
    })
     })
  observeEvent(input$remove_signal, {
    tryCatch({
    reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data_check[-input$roi_profiles_rows_selected,]
    reactiveprogramdata$list=setdiff(reactiveprogramdata$list,input$roi_profiles_rows_selected)
    }, error = function(e) {
      print('Error. Please explain the issue in the Github page.')
    })
     })
  observeEvent(input$save_changes, {
    tryCatch({
      if (any(duplicated(reactiveprogramdata$ROI_data_check[,4:5])==T)) {
        print("Revise duplicated signal IDs.")
        return(NULL)
      }

      trel=which(order(reactiveprogramdata$ROI_data_check[,1]) %in% seq(length(reactiveprogramdata$list))==T)
      sA <- apply(reactiveprogramdata$ROI_data[reactiveprogramdata$list,4:5],1,paste,collapse=' ')
      sB <- apply(reactiveprogramdata$ROI_data_check[seq(length(reactiveprogramdata$list)),4:5],1,paste,collapse=' ')
      knoc=reactiveprogramdata$ROI_data[which((sB %in% sA)==F),4:5]
      knoc2=reactiveprogramdata$ROI_data_check[which((sA %in% sB)==F),4:5]
      reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data_check[order(reactiveprogramdata$ROI_data_check[,1]),]


      new_validation_data=rep(list(matrix(NA,nrow(reactiveprogramdata$final_output$signal_area_ratio),nrow(reactiveprogramdata$ROI_data_check),dimnames=list(reactiveprogramdata$imported_data$Experiments,paste(reactiveprogramdata$ROI_data_check[,4],reactiveprogramdata$ROI_data_check[,5],sep='_')))), length(reactiveprogramdata$validation_data$alarm_matrix))
      names(new_validation_data)=names(reactiveprogramdata$validation_data$alarm_matrix)
      new_shown_matrix=new_fitting_error=new_intensity=new_signal_area_ratio=new_shift=new_width=new_Area=matrix(NA,nrow(reactiveprogramdata$final_output$signal_area_ratio),nrow(reactiveprogramdata$ROI_data_check),dimnames=list(reactiveprogramdata$imported_data$Experiments,make.names(paste(reactiveprogramdata$ROI_data_check[,4],reactiveprogramdata$ROI_data_check[,5],sep='_'))))
      new_useful_data=reactiveprogramdata$reproducibility_data
      for (i in 1:length(new_useful_data)) new_useful_data[[i]]=vector("list", nrow(reactiveprogramdata$ROI_data_check))
      new_fitting_error[,trel]=reactiveprogramdata$final_output$fitting_error[,reactiveprogramdata$list]
      new_intensity[,trel]=reactiveprogramdata$final_output$intensity[,reactiveprogramdata$list]
      new_signal_area_ratio[,trel]=reactiveprogramdata$final_output$signal_area_ratio[,reactiveprogramdata$list]
      new_shift[,trel]=reactiveprogramdata$final_output$chemical_shift[,reactiveprogramdata$list]
      new_width[,trel]=reactiveprogramdata$final_output$half_bandwidth[,reactiveprogramdata$list]
      new_Area[,trel]=reactiveprogramdata$final_output$quantification[,reactiveprogramdata$list]
      new_shown_matrix[,trel]=as.matrix(reactiveprogramdata$validation_data$shown_matrix[,reactiveprogramdata$list])
      for (j in 1:length(new_validation_data)) new_validation_data[[j]][,trel]=as.matrix(reactiveprogramdata$validation_data$alarm_matrix[[j]][,reactiveprogramdata$list])
      for (j in 1:length(new_useful_data)) new_useful_data[[j]][trel]=reactiveprogramdata$reproducibility_data[[j]][reactiveprogramdata$list]
      # for (j in 1:length(new_useful_data)) new_useful_data[[j]][setdiff(seq(nrow(reactiveprogramdata$ROI_data_check)),trel)]=list(Ydata=NA,Xdata=NA,ROI_profile=NA,program_parameters=NA,plot_data=NA,FeaturesMatrix=NA,signals_parameters=NA,results_to_save=NA,error1=1000000)

      inde=c()
      for (j in seq(nrow(knoc))) {
        for (k in 1) {
          for (l in trel) {
            dsd=any(apply(new_useful_data[[k]][[l]]$ROI_profile[,4:5],1,function(x)identical(paste(x,collapse='_'),paste(knoc[j,],collapse='_')))  ==T)
            dsd2=length(grep(paste(knoc[j,],collapse='_'),rownames(new_useful_data[[k]][[l]]$plot_data))>0)==T
            if (dsd==T||dsd2==T) inde=c(inde,l)
          }
        }
      }
      for (j in seq(nrow(knoc))) {
        for (k in 1:length(new_useful_data)) {
          for (l in seq_along(inde)) {
            new_useful_data[[k]][[inde[l]]]$ROI_profile[apply(new_useful_data[[k]][[inde[l]]]$ROI_profile[,4:5],1,function(x)identical(paste(x,collapse='_'),paste(knoc[j,],collapse='_'))),4:5]=knoc2[j,]
            if (!is.null(new_useful_data[[k]][[inde[l]]]$plot_data)) {
            rownames(new_useful_data[[k]][[inde[l]]]$plot_data)=gsub(paste(knoc[j,],collapse='_'),paste(knoc2[j,],collapse='_'),rownames(new_useful_data[[k]][[inde[l]]]$plot_data))
            }
            }
        }
      }



      reactiveprogramdata$final_output$fitting_error=new_fitting_error
      reactiveprogramdata$final_output$intensity=new_intensity
      reactiveprogramdata$final_output$signal_area_ratio=new_signal_area_ratio
      reactiveprogramdata$final_output$chemical_shift=new_shift
      reactiveprogramdata$final_output$half_bandwidth=new_width
      reactiveprogramdata$final_output$quantification=new_Area
      reactiveprogramdata$reproducibility_data=new_useful_data
      reactiveprogramdata$ROI_data=reactiveprogramdata$ROI_data_check
      reactiveprogramdata$list=seq(nrow(reactiveprogramdata$ROI_data))
      reactiveprogramdata$imported_data$signals_names=paste(reactiveprogramdata$ROI_data[,4],reactiveprogramdata$ROI_data[,5],sep='_')
      reactiveprogramdata$validation_data$alarm_matrix=new_validation_data
      reactiveprogramdata$validation_data$shown_matrix=new_shown_matrix

      dummy = which(is.na(reactiveprogramdata$ROI_data[, 1]))
      if (length(dummy)==0) dummy=dim(reactiveprogramdata$ROI_data)[1]+1
      lal=which(duplicated(reactiveprogramdata$ROI_data[-dummy,1:2])==FALSE)
      ROI_separator = cbind(lal, c(lal[-1] - 1, dim(reactiveprogramdata$ROI_data[-dummy,])[1]))
      ROI_names=paste(reactiveprogramdata$ROI_data[ROI_separator[, 1],1],reactiveprogramdata$ROI_data[ROI_separator[, 1],2])
      reactiveprogramdata$select_options=1:length(ROI_names)
      names(reactiveprogramdata$select_options)=ROI_names
      print('ROI profiles modified.')

    updateSelectInput(session, "select",
      choices = reactiveprogramdata$select_options,selected = input$x1_rows_selected
    )
    }, error = function(e) {
      print('Problem during update of data. Please explain the issue in the Github page.')
})
  })

  #ROI Profiles table
  output$roi_profiles  = DT::renderDataTable(reactiveprogramdata$ROI_data_check, selection = 'multiple', rownames = FALSE,options=list(pageLength = 100,
                                             lengthMenu = c(50,100)),editable=T)

  proxy_roi_profiles = dataTableProxy('roi_profiles')

  observeEvent(input$roi_profiles_cell_edit, {
    info = input$roi_profiles_cell_edit
    i = info$row
    j = info$col + 1
    v = info$value
    reactiveprogramdata$ROI_data_check[i, j] <<- DT:::coerceValue(v, reactiveprogramdata$ROI_data_check[i, j])
    replaceData(proxy_roi_profiles, reactiveprogramdata$ROI_data_check, resetPaging = FALSE, rownames = FALSE)
  })
#
#   ## FIFTH TAB REACTIVE OUTPUTS
#
#   #Boxplot plot
#   output$plot_p_value_2 <- renderPlotly({
#     if(all(is.na(reactiveprogramdata$final_output$quantification))) return()
#     tryCatch({type_analysis_plot(reactiveprogramdata$final_output$quantification,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,type='boxplot')
#     }, error = function(e) {
#       print('Error. Please explain the issue in the Github page.')
#       return(NULL)
#     })
#   })
#   #PCA plot
#   output$pcascores <- renderPlotly({
#     if(all(is.na(reactiveprogramdata$final_output$quantification))) return()
#
#     tryCatch({type_analysis_plot(reactiveprogramdata$final_output$quantification,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,type='pca')
#   }, error = function(e) {
#     print('Generation of Regions of Interest not possible. Please explain the issue in the Github page.')
#     return(NULL)
#   })
#   })

  ## SIXTH TAB REACTIVE OUTPUTS

# STOCSY generation
  observeEvent(input$stocsy, {
  tryCatch({
    left_ppm <- renderText({ input$left_ppm })
    right_ppm <- renderText({ input$right_ppm })
      output$stocsy_plot=renderPlotly({
        if (input$stocsy==1) {
          reactiveprogramdata$clusterplot
          } else {
        identification_tool(reactiveprogramdata$imported_data$dataset,reactiveprogramdata$imported_data$ppm,c(input$left_ppm,input$right_ppm),input$correlation_method)
          }
        })
  },error=function(e) {return(NULL) })
    })



  #Dengrogran heatmaps for quantification and chemical $chemical_shift
  tryCatch(output$dendheatmapareadata <- renderPlotly({
    if(all(is.na(reactiveprogramdata$final_output$quantification))) return()

    tryCatch({suppressWarnings(suppressMessages(type_analysis_plot(reactiveprogramdata$final_output$quantification,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,type='dendrogram_heatmap')))
  }, error = function(e) {
    print('Error. Please explain the issue in the Github page.')
    return(NULL)
  })
  }))

  tryCatch(output$dendheatmapshiftdata <- renderPlotly({
    if(all(is.na(reactiveprogramdata$final_output$quantification))) return()

    tryCatch({suppressWarnings(suppressMessages(type_analysis_plot(reactiveprogramdata$final_output$chemical_shift,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,type='dendrogram_heatmap')))
  }, error = function(e) {
    print('Error. Please explain the issue in the Github page.')
    return(NULL)
  })
}))

  roifunc <- function(ROI_data,Metadata,Experiments) {
    dummy = which(is.na(ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(ROI_data)[1]+1
    lal=which(duplicated(ROI_data[-dummy,1:2])==F)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(ROI_data[-dummy,])[1]))

    ROI_names=paste(ROI_data[ROI_separator[, 1],1],ROI_data[ROI_separator[, 1],2])
    select_options=seq_along(ROI_names)
    names(select_options)=ROI_names
    mm=matrix(NA,2,dim(Metadata)[2])
    colnames(mm)=colnames(Metadata)
    spectra=cbind(c('Exemplars','Median Spectrum per group',Experiments),rbind(mm,Metadata))
    colnames(spectra)=c('spectrum',colnames(mm))
    dummy=list(select_options=select_options,spectra=spectra)
    return(dummy)
  }

}



