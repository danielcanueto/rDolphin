server = function(input, output,session) {

  #Increase of maximum memory size that can be uploaded
  options(shiny.maxRequestSize=1000*1024^2)
  options(warn =-1)
  #Setting of reactive parameters in the Shiny GUI
  reactiveROItestingdata <- reactiveValues(signpar = matrix(NA,2,7,dimnames = list(seq(2),c("intensity",	"chemical shift",	"half bandwidth",	"gaussian %",	"J coupling",	"multiplicities",	"roof effect"))),qualitypar = matrix(NA,2,3,dimnames=list(seq(2),c('Quantification','fitting_error','signal/total area ratio'))))
  reactivequantdata <- reactiveValues(method1=NA)
 reactiveprogramdata <- reactiveValues(ROIdata_subset=NA,ind=NA,beginning=FALSE,dataset=NA,final_output=list(),useful_data=list(),imported_data=NA,p_value_final=NA,ROI_data=NA,ROI_data_check=NA,info=c(),select_options=NA,new_roi_profile=NA,p=NA,bgColScales=NA,autorun_plot=NA,ROI_names=NA,clusterplot=NA,medianplot=NA,jres_plot=NA)

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

    session$onSessionEnded(stopApp)



	#Imported data is loaded to 'dummy'. Only after the check that the parameters are correct, they are stored in 'reactiveprogramdata'.
	# dummy=list(imported_data=NA,autorun_plot=NA,select_options=NA,spectra=NA,clusterplot=NA,medianplot=NA,beginning=FALSE,jres_plot=NA)
	# Import of data
	# dummy = tryCatch({helperimport(reactiveprogramdata$inFile$datapath,dummy)
	  # }, error = function(e) {
	# print('Error. Please explain the issue in the Github page if necessary.')
	# return(dummy)
	# })
	tryCatch({
	reactiveprogramdata$imported_data = tryCatch({suppressWarnings(import_data(reactiveprogramdata$inFile$datapath))}, error = function(e) {
	  print('Import of data not posible with current input')
	  return(NULL)
	})
	reset("file1")

  reactiveprogramdata$final_output=reactiveprogramdata$imported_data$final_output
	reactiveprogramdata$useful_data=reactiveprogramdata$imported_data$useful_data
	reactiveprogramdata$ROI_data=reactiveprogramdata$ROI_data_check=reactiveprogramdata$imported_data$ROI_data
	reactiveprogramdata$imported_data$final_output=reactiveprogramdata$imported_data$useful_data=reactiveprogramdata$imported_data$ROI_data=NULL
	colnames(reactiveprogramdata$ROI_data)=c("ROI left edge","ROI right edge","Quantification Mode","Metabolite","Quantification Signal","Chemical shift","Chemical shift tolerance","Half bandwidth","Multiplicity","J coupling","Roof effect","Relative intensity")

	dummy=tryCatch({profile_model_spectrum(reactiveprogramdata$imported_data,reactiveprogramdata$ROI_data)}, error = function(e) {
	print('Automatic quantification of model spectrum not possible.')
	})
	reactiveprogramdata$autorun_plot=dummy$p
	reactiveprogramdata$total_signals_parameters=dummy$total_signals_parameters
	# dummy$indicators=dummy2$indicators

	print('Generating additional information...')

	#plots of representative spectra and median spectra per group to help setting the right ROI parameters

    reactiveprogramdata$clusterplot=tryCatch({clustspectraplot(reactiveprogramdata$imported_data)  }, error = function(e) {
	print('Generation of subsets or representative spectra not possible.')
	})
    reactiveprogramdata$medianplot=tryCatch({medianplot(reactiveprogramdata$imported_data)
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

	# if (dummy$beginning==TRUE) {
	  # plo=names(sapply(dummy, names))
	  # for (i in 1:length(plo)) reactiveprogramdata[[plo[i]]]=dummy[[plo[i]]]

	 #Variables that can change during the use of the GUI are separated from 'imported_data'.
	output$moreControls <- renderUI({
	  if (reactiveprogramdata$beginning==T) selectInput("select",label=NULL,choices = reactiveprogramdata$select_options,selected = 1)
	})
	shinyjs::show('autorun')
	shinyjs::show('alignment')
	shinyjs::show('model_spectrum')
	#When the session is prepared, the tabs and some inputs become active
	print('Done!')

    # updateSelectInput(session, "select",choices = reactiveprogramdata$select_options,selected = 1)
    updateSelectInput(session, "select_validation",selected = 1)
    # session$sendCustomMessage('activeNavs', 'Individual Quantification')
    # session$sendCustomMessage('activeNavs', 'Quantification Validation')
    # session$sendCustomMessage('activeNavs', 'Uni and multivariate analysis')
    # session$sendCustomMessage('activeNavs', 'ROI Profiles')
    # session$sendCustomMessage('activeNavs', 'STOCSY and dendrogram heatmaps')
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
      rm(savedreactivedata)
reset("file2")
}, error = function(e) {
      print('Not possible to load the session. Please revise your choice.')
      return(NULL)
    })
    colnames(reactiveprogramdata$ROI_data)=c("ROI left edge","ROI right edge","Quantification Mode","Metabolite","Quantification Signal","Chemical shift","Chemical shift tolerance","Half bandwidth","Multiplicity","J coupling","Roof effect","Relative intensity")

    reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data
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
	reactiveprogramdata$beginning =TRUE
	shinyjs::show('autorun')
	shinyjs::show('alignment')
	shinyjs::show('model_spectrum')
	p=plot_ly(x=reactiveprogramdata$imported_data$ppm,y=reactiveprogramdata$dataset)
	# if (is.null(reactiveprogramdata$autorun_plot$dependencies))
	  reactiveprogramdata$autorun_plot$dependencies=reactiveprogramdata$clusterplot$dependencies=reactiveprogramdata$medianplot$dependencies=p$dependencies
# 	  reactiveprogramdata$autorun_plot$dependencies[[1]]=reactiveprogramdata$clusterplot$dependencies[[1]]=reactiveprogramdata$medianplot$dependencies[[1]]=list(version="1.11.3",src=list(file=file.path(system.file(package = "crosstalk"),"lib","jquery")))
#
#
# 	   reactiveprogramdata$autorun_plot$dependencies[[2]]=reactiveprogramdata$clusterplot$dependencies[[2]]=reactiveprogramdata$medianplot$dependencies[[2]]=list(version="1.0.0",src=list(file=file.path(system.file(package = "crosstalk"),"www")))
#
# 	  }
# 	   if (is.null(reactiveprogramdata$autorun_plot$dependencies[[3]])) {
# 	     reactiveprogramdata$autorun_plot$dependencies[[3]]=reactiveprogramdata$clusterplot$dependencies[[3]]=reactiveprogramdata$medianplot$dependencies[[3]]=list(version="0.1",src=list(file=file.path(system.file(package = "plotly"),"htmlwidgets/lib/typedarray")))
# }
	 	output$moreControls <- renderUI({
	  if (reactiveprogramdata$beginning==T)  {
	    print("Done!")

	    selectInput("select",label=NULL,choices = reactiveprogramdata$select_options,selected = 1)
	  }
	})
	  #When the session is loaded the tabs and some inputs become active
    updateSelectInput(session, "select_validation",selected = 1)
    # session$sendCustomMessage('activeNavs', 'Individual Quantification')
    # session$sendCustomMessage('activeNavs', 'Quantification Validation')
    # session$sendCustomMessage('activeNavs', 'Uni and multivariate analysis')
    # session$sendCustomMessage('activeNavs', 'ROI Profiles')
    # session$sendCustomMessage('activeNavs', 'STOCSY and dendrogram heatmaps')
    # updateTabsetPanel(session, "mynavlist",selected = "tab3")

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
    write_plots(input$caption,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,reactiveprogramdata$useful_data)},
    error= function(e) {       print('Not possible to generate the plot folder. Please check that you have permissions for the path specified.')

})})



  #Appearance of autorun and aligment buttons only after beginning or loading session
  # output$varselect <- renderUI({
  #   if(reactiveprogramdata$beginning==FALSE){return()}
  #   actionButton('autorun', 'Autorun all spectra', class = "btn-primary")
  # })
  # output$align_button <- renderUI({
  #   if(reactiveprogramdata$beginning==FALSE){return()}
  #   actionButton('alignment', 'Alignment of signals')
  # })
  # output$model_button <- renderUI({
  #   if(reactiveprogramdata$beginning==FALSE){return()}
  #   actionButton('model_spectrum', 'Profile model spectrum again')
  # })
  output$sp = DT::renderDataTable(
    reactiveprogramdata$total_signals_parameters , selection = list(selected = NULL),server = TRUE)

#Removed by now because it is unstable.
 # output$peak_analysis <- renderUI({
  #   if(reactiveprogramdata$beginning==FALSE){return()}
  #   actionButton('peak_analysis', 'Peak analysis')
  # })

   #Automatic quantification of all ROIs in all spectra
  observeEvent(input$autorun, {
    tryCatch({
    quantification_variables = autorun(reactiveprogramdata$imported_data, reactiveprogramdata$final_output,reactiveprogramdata$useful_data,reactiveprogramdata$ROI_data)
    reactiveprogramdata$final_output=quantification_variables$final_output
    reactiveprogramdata$useful_data=quantification_variables$useful_data
    },
    error = function(e) {
      print('Error. Please explain the issue in the Github page if necessary.')
      quantification_variables=NA
      return(quantification_variables)
    })
  })
  #Alignment of signals
  tryCatch({observeEvent(input$alignment, {
    reactiveprogramdata$imported_data$dataset= alignment(reactiveprogramdata$imported_data$dataset,reactiveprogramdata$imported_data$buck_step)
      reactiveprogramdata$clusterplot=clustspectraplot(reactiveprogramdata$imported_data)
      reactiveprogramdata$medianplot=medianplot(reactiveprogramdata$imported_data)
      }	)}, error = function(e) {
	print('Error during alignment. Please explain the issue in the Github page if necessary.')
	return(NULL)

  })
   #Automatic quantification of all ROIs in all spectra
  observeEvent(input$model_spectrum, {
    tryCatch({
    dummy=profile_model_spectrum(reactiveprogramdata$imported_data,reactiveprogramdata$ROI_data)
    reactiveprogramdata$autorun_plot=dummy$p
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
  output$autorun_plot <- renderPlotly({
    if (reactiveprogramdata$beginning==FALSE) return()
    reactiveprogramdata$autorun_plot
  })



  ## SECOND TAB REACTIVE OUTPUTS

  #Selection of ROI


  # ROI parameters are loaded when a ROI is selected
  observeEvent(input$select, {
    if (reactiveprogramdata$beginning==FALSE) return()
	#Splitting of ROI data into individual ROIs to be quantified
    tryCatch({
      reset("fit_selection_cell_clicked")

    dummy = which(is.na(reactiveprogramdata$ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(reactiveprogramdata$ROI_data)[1]+1
    lal=which(duplicated(reactiveprogramdata$ROI_data[-dummy,1:2])==FALSE)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(reactiveprogramdata$ROI_data[-dummy,])[1]))
      reactiveprogramdata$ROIdata_subset=reactiveprogramdata$ROI_data[ROI_separator[as.numeric(input$select), 1]:ROI_separator[as.numeric(input$select), 2],]
    reactiveROItestingdata$ROIpar <- reactiveprogramdata$ROIdata_subset
    reactiveROItestingdata$signpar <- rbind(rep(NA,7),rep(NA,7))
    colnames(reactiveROItestingdata$signpar)=c("intensity",	"chemical shift",	"half bandwidth",	"gaussian %",	"J coupling",	"multiplicities",	"roof effect")
    reactiveROItestingdata$qualitypar <- rbind(rep(NA,3),rep(NA,3))
    rownames(reactiveROItestingdata$qualitypar)=NULL
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','fitting_error','signal/total area ratio')

    # Plot is prepared
    ROI_limits=c(reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,1]))],reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,2]))])
    ind=ifelse(is.null(input$x1_rows_selected),1,input$x1_rows_selected)
    dummy=type_plot(reactiveprogramdata$imported_data,ROI_limits,ind,reactiveprogramdata$medianplot,reactiveprogramdata$clusterplot)
    if (!is.null(dummy)) reactiveprogramdata$plot=dummy
    }, error = function(e) {
      print('Error. Please explain the issue on Github website.')
      return()
    })
	#Setting of reactivity to edition of parameters prepared. Probably improvable, but it is a delicate matter. Now that I can debug it easily, it can be 'cleaned' in the future.
    reactiveprogramdata$change=reactiveprogramdata$change2=1
    reactiveprogramdata$stop=reactiveprogramdata$stop2=0

    reactiveprogramdata$info=c()

	  # Tables of ROI paramters ('ROIpar'), calculated deconvolution parameters ('signpar') and indicators of quantification ('qualitypar') are reset
	  resetInput(session, "ROIdata_edit")
    resetInput(session, "directedition_edit")

	#Analysis of ROI edition. If edition is not correct (for example, there are characters in a numeric input), the edition is rejected and shown with red colour. If correct, the change is accepted with green colour.
	#TODO: it seems sometiems the edition fails if the change was too quick. Revise possible ways to control it.
    output$ROIdata <- renderD3tf({
      tableProps <- list(
        btn_reset = FALSE,
        sort = TRUE,
        sort_config = list(
          sort_types = c("String", rep("Number", ncol(reactiveprogramdata$ROIdata_subset)))
        ))

      d3tf(reactiveROItestingdata$ROIpar,
           tableProps = tableProps,
           enableTf = FALSE,
           edit=TRUE,
           tableStyle = "table table-bordered")

    })
    observe({
      # if(is.null(input$fit_selection_cell_clicked)) {


      edit <- input$ROIdata_edit
      isolate({
        if (is.null(edit)) return()
        id <- edit$id
        row <- as.integer(edit$row)
        col <- as.integer(edit$col)
        val <- edit$val

        if(col == 0) {
          oldval <- rownames(reactiveROItestingdata$ROIpar)[row]
        } else if (col %in% c(1:2,6:11)){
          # numeric columns
          if(is.na(suppressWarnings(as.numeric(val)))) {
            oldval <- reactiveROItestingdata$ROIpar[row, col]
            rejectEdit(session, tbl = "ROIdata", row = row, col = col, id = id, value = oldval)
            # reactiveprogramdata$roi=0
            return(NULL)
          }
        } else if (col %in% c(3)) {
          if(!val %in% c('Clean Sum','Baseline Sum','Clean Fitting','Baseline Fitting')) {
            oldval <- reactiveROItestingdata$ROIpar[row, col]
            rejectEdit(session, tbl = "ROIdata", row = row, col = col, id = id, value = oldval)
            return(NULL)
          }
        }
        # if (reactiveprogramdata$change==1){
        #   reactiveprogramdata$change=0
        #   reactiveprogramdata$stop=1
        #   disableEdit(session, "ROIdata", c(1:11))
        # } else{
          if(col == 0) {
          } else if (col %in% c(1:2,5:11)) {
            reactiveROItestingdata$ROIpar[row, col] <- as.numeric(val)
            # reactiveprogramdata$roi=1
          } else if (col %in% c(3,4)) {
            reactiveROItestingdata$ROIpar[row, col] <- val
          }
          confirmEdit(session, tbl = "ROIdata", row = row, col = col, id = id, value = val)
        # }
      })
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
    # reactiveprogramdata$change=reactiveprogramdata$change2=1
    reactiveprogramdata$stop=reactiveprogramdata$stop2=0


    #reactivequantdata$stop3=0
        #reactiveprogramdata$roi=NULL
    reactiveprogramdata$info=c()
	resetInput(session, "directedition_edit")
    # reactiveROItestingdata$ROIpar <- reactiveprogramdata$ROIdata_subset
    reactiveROItestingdata$signpar <- rbind(rep(NA,7),rep(NA,7))
    colnames(reactiveROItestingdata$signpar)=c("intensity",	"shift",	"half_band_width",	"gaussian",	"J_coupling",	"multiplicities",	"roof_effect")
    reactiveROItestingdata$qualitypar <- rbind(rep(NA,3),rep(NA,3))
    rownames(reactiveROItestingdata$qualitypar)=NULL

    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','fitting_error','signal/total area ratio')
    ROI_limits=c(reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,1]))],reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,2]))])

	#Generation of plot
    dummy=type_plot(reactiveprogramdata$imported_data,ROI_limits,input$x1_rows_selected,reactiveprogramdata$medianplot,reactiveprogramdata$clusterplot)
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
    reactivequantdata$method1 <- tryCatch({not_automatic_quant(reactiveprogramdata$imported_data, reactiveprogramdata$final_output, reactiveprogramdata$ind,reactiveROItestingdata$ROIpar,reactiveprogramdata$useful_data,interface=TRUE)}, warning = function(w) {},error=function(e) {
      print("There was a problem.")
      return(NULL)
      })

	#Update of tables of tab
    if (!is.null(reactivequantdata$method1$results_to_save)) {
      reactiveprogramdata$plot=reactivequantdata$method1$p
      #reactivequantdata$stop3=1
      reactiveROItestingdata$qualitypar=cbind(reactivequantdata$method1$results_to_save$quantification,reactivequantdata$method1$results_to_save$fitting_error,reactivequantdata$method1$results_to_save$signal_area_ratio)
      colnames(reactiveROItestingdata$qualitypar)=c('Quantification','Fitting Error','Signal/total area ratio')
      ind=which(reactiveROItestingdata$ROI_par[,5]==1)+3

      rownames(reactiveROItestingdata$qualitypar)=rownames(reactivequantdata$method1$plot_data)[ind]
      # ind=which(rownames(reactiveROItestingdata$qualitypar)=='additional signal')

      # reactiveprogramdata$bgColScales = rep(c("","info"),times=c(length(rownames(reactiveROItestingdata$qualitypar))-length(ind),length(ind)))
      if (!is.null(reactivequantdata$method1$signals_parameters)) reactiveROItestingdata$signpar <- t(reactivequantdata$method1$signals_parameters)
      reactiveprogramdata$stop=0
         # reactiveprogramdata$roi=1
    }
    },error=function(e) {
      print("Problem. Please explain the issue in the Github page")
    })
  })



  #Quantification of all spectra in the ROI:
  observeEvent(input$autorun_signal, {
    tryCatch({
    dummy <- not_automatic_quant(reactiveprogramdata$imported_data, reactiveprogramdata$final_output, seq(nrow(reactiveprogramdata$imported_data$dataset)),reactiveROItestingdata$ROIpar,reactiveprogramdata$useful_data,interface=TRUE)
    reactiveprogramdata$final_output=dummy$final_output
    reactiveprogramdata$useful_data=dummy$useful_data

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

    dummy=save_roi_testing(reactivequantdata$method1,reactiveprogramdata$imported_data, reactiveprogramdata$final_output,reactiveprogramdata$useful_data)
    reactiveprogramdata$final_output=dummy$final_output
    reactiveprogramdata$useful_data=dummy$useful_data
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

    reactiveprogramdata$ROI_data[ROI_separator[ind, 1]:ROI_separator[ind, 2],]=reactiveprogramdata$ROI_data_check[ROI_separator[ind, 1]:ROI_separator[ind, 2],]=reactiveROItestingdata$ROIpar
    ROI_names=paste(reactiveprogramdata$ROI_data[ROI_separator[, 1],1],reactiveprogramdata$ROI_data[ROI_separator[, 1],2])
    names(reactiveprogramdata$select_options)=ROI_names
    },error=function(e) {print('Error. Please explain the issue on the Github website')})
    })



  #Spectra table.
    output$x1 = tryCatch({DT::renderDataTable(reactiveprogramdata$spectra , selection = list(mode = 'multiple', selected = 1),server = TRUE)},error=function(e){})
    # proxy = dataTableProxy('x1')



	#Plotly figure where to analyze peak shape fitting
  output$plot <- renderPlotly({
    if (reactiveprogramdata$beginning==FALSE) return()
     reactiveprogramdata$plot
  })



  #Table where to analyze quantifications
  output$qualitypar <- renderD3tf({
    tryCatch({
    tableProps =
    d3tf(reactiveROItestingdata$qualitypar,
      tableProps = list(btn_reset = TRUE),
      enableTf = FALSE,
      edit=FALSE,
      showRowNames = TRUE,
      tableStyle = "table table-bordered")
    },error=function(e) {return(NULL) })

  })

  #Repository table
  observe({
    suppressWarnings(
    if (!is.na(reactiveprogramdata$imported_data))  {
  output$repository = DT::renderDataTable(
    reactiveprogramdata$imported_data$repository[which(reactiveprogramdata$imported_data$repository[,3]>reactiveROItestingdata$ROIpar[1,2]&reactiveprogramdata$imported_data$repository[,3]<reactiveROItestingdata$ROIpar[1,1]),] , server = TRUE)
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
  output$directedition <- renderD3tf({
    if(reactiveprogramdata$beginning==F) return(NULL)
    observe({
      if(is.null(input$directedition_edit)||is.null(reactiveprogramdata$stop2)|| (reactiveprogramdata$stop2==1)) {
        reactiveprogramdata$change2=0
        return(NULL)
      }
      edit <- input$directedition_edit
      isolate({
        id <- edit$id
        row <- as.integer(edit$row)
        col <- as.integer(edit$col)
        val <- as.numeric(edit$val)
        if(col == 0) {
          oldval <- rownames(reactiveROItestingdata$signpar)[row]
        } else if (col %in% c(1:7)){
          if(is.na(suppressWarnings(as.numeric(val)))) {
            oldval <- reactiveROItestingdata$signpar[row, col]
            rejectEdit(session, tbl = "directedition_edit", row = row, col = col, id = id, value = oldval)
            # reactiveprogramdata$roi=0
            return(NULL)
          }
        }
        if (reactiveprogramdata$change2==1){
          reactiveprogramdata$change2=0
          reactiveprogramdata$stop2=1
        } else {
          reactiveROItestingdata$signpar[row, col] <- val
          confirmEdit(session, tbl = "directedition_edit", row = row, col = col, id = id, value = round(val,4))
        }
      })
    })
    d3tf(reactiveROItestingdata$signpar,
         tableProps = list(btn_reset = TRUE),
         enableTf = F,
         edit=TRUE,
         showRowNames = TRUE,
         tableStyle = "table table-bordered")
  })

  #Quantification after direct edition of paramters
  observeEvent(input$direct_edition, {
    if (all(is.na(reactiveROItestingdata$signpar))) {
      print("You can only perform direct edition of line shape fitted ROIs")
      return(NULL)
    }
    tryCatch({
    reactivequantdata$method1=signals_int(reactiveprogramdata$imported_data, reactiveprogramdata$final_output,reactiveprogramdata$ind,reactiveROItestingdata$signpar,reactiveROItestingdata$ROIpar)
    reactiveprogramdata$plot=reactivequantdata$method1$p
    #reactivequantdata$stop3=1
    reactiveROItestingdata$qualitypar=cbind(reactivequantdata$method1$results_to_save$quantification,reactivequantdata$method1$results_to_save$fitting_error,reactivequantdata$method1$results_to_save$signal_area_ratio)
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','Fitting Error','Signal/total area ratio')
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
  observe({
    if (reactiveprogramdata$beginning==FALSE) return()
    tryCatch({
    validation_data=validation(reactiveprogramdata$final_output,input$select_validation,reactiveprogramdata$ROI_data,reactiveprogramdata$imported_data$Metadata)
    output$fit_selection = DT::renderDataTable({ datatable(round(validation_data$alarmmatrix,4),selection = list(mode = 'single', target = 'cell')) %>% formatStyle(colnames(validation_data$alarmmatrix), backgroundColor = styleInterval(validation_data$brks, validation_data$clrs))
    })},error=function(e) {
      print("Not enough data to model it.")
    })
  })

  #Loading of quantification
  observeEvent(input$fit_selection_cell_clicked, {
if (length(input$fit_selection_cell_clicked)<1) return()
    #Checks and setting of parameters
    tryCatch({reactiveprogramdata$info=input$fit_selection_cell_clicked
    reactiveprogramdata$ind=reactiveprogramdata$info$row
    reactiveprogramdata$change=reactiveprogramdata$change2=1
    reactiveprogramdata$stop=reactiveprogramdata$stop2=0
    # proxy %>% selectRows(as.numeric(input$fit_selection_cell_clicked$col+2))
    #if (length(reactiveprogramdata$info$row)>0) reactivequantdata$stop3=1
    resetInput(session, "ROIdata_edit")
    resetInput(session, "directedition_edit")
    updateSelectInput(session, "select",selected = NULL)

    # if (length(reactiveprogramdata$info$row)!=1) return(NULL)

    dummy=load_quantification(reactiveprogramdata$useful_data,reactiveprogramdata$imported_data,reactiveprogramdata$final_output,reactiveprogramdata$info,reactiveprogramdata$ROI_data)

      reactiveprogramdata$plot=dummy$plot
      reactiveROItestingdata$signpar=dummy$signpar
      reactiveROItestingdata$ROIpar=dummy$ROIpar
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
            reactiveprogramdata$medianplot
            }
        })
  },error=function(e) {print('Error. Please explain the issue on the Github website')})


    })

  #Add and remove signals and save changes
  observeEvent(input$add_hmdb_signal, {
    tryCatch({
    dummy=reactiveprogramdata$imported_data$repository[input$repository2_rows_selected,]
    dummy=c(dummy[,3]+0.02,dummy[,3]-0.02,'Baseline Fitting',dummy[,1],1,dummy[,3],0.005,median(reactiveprogramdata$ROI_data_check[,8]),dummy[,4],dummy[,5],0,dummy[,6])
    if (dummy[9]=='d') {
      dummy[9]=2
    } else if (dummy[9]=='t') {
      dummy[9]=3
    } else {
      dummy[9]=1
    }

    if (is.na(as.numeric(dummy[10])))  dummy[10]=0
    dummy=as.list(dummy)
    dummy[-c(3,4)]=as.numeric(dummy[-c(3,4)])

    reactiveprogramdata$ROI_data_check=rbind(dummy,reactiveprogramdata$ROI_data_check)
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
    reactiveprogramdata$ROI_data_check=rbind(rep(NA,ncol(reactiveprogramdata$ROI_data_check)),reactiveprogramdata$ROI_data_check)
    }, error = function(e) {
      print('Error. Please explain the issue in the Github page.')
    })
     })
  observeEvent(input$remove_signal, {
    tryCatch({
    reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data_check[-input$roi_profiles_select,]
    resetInput(session, "roi_profiles_edit")
    }, error = function(e) {
      print('Error. Please explain the issue in the Github page.')
    })
     })
  observeEvent(input$save_changes, {
    tryCatch({
      reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data_check[!duplicated(reactiveprogramdata$ROI_data_check[,4:5]),]

      reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data_check[sort(reactiveprogramdata$ROI_data_check[,1],index.return=TRUE)$ix,]
    new_fitting_error=new_intensity=new_signal_area_ratio=new_shift=new_width=new_Area=matrix(NA,nrow(reactiveprogramdata$final_output$signal_area_ratio),nrow(reactiveprogramdata$ROI_data_check),dimnames=list(reactiveprogramdata$imported_data$Experiments,paste(reactiveprogramdata$ROI_data_check[,4],reactiveprogramdata$ROI_data_check[,5],sep='_')))
    new_signals_codes=new_signals_names=rep(NA,nrow(reactiveprogramdata$ROI_data_check))
    new_useful_data=reactiveprogramdata$useful_data
    for (i in 1:length(new_useful_data)) new_useful_data[[i]]=vector("list", nrow(reactiveprogramdata$ROI_data_check))
    for (i in 1:nrow(reactiveprogramdata$ROI_data_check)) {
      ind=which(reactiveprogramdata$ROI_data[,4]==reactiveprogramdata$ROI_data_check[i,4]&reactiveprogramdata$ROI_data[,5]==reactiveprogramdata$ROI_data_check[i,5])
      if (length(ind)>0) {
        new_fitting_error[,i]=reactiveprogramdata$final_output$fitting_error[,ind]
        new_intensity[,i]=reactiveprogramdata$final_output$intensity[,ind]
        new_signal_area_ratio[,i]=reactiveprogramdata$final_output$signal_area_ratio[,ind]
        new_shift[,i]=reactiveprogramdata$final_output$shift[,ind]
        new_width[,i]=reactiveprogramdata$final_output$half_band_width[,ind]
        new_Area[,i]=reactiveprogramdata$final_output$quantification[,ind]
        for (j in 1:length(new_useful_data)) new_useful_data[[j]][[i]]=reactiveprogramdata$useful_data[[j]][[ind]]
      }
    }
    reactiveprogramdata$final_output$fitting_error=new_fitting_error
    reactiveprogramdata$final_output$intensity=new_intensity
    reactiveprogramdata$final_output$signal_area_ratio=new_signal_area_ratio
    reactiveprogramdata$final_output$shift=new_shift
    reactiveprogramdata$final_output$half_band_width=new_width
    reactiveprogramdata$final_output$quantification=new_Area
    reactiveprogramdata$useful_data=new_useful_data
    reactiveprogramdata$ROI_data=reactiveprogramdata$ROI_data_check
    reactiveprogramdata$imported_data$signals_codes=seq(nrow(reactiveprogramdata$ROI_data))
    reactiveprogramdata$imported_data$signals_names=paste(reactiveprogramdata$ROI_data[,4],reactiveprogramdata$ROI_data[,5],sep='_')

    dummy = which(is.na(reactiveprogramdata$ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(reactiveprogramdata$ROI_data)[1]+1
    lal=which(duplicated(reactiveprogramdata$ROI_data[-dummy,1:2])==FALSE)
    ROI_separator = cbind(lal, c(lal[-1] - 1, dim(reactiveprogramdata$ROI_data[-dummy,])[1]))
    ROI_names=paste(reactiveprogramdata$ROI_data[ROI_separator[, 1],1],reactiveprogramdata$ROI_data[ROI_separator[, 1],2])
    reactiveprogramdata$select_options=1:length(ROI_names)
    names(reactiveprogramdata$select_options)=ROI_names
    updateSelectInput(session, "select",
      choices = reactiveprogramdata$select_options,selected = input$x1_rows_selected
    )
    }, error = function(e) {
      print('Problem during update of data. Please explain the issue in the Github page.')
})
  })

  #ROI Profiles table
  output$roi_profiles <- renderD3tf({
    observe({
      edit <- input$roi_profiles_edit
      if (!is.null(edit)) {
        isolate({
          id <- edit$id
          row <- as.integer(edit$row)
          col <- as.integer(edit$col)
          val <- edit$val

          if(col == 0) {
            oldval <- rownames(reactiveprogramdata$ROI_data_check)[row]
          } else if (col %in% c(1:2,5:12)){
            # numeric columns
            if(is.na(suppressWarnings(as.numeric(val)))) {
              oldval <- reactiveprogramdata$ROI_data_check[row, col]
              rejectEdit(session, tbl = "roi_profiles", row = row, col = col, id = id, value = oldval)
              #    reactiveprogramdata$roi=0
              return(NULL)
            }
          } else if (col %in% c(3,4)) {
            if(is.na(suppressWarnings(val))) {
              oldval <- reactiveprogramdata$ROI_data_check[row, col]
              rejectEdit(session, tbl = "roi_profiles", row = row, col = col, id = id, value = oldval)
              return(NULL)
            }
          }

          if (reactiveprogramdata$change==1){
            reactiveprogramdata$change=0
            reactiveprogramdata$stop=1
            disableEdit(session, "roi_profiles", c(1:12))
          } else{
            if(col == 0) {
            } else if (col %in% c(1:2,5:12)) {
              reactiveprogramdata$ROI_data_check[row, col] <- as.numeric(val)
              #    reactiveprogramdata$roi=1
            } else if (col %in% c(3,4)) {
              reactiveprogramdata$ROI_data_check[row, col] <- val
            }
            confirmEdit(session, tbl = "roi_profiles", row = row, col = col, id = id, value = val)
          }
        })
      }
    })

    d3tf(reactiveprogramdata$ROI_data_check,
      tableProps = list(btn_reset = FALSE),
      enableTf = FALSE,
      edit=TRUE,
      selectableRows = "single",
      tableStyle = "table table-bordered")

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
        STOCSY(reactiveprogramdata$imported_data$dataset,reactiveprogramdata$imported_data$ppm,c(input$left_ppm,input$right_ppm),input$correlation_method)
          }
        })
  },error=function(e) {return(NULL) })
    })



  #Dengrogran heatmaps for quantification and chemical shift
  tryCatch(output$dendheatmapareadata <- renderPlotly({
    if(all(is.na(reactiveprogramdata$final_output$quantification))) return()

    tryCatch({type_analysis_plot(reactiveprogramdata$final_output$quantification,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,type='dendrogram_heatmap')
  }, error = function(e) {
    print('Error. Please explain the issue in the Github page.')
    return(NULL)
  })
  }))

  tryCatch(output$dendheatmapshiftdata <- renderPlotly({
    if(all(is.na(reactiveprogramdata$final_output$quantification))) return()

    tryCatch({type_analysis_plot(reactiveprogramdata$final_output$shift,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,type='dendrogram_heatmap')
  }, error = function(e) {
    print('Error. Please explain the issue in the Github page.')
    return(NULL)
  })
}))



}
