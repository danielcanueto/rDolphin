

server = function(input, output,session) {
  
  #Increase of maximum memory size that can be uploaded
  options(shiny.maxRequestSize=1000*1024^2)

  #Setting of reactive parameters in the Shiny GUI
  reactiveROItestingdata <- reactiveValues()
  reactivequantdata <- reactiveValues(method1=NA, method2 = NA)
  #reactivequantdata <- reactiveValues(method2=NA, method1 = NA,stop3=0)

  reactiveprogramdata <- reactiveValues(ROIdata_subset=NA,ind=NA,beginning=F,dataset=NA,final_output=list(),useful_data=list(),imported_data=NA,p_value_final=NA,ROI_data=NA,ROI_data_check=NA,info=c(),ROI_separator=NA,select_options=NA,new_roi_profile=NA,p=NA,bgColScales=NA,autorun_plot=NA,ROI_names=NA,clusterplot=NA,medianplot=NA,jres_plot=NA)
  
  ## FIRST TAB REACTIVE OUTPUTS

  #Read of input provided by user
  tryCatch(observeEvent(input$file1, {
    reactiveprogramdata$inFile <- input$file1
    if (is.null(reactiveprogramdata$inFile)) {
      return(NULL)
    }

	#Imported data is loaded to 'dummy'. Only after the check that the parameters are correct, they are stored in 'reactiveprogramdata'.
	dummy=list(imported_data=NA,autorun_plot=NA,select_options=NA,spectra=NA,clusterplot=NA,medianplot=NA,beginning=F,jres_plot=NA)
	#Import of data
	dummy = tryCatch({helperimport(reactiveprogramdata$inFile$datapath,dummy)
	  }, error = function(e) {
	print('Error. Please explain the issue in the Github page if necessary.')
	return(dummy)
	})

	if (dummy$beginning==T) {
	  plo=names(sapply(dummy, names))
	  for (i in 1:length(plo)) reactiveprogramdata[[plo[i]]]=dummy[[plo[i]]]

	 #Variables that can change during the use of the GUI are separated from 'imported_data'.  
    reactiveprogramdata$final_output=dummy$imported_data$final_output
    reactiveprogramdata$useful_data=dummy$imported_data$useful_data
    reactiveprogramdata$ROI_separator=dummy$imported_data$ROI_separator
    reactiveprogramdata$ROI_data=reactiveprogramdata$ROI_data_check=dummy$imported_data$ROI_data
    reactiveprogramdata$imported_data$final_output=reactiveprogramdata$imported_data$useful_data=reactiveprogramdata$imported_data$ROI_separator=reactiveprogramdata$imported_data$ROI_data=NULL
    
	#When the session is prepared, the tabs and some inputs become active
	output$sp = DT::renderDataTable(
      dummy$total_signals_parameters , selection = list(selected = NULL),server = T)
    updateSelectInput(session, "select",choices = reactiveprogramdata$select_options,selected = 1)
    updateSelectInput(session, "select_validation",selected = 1)
    session$sendCustomMessage('activeNavs', 'Individual Quantification')
    session$sendCustomMessage('activeNavs', 'Quantification Validation')
    session$sendCustomMessage('activeNavs', 'Uni and multivariate analysis')
    session$sendCustomMessage('activeNavs', 'ROI Profiles')
    session$sendCustomMessage('activeNavs', 'STOCSY and dendrogram heatmaps')
	 }}))

	 
	 
  #Loading of previous session
  
    #Read of input provided by user
  observeEvent(input$file2, {
    reactiveprogramdata$inFile2 <- input$file2
    if (is.null(reactiveprogramdata$inFile2))
      return(NULL)

    #Session is loaded in 'savedreactiveddata' variable and passed to 'reactiveprogramdata'.
    tryCatch({load(reactiveprogramdata$inFile2$datapath)
      print("Uploading saved session.")
      plo=names(sapply(savedreactivedata, names))
      for (i in 1:length(plo)) reactiveprogramdata[[plo[i]]]=savedreactivedata[[plo[i]]]
      rm(savedreactivedata)
    }, error = function(e) {
      print('Not possible to load the session. Please revise your choice.')
      return(NULL)
    })

    
    #Names of ROIS are prepared
	dummy=NULL
	dummy=tryCatch({roifunc(reactiveprogramdata$ROI_data,reactiveprogramdata$ROI_separator,reactiveprogramdata$imported_data$Metadata,reactiveprogramdata$imported_data$Experiments)
  }, error = function(e) {
	print('Generation of Regions of Interest not possible. Please explain the issue in the Github page.')
	return(NULL)
	})
	if (!is.null(dummy)) {
	reactiveprogramdata$select_options=dummy$select_options
	reactiveprogramdata$spectra=dummy$spectra
	reactiveprogramdata$beginning =T
	print("Done!")

    #When the session is loaded the tabs and some inputs become active
    updateSelectInput(session, "select",choices = reactiveprogramdata$select_options,selected = 1)
    updateSelectInput(session, "select_validation",selected = 1)
    session$sendCustomMessage('activeNavs', 'Individual Quantification')
    session$sendCustomMessage('activeNavs', 'Quantification Validation')
    session$sendCustomMessage('activeNavs', 'Uni and multivariate analysis')
    session$sendCustomMessage('activeNavs', 'ROI Profiles')
    session$sendCustomMessage('activeNavs', 'STOCSY and dendrogram heatmaps')
	}
  })


  

  #Choice and storage of data associated to session
  
  observe({
    volumes <- c("UserFolder"="C:/")
    shinyFileSave(input, "save", roots=volumes, session=session)
    fileinfo <- parseSavePath(volumes, input$save)
    savedreactivedata=isolate(reactiveValuesToList(reactiveprogramdata))
    if (nrow(fileinfo) > 0) {
      save(savedreactivedata, file=as.character(fileinfo$datapath))
      export_path=paste(substr(as.character(fileinfo$datapath),1,(nchar(as.character(fileinfo$datapath))-6)),'_associated_data',sep='')
      tryCatch(write_info(export_path, reactiveprogramdata$final_output, reactiveprogramdata$ROI_data),error= function(e) print('Not possible to overwrite open files'))
    }
  })


  
  #Choice of folder to save plots
  
  folderInput1 <- reactive({
    volumes = c("UserFolder"="C:/")
    shinyDirChoose(input, 'folder', roots = volumes, session = session,
      restrictions = system.file(package = 'base'))
    return(parseDirPath(volumes, input$folder))
  })
  observe({
    if (length(folderInput1()) > 0) tryCatch(write_plots(folderInput1(),reactiveprogramdata$final_output,reactiveprogramdata$imported_data,reactiveprogramdata$useful_data),error= function(e) print('Not possible to overwrite open files'))
  })

  #Load of quantifications of previous session to combine with current session, to avoid repeating already performed quantifications. UNSTABLE!!!!
  # observeEvent(input$file3, {
  #   reactiveprogramdata$inFile2 <- input$file2
  #   if (is.null(reactiveprogramdata$inFile2))
  #     return(NULL)
  #   load(reactiveprogramdata$inFile2$datapath)
  #   plo=names(sapply(savedreactivedata, names))
  #   for (i in 1:length(plo)) {
  #     added_data[[plo[i]]]=savedreactivedata[plo[i]]
  #   }
  #   ind=which(reactiveprogramdata$imported_data$Experiments %in% added_data$imported_data$Experiments==T)
  #   ind2=which(reactiveprogramdata$imported_data$signals_names %in% added_data$imported_data$signals_names==T)
  #   ind3=which(added_data$imported_data$Experiments %in% reactiveprogramdata$imported_data$Experiments==T)
  #   ind4=which(added_data$imported_data$signals_names %in% reactiveprogramdata$imported_data$signals_names==T)
  #
  #   for (i in 1:length(reactiveprogramdata$final_output)) {
  #     reactiveprogramdata$final_output[[i]][ind,ind2]=added_data$final_output[[i]][ind3,ind4]
  #   }
  #   for (i in 1:length(reactiveprogramdata$useful_data)) {
  #     for (j in 1:length(reactiveprogramdata$final_output[[j]])) {
  #
  #       reactiveprogramdata$useful_data[[ind[i]]][[ind2[j]]]=added_data$useful_data[[ind3[i]]][[ind4[j]]]
  #     }}
  # })

  
  
  #Appearance of autorun and aligment buttons only after beginning or loading session
  output$varselect <- renderUI({
    if(reactiveprogramdata$beginning==F){return()}
    actionButton('autorun', 'Autorun all spectra')
  })
  output$align_button <- renderUI({
    if(reactiveprogramdata$beginning==F){return()}
    actionButton('alignment', 'Alignment of signals')
  })
 

#Removed by now because it is unstable.
 # output$peak_analysis <- renderUI({
  #   if(reactiveprogramdata$beginning==F){return()}
  #   actionButton('peak_analysis', 'Peak analysis')
  # })

   #Automatic quantification of all ROIs in all spectra
  tryCatch({observeEvent(input$autorun, {
    quantification_variables = autorun(reactiveprogramdata$imported_data, reactiveprogramdata$final_output,reactiveprogramdata$useful_data,reactiveprogramdata$ROI_data,reactiveprogramdata$ROI_separator)
    reactiveprogramdata$final_output=quantification_variables$final_output
    reactiveprogramdata$useful_data=quantification_variables$useful_data
    })},
    error = function(e) {
      print('Error. Please explain the issue in the Github page if necessary.')
      quantification_variables=NA
      return(quantification_variables)
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
    if (reactiveprogramdata$beginning==F) return()
    reactiveprogramdata$autorun_plot
  })

  
  
  ## SECOND TAB REACTIVE OUTPUTS

  #Selection of ROI
  
  
  # ROI parameters are loaded when a ROI is selected
  tryCatch(observeEvent(input$select, {
    if (reactiveprogramdata$beginning==F) return()
    if (reactiveprogramdata$beginning ==T) {
      reactiveprogramdata$ROIdata_subset=reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[as.numeric(input$select), 1]:reactiveprogramdata$ROI_separator[as.numeric(input$select), 2],]
    }

	#Setting of reactivity to edition of parameters prepared. Probably improvable, but it is a delicate matter. Now that I can debug it easily, it can be 'cleaned' in the future.
    reactiveprogramdata$change=reactiveprogramdata$change2=1
    reactiveprogramdata$stop=reactiveprogramdata$stop2=0
    
    
    #reactivequantdata$stop3=0
    # reactiveprogramdata$roi=NULL
    reactiveprogramdata$info=c()
   
	  # Tables of ROI paramters ('ROIpar'), calculated deconvolution parameters ('signpar') and indicators of quantification ('qualitypar') are reset
	  resetInput(session, "ROIdata_edit")
    resetInput(session, "directedition_edit")
    reactiveROItestingdata$ROIpar <- reactiveprogramdata$ROIdata_subset
    reactiveROItestingdata$signpar <- rbind(rep(NA,7),rep(NA,7))
    colnames(reactiveROItestingdata$signpar)=c("intensity",	"chemical shift",	"half bandwidth",	"gaussian %",	"J coupling",	"multiplicities",	"roof effect")
    reactiveROItestingdata$qualitypar <- rbind(rep(NA,3),rep(NA,3))
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','fitting_error','signal/total area ratio')
    
	# Plot is prepared
	ROI_limits=c(reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,1]))],reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,2]))])
    ind=ifelse(is.null(input$x1_rows_selected),1,input$x1_rows_selected)
    dummy=type_plot(reactiveprogramdata$imported_data,ROI_limits,ind,reactiveprogramdata$medianplot,reactiveprogramdata$clusterplot)
    if (!is.null(dummy)) reactiveprogramdata$plot=dummy

	#Analysis of ROI edition. If edition is not correct (for example, there are characters in a numeric input), the edition is rejected and shown with red colour. If correct, the change is accepted with green colour.
	#TODO: it seems sometiems the edition fails if the change was too quick. Revise possible ways to control it.
	output$ROIdata <- renderD3tf({
      tableProps <- list(
        btn_reset = F,
        sort = TRUE,
        sort_config = list(
          sort_types = c("String", rep("Number", ncol(reactiveprogramdata$ROIdata_subset)))
        )
      )
      observe({
        if(is.null(input$ROIdata_edit)|(reactiveprogramdata$stop==1)) {
          reactiveprogramdata$change=0
          return(NULL)
        }
        edit <- input$ROIdata_edit
        isolate({
          id <- edit$id
          row <- as.integer(edit$row)
          col <- as.integer(edit$col)
          val <- edit$val

          if(col == 0) {
            oldval <- rownames(reactiveROItestingdata$ROIpar)[row]
          } else if (col %in% c(1:2,5:11)){
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
          if (reactiveprogramdata$change==1){
            reactiveprogramdata$change=0
            reactiveprogramdata$stop=1
            disableEdit(session, "ROIdata", c(1:11))
          } else{
            if(col == 0) {
            } else if (col %in% c(1:2,5:11)) {
              reactiveROItestingdata$ROIpar[row, col] <- as.numeric(val)
                 # reactiveprogramdata$roi=1
            } else if (col %in% c(3)) {
              reactiveROItestingdata$ROIpar[row, col] <- val
            }
            confirmEdit(session, tbl = "ROIdata", row = row, col = col, id = id, value = val)
          }
        })
      })

      d3tf(reactiveROItestingdata$ROIpar,
        tableProps = tableProps,
        enableTf = F,
        edit=TRUE,
        tableStyle = "table table-bordered")
    })}))

  
  
  #Selection of spectra, or of cluster or median plots
  tryCatch(observeEvent(input$x1_rows_selected, {
    if (reactiveprogramdata$beginning==F) return()
    if (reactiveprogramdata$beginning ==T) {
      reactiveprogramdata$ROIdata_subset=reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[as.numeric(input$select), 1]:reactiveprogramdata$ROI_separator[as.numeric(input$select), 2],]
    }
	#Reset of parameters
    reactiveprogramdata$change=reactiveprogramdata$change2=1
    reactiveprogramdata$stop=reactiveprogramdata$stop2=0
    
    
    #reactivequantdata$stop3=0
        #reactiveprogramdata$roi=NULL
    reactiveprogramdata$info=c()
	resetInput(session, "directedition_edit")
    reactiveROItestingdata$ROIpar <- reactiveprogramdata$ROIdata_subset
    reactiveROItestingdata$signpar <- rbind(rep(NA,7),rep(NA,7))
    colnames(reactiveROItestingdata$signpar)=c("intensity",	"shift",	"half_band_width",	"gaussian",	"J_coupling",	"multiplicities",	"roof_effect")
    reactiveROItestingdata$qualitypar <- rbind(rep(NA,3),rep(NA,3))
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','fitting_error','signal/total area ratio')
    ROI_limits=c(reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,1]))],reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,2]))])

	#Generation of plot
    dummy=type_plot(reactiveprogramdata$imported_data,ROI_limits,input$x1_rows_selected,reactiveprogramdata$medianplot,reactiveprogramdata$clusterplot)
    if (!is.null(dummy)) reactiveprogramdata$plot=dummy
  }))

  #Individual quantification
  tryCatch(observeEvent(input$action, {
	#Only when there is a spectrum selected and with number >2 (to avoid cluster and median plot) the quantifications is performed, with prior adaptation to row of dataset
    if(length(reactiveprogramdata$info)==0) reactiveprogramdata$ind=input$x1_rows_selected-2
    if (length(reactiveprogramdata$ind)!=1|reactiveprogramdata$ind<1) {
      print('Select one valid spectrum')
      return(NULL)
    }

	#The automatic quantification
    reactivequantdata$method1 <- tryCatch({not_automatic_quant(reactiveprogramdata$imported_data, reactiveprogramdata$final_output, reactiveprogramdata$ind,reactiveROItestingdata$ROIpar,reactiveprogramdata$useful_data,interface=T)},error=function(e) {
      print("There was a problem.")
      return(NULL)
      })
    
	#Update of tables of tab
    if (!is.null(reactivequantdata$method1$signals_parameters)) {
      reactiveprogramdata$plot=reactivequantdata$method1$p
      #reactivequantdata$stop3=1
      reactiveROItestingdata$qualitypar=cbind(reactivequantdata$method1$results_to_save$Area,reactivequantdata$method1$results_to_save$fitting_error,reactivequantdata$method1$results_to_save$signal_area_ratio)
      colnames(reactiveROItestingdata$qualitypar)=c('Quantification','Fitting Error','Signal/total area ratio')
      rownames(reactiveROItestingdata$qualitypar)=rownames(reactivequantdata$method1$plot_data)[4:(3+nrow(reactiveROItestingdata$qualitypar))]
      ind=which(rownames(reactiveROItestingdata$qualitypar)=='additional signal')
      reactiveprogramdata$bgColScales = rep(c("","info"),times=c(length(rownames(reactiveROItestingdata$qualitypar))-length(ind),length(ind)))
      reactiveROItestingdata$signpar <- t(reactivequantdata$method1$signals_parameters)
      reactiveprogramdata$stop=0
         # reactiveprogramdata$roi=1
    }
  }))

  
  
  #Quantification of all spectra in the ROI:
  tryCatch(observeEvent(input$autorun_signal, {
    is_autorun='Y'
    dummy <- not_automatic_quant(reactiveprogramdata$imported_data, reactiveprogramdata$final_output, seq(nrow(reactiveprogramdata$imported_data$dataset)),reactiveROItestingdata$ROIpar,reactiveprogramdata$useful_data,interface=T)
    reactiveprogramdata$final_output=dummy$final_output
    reactiveprogramdata$useful_data=dummy$useful_data
  }))


  
  #Remove quantification or save quantification or ROI profile edited
  tryCatch(observeEvent(input$remove_q, {
    if (!is.null(reactiveprogramdata$imported_data$signals_names[reactiveprogramdata$info$col])) {
      ind=which(reactiveprogramdata$ROI_data[,4]==reactiveprogramdata$imported_data$signals_names[reactiveprogramdata$info$col])
    } else {
      ind=as.numeric(input$select)
    }
    reactiveprogramdata$final_output <- remove_quant(reactiveprogramdata$info,reactiveprogramdata$imported_data, reactiveprogramdata$final_output)
  }))
  
  tryCatch(observeEvent(input$save_results, {
    if (is.null(reactivequantdata$method1$Ydata)) {return(NULL)
	}
    dummy=save_roi_testing(reactivequantdata$method1,reactiveprogramdata$imported_data, reactiveprogramdata$final_output,reactiveprogramdata$useful_data)
    reactiveprogramdata$final_output=dummy$final_output
    reactiveprogramdata$useful_data=dummy$useful_data
  }))
  
  
  
    #Save edition of ROI profile
  tryCatch(observeEvent(input$save_profile, {
    if (length(reactiveprogramdata$info$col)>0) {
      ind=which(reactiveprogramdata$ROI_separator[,2]-reactiveprogramdata$info$col>=0)[1]
    } else {
      ind=as.numeric(input$select)
    }
    reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[ind, 1]:reactiveprogramdata$ROI_separator[ind, 2],]=reactiveprogramdata$ROI_data_check[reactiveprogramdata$ROI_separator[ind, 1]:reactiveprogramdata$ROI_separator[ind, 2],]=reactiveROItestingdata$ROIpar
    ROI_names=paste(reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],1],reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],2])
    names(reactiveprogramdata$select_options)=ROI_names
  }))


  
  #Spectra table. 
  output$x1 = tryCatch(DT::renderDataTable(reactiveprogramdata$spectra , selection = list(mode = 'multiple', selected = 1),server = T))


  
	#Plotly figure where to analyze peak shape fitting
  tryCatch(output$plot <- renderPlotly({
    if (reactiveprogramdata$beginning==F | is.null(input$x1_rows_selected)) return()
      print(reactiveprogramdata$plot)
  }))

  
  
  #Table where to analyze quantifications
  tryCatch(output$qualitypar <- renderD3tf({
    tableProps =
    d3tf(reactiveROItestingdata$qualitypar,
      tableProps = list(btn_reset = TRUE),
      enableTf = F,
      edit=F,
      showRowNames = TRUE,
      tableStyle = "table table-bordered")
  }))

  #Repository table
  tryCatch(observe({
    if (!is.na(reactiveprogramdata$imported_data))  {
  output$repository = DT::renderDataTable(
    reactiveprogramdata$imported_data$repository[which(reactiveprogramdata$imported_data$repository[,3]>reactiveROItestingdata$ROIpar[1,2]&reactiveprogramdata$imported_data$repository[,3]<reactiveROItestingdata$ROIpar[1,1]),] , server = T)
    }
  }))

  #Direct edition of parameters before quantification
  tryCatch(output$directedition <- renderD3tf({
    observe({
      if(is.null(input$directedition_edit)|(reactiveprogramdata$stop2==1)) {
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
      rowStyles = reactiveprogramdata$bgColScales,
      tableStyle = "table table-bordered")
  }))

  #Quantification after direct edition of paramters
  tryCatch(observeEvent(input$direct_edition, {
    reactivequantdata$method1=reactivequantdata$method2=signals_int(reactiveprogramdata$imported_data, reactiveprogramdata$final_output,reactiveprogramdata$ind,reactiveROItestingdata$signpar,reactiveROItestingdata$ROIpar)
    reactiveROItestingdata$qualitypar=cbind(reactivequantdata$method2$results_to_save$Area,reactivequantdata$method2$results_to_save$fitting_error,reactivequantdata$method2$results_to_save$signal_area_ratio)
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','fitting_error','signal/total area ratio')
    rownames(reactiveROItestingdata$qualitypar)=rownames(reactivequantdata$method2$plot_data)[-c(1, 2, 3)]
  }))

  
  
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
  tryCatch(observe({
    if (reactiveprogramdata$beginning==F) return()
    validation_data=validation(reactiveprogramdata$final_output,input$select_validation,reactiveprogramdata$ROI_data,reactiveprogramdata$imported_data$Metadata)
    output$fit_selection = DT::renderDataTable({ datatable(round(validation_data$alarmmatrix,4),selection = list(mode = 'single', target = 'cell')) %>% formatStyle(colnames(validation_data$alarmmatrix), backgroundColor = styleInterval(validation_data$brks, validation_data$clrs))
    })
  }))

  #Loading of quantification
  tryCatch(observeEvent(input$fit_selection_cell_clicked, {
if (length(input$fit_selection_cell_clicked)<1) return()
    #Checks and setting of parameters
    reactiveprogramdata$info=input$fit_selection_cell_clicked
    reactiveprogramdata$ind=reactiveprogramdata$info$row
    reactivequantdata$method2=NULL
    reactiveprogramdata$change=reactiveprogramdata$change2=1
    reactiveprogramdata$stop=reactiveprogramdata$stop2=0

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
    updateTabsetPanel(session, "mynavlist",selected = "Individual Quantification")
  }))


  ## FOURTH TAB REACTIVE OUTPUTS

  #Add and remove signals and save changes
  tryCatch(observeEvent(input$add_signal, {
    reactiveprogramdata$ROI_data_check=rbind(reactiveprogramdata$ROI_data_check,rep(NA,ncol(reactiveprogramdata$ROI_data_check)))
  }))
  tryCatch(observeEvent(input$remove_signal, {
    reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data_check[-input$roi_profiles_select,]
    resetInput(session, "roi_profiles_edit")
  }))
  tryCatch(observeEvent(input$save_changes, {
    reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data_check[sort(reactiveprogramdata$ROI_data_check[,1],index.return=T)$ix,]
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
        new_Area[,i]=reactiveprogramdata$final_output$Area[,ind]
        for (j in 1:length(new_useful_data)) new_useful_data[[j]][[i]]=reactiveprogramdata$useful_data[[j]][[ind]]
      }
    }
    reactiveprogramdata$final_output$fitting_error=new_fitting_error
    reactiveprogramdata$final_output$intensity=new_intensity
    reactiveprogramdata$final_output$signal_area_ratio=new_signal_area_ratio
    reactiveprogramdata$final_output$shift=new_shift
    reactiveprogramdata$final_output$half_band_width=new_width
    reactiveprogramdata$final_output$Area=new_Area
    reactiveprogramdata$useful_data=new_useful_data
    reactiveprogramdata$ROI_data=reactiveprogramdata$ROI_data_check
    reactiveprogramdata$imported_data$signals_codes=seq(nrow(reactiveprogramdata$ROI_data))
    reactiveprogramdata$imported_data$signals_names=paste(reactiveprogramdata$ROI_data[,4],reactiveprogramdata$ROI_data[,5],sep='_')

    dummy = which(is.na(reactiveprogramdata$ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(reactiveprogramdata$ROI_data)[1]+1
    lal=which(duplicated(reactiveprogramdata$ROI_data[-dummy,1:2])==F)
    reactiveprogramdata$ROI_separator = cbind(lal, c(lal[-1] - 1, dim(reactiveprogramdata$ROI_data[-dummy,])[1]))
    ROI_names=paste(reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],1],reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],2])
    reactiveprogramdata$select_options=1:length(ROI_names)
    names(reactiveprogramdata$select_options)=ROI_names
    updateSelectInput(session, "select",
      choices = reactiveprogramdata$select_options,selected = 1
    )

  }))

  #ROI Profiles table
  tryCatch(output$roi_profiles <- renderD3tf({
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
          } else if (col %in% c(1:2,5:11)){
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
            disableEdit(session, "roi_profiles", c(1:11))
          } else{
            if(col == 0) {
            } else if (col %in% c(1:2,5:11)) {
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
      tableProps = list(btn_reset = F),
      enableTf = F,
      edit=TRUE,
      selectableRows = "multi",
      tableStyle = "table table-bordered")

  }))

  ## FIFTH TAB REACTIVE OUTPUTS

  #Boxplot plot
  tryCatch(output$plot_p_value_2 <- renderPlotly({
    if(all(is.na(reactiveprogramdata$final_output$Area))) return()
    tryCatch({type_analysis_plot(reactiveprogramdata$final_output$Area,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,reactiveprogramdata$ROI_data,type='boxplot')
    }, error = function(e) {
      print('Generation of Regions of Interest not possible. Please explain the issue in the Github page.')
      return(NULL)
    })
  }))
  #PCA plot
  tryCatch(output$pcascores <- renderPlotly({
    if(all(is.na(reactiveprogramdata$final_output$Area))) return()

    tryCatch({type_analysis_plot(reactiveprogramdata$final_output$Area,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,reactiveprogramdata$ROI_data,type='pca')
  }, error = function(e) {
    print('Generation of Regions of Interest not possible. Please explain the issue in the Github page.')
    return(NULL)
  })
  }))

  ## SIXTH TAB REACTIVE OUTPUTS

# STOCSY generation
  tryCatch({observeEvent(input$stocsy, {
    left_ppm <- renderText({ input$left_ppm })
    right_ppm <- renderText({ input$right_ppm })
      output$stocsy_plot=renderPlotly({
        if (input$stocsy==1) {
          reactiveprogramdata$clusterplot
        } else {
        STOCSY(reactiveprogramdata$imported_data$dataset,reactiveprogramdata$imported_data$ppm,c(input$left_ppm,input$right_ppm),input$correlation_method)
          }
        })
    })
  },
    error = function(e) {
      print('Error. Please explain the issue in the Github page if necessary.')
    })


  #Dengrogran heatmaps for quantification and chemical shift
  tryCatch(output$dendheatmapareadata <- renderPlotly({
    if(all(is.na(reactiveprogramdata$final_output$Area))) return()

    tryCatch({type_analysis_plot(reactiveprogramdata$final_output$Area,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,reactiveprogramdata$ROI_data,type='dendrogram_heatmap')
  }, error = function(e) {
    print('Generation of Regions of Interest not possible. Please explain the issue in the Github page.')
    return(NULL)
  })
  }))

  tryCatch(output$dendheatmapshiftdata <- renderPlotly({
    if(all(is.na(reactiveprogramdata$final_output$Area))) return()

    tryCatch({type_analysis_plot(reactiveprogramdata$final_output$shift,reactiveprogramdata$final_output,reactiveprogramdata$imported_data,reactiveprogramdata$ROI_data,type='dendrogram_heatmap')
  }, error = function(e) {
    print('Generation of Regions of Interest not possible. Please explain the issue in the Github page.')
    return(NULL)
  })
}))



}
