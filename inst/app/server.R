

server = function(input, output,session) {
  #Increase of maximum memory size that can be uploaded and setting of margins for plots
  options(shiny.maxRequestSize=1000*1024^2)

  #Setting of reactive parameters in the program
  reactiveROItestingdata <- reactiveValues()
  reactivequantdata <- reactiveValues(method2=NULL, method1 = NULL,stop3=0)
  reactiveprogramdata <- reactiveValues(ROIdata_subset=NULL,ind=NULL,beginning=F,dataset=NULL,finaloutput=NULL,useful_data=list(),imported_data=NULL,p_value_final=NULL,ROI_data=NULL,ROI_data_check=NULL,info=c(),ROI_separator=NULL,select_options=NULL,new_roi_profile=NULL,p=NULL,bgColScales=NULL,autorun_plot=NULL,ROI_names=NULL,clusterplot=NULL,medianplot=NULL)

  ## FIRST TAB REACTIVE OUTPUTS

  #Import of data and inputs provided by user
  observeEvent(input$file1, {
    reactiveprogramdata$inFile <- input$file1

    if (is.null(reactiveprogramdata$inFile))
      return(NULL)
    #Import of data
    reactiveprogramdata$imported_data = import_data(reactiveprogramdata$inFile$datapath)


    #j-res plot is prepared
    reactiveprogramdata$jres_path=reactiveprogramdata$imported_data$jres_path
    if (reactiveprogramdata$jres_path!='')
      output$jres_plot <- try(renderPlotly({
        fhs(reactiveprogramdata$jres_path)
      }))

    # reactiveprogramdata$originaldataset=imported_data$dataset
    reactiveprogramdata$beginning =T

    #plot of quantification in model spectrum with current roi profiles
    reactiveprogramdata$autorun_plot=autorun_model_spectrum(reactiveprogramdata$imported_data)
    reactiveprogramdata$finaloutput=reactiveprogramdata$imported_data$finaloutput
    reactiveprogramdata$useful_data=reactiveprogramdata$imported_data$useful_data
    reactiveprogramdata$ROI_separator=reactiveprogramdata$imported_data$ROI_separator
    reactiveprogramdata$ROI_data=reactiveprogramdata$ROI_data_check=reactiveprogramdata$imported_data$ROI_data
    reactiveprogramdata$imported_data$finaloutput=reactiveprogramdata$imported_data$useful_data=reactiveprogramdata$imported_data$ROI_separator=reactiveprogramdata$imported_data$ROI_data=NULL

    #plots of representative spectra and median spectra per group to help setting the right ROI parameters
    reactiveprogramdata$clusterplot=clustspectraplot(reactiveprogramdata$imported_data)
    reactiveprogramdata$medianplot=medianplot(reactiveprogramdata$imported_data)


    #Subsetting of ROIs is prepared


    #Names of ROIS and cluster and median spectra are prepared
    ROI_names=paste(reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],1],reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],2])
    reactiveprogramdata$select_options=seq_along(ROI_names)
    names(reactiveprogramdata$select_options)=ROI_names
    mm=matrix(NA,2,dim(reactiveprogramdata$imported_data$Metadata)[2])
    colnames(mm)=colnames(reactiveprogramdata$imported_data$Metadata)
    spectra=cbind(c('Exemplars','Median Spectrum per group',reactiveprogramdata$imported_data$Experiments),rbind(mm,reactiveprogramdata$imported_data$Metadata))
    colnames(spectra)=c('spectrum',colnames(mm))
    output$x1 = DT::renderDataTable(
      spectra , selection = list(mode = 'multiple', selected = 1),server = T)

    #When the session is prepared the tabs and some inputs become active
    updateSelectInput(session, "select",
      choices = reactiveprogramdata$select_options,selected = 1
    )
    updateSelectInput(session, "select_validation",selected = 1)
    session$sendCustomMessage('activeNavs', 'Individual Quantification')
    session$sendCustomMessage('activeNavs', 'Quantification Validation')
    session$sendCustomMessage('activeNavs', 'Uni and multivariate analysis')
    session$sendCustomMessage('activeNavs', 'ROI Profiles')
    session$sendCustomMessage('activeNavs', 'Dendrogram heatmaps')

  })

  #Loading of previous session
  observeEvent(input$file2, {
    reactiveprogramdata$inFile2 <- input$file2
    if (is.null(reactiveprogramdata$inFile2))
      return(NULL)
    load(reactiveprogramdata$inFile2$datapath)

    #Session is loaded in 'savedreactiveddata' variable and passed to the variable that collects reactive data
    plo=names(sapply(savedreactivedata, names))
    for (i in 1:length(plo)) {
      reactiveprogramdata[[plo[i]]]=savedreactivedata[plo[i]]
    }
    reactiveprogramdata$finaloutput=savedreactivedata$finaloutput
    reactiveprogramdata$jres_path=savedreactivedata$jres_path
    reactiveprogramdata$imported_data$repository=savedreactivedata$repository
    reactiveprogramdata$imported_data=savedreactivedata$imported_data
    reactiveprogramdata$useful_data=savedreactivedata$useful_data
    # reactiveprogramdata$originaldataset=savedreactivedata$originaldataset
    reactiveprogramdata$p_value_final=savedreactivedata$p_value_final
    # reactiveprogramdata$ROI_data=reactiveprogramdata$ROI_data_check=savedreactivedata$ROI_data
    reactiveprogramdata$ROI_separator=savedreactivedata$ROI_separator
    reactiveprogramdata$autorun_plot=savedreactivedata$autorun_plot
    reactiveprogramdata$clusterplot=savedreactivedata$clusterplot
    reactiveprogramdata$medianplot=savedreactivedata$medianplot

    rm(savedreactivedata)

    #Names of ROIS and cluster and median spectra are prepared
    ROI_names=paste(reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],1],reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],2])
    reactiveprogramdata$select_options=1:length(ROI_names)
    names(reactiveprogramdata$select_options)=ROI_names
    mm=matrix(NA,2,dim(reactiveprogramdata$imported_data$Metadata)[2])
    colnames(mm)=colnames(reactiveprogramdata$imported_data$Metadata)
    spectra=cbind(c('Exemplars','Median Spectrum per group',reactiveprogramdata$imported_data$Experiments),rbind(mm,reactiveprogramdata$imported_data$Metadata))
    colnames(spectra)=c('spectrum',colnames(mm))
    output$x1 = DT::renderDataTable(
      spectra , selection = list(mode = 'multiple', selected = 1),server = T)
    reactiveprogramdata$beginning =T


    #j-res spectrum is prepared
    if (reactiveprogramdata$jres_path!='')
      output$jres_plot <- try(renderPlotly({
        fhs(reactiveprogramdata$jres_path)
      }))

    #When the session is loaded the tabs and some inputs become active
    updateSelectInput(session, "select",
      choices = reactiveprogramdata$select_options,selected = 1
    )
    updateSelectInput(session, "select_validation",selected = 1)
    session$sendCustomMessage('activeNavs', 'Individual Quantification')
    session$sendCustomMessage('activeNavs', 'Quantification Validation')
    session$sendCustomMessage('activeNavs', 'Uni and multivariate analysis')
    session$sendCustomMessage('activeNavs', 'ROI Profiles')
    session$sendCustomMessage('activeNavs', 'Dendrogram heatmaps')
  })

  #Choice of folder to save plots
  folderInput1 <- reactive({
    volumes = c("UserFolder"="C:/")

    shinyDirChoose(input, 'folder', roots = volumes, session = session,
      restrictions = system.file(package = 'base'))
    return(parseDirPath(volumes, input$folder))
  })
  #If a choice of plots folder has been set
  observe({
    if (length(folderInput1()) > 0) tryCatch(write_plots(folderInput1(),reactiveprogramdata$finaloutput,reactiveprogramdata$imported_data,reactiveprogramdata$useful_data),error= function(e) print('Not possible to overwrite open files'))

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
      tryCatch(write_info(export_path, reactiveprogramdata$finaloutput, reactiveprogramdata$ROI_data),error= function(e) print('NOt possible to overwrite open files'))
    }
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
  #   for (i in 1:length(reactiveprogramdata$finaloutput)) {
  #     reactiveprogramdata$finaloutput[[i]][ind,ind2]=added_data$finaloutput[[i]][ind3,ind4]
  #   }
  #   for (i in 1:length(reactiveprogramdata$useful_data)) {
  #     for (j in 1:length(reactiveprogramdata$finaloutput[[j]])) {
  #
  #       reactiveprogramdata$useful_data[[ind[i]]][[ind2[j]]]=added_data$useful_data[[ind3[i]]][[ind4[j]]]
  #     }}
  # })

  #Appearance of buttons only after beginnig or laoding session
  output$varselect <- renderUI({
    if(reactiveprogramdata$beginning==F){return()}
    actionButton('autorun', 'Autorun all spectra')
  })
  output$align_button <- renderUI({
    if(reactiveprogramdata$beginning==F){return()}
    actionButton('alignment', 'Alignment of signals')
  })
  # output$peak_analysis <- renderUI({
  #   if(reactiveprogramdata$beginning==F){return()}
  #   actionButton('peak_analysis', 'Peak analysis')
  # })
  #Automatic quantification of all ROIs in all spectra
  observeEvent(input$autorun, {

    quantification_variables = autorun(reactiveprogramdata$imported_data, reactiveprogramdata$finaloutput,reactiveprogramdata$useful_data,reactiveprogramdata$ROI_data,reactiveprogramdata$ROI_separator)
    reactiveprogramdata$finaloutput=quantification_variables$finaloutput
    reactiveprogramdata$useful_data=quantification_variables$useful_data


  })

  #Alignment of signals
  observeEvent(input$alignment, {

    reactiveprogramdata$imported_data$dataset=alignment(reactiveprogramdata$imported_data$dataset,reactiveprogramdata$imported_data$buck_step)
    reactiveprogramdata$alignment_check=1
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

  #Selection of ROI. ROI paramters have to be loaded and setting of reactivity prepared
  observeEvent(input$select, {
    if (reactiveprogramdata$beginning==F) return()
    if (reactiveprogramdata$beginning ==T) {
      reactiveprogramdata$ROIdata_subset=reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[as.numeric(input$select), 1]:reactiveprogramdata$ROI_separator[as.numeric(input$select), 2],]
    }

    reactiveprogramdata$change=1
    reactiveprogramdata$stop=0
    reactiveprogramdata$change2=1
    reactiveprogramdata$stop2=0
    reactivequantdata$stop3=0
    reactiveprogramdata$roi=NULL
    reactiveprogramdata$info=c()


    resetInput(session, "ROIdata_edit")
    resetInput(session, "directedition_edit")
    reactiveROItestingdata$ROIpar <- reactiveprogramdata$ROIdata_subset
    reactiveROItestingdata$signpar <- rbind(rep(NA,7),rep(NA,7))
    colnames(reactiveROItestingdata$signpar)=c("intensity",	"shift",	"half_band_width",	"gaussian",	"J_coupling",	"multiplicities",	"roof_effect")
    reactiveROItestingdata$qualitypar <- rbind(rep(NA,3),rep(NA,3))

    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','fitting_error','signal/total spectrum ratio')
    ROI_limits=c(reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,1]))],reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,2]))])
    ind=ifelse(is.null(input$x1_rows_selected),1,input$x1_rows_selected)
    dummy=type_plot(reactiveprogramdata$imported_data,ROI_limits,ind,reactiveprogramdata$medianplot,reactiveprogramdata$clusterplot)
    if (!is.null(dummy)) reactiveprogramdata$plot=dummy

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
            reactiveprogramdata$roi=0
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
              reactiveprogramdata$roi=1
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

    })


  })

  #Selection of spectra, or of cluster or median plots
  observeEvent(input$x1_rows_selected, {
    if (reactiveprogramdata$beginning==F) return()
    if (reactiveprogramdata$beginning ==T) {
      reactiveprogramdata$ROIdata_subset=reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[as.numeric(input$select), 1]:reactiveprogramdata$ROI_separator[as.numeric(input$select), 2],]
    }

    reactiveprogramdata$change=1
    reactiveprogramdata$stop=0
    reactiveprogramdata$change2=1
    reactiveprogramdata$stop2=0
    reactivequantdata$stop3=0
    reactiveprogramdata$roi=NULL
    reactiveprogramdata$info=c()


    resetInput(session, "directedition_edit")

    reactiveROItestingdata$ROIpar <- reactiveprogramdata$ROIdata_subset
    reactiveROItestingdata$signpar <- rbind(rep(NA,7),rep(NA,7))
    colnames(reactiveROItestingdata$signpar)=c("intensity",	"shift",	"half_band_width",	"gaussian",	"J_coupling",	"multiplicities",	"roof_effect")
    reactiveROItestingdata$qualitypar <- rbind(rep(NA,3),rep(NA,3))
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','fitting_error','signal/total spectrum ratio')

    ROI_limits=c(reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,1]))],reactiveprogramdata$imported_data$ppm[which.min(abs(reactiveprogramdata$imported_data$ppm-reactiveprogramdata$ROIdata_subset[1,2]))])

    dummy=type_plot(reactiveprogramdata$imported_data,ROI_limits,input$x1_rows_selected,reactiveprogramdata$medianplot,reactiveprogramdata$clusterplot)
    if (!is.null(dummy)) reactiveprogramdata$plot=dummy


  })

  #Individual quantification
  observeEvent(input$action, {
    print(reactiveprogramdata$info)
    print(input$x1_rows_selected)
    print(reactiveROItestingdata$ROIpar)
    is_autorun='N'
    if(length(reactiveprogramdata$info)==0) reactiveprogramdata$ind=input$x1_rows_selected-2
    if (length(reactiveprogramdata$ind)!=1) {
      print('Select one valid spectrum')
      return(NULL)
    }
    print(reactiveprogramdata$ind)
    print(dim(reactiveprogramdata$imported_data$dataset))
    # print(str(reactiveprogramdata$useful_data))
    print(dim(reactiveprogramdata$finaloutput$fitting_error))


    reactivequantdata$method1 <- not_automatic_quant(reactiveprogramdata$imported_data, reactiveprogramdata$finaloutput, reactiveprogramdata$ind,reactiveROItestingdata$ROIpar,reactiveprogramdata$useful_data)
    # if ( is.null(reactivequantdata$method1$Ydata)) {
    #   print('Quantification probably worse than current one. Quantification not changed')
    #   return()
    # }
    reactiveprogramdata$plot=reactivequantdata$method1$p
    reactivequantdata$stop3=1
    reactiveROItestingdata$qualitypar=cbind(reactivequantdata$method1$results_to_save$Area,reactivequantdata$method1$results_to_save$fitting_error,reactivequantdata$method1$results_to_save$signal_area_ratio,interface=='T')
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','Fitting Error','Signal/total area ratio')
    rownames(reactiveROItestingdata$qualitypar)=rownames(reactivequantdata$method1$plot_data)[4:(3+nrow(reactiveROItestingdata$qualitypar))]

    if (!is.null(reactivequantdata$method1$signals_parameters)) {
      # reactiveprogramdata$bgColScales = c(rep("", dim(reactivequantdata$method1$signals_parameters)[1]), rep("info", dim(reactivequantdata$method1$signals_parameters_2)[1]-dim(reactivequantdata$method1$signals_parameters)[1]))
      ind=which(rownames(reactiveROItestingdata$qualitypar)=='additional signal')
      reactiveprogramdata$bgColScales = rep(c("","info"),times=c(length(rownames(reactiveROItestingdata$qualitypar))-length(ind),length(ind)))
      reactiveROItestingdata$signpar <- t(reactivequantdata$method1$signals_parameters)
      reactiveprogramdata$stop=0
      reactiveprogramdata$roi=1
    }
  })

  #Qunatification of all spectra in the ROI:
  observeEvent(input$autorun_signal, {
    is_autorun='Y'
    reactivequantdata$chor <- not_automatic_quant(reactiveprogramdata$imported_data, reactiveprogramdata$finaloutput, seq(nrow(reactiveprogramdata$imported_data)),reactiveROItestingdata$ROIpar,reactiveprogramdata$useful_data,interface=='T')
    reactiveprogramdata$finaloutput=reactivequantdata$chor$finaloutput
    reactiveprogramdata$useful_data=reactivequantdata$chor$useful_data
  })


  #Remove quantification or save qunatification or ROI profile edited
  observeEvent(input$remove_q, {
    if (!is.null(reactiveprogramdata$imported_data$signals_names[reactiveprogramdata$info$col])) {
      ind=which(reactiveprogramdata$ROI_data[,4]==reactiveprogramdata$imported_data$signals_names[reactiveprogramdata$info$col])
    } else {
      ind=as.numeric(input$select)
    }

    reactiveprogramdata$finaloutput <- remove_quant(reactiveprogramdata$info,reactiveprogramdata$imported_data, reactiveprogramdata$finaloutput)
  })
  observeEvent(input$save_results, {
    if (is.null(reactivequantdata$method1$Ydata)) {
      return(NULL)
    }

    dummy=save_roi_testing(reactivequantdata$method1,reactiveprogramdata$imported_data, reactiveprogramdata$finaloutput,reactiveprogramdata$useful_data)
    reactiveprogramdata$finaloutput=dummy$finaloutput
    reactiveprogramdata$useful_data=dummy$useful_data


  })
  observeEvent(input$save_profile, {
    if (length(reactiveprogramdata$info$col)>0) {
      ind=which(reactiveprogramdata$ROI_separator[,2]-reactiveprogramdata$info$col>=0)[1]
    } else {
      ind=as.numeric(input$select)
    }
    reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[ind, 1]:reactiveprogramdata$ROI_separator[ind, 2],]=reactiveROItestingdata$ROIpar
    ROI_names=paste(reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],1],reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],2])
    names(reactiveprogramdata$select_options)=ROI_names

  })

  output$qualitypar <- renderD3tf({
    d3tf(reactiveROItestingdata$qualitypar,
      tableProps = list(btn_reset = TRUE,sort_config = list(sort_types = c("String", rep("Number", ncol(reactiveROItestingdata$qualitypar))))),
      enableTf = F,
      edit=F,
      showRowNames = TRUE,
      tableStyle = "table table-bordered")
  })


  output$plot <- renderPlotly({
    if (reactiveprogramdata$beginning==F | is.null(input$x1_rows_selected)) return()
      print(reactiveprogramdata$plot)
  })
  #Repository table
  output$repository = DT::renderDataTable(
    reactiveprogramdata$imported_data$repository[which(reactiveprogramdata$imported_data$repository[,3]>reactiveROItestingdata$ROIpar[1,2]&reactiveprogramdata$imported_data$repository[,3]<reactiveROItestingdata$ROIpar[1,1]),] , server = T)
  proxy2 = dataTableProxy('repository')
  observe({
    replaceData(proxy2,  reactiveprogramdata$imported_data$repository[which(reactiveprogramdata$imported_data$repository[,3]>reactiveROItestingdata$ROIpar[1,2]&reactiveprogramdata$imported_data$repository[,3]<reactiveROItestingdata$ROIpar[1,1]),] )
  })

  #Direct edition of parameters before quantification
  output$directedition <- renderD3tf({

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
            reactiveprogramdata$roi=0
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
      tableProps = list(btn_reset = TRUE,sort_config = list(sort_types = c("String", rep("Number", ncol(reactiveROItestingdata$signpar))))),
      enableTf = F,
      edit=TRUE,
      rowStyles = reactiveprogramdata$bgColScales,
      tableStyle = "table table-bordered")

  })

  #Quantification after direct edition of paramters
  observeEvent(input$direct_edition, {
    reactivequantdata$method1=reactivequantdata$method2=signals_int(reactiveprogramdata$imported_data, reactiveprogramdata$finaloutput,reactiveprogramdata$ind,reactiveROItestingdata$signpar,reactiveROItestingdata$ROIpar)
    reactiveROItestingdata$qualitypar=cbind(reactivequantdata$method2$results_to_save$Area,reactivequantdata$method2$results_to_save$fitting_error,reactivequantdata$method2$results_to_save$signal_area_ratio)
    colnames(reactiveROItestingdata$qualitypar)=c('Quantification','fitting_error','signal/total spectrum ratio')
    rownames(reactiveROItestingdata$qualitypar)=rownames(reactivequantdata$method2$plot_data)[-c(1, 2, 3)]
  })


  ## THIRD TAB REACTIVE OUTPUTS

  #Creation of table to check quantifications with the parameter chosen by the user
  observe({
    if (reactiveprogramdata$beginning==F) return()
    validation_data=validation(reactiveprogramdata$finaloutput,reactiveprogramdata$autourn_data$program_parameters,input$select_validation,reactiveprogramdata$ROI_data,reactiveprogramdata$imported_data$Metadata)
    output$fit_selection = DT::renderDataTable({ datatable(round(validation_data$alarmmatrix,4),selection = list(mode = 'single', target = 'cell')) %>% formatStyle(colnames(validation_data$alarmmatrix), backgroundColor = styleInterval(validation_data$brks, validation_data$clrs))
    })
  })

  #Loading of quantification
  observeEvent(input$fit_selection_cell_clicked, {
if (length(input$fit_selection_cell_clicked)<1) return()
    #Checks and setting of parameters
    reactiveprogramdata$info=input$fit_selection_cell_clicked
    reactiveprogramdata$ind=reactiveprogramdata$info$row
    reactivequantdata$method2=NULL
    reactiveprogramdata$change=reactiveprogramdata$change2=1
    reactiveprogramdata$stop=reactiveprogramdata$stop2=0

    if (length(reactiveprogramdata$info$row)>0) reactivequantdata$stop3=1
    resetInput(session, "ROIdata_edit")
    resetInput(session, "directedition_edit")
    updateSelectInput(session, "select",selected = NULL)

    # if (length(reactiveprogramdata$info$row)!=1) return(NULL)

    #Loading of necessary paramteers
    # dummy=tryCatch({load_quantification(reactiveprogramdata$useful_data,reactiveprogramdata$imported_data,reactiveprogramdata$finaloutput,c(reactiveprogramdata$info$row,reactiveprogramdata$info$col),reactiveprogramdata$ROI_separator)},error= function(e) {
    #   print('Improper input parameters')
    #   return()
    # })
    dummy=load_quantification(reactiveprogramdata$useful_data,reactiveprogramdata$imported_data,reactiveprogramdata$finaloutput,reactiveprogramdata$info,reactiveprogramdata$ROI_separator)

      reactiveprogramdata$plot=dummy$plot
      reactiveROItestingdata$signpar=dummy$signpar
      reactiveROItestingdata$ROIpar=dummy$ROIpar
      reactiveROItestingdata$qualitypar=dummy$qualitypar

    #Redirect to quantification tab
    updateTabsetPanel(session, "mynavlist",selected = "Individual Quantification")

  })


  ## FOURTH TAB REACTIVE OUTPUTS

  #Add and remove signals and save changes
  observeEvent(input$add_signal, {
    reactiveprogramdata$ROI_data_check=rbind(reactiveprogramdata$ROI_data_check,rep(NA,ncol(reactiveprogramdata$ROI_data_check)))
  })
  observeEvent(input$remove_signal, {
    reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data_check[-input$roi_profiles_select,]
    reactiveprogramdata$imported_data$signals_names=reactiveprogramdata$imported_data$signals_names[-input$roi_profiles_select]
    resetInput(session, "roi_profiles_edit")
  })
  observeEvent(input$save_changes, {
    reactiveprogramdata$ROI_data_check=reactiveprogramdata$ROI_data_check[sort(reactiveprogramdata$ROI_data_check[,1],index.return=T)$ix,]
    new_fitting_error=new_intensity=new_signal_area_ratio=new_shift=new_width=new_Area=matrix(NA,nrow(reactiveprogramdata$finaloutput$signal_area_ratio),nrow(reactiveprogramdata$ROI_data_check),dimnames=list(reactiveprogramdata$imported_data$Experiments,paste(reactiveprogramdata$ROI_data_check[,4],reactiveprogramdata$ROI_data_check[,5],sep='_')))
    new_signals_codes=new_signals_names=rep(NA,nrow(reactiveprogramdata$ROI_data_check))
    new_useful_data=reactiveprogramdata$useful_data
    for (i in 1:length(new_useful_data)) new_useful_data[[i]]=vector("list", nrow(reactiveprogramdata$ROI_data_check))
    for (i in 1:nrow(reactiveprogramdata$ROI_data_check)) {
      ind=which(reactiveprogramdata$ROI_data[,4]==reactiveprogramdata$ROI_data_check[i,4]&reactiveprogramdata$ROI_data[,5]==reactiveprogramdata$ROI_data_check[i,5])
      if (length(ind)>0) {
        new_fitting_error[,i]=reactiveprogramdata$finaloutput$fitting_error[,ind]
        new_intensity[,i]=reactiveprogramdata$finaloutput$intensity[,ind]
        new_signal_area_ratio[,i]=reactiveprogramdata$finaloutput$signal_area_ratio[,ind]
        new_shift[,i]=reactiveprogramdata$finaloutput$shift[,ind]
        new_width[,i]=reactiveprogramdata$finaloutput$half_band_width[,ind]
        new_Area[,i]=reactiveprogramdata$finaloutput$Area[,ind]
        new_signals_codes[i]=reactiveprogramdata$imported_data$signals_codes[ind]
        new_signals_names[i]=reactiveprogramdata$imported_data$signals_names[ind]
        for (j in 1:length(new_useful_data)) new_useful_data[[j]][[i]]=reactiveprogramdata$useful_data[[j]][[ind]]
      }
    }
    reactiveprogramdata$finaloutput$fitting_error=new_fitting_error
    reactiveprogramdata$finaloutput$intensity=new_intensity
    reactiveprogramdata$finaloutput$signal_area_ratio=new_signal_area_ratio
    reactiveprogramdata$finaloutput$shift=new_shift
    reactiveprogramdata$finaloutput$half_band_width=new_width
    reactiveprogramdata$finaloutput$Area=new_Area
    reactiveprogramdata$imported_data$signals_codes=new_signals_codes
    reactiveprogramdata$imported_data$signals_names=new_signals_names
    reactiveprogramdata$useful_data=new_useful_data
    reactiveprogramdata$ROI_data=reactiveprogramdata$ROI_data_check

    dummy = which(is.na(reactiveprogramdata$ROI_data[, 1]))
    if (length(dummy)==0) dummy=dim(reactiveprogramdata$ROI_data)[1]+1
    lal=which(duplicated(reactiveprogramdata$ROI_data[-dummy,1:2])==F)
    reactiveprogramdata$ROI_separator = cbind(lal, c(lal[-1] - 1, dim(reactiveprogramdata$ROI_data[-dummy,])[1]))
    ROI_names=paste(reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],1],reactiveprogramdata$ROI_data[reactiveprogramdata$ROI_separator[, 1],2])
    reactiveprogramdata$select_options=1:length(ROI_names)
    names(reactiveprogramdata$select_options)=ROI_names
    mm=matrix(NA,2,dim(reactiveprogramdata$imported_data$Metadata)[2])
    colnames(mm)=colnames(reactiveprogramdata$imported_data$Metadata)
    spectra=cbind(c('Exemplars','Median Spectrum per group',reactiveprogramdata$imported_data$Experiments),rbind(mm,reactiveprogramdata$imported_data$Metadata))
    colnames(spectra)=c('spectrum',colnames(mm))
    output$x1 = DT::renderDataTable(
      spectra , selection = list(mode = 'multiple', selected = 1),server = T)


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
          } else if (col %in% c(1:2,5:11)){
            # numeric columns
            if(is.na(suppressWarnings(as.numeric(val)))) {
              oldval <- reactiveprogramdata$ROI_data_check[row, col]
              rejectEdit(session, tbl = "roi_profiles", row = row, col = col, id = id, value = oldval)
              reactiveprogramdata$roi=0
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
              reactiveprogramdata$roi=1
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

  })

  ## FIFTH TAB REACTIVE OUTPUTS

  #Boxplot plot
  output$plot_p_value_2 <- renderPlotly({
    type_analysis_plot(reactiveprogramdata$finaloutput$Area,reactiveprogramdata$finaloutput,reactiveprogramdata$imported_data,type='boxplot')
  })
  #PCA plot
  output$pcascores <- renderPlotly({
    type_analysis_plot(reactiveprogramdata$finaloutput$Area,reactiveprogramdata$finaloutput,reactiveprogramdata$imported_data,type='pca')
  })

  ## SIXTH TAB REACTIVE OUTPUTS
  #Dengrogran heatmaps for quantification and chemical shift
  output$dendheatmapareadata <- renderPlotly({
    type_analysis_plot(reactiveprogramdata$finaloutput$Area,reactiveprogramdata$finaloutput,reactiveprogramdata$imported_data,type='dendrogram_heatmap')
  })

  output$dendheatmapshiftdata <- renderPlotly({
    type_analysis_plot(reactiveprogramdata$finaloutput$shift,reactiveprogramdata$finaloutput,reactiveprogramdata$imported_data,type='dendrogram_heatmap')
  })



}
