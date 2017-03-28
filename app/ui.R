library(shiny)
library(shinyFiles)
library(plotly)
library(DT)
library(D3TableFilter)


shinyUI(fluidPage(
  #Tabs are disabled until data is not laoded
  tags$head(tags$script("
    window.onload = function() {
    $('#mynavlist a:contains(\"Individual Quantification\")').parent().addClass('disabled')
    $('#mynavlist a:contains(\"Quantification Validation\")').parent().addClass('disabled')
    $('#mynavlist a:contains(\"Uni and multivariate analysis\")').parent().addClass('disabled')
    $('#mynavlist a:contains(\"ROI Profiles\")').parent().addClass('disabled')
    $('#mynavlist a:contains(\"Dendrogram heatmaps\")').parent().addClass('disabled')
    }
    Shiny.addCustomMessageHandler('activeNavs', function(nav_label) {
    $('#mynavlist a:contains(\"' + nav_label + '\")').parent().removeClass('disabled')
    })
    ")),

  titlePanel("Dolphin Demo"),
  #First tab
  tabsetPanel(selected="Data Upload and Processing", id='mynavlist',
    tabPanel("Data Upload and Processing",
      sidebarLayout(
        sidebarPanel(
          fileInput("file1", "Load the parameters file. It will automatically do an autorun of the ROI Profiles in the model spectrum. Then click Autorun if you wanna use these profiles for all spectra",
            accept = c("text/csv")
          ),
          fileInput("file2", "Reanudate a saved session",
            accept = c("text/RData")),
          shinySaveButton("save", "Save session", "Save session as ...", filetype=list(RData="RData")),
          shinyDirButton('folder', 'Save quantification plots', 'Please select a folder', FALSE)
          # ,
          # fileInput("file3", "Combine data of other sessions",
          #   accept = c("text/RData"))



        ),
        mainPanel(
          div(style="display:inline-block",uiOutput('varselect')),
          div(style="display:inline-block",uiOutput('align_button')),
          # div(style="display:inline-block",uiOutput('peak_analysis')),
          fluidRow(column(width = 12, h4("You can watch how the signals have been quantified in the spectrum model and, at the same time, an univariate analysis of every bin in the spectrum, according to the metadata given by the user.The idea is that you can analyze other parts of the spectrum with significant differences and add a ROI profile through the 'Profiles' tab."))),
          plotlyOutput("autorun_plot"),
          div(dataTableOutput("sp"), style = "font-size:80%"),
          div(dataTableOutput("indicators"), style = "font-size:80%")


        ))),

    #Second tab
    tabPanel("Individual Quantification",
      fluidRow(column(width = 12, h4("Here you can change the quantifications and save them or remove them if they can't be well quantified. You can also save the edited profile of the ROI"))),
      sidebarLayout(

        sidebarPanel(

          actionButton("save_results", label = "Save Quantification"),
          actionButton("save_profile", label = "Save profile"),
          actionButton("autorun_signal", label = "Autorun of the signal"),
          actionButton("remove_q", label = "Remove quantification"),

          actionButton("action", label = "Quantification (without saving!)"),
          fluidRow(column(width = 12, h4("Select ROI"))),
          selectInput("select",label=NULL,choices=""),
          fluidRow(column(width = 12, h4("Select spectrum"))),
          div(dataTableOutput('x1'), style = "font-size:80%"),
          width=3
        ),


        mainPanel(
          plotlyOutput("plot",height = "250px"),
          div(d3tfOutput('qualitypar',width = "100%", height = "auto"), style = "font-size:80%"),

          fluidRow(column(width = 12, h4("You can edit the ROI Profile and quantify it"))),

          div(d3tfOutput('ROIdata',width = "100%", height = "auto"), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("Here you can see the signals in the HMDB Repository located at the same zone of the spectrum, selected by biofluid"))),

          div(dataTableOutput("repository"), style = "font-size:80%"),

          fluidRow(column(width = 12, h4("You can directly edit the signals parameters if you are not satisfied with the calculated parameters."))),
          actionButton("direct_edition", label = "Direct edition"),

          div(d3tfOutput('directedition',width = "100%", height = "auto"), style = "font-size:80%"),


          fluidRow(column(width = 12, h4("You can watch the uploaded 2D file here"))),

          plotlyOutput("jres_plot",height='250px')

        )
      )
    ),

    #Third tab
    tabPanel("Quantification Validation",
      fluidRow(column(width = 12, h4("Here you some indicators of quality for every quantification. Press one cell to analyze the quantification."))),
      selectInput("select_validation",label=NULL,choices=c('Fitting Error'=1,'Signal/total spectrum ratio'=2,'Shift'=3,'Halfwidth'=4,'Outliers'=5,'Relative Intensity'=6),selected=NULL),
      div(dataTableOutput("fit_selection"), style = "font-size:80%")
    ),

    #Fourth tab
    tabPanel("ROI Profiles",
      actionButton("add_signal", label = "Add signal"),actionButton("remove_signal", label = "Remove signals"),actionButton("save_changes", label = "Save changes"),
      fluidRow(column(width = 12, h4("Here you have the ROI profiles")),
        div(d3tfOutput('roi_profiles',width = "100%", height = "auto"), style = "font-size:80%")
      )),

    #Fifth tab
    tabPanel("Uni and multivariate analysis",
      fluidRow(column(width = 12, h4("Here you have boxplots for every quantified signal, with p values on the x axis"))),
      plotlyOutput(outputId = "plot_p_value_2"),
      fluidRow(column(width = 12, h4("PCA with loadings and scores"))),
      plotlyOutput(outputId = "pcascores"))

    #Sixth tab
    ,tabPanel("Dendrogram heatmaps",
      fluidRow(column(width = 12, h4("Here you have the dendrogram heatmap of quantification, so you can analyze relationships between spectra and between signals "))),
      plotlyOutput(outputId = "dendheatmapareadata"),
      fluidRow(column(width = 12, h4("Here you have the dendrogram heatmap of chemical shift, so you can analyze relationships between spectra and between signals"))),
      plotlyOutput(outputId = "dendheatmapshiftdata"))
  )))

