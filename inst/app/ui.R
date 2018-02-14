suppressMessages(library(shiny))
suppressMessages(library(plotly))
suppressMessages(library(DT))
suppressMessages(library(shinyjs))
css <- "
.shiny-output-error { visibility: hidden; }
.shiny-output-error:before {visibility: hidden; }"
shinyUI(fluidPage(
  tags$style(type="text/css", css),
  shinyjs::useShinyjs(),
  titlePanel("rDolphin GUI"),
  #First tab
  tabsetPanel(selected="tab1", id='mynavlist',
       tabPanel("Main Tab",value = "tab1",
                       fluidRow(column(width = 12, h4("Here you can import spectra datasets and save and load profiling sessions. After importing data, profiling results on a model spectrum will appear on the right panel."))),
                                sidebarLayout(
                                  sidebarPanel(
                                    fileInput("file1", "Load a file of parameters of the dataset to profile.",
                                              accept = c("text/csv")
                                    ),
                                    fileInput("file2", "Resume a saved profiling session.",
                                              accept = c("text/RData")),
                                    actionButton("save", "Save a profiling session"),
                                    actionButton('folder', 'Save quantification plots'),
                                    br(),
                                    br(),

                                    textInput("caption", "Specify the path of the session or of the quantification figures to save ", '')
                                    # ,
                                    # fileInput("file3", "Combine data of other sessions",
                                    #   accept = c("text/RData"))

                                  ),
                                  mainPanel(
                                    shinyjs::hidden(actionButton('automatic_profiling', 'Autorun all spectra', class = "btn-primary")),
                                    shinyjs::hidden(actionButton('alignment', 'Alignment of signals')),
                                    shinyjs::hidden(actionButton('model_spectrum', 'Profile model spectrum again')),

                                    fluidRow(column(width = 12, h4())),
                                    plotlyOutput("automatic_profiling_plot"),
                                    div(dataTableOutput("sp"), style = "font-size:80%")
  ))),
                       #Second tab
tabPanel("ROI Profiles",value = "tab2",
  fluidRow(column(width = 12, h4("Here you can visually analyze the dataset characteristic traits."))),
           selectInput("roi_profile_option",label="Select a possibility",choices=c('Exemplars'=1,'Median spectrum for each kind of sample'=2),selected=1),
           plotlyOutput(outputId = "roi_profiles_plot"),
           fluidRow(column(width = 12, h4("Here you have a HMDB repository to help with the identification of signals and the choice of ROI parameters."))),
          div(dataTableOutput("repository2"), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("Here you have the current ROI profiles. They can be edited to optimize the profiling."))),
                   actionButton("add_hmdb_signal", label = "Add signal from repository"),actionButton("open_hmdb_url", label = "Open signal HMDB website"),actionButton("add_signal", label = "Add signal"),actionButton("remove_signal", label = "Remove signal"),actionButton("save_changes", label = "Confirm changes"),
  div(dataTableOutput("roi_profiles"), style = "font-size:80%")
          ),
#Third tab
    tabPanel("Individual Quantification",value = "tab3",
      fluidRow(column(width = 12, h4("Here you can supervise the ROI profiles to edit them before the automatic profiling. Here you can also see loaded quantifications on 'Quantifiction validation' tab and optimize them if necessary."))),
      sidebarLayout(

        sidebarPanel(
          actionButton("action", label = "Check quantification"),
          actionButton("save_results", label = "Save quantification"),
          actionButton("remove_q", label = "Remove quantification"),
          actionButton("automatic_profiling_signal", label = "ROI automatic_profiling"),
          actionButton("save_profile", label = "Save ROI profile"),
          fluidRow(column(width = 12, h4("Select ROI"))),
          uiOutput("moreControls"),
          fluidRow(column(width = 12, h4("Select spectrum"))),
          div(dataTableOutput('x1'), style = "font-size:80%")
        ),
       mainPanel(
          plotlyOutput("plot",height = "250px"),
		  fluidRow(column(width = 12, h4("Here you have some indicators of the quantification."))),
          div(dataTableOutput('qualitypar',width = "100%", height = "auto"), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("You can edit the ROI Profile and quantify it"))),
          div(dataTableOutput('ROIdata'), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("Here you can see the signals in the HMDB Repository located at the same zone of the spectrum, selected by biofluid."))),
          div(dataTableOutput("repository"), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("You can directly edit the signals parameters if you are not satisfied with the calculated parameters."))),
          actionButton("direct_edition", label = "Direct edition"),
          div(dataTableOutput('directedition'), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("If you have selected to upload a 2D file, you can watch it here."))),
          div(style="display:inline-block",uiOutput('jres_plot'))

        )
      )
    ),

    #Fourth tab
    tabPanel("Quantification Validation",value = "tab4",
      fluidRow(column(width = 12, h4("Here you have some indicators of quality for every quantification. Press one cell to analyze the quantification."))),
      selectInput("select_validation",label=NULL,choices=c('Choose a method'=0,'Fitting Error'=1,'Signal/total area ratio'=2,'Shift'=3,'Halfwidth'=4,'Intensity'=5),selected=0),
      div(dataTableOutput("fit_selection"), style = "font-size:80%")
    ),



    # #Fifth tab
    # tabPanel("Uni and multivariate analysis",value = "tab5",
    #   fluidRow(column(width = 12, h4("Here you have boxplots for every quantified signal, with p values on the x axis labels."))),
    #   plotlyOutput(outputId = "plot_p_value_2"),
    #   fluidRow(column(width = 12, h4("PCA with loadings and scores"))),
    #   plotlyOutput(outputId = "pcascores"))

    #Sixth tab
    tabPanel("Identification tools",value = "tab6",
      fluidRow(column(width = 12, h4("Here you can perform spectroscopic methods to identify unknown signals."))),
      div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("left_ppm", "Left edge of region", NA)),
          div(style="display: inline-block;vertical-align:top; width: 150px;",numericInput("right_ppm", "Right edge of region", NA)),
              div(style="display: inline-block;vertical-align:top; width: 150px;",selectInput("correlation_method",label="Select method",choices=c('STOCSY Pearson'='pearson','STOCSY Spearman'='spearman', 'RANSY'='ransy'),selected='pearson')),
                  div(style="display: inline-block;vertical-align:top; width: 150px;",selectInput("stocsy",label="Select a possibility",choices=c('Exemplars'=1,'Run spectroscopic method'=2),selected=1)),
      plotlyOutput(outputId = "stocsy_plot"),
      fluidRow(column(width = 12, h4("Here you have the dendrogram heatmap of quantification, so you can analyze relationships between spectra and between signals."))),
      plotlyOutput(outputId = "dendheatmapareadata",height="1100px"),
      fluidRow(column(width = 12, h4("Here you have the dendrogram heatmap of chemical $chemical_shift, so you can analyze relationships between spectra and between signals."))),
      plotlyOutput(outputId = "dendheatmapshiftdata"))
  )))

