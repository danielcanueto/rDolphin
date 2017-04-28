library(shiny)
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
    $('#mynavlist a:contains(\"STOCSY and dendrogram heatmaps\")').parent().addClass('disabled')
    }
    Shiny.addCustomMessageHandler('activeNavs', function(nav_label) {
    $('#mynavlist a:contains(\"' + nav_label + '\")').parent().removeClass('disabled')
    })
    ")),

  titlePanel("rDolphin GUI"),
  #First tab
  tabsetPanel(selected="Data Upload", id='mynavlist',
       tabPanel("Data Upload",
                       fluidRow(column(width = 12, h4("Here you can import spectra datasets and save and load profiling sessions. After importing data, profiling results on a model spectrum will appear on the right panel."))),
                                sidebarLayout(
                                  sidebarPanel(
                                    fileInput("file1", "Load a file of parameters of the dataset to profile.",
                                              accept = c("text/csv")
                                    ),
                                    fileInput("file2", "Reanudate a saved profiling session.",
                                              accept = c("text/RData")),
                                    actionButton("save", "Save a profiling session"),
                                    actionButton('folder', 'Save quantification plots'),
                                    br(),
                                    br(),

                                    textInput("caption", "Specify the path of the session to save (e.g. C:/session.RData) or of the folder where to generate the plots folder (e.g. C:/session)", '')
                                    # ,
                                    # fileInput("file3", "Combine data of other sessions",
                                    #   accept = c("text/RData"))



                                  ),
                                  mainPanel(
                                    div(style="display:inline-block",uiOutput('varselect')),
                                    div(style="display:inline-block",uiOutput('align_button')),
                                    div(style="display:inline-block",uiOutput('model_button')),
                                    fluidRow(column(width = 12, h4())),
                                    plotlyOutput("autorun_plot"),
                                    div(dataTableOutput("sp"), style = "font-size:80%")
                                    # div(dataTableOutput("indicators"), style = "font-size:80%")


                                  ))),
                       #Second tab
                       tabPanel("ROI Profiles",
                                fluidRow(column(width = 12, h4("Here you can visually analyze the dataset characteristic traits."))),

                                         selectInput("roi_profile_option",label="Select a possibility",choices=c('Exemplars'=1,'Median spectrum for each kind of sample'=2),selected=1),
                                         plotlyOutput(outputId = "roi_profiles_plot"),
                                         fluidRow(column(width = 12, h4("Here you have a HMDB repository to help with the identification of signals and the choice of ROI parameters."))),

                                                  div(dataTableOutput("repository2"), style = "font-size:80%"),
                                                  fluidRow(column(width = 12, h4("Here you have the current ROI profiles. They can be edited to optimize the profiling."))),

                                                           actionButton("add_hmdb_signal", label = "Add signal from repository"),actionButton("open_hmdb_url", label = "Open signal HMDB website"),actionButton("add_signal", label = "Add signal"),actionButton("remove_signal", label = "Remove signal"),actionButton("save_changes", label = "Save changes"),

                                                           div(d3tfOutput('roi_profiles',width = "100%", height = "auto"), style = "font-size:80%")
                                                  ),
    #Third tab
    tabPanel("Individual Quantification",
      fluidRow(column(width = 12, h4("Here you can supervise the ROI profiles to edit them before the automatic profiling. Here you can also see loaded quantifications on 'Quantifiction validation' tab and optimize them if necessary."))),
      sidebarLayout(

        sidebarPanel(
          actionButton("action", label = "Check quantification"),
          actionButton("save_results", label = "Save quantification"),
          actionButton("remove_q", label = "Remove quantification"),
          actionButton("autorun_signal", label = "ROI autorun"),
          actionButton("save_profile", label = "Save ROI profile"),
          fluidRow(column(width = 12, h4("Select ROI"))),
          selectInput("select",label=NULL,choices=""),
          fluidRow(column(width = 12, h4("Select spectrum"))),
          div(dataTableOutput('x1'), style = "font-size:80%")
        ),


        mainPanel(
          plotlyOutput("plot",height = "250px"),
		  fluidRow(column(width = 12, h4("Here you have some indicators of the quantification."))),
          div(d3tfOutput('qualitypar',width = "100%", height = "auto"), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("You can edit the ROI Profile and quantify it"))),
          div(d3tfOutput('ROIdata',width = "100%", height = "auto"), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("Here you can see the signals in the HMDB Repository located at the same zone of the spectrum, selected by biofluid."))),
          div(dataTableOutput("repository"), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("You can directly edit the signals parameters if you are not satisfied with the calculated parameters."))),
          actionButton("direct_edition", label = "Direct edition"),
          div(d3tfOutput('directedition',width = "100%", height = "auto"), style = "font-size:80%"),
          fluidRow(column(width = 12, h4("If you have selected to upload a 2D file, you can watch it here."))),
          div(style="display:inline-block",uiOutput('jres_plot'))

        )
      )
    ),

    #Fourth tab
    tabPanel("Quantification Validation",
      fluidRow(column(width = 12, h4("Here you have some indicators of quality for every quantification. Press one cell to analyze the quantification."))),
      selectInput("select_validation",label=NULL,choices=c('Fitting Error'=1,'Signal/total area ratio'=2,'Shift'=3,'Halfwidth'=4,'Outliers'=5,'Relative Intensity'=6),selected=NULL),
      div(dataTableOutput("fit_selection"), style = "font-size:80%")
    ),



    #Fifth tab
    tabPanel("Uni and multivariate analysis",
      fluidRow(column(width = 12, h4("Here you have boxplots for every quantified signal, with p values on the x axis labels."))),
      plotlyOutput(outputId = "plot_p_value_2"),
      fluidRow(column(width = 12, h4("PCA with loadings and scores"))),
      plotlyOutput(outputId = "pcascores"))

    #Sixth tab
    ,tabPanel("STOCSY and dendrogram heatmaps",
      fluidRow(column(width = 12, h4("Here you can perform STOCSY to identify unknown signals."))),
      numericInput("left_ppm", "Left edge of region", NA),numericInput("right_ppm", "Right edge of region", NA),selectInput("correlation_method",label=NULL,choices=c('Pearson'='pearson','Spearman'='spearman'),selected='pearson'),selectInput("stocsy",label="Select a possibility",choices=c('Exemplars'=1,'Run STOCSY'=2),selected=1),
      plotlyOutput(outputId = "stocsy_plot"),
      fluidRow(column(width = 12, h4("Here you have the dendrogram heatmap of quantification, so you can analyze relationships between spectra and between signals."))),
      plotlyOutput(outputId = "dendheatmapareadata"),
      fluidRow(column(width = 12, h4("Here you have the dendrogram heatmap of chemical shift, so you can analyze relationships between spectra and between signals."))),
      plotlyOutput(outputId = "dendheatmapshiftdata"))
  )))

