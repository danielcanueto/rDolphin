# Dolphin

## Please run these commands in the R console to install the package:

1. `install.packages(c("devtools","htmltools","shiny","spam","maps","gtable","scales","dendextend","base64enc","dplyr","purrr","tidyr","viridisLite","seriation","randomForest","itertools","fit.models","rrcov","RJSONIO","DEoptimR"))`               (these are dependencies from packages called by Dolphin not correctly installed by install_github) 
2. `library(devtools)`           (loads package devtools to use the install_github function)
3. `devtools::install_github("ThomasSiegmund/D3TableFilter")`             (installs D3TableFilter) 
4. `source("https://bioconductor.org/biocLite.R"); biocLite("MassSpecWavelet")`            (installs MassSpecWavelet) 
5. `devtools::install_github("danielcanueto/Dolphin", build_vignettes = TRUE)`           (installs Dolphin)
6. `library(Dolphin)`          (loads Dolphin)

If, at some moment, the installation fails, probably the reason is a missing dependency. Please read on the console for "there is no package called" messages, install the package required with `install.packages` and run `devtools::install_github("danielcanueto/Dolphin")` again.


## Introduction to the package:

For a quick look to the capabilitites of the package, examples of sessions of profiling of MTBLS1 and MTBLS242 Metabolights data sets are available on the next links:

MTBLS1: https://www.dropbox.com/s/8onmpp9hsngtgjc/MTBLS1_example.RData?dl=0
MTBLS242: https://www.dropbox.com/s/zlk2vo0e0iu1z0f/MTBLS242_example.RData?dl=0

For example, load e.g. the MTBLS242 profiling session, run “Dolphin_GUI()” to open the GUI and click on the "Reanudate a saved session" button on the left panel.  


A vignette is available for the use of functions outside the GUI and can be run with this command `browseVignettes("Dolphin")`. A introductory tutorial for the use of Dolphin inside the GUI is available on this link: https://docs.google.com/document/d/12OJHLgs4dbboBQJHUP0YExfIeZQV8l7n4HQmkWZBWwo/edit

