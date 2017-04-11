# Dolphin

## Please run these commands in the R console to install the package:

1. `install.packages(c("devtools","htmltools","shiny","spam","maps","gtable","scales","dendextend","base64enc","dplyr","purrr","tidyr","viridisLite","seriation","randomForest","itertools","csvy","data.table","haven","openxlsx","readODS","readxl","urltools","fit.models","rrcov","RJSONIO","DEoptimR"))`               (these are dependencies from packages called by Dolphin not correctly installed by install_github) 
2. `library(devtools)`           (loads package devtools to use the install_github function)
3. `devtools::install_github("ThomasSiegmund/D3TableFilter")`             (installs D3TableFilter) 
4. `source("https://bioconductor.org/biocLite.R"); biocLite("MassSpecWavelet")`            (installs MassSpecWavelet) 
5. `devtools::install_github("danielcanueto/Dolphin")`           (installs Dolphin)
6. `library(Dolphin)`          (loads Dolphin)

If, at some moment, the installation fails, probably the reason is a missing dependency. Please read on the console for "there is no package called" messages, install the package required with `install.packages` and run `devtools::install_github("danielcanueto/Dolphin")` again.


A vignette is available for the use of functions outside the GUI. `Dolphin_GUI()` opens the GUI on the browser.
