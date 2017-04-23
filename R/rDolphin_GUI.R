#' rDolphin Shiny GUI
#'
#' @return rDolphin Shiny GUI
#' @export rDolphin_GUI
#' @import shiny
#' @import shinyFiles

rDolphin_GUI <- function() {
  runApp(system.file("app", package = "rDolphin"), launch.browser = TRUE)
}
