#' rDolphin Shiny GUI
#'
#' @return rDolphin Shiny GUI
#' @export Dolphin_GUI
#' @import shiny
#' @import shinyFiles
#'
#' @examples Dolphin_GUI()
#'
Dolphin_GUI <- function() {
  runApp(system.file("app", package = "rDolphin"), launch.browser = TRUE)
}
