#' Dolphin Shiny GUI
#'
#' @return Dolphin Shiny GUI
#' @export Dolphin_GUI
#' @import shiny
#' @import shinyFiles
#'
#' @examples Dolphin_GUI()
#'
Dolphin_GUI <- function() {
  runApp(system.file("app", package = "Dolphin"), launch.browser = TRUE)
}
