#' rDolphin Shiny GUI
#'
#' @return rDolphin Shiny GUI
#' @export rDolphin_GUI
#' @import shiny

rDolphin_GUI <- function() {
  runApp(system.file("app", package = "rDolphin"), launch.browser = TRUE)
}
