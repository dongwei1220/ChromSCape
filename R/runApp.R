#' Launch ChromSCape
#'
#' @export
#'
#' @examples ChromScape::launchApp()
#' 
#' @import shiny
#'  
launchApp <- function() {
    shinyApp(ui = shinyAppUI, server = shinyAppServer)
}
