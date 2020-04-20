#' Launch ChromSCape
#'
#' @export
#'
#' @examples ChromScape::launchApp()
#' 
#' @import shiny
#' @import shinyjs
#'  
launchApp <- function() {
  #Include js.cookie.js
  # shiny::includeScript(file.path(system.file(package="ChromSCape"),"js.cookie.js"))
  # print(file.path(system.file(package="ChromSCape"),"www","js.cookie.js"))
  # print(file.exists(file.path(system.file(package="ChromSCape"),"js.cookie.js")))
  # #Load shinyJS added functions
  # print(file.path(system.file(package="ChromSCape"),"www","shiny_js_functions.js"))
  # print(file.exists(file.path(system.file(package="ChromSCape"),"www","shiny_js_functions.js")))
  # shinyjs::extendShinyjs(script = file.path(system.file(package="ChromSCape"),"www","shiny_js_functions.js"),
  #                        functions = c("init_directory","save_cookie","disableTab","enableTab"))
  print("A")
  # shinyApp(ui = shinyAppUI, server = shinyAppServer)
  shiny::runApp(system.file(package="ChromSCape"))
}
