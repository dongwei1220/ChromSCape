.onAttach <- function(libname, pkgname) {
  library(shiny)
  library(shinyjs)
  #Inclue js.cookie.js
  shiny::includeScript(file.path(system.file(package="ChromSCape"),"js.cookie.js"))
  print(file.path(system.file(package="ChromSCape"),"js.cookie.js"))
  print(file.exists(file.path(system.file(package="ChromSCape"),"js.cookie.js"))))
  #Load shinyJS added functions
  shinyjs::extendShinyjs(script = file.path(system.file(package="ChromSCape"),"shiny_js_functions.js"),
                         functions = c("init_directory","save_cookie","disableTab","enableTab"))
  shiny::addResourcePath('www',
                         system.file('www',
                                     package = 'ChromSCape'))
}