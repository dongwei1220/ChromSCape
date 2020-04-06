# 0. Importing packages ####
library(scater)
library(scran)
library(tibble)
library(dplyr)
library(stringr)
library(irlba)
library(reshape2)
library(Rtsne)
library(ConsensusClusterPlus)
library(tidyr)
library(IRanges)
library(GenomicRanges)
library(splitstackshape)
library(rlist)

# Shiny
library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)
library(shinyDirectoryInput)

# Graphics
library(RColorBrewer)
library(colorRamps)
library(colourpicker)
library(kableExtra)
library(knitr)
library(viridis)
library(ggplot2)
library(gplots)
library(png)
library(grid)
library(gridExtra)
library(DT)

# global variables options(scipen = 999)
jscode <- "shinyjs.closeWindow = function() { window.close(); }"

gg_fill_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


