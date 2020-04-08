## Authors : Pacôme Prompsy Title : Wrappers & function to create variety of
## plot to uncover heterogeneity in single cell Dataset

#Wrapper for plotting distribution of signal
plot_distribution_scExp <- function(scExp, raw = T, log10 = F, 
                                    pseudo_counts = 1, bins = 150)
  {
  
  stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(pseudo_counts))
  if( ! raw %in% c(T,F) | ! log10 %in% c(T,F) ) stop(
    "ChromSCape::plot_distribution_scExp - raw and log10 must be true or false.")
  if( raw ==F && !("normcounts" %in% assayNames(scExp))) stop(
    "ChromSCape::plot_distribution_scExp - If raw is false, normcounts must not be empty - run normalize_scExp first.")
  
  if(raw) cell_cov_df = data.frame(coverage=Matrix::colSums(counts(scExp))) else
    cell_cov_df = data.frame(coverage=Matrix::colSums(normcounts(scExp)))
  
  if(log10) cell_cov_df$coverage = log10(cell_cov_df$coverage + pseudo_counts)
  
  ggplot(cell_cov_df, aes(x=coverage)) +
    geom_histogram(color="black", fill="steelblue", bins=bins) +
    labs(x="read coverage per cell") + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
    panel.background=element_blank(), axis.line=element_line(colour="black"),
    panel.border=element_rect(colour="black", fill=NA))
}


#Function to add color to each cell & cell features
# e.g. sample_id, batch_id, cell_count, chromatinGroup...

color_scExp <- function(scExp, annotCol=c("sample_id","total_counts"),
                        color_by ="sample_id"){
  annot = colData(scExp)
  annotCol <-  c("sample_id","total_counts")
  annotColors=NULL
  
  #anocol object
  anocol <- geco.annotToCol2(annotS=annot[, annotCol], annotT=annot, plotLegend=F, categCol=NULL)
  annot[,paste0(annotCol,"_color")] = as.data.frame(anocol)[,annotCol]


  output$color_by <- renderUI({selectInput("color_by", "Color by", choices=c('sample_id', 'total_counts')) })
  input$color_by = "sample_id"
  
  #levels_selected
  levels_selected = annot[,input$color_by] %>% unique() %>% as.vector()
  
  #Color picker
  colsModif <- annot[,c(input$color_by,paste0(input$color_by,"_color"))] %>% unique()
  lapply(seq_along(levels_selected), function(i) {
    colourpicker::colourInput(inputId=paste0("color_", levels_selected[i]),
                              label=paste0("Choose colour for ", levels_selected[i]),
                              value=colsModif[i,paste0(input$color_by,"_color")], returnName=TRUE) ## Add ", palette = "limited"" to get selectable color panel       
  })

  #Update colorbox - colorpickers
  cols <- gg_fill_hue(length(levels_selected))
  lapply(seq_along(levels_selected), function(i) {
    do.call(what="updateColourInput", args=list(session=session, inputId=paste0("color_", levels_selected[i]), value=cols[i]))
  })
  
  #Update anocol
  anocol <- geco.annotToCol2(annotS=annot[, annotCol], annotT=annot, plotLegend=F, categCol=NULL)
  annot[,paste0(annotCol,"_color")] = as.data.frame(anocol)[,annotCol]
  
  if(!is.null(annotColors)){
    lapply(colnames(ac), function(col){ if(col != "total_counts"){
      ac[, col] <<- as.character(annotColors[which(annotColors$Sample %in% rownames(ac)),
                                             paste0(col, '_Color')]) } })
  }
  anocol = ac
}

#Internal function to get color datafram from list of colourpicker::colorInput
get_color_dataframe_from_input <- function(input, levels_selected, color_by)
{
  stopifnot(!is.null(input), is.character(levels_selected), is.character(color_by))

  #Get colors
  color_list <- paste0("list(",paste0(levels_selected," = input$color_",
                                      levels_selected, collapse = ", "), ")")
  color_list <- eval(parse(text = color_list))
  # Transform into dataframe with right column names
  color_df = as.matrix(color_list) %>% as.data.frame %>% rownames_to_column(color_by)
  colnames(color_df)[2] = paste0(color_by,"_color")
  return(color_df)
}

#Wrapper for plotting PCA & TSNE
plot_reduced_dim_scExp <- function(){
  plot_df <- as.data.frame(cbind(pca(), annot()))
  p <- ggplot(plot_df, aes_string(x=input$pc_select_x, y=input$pc_select_y)) + geom_point(alpha=0.6, aes(color=annot()[, input$color_by])) +
    labs(color=input$color_by) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.background=element_blank(), axis.line=element_line(colour="black"),
          panel.border=element_rect(colour="black", fill=NA))
  if(input$pca_anno_2D != 'none'){ p <- p + geom_text(aes(label=annot()[, input$pca_anno_2D])) }
  if(input$color_by == 'total_counts'){
    p <- p + scale_color_gradientn(colours = matlab.like(100))
  }else{
    cols <- paste0("c(", paste0("input$color_", sort(levels_selected), collapse = ", "), ")")
    cols <- eval(parse(text = cols))
    req(cols)
    lev <- sort(unique(levels_selected))
    names(cols) <- lev
    p <- p + scale_color_manual(values = cols)
  }
  unlocked$list$pca=T
  p
}

# Wrapper for geco.hclustAnnotHeatmapPlot Heatmap
plot_heatmap_scExp <- function(){
  #make a variety of sanity check 
  # same rownames colnames cor mat
  # got hierarchical clustering or not ?
  # same order of rowData (coloring ) -> implement later 
}
