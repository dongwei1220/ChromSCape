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
# e.g. sample_id, batch_id, cell_count, chromatin_group...

colors_scExp <- function(scExp, annotCol = "sample_id",
                        color_by = "sample_id", color_df = NULL){
  
  stopifnot(is(scExp,"SingleCellExperiment"), is.character(annotCol), is.character(color_by))
  
  # initialization
  annot = as.data.frame(colData(scExp))
  anocol <- geco.annotToCol2(annotS = annot[,annotCol,drop = F], annotT = annot, plotLegend=F, categCol=NULL)
  colData(scExp)[,paste0(annotCol,"_color")] = as.data.frame(anocol,stringsAsFactors = F)[,annotCol] # factor or not ?
  
  if(!is.null(color_df)) { #add custom colors
    if(! color_by %in% colnames(color_df)) 
      stop("ChromSCape::color_scExp - color_by must be present in colnames of color_df is not null.")

    colData(scExp)[,paste0(annotCol,"_color")] = color_df[match(colData(scExp)[,color_by],
                                                                color_df[,color_by]), paste0(color_by,"_color"), drop=F] 
  }
  return(scExp)
}

#Internal function to get color dataframefrom list of colourpicker::colorInput
get_color_dataframe_from_input <- function(input, levels_selected, color_by =c("sample_id","total_counts"),
                                           input_id_prefix = "color_")
{
  stopifnot(!is.null(input), is.character(levels_selected), is.character(color_by))

  #Get colors
  color_list <- paste0("list(",paste0(levels_selected," = ",input_id_prefix,
                                      levels_selected, collapse = ", "), ")")
  print(color_list)
  color_list <- eval(parse(text = color_list))
  # Transform into dataframe with right column names
  print(color_list)
  color_df = as.matrix(color_list) %>% as.data.frame(stringsAsFactors = F) %>% rownames_to_column(color_by)
  print(color_df)
  color_df[,2] = as.character(color_df[,2])
  colnames(color_df)[2] = paste0(color_by,"_color")
  
  return(color_df)
}

#Wrapper for plotting PCA & TSNE
plot_reduced_dim_scExp <- function(scExp, color_by = "sample_id", reduced_dim = c("PCA","TSNE"),
                                   select_x = "PC1", select_y = "PC2", annot_label = "none")
  {
  
  stopifnot(is(scExp,"SingleCellExperiment"), is.character(color_by), 
            is.character(reduced_dim),  is.character(select_x), is.character(select_y),
            is.character(annot_label))
  
  if(! reduced_dim[1] %in% reducedDimNames(scExp)) 
    stop(paste0("ChromSCape::plot_reduced_dim_scExp - ",reduced_dim[1],
                " is not present in object, please run normalize_scExp first."))
  
  if(! color_by %in% colnames(colData(scExp))) 
    stop("ChromSCape::plot_reduced_dim_scExp - color_by must be present in colnames of colData(scExp).")
  
  if(! select_x %in% colnames(reducedDim(scExp,reduced_dim[1])) ) 
    stop("ChromSCape::plot_reduced_dim_scExp - select_x must be present in colnames of PCA of scExp.")
  
  if(! select_y %in% colnames(reducedDim(scExp,reduced_dim[1])) ) 
    stop("ChromSCape::plot_reduced_dim_scExp - select_y must be present in colnames of PCA of scExp.")
  
  plot_df = as.data.frame(cbind(reducedDim(scExp,reduced_dim[1]),colData(scExp)) )
 
  p <- ggplot(plot_df, aes_string(x = select_x, y = select_y)) + 
    geom_point(alpha=0.6, aes(color = colData(scExp)[, color_by])) +
    labs(color = color_by) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour="black"),
          panel.border = element_rect(colour="black", fill=NA))
  
  if(annot_label != 'none'){
    p <- p + geom_text(aes(label=colData(scExp)[, annot_label])) 
  }
  if(color_by == 'total_counts'){
    p <- p + scale_color_gradientn(colours = matlab.like(100))
  } else {
    p <- p + scale_color_manual(values = unique(as.character(colData(scExp)[,paste0(color_by,"_color")])) )
  }
  return(p)
}

# Wrapper for geco.hclustAnnotHeatmapPlot Heatmap
plot_heatmap_scExp <- function(scExp, name_hc = "hc_cor", corColors = colorRampPalette(c("royalblue","white","indianred1"))(256)){
  
  #make a variety of sanity check
  stopifnot(is(scExp,"SingleCellExperiment"))
  if( ! "Cor" %in% reducedDimNames(scExp)) 
    stop("ChromSCape::plot_heatmap_scExp - No correlation, run correlation_and_hierarchical_clust_scExp before filtering.")
  
  if( ! name_hc %in% names(scExp@metadata))
    stop("ChromSCape::plot_heatmap_scExp - No dendrogram, run correlation_and_hierarchical_clust_scExp before filtering.")
  
  if( length(scExp@metadata[[name_hc]]$order) != ncol(scExp))
    stop("ChromSCape::plot_heatmap_scExp - Dendrogram has different number of cells than dataset.")
  
  anocol = as.matrix(colData(scExp)[,grep("_color",colnames(colData(scExp)))])
                    
  geco.hclustAnnotHeatmapPlot(
    x=reducedDim(scExp,"Cor")[scExp@metadata[[name_hc]]$order,scExp@metadata[[name_hc]]$order],
    hc=scExp@metadata[[name_hc]],
    hmColors=corColors,
    anocol=anocol[scExp@metadata[[name_hc]]$order,],
    xpos=c(0.15, 0.9, 0.164, 0.885),
    ypos=c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
    dendro.cex=0.04,
    xlab.cex=0.8,
    hmRowNames=FALSE)
}


# Differential barplot
plot_differential_summary_scExp <- function(scExp_cf){
  
  #make a variety of sanity check
  stopifnot(is(scExp_cf,"SingleCellExperiment"))
  if( is.null(scExp_cf@metadata$diff)) 
    stop("ChromSCape::differential_barplot_scExp - No DA, please run differential_analysis_scExp first.")
  
  summary = scExp_cf@metadata$diff$summary
  myylim <- range(c(summary["over",], -summary["under", ]))
  barplot(summary["over",], col="red", las = 1, ylim = myylim, main="Differentially bound regions",
          ylab="Number of regions", axes = F)
  barplot(-summary["under", ], col="forestgreen", ylim = myylim, add = T, axes = F, names.arg="")
  z <- axis(2, pos=-10)
  axis(2, at = z, labels = abs(z), las = 1)
}

# Differential H1 distribution plot
plot_differential_H1_scExp <- function(scExp_cf, chromatin_group = "C1"){
  
  #make a variety of sanity check
  stopifnot(is(scExp_cf,"SingleCellExperiment"), is.character(chromatin_group))
  if( is.null(scExp_cf@metadata$diff)) 
    stop("ChromSCape::differential_H1_plot_scExp - No DA, please run differential_analysis_scExp first.")
  
  if( ! chromatin_group %in% scExp_cf@metadata$diff$groups) 
    stop("ChromSCape::differential_H1_plot_scExp - Chromatin group specified doesn't correspond to differential 
         analysis, please rerun run_differential_analysis_scExp first with correct parameters.")
  
  res = scExp_cf@metadata$diff$res
  
  tmp <- geco.H1proportion(res[, paste("pval", chromatin_group, sep=".")])
  hist(res[, paste("pval", chromatin_group, sep=".")], breaks = seq(0, 1, by = 0.05), xlab="P-value",
       ylab="Frequency", main = paste(chromatin_group, "vs the rest", "\n", "H1 proportion:", round(tmp, 3)))
}

# Differential H1 distribution plot
plot_differential_volcano_scExp <- function(scExp_cf, chromatin_group = "C1", cdiff.th = 1, qval.th = 0.01){
  
  #make a variety of sanity check
  stopifnot(is(scExp_cf,"SingleCellExperiment"), is.character(chromatin_group),
            is.numeric(qval.th), is.numeric(cdiff.th))
  
  if( is.null(scExp_cf@metadata$diff)) 
    stop("ChromSCape::differential_volcano_plot_scExp - No DA, please run differential_analysis_scExp first.")
  
  if( ! chromatin_group %in% scExp_cf@metadata$diff$groups) 
    stop("ChromSCape::differential_volcano_plot_scExp - Chromatin group specified doesn't correspond to differential 
         analysis, please rerun differential_analysis_scExp first with correct parameters.")
  
  res = scExp_cf@metadata$diff$res
  summary = scExp_cf@metadata$diff$summary
  
  mycol <- rep("black", nrow(res))
  mycol[which(res[, paste("qval", chromatin_group, sep=".")] <= qval.th & res[, paste("cdiff", chromatin_group, sep=".")] > cdiff.th)] <- "red"
  mycol[which(res[, paste("qval", chromatin_group, sep=".")] <= qval.th & res[, paste("cdiff", chromatin_group, sep=".")] < -cdiff.th)] <- "forestgreen"
  plot(res[, paste("cdiff", chromatin_group, sep=".")], -log10(res[, paste("qval", chromatin_group, sep=".")]),
       col = mycol, cex = 0.7, pch = 16, xlab="count difference", ylab="-log10(adjusted p-value)", las = 1,
       main = paste(chromatin_group, "vs the rest","\n", summary["over",chromatin_group], "enriched,", summary["under", chromatin_group], "depleted"))
  abline(v = cdiff.th, lty = 2)
  abline(h=-log10(qval.th), lty = 2)
  abline(v=-cdiff.th, lty = 2)
  
}


