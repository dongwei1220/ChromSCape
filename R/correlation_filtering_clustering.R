## Authors : Pacôme Prompsy Title : Wrappers & functions to filter and cluser
## single cell data based on correlation between cells Description : Wrappers &
## functions to filter and cluser single cell data based on correlation between
## cells

correlation_and_hierarchical_clust_scExp <- function(scExp)
  {
  
  stopifnot(is(scExp, "SingleCellExperiment"))
  
  if(is.null(reducedDim(scExp,"PCA"))) stop("ChromSCape::correlation_and_hierarchical_clust_scExp - Run PCA on the object before correlation") 
    
  pca = reducedDim(scExp,"PCA")
  pca_t <- Matrix::t(pca)
  
  hc_cor = hclust(as.dist(1 - cor(pca_t)), method="ward.D") 
  hc_cor$labels = rep("",length(hc_cor$labels))
  
  scExp@metadata$hc_cor = hc_cor
  reducedDim(scExp,"Cor") <- cor(pca_t)
  
  return(scExp)
}

filter_correlated_cell_scExp <- function(scExp, random_iter = 50, corr_threshold = 99, percent_correlation = 1){
  
  stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(random_iter),
            is.numeric(corr_threshold), is.numeric(percent_correlation))
  
  if(is.null(reducedDim(scExp,"Cor"))) stop(
    "ChromSCape::correlation_and_hierarchical_clust_scExp - No correlation, run correlation before filtering.")
  if(is.null(reducedDim(scExp,"PCA"))) stop(
    "ChromSCape::correlation_and_hierarchical_clust_scExp - No PCA, run reduced_dim before filtering.")
  
  pca_t = Matrix::t(reducedDim(scExp,"PCA"))
  correlation_values <- vector(length=random_iter)
  corChIP <- reducedDim(scExp,"Cor")
  limitC <- 0
  
  set.seed(47)
  for(i in 1:random_iter){
    random_mat <-  matrix(sample(pca_t), nrow=dim(pca_t)[1])
    threshold <- quantile(cor(random_mat), probs=seq(0,1,0.01))
    limitC <-  threshold[corr_threshold+1]
    correlation_values[i] = limitC
  }
  
  limitC_mean = mean(correlation_values,na.rm = T)
  
  selection_cor_filtered <- (apply(corChIP, 1, function(x) length(which(x > limitC_mean)))
              ) > (percent_correlation*0.01)*dim(corChIP)[1]
  
  scExp <- scExp[,selection_cor_filtered]
  
  hc_cor_cor_filtered <- hclust(as.dist(1 - reducedDim(scExp,"Cor")), method="ward.D")
  hc_cor_cor_filtered$labels = rep("",length(hc_cor_cor_filtered$labels))
  cor_cor_filtered <- Matrix::t(reducedDim(scExp,"PCA"))[, hc_cor_cor_filtered$order]
  
 
  scExp@metadata$hc_cor = hc_cor_cor_filtered
  scExp()@metadata$limitC = limitC_mean
  return(scExp)
}


num_cell_before_cor_filt <- function(annot, annotColors){
  
  table <- as.data.frame(table(annot$sample_id))
  colnames(table) = c("Sample","#Cells")
  rownames(table) = NULL
  
  #Retrieve sample colors from user specified colors & add to table
  colors = unique(annotColors[,c("sample_id","sample_id_Color")])
  colors = as.vector(as.character(left_join(table,colors,by=c("Sample"="sample_id"))[,"sample_id_Color"]))
  colors = c(col2hex(colors),"")
  
  
  table[,1] = as.character(table[,1])
  table = table %>% bind_rows(., tibble(Sample="",`#Cells`=sum(table[,-1])) )
  
  table %>% mutate(Sample=cell_spec(Sample, color="white", bold=T, background=colors)) %>%
    kable(escape=F, align="c") %>% kable_styling(c("striped", "condensed"), full_width = T) %>%
    group_rows("Total cell count", dim(table)[1], dim(table)[1])
}

num_cell_after_cor_filt <- function(annot,annot_sel){
  
  table <- as.data.frame(table(annot$sample_id))
  table_filtered <- as.data.frame(table(annot_sel$sample_id))
  colnames(table) = c("Sample","#Cells Before Filtering")
  rownames(table) = NULL 
  colnames(table_filtered) = c("Sample","#Cells After Filtering")
  rownames(table_filtered) = NULL 
  
  #Retrieve sample colors from user specified colors & add to table
  colors = unique(annotColors[,c("sample_id","sample_id_Color")])
  colors = as.vector(as.character(left_join(table,colors,by=c("Sample"="sample_id"))[,"sample_id_Color"]))
  colors = c(col2hex(colors),"")
  
  table_both = left_join(table,table_filtered, by=c("Sample"))
  table_both[,1] = as.character(table_both[,1])
  table_both = table_both %>% bind_rows(., tibble(Sample="",`#Cells Before Filtering`=sum(table_both[,2]),`#Cells After Filtering`=sum(table_both[,3]) ) )
  
  table_both %>% mutate(Sample=cell_spec(Sample, color="white", bold=T, background=colors)) %>%
    kable(escape=F, align="c") %>% kable_styling(c("striped", "condensed"), full_width = T) %>% 
    group_rows("Total cell count", dim(table_both)[1], dim(table_both)[1])
}