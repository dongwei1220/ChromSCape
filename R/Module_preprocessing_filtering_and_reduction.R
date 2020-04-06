moduleFiltering_and_ReductionUI <- function(id, label = "Filtering_and_Reduction module") {
    ns <- NS(id)
}

moduleFiltering_and_Reduction <- function(input, output, session, raw_dataset_name, min_cov_cell, percentMin, quant_removal, datamatrix, 
    annot_raw, data_folder, annotation_id, exclude_regions) {
    withProgress(message = "Processing data set...", value = 0, {
        
        batch_string <- "uncorrected"
        
        ############################################################### 1. Data loading
        
        incProgress(amount = 0.1, detail = paste("Loading raw data..."))
        
        scExp = create_scExp(datamatrix(), annot_raw(), removeZeroCells = T, removeZeroFeatures = T)
        
        ############################################################### 2. Filtering & Window selection
        incProgress(amount = 0.2, detail = paste("Filtering dataset..."))
        
        scExp = filter_scExp(scExp, min_cov_cell = min_cov_cell(), quant_removal = quant_removal(), percentMin = percentMin())
        
        # Filtering based on exclude-regions from bed file, if provided
        if (!is.null(exclude_regions())) {
            exclude_features_scExp(scExp, exclude_regions(), by = "region")
        }
        
        ############################################################### 3. Normalizing
        incProgress(amount = 0.1, detail = paste("Normalization..."))
        
        scExp = normalize_scExp(scExp, type = "RPKM")
        
        ############################################################### 4. Feature annotation
        
        incProgress(amount = 0.3, detail = paste("Feature annotation..."))
        
        reference_annotation = read.table(file.path("annotation", annotation_id(), "Gencode_TSS_pc_lincRNA_antisense.bed"), col.names = c("chr", 
            "start", "end", "Gene"))
        
        scExp = feature_annotation_scExp(scExp, reference_annotation = reference_annotation)
        
        # ############################################################### # 5. Batch correction
        # ############################################################### pca <- NULL batches <- list() if(doBatchCorr()){
        # incProgress(amount=0.2, detail=paste('Performing batch correction and dimensionality reduction - batch Corr'))
        # annot$batch_name <- 'unspecified' for(i in 1:num_batches()){ for(s_id in batch_sels()[[i]]){ annot[annot$sample_id==s_id,
        # 'batch_name'] <- batch_names()[i] } } adj_annot <- data.frame() # annot must be reordered based on batches b_names <-
        # unique(annot$batch_name) for(i in 1:length(b_names)){ b_name <- b_names[i] batches[[i]] <- mat[, annot$batch_name==b_name]
        # adj_annot <- rbind(adj_annot, annot[annot$batch_name==b_name, ]) } mnn.out <- do.call(fastMNN, c(batches, list(k=20, d=50,
        # approximate=TRUE, auto.order=TRUE, cos.norm=FALSE))) pca <- mnn.out$corrected colnames(pca) <- paste0('PC', 1:dim(pca)[2])
        # annot <- adj_annot } Extract the rotation to do the cross product out.2 <- fastMNN(pcs[[1]], pcs[[2]], pc.input=TRUE)
        # all.equal(head(out,-1), out.2) # should be TRUE (no rotation) # Obtaining corrected expression values for genes 1 and 10.
        # cor.exp <- tcrossprod(out$rotation[c(1,10),], out$corrected) dim(cor.exp)
        
        ############################################################### 5. Dimensionality Reduction
        
        # batch_ids <- unique(annot$batch_id) ## !!!!!!! Put Batch_ids = 1 for testing purpose !!!!!!  batch_ids = 1
        # ###################################################'###
        
        # if(!doBatchCorr()){ print('No Batch Correction') incProgress(amount=0.2, detail=paste('Performing dimensionality reduction- NO
        # batch Corr')) # if(length(batch_ids) > 1){ # print('Detecting different samples') # print(length(unique(annot$batch_id))) #
        # for(i in batch_ids){ # batches[[i]] <- mat[, annot$batch_id==i] # } # }else{ # simulate two batches by splitting data rand <-
        # runif(dim(mat)[2]) > 0.5 batches[[1]] <- mat[, rand] batches[[2]] <- mat[, !rand] # } pca_batches <- do.call(multiBatchPCA,
        # batches) pca <- do.call(rbind, pca_batches) colnames(pca) <- paste0('PC', 1:50) pca <- pca[rownames(annot),]
        
        
        # Original PCA
        print("Running Dimensionality Reduction...")
        
        incProgress(amount = 0.3, detail = paste("Performing Dimensionality Reduction..."))
        
        scExp = reduce_dims_scExp(scExp, verbose = F)
        
        ############################################################### 7. Save data
        
        save(scExp, file = file.path(data_folder(), "datasets", raw_dataset_name(), "reduced_data", paste0(paste(raw_dataset_name(), 
            min_cov_cell(), (percentMin() * 100), quant_removal(), batch_string, sep = "_"), ".RData")))
    })
}
