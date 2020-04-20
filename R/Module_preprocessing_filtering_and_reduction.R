moduleFiltering_and_ReductionUI <- function(id, label = "Filtering_and_Reduction module") {
    ns <- NS(id)
}

moduleFiltering_and_Reduction <- function(input, output, session, raw_dataset_name, min_cov_cell, percentMin, quant_removal, datamatrix, 
                                          annot_raw, data_folder, annotation_id, exclude_regions, annotCol, doBatchCorr,
                                          num_batches, batch_names, batch_sels) {
    withProgress(message = "Processing data set...", value = 0, {
        
        
        ### 1. Data loading ###
        batch_string <- if (doBatchCorr()) "batchCorrected" else "uncorrected"
        incProgress(amount = 0.1, detail = paste("Loading raw data..."))
        
        scExp = create_scExp(datamatrix(), annot_raw(), remove_zero_cells = T, remove_zero_features = T)
        
        ### 2. Filtering & Window selection ###
        
        incProgress(amount = 0.2, detail = paste("Filtering dataset..."))
        
        scExp = filter_scExp(scExp, min_cov_cell = min_cov_cell(), quant_removal = quant_removal(), percentMin = percentMin())
        
        # Filtering based on exclude-regions from bed file, if provided
        if (!is.null(exclude_regions())) {
            scExp = exclude_features_scExp(scExp, exclude_regions(), by = "region")
        }
        
        ### 3. Normalizing ###
        
        incProgress(amount = 0.1, detail = paste("Normalization..."))
        
        scExp = normalize_scExp(scExp, type = "CPM")
        
        ### 4. Feature annotation ###
        
        incProgress(amount = 0.3, detail = paste("Feature annotation..."))
        
        scExp = feature_annotation_scExp(scExp, ref = annotation_id())
        
        ### 6. Dimensionality Reduction ###
        
        # pca <- NULL batches <- list() if(doBatchCorr()){
        # incProgress(amount=0.2, detail=paste('Performing batch correction and dimensionality reduction - batch Corr'))
        # annot$batch_name <- 'unspecified' for(i in 1:num_batches()){ for(s_id in batch_sels()[[i]]){ annot[annot$sample_id==s_id,
        # 'batch_name'] <- batch_names()[i] } } adj_annot <- data.frame() # annot must be reordered based on batches b_names <-
        # unique(annot$batch_name) for(i in 1:length(b_names)){ b_name <- b_names[i] batches[[i]] <- mat[, annot$batch_name==b_name]
        # adj_annot <- rbind(adj_annot, annot[annot$batch_name==b_name, ]) } mnn.out <- do.call(fastMNN, c(batches, list(k=20, d=50,
        # approximate=TRUE, auto.order=TRUE, cos.norm=FALSE))) pca <- mnn.out$corrected colnames(pca) <- paste0('PC', 1:dim(pca)[2])
        # annot <- adj_annot } Extract the rotation to do the cross product out.2 <- fastMNN(pcs[[1]], pcs[[2]], pc.input=TRUE)
        # all.equal(head(out,-1), out.2) # should be TRUE (no rotation) # Obtaining corrected expression values for genes 1 and 10.
        # cor.exp <- tcrossprod(out$rotation[c(1,10),], out$corrected) dim(cor.exp)
        
        # if(!doBatchCorr()){ print('No Batch Correction') incProgress(amount=0.2, detail=paste('Performing dimensionality reduction- NO
        # batch Corr')) # if(length(batch_ids) > 1){ # print('Detecting different samples') # print(length(unique(annot$batch_id))) #
        # for(i in batch_ids){ # batches[[i]] <- mat[, annot$batch_id==i] # } # }else{ # simulate two batches by splitting data rand <-
        # runif(dim(mat)[2]) > 0.5 batches[[1]] <- mat[, rand] batches[[2]] <- mat[, !rand] # } pca_batches <- do.call(multiBatchPCA,
        # batches) pca <- do.call(rbind, pca_batches) colnames(pca) <- paste0('PC', 1:50) pca <- pca[rownames(annot),]
        
        # Original PCA
        print("Running Dimensionality Reduction...")
        
        incProgress(amount = 0.3, detail = paste("Performing Dimensionality Reduction..."))
        
        scExp = reduce_dims_scExp(scExp, verbose = F)
        
        ### 7. Add default colors ###
        
        scExp = colors_scExp(scExp, annotCol()) # add colors 
        
        ### 8. Save data ###
        save(scExp, file = file.path(data_folder(), "datasets", raw_dataset_name(), "QC_filtering", paste0(paste(raw_dataset_name(), 
                                                                                                                 min_cov_cell(), percentMin(), quant_removal(), batch_string(), sep = "_"), ".RData")))
        
        print("Filtering & Reduction done !")
    })
}
