context("Testing preprocessing, filtering & reduction functions")



# Functions for testing purposes
create_scDataset_raw <- function(cells=300,features=600,
                                 featureType = c("window","peak","gene"),
                                sparse=T, nsamp=4, ref = "hg38",batch_id = rep(1,nsamp)) {
  
  stopifnot(featureType %in% c("window","peak","gene"), ref %in% c("mm10","hg38"),
            nsamp >= 1, cells >= nsamp, features >=1, length(batch_id) == nsamp)
  
  stopifnot()
  
  set.seed(47)
  
  # Create cell names
  cell_counts = sapply(split( 1:cells , sample(nsamp, cells , repl = TRUE) ), length)
  cell_names = sample = batches = list()
  for(i in 1:length(cell_counts)) {
    cell_names[[i]] = paste0("sample_",i,"_c", 1:cell_counts[i])
    sample[[i]] = rep(paste0("sample_",i),cell_counts[i])
    batches[[i]] = rep(batch_id[i],cell_counts[i])
  }
  cell_names = as.character(unlist(cell_names))
  sample = as.character(unlist(sample))
  batches = as.numeric(unlist(batches))
  
  # Create feature names
  eval(parse(text = paste0("chr <- ChromSCape::",ref,".chromosomes"))) #load species chromsizes
  chr = GRanges(chr)
  
  if(featureType[1] == "window") {
    chr_ranges = unlist(tileGenome(setNames(width(chr),seqnames(chr)),
                                      ntile = features))[1:features] # ~constant window size
    features_names = paste(as.data.frame(chr_ranges)$seqnames,
                           as.data.frame(chr_ranges)$start,
                           as.data.frame(chr_ranges)$end, sep="_")
  }
  if(featureType[1] == "peak") {
    size_peaks = c(1000,2500,7999,10000,150000,10^6) #Different size of peaks
    peaks = sapply(split( 1:features , sample(length(size_peaks), features , repl = TRUE) ), length)
    chr_ranges_list = GRangesList()
    for(i in 1:length(peaks)){
      chr_ranges = unlist(tileGenome(setNames(width(chr),seqnames(chr)),
                                     tilewidth = size_peaks[i], cut.last.tile.in.chrom = F))
      chr_ranges_list[[i]] = chr_ranges[sample(1:length(chr_ranges),size = peaks[i]),]
    }
    chr_ranges = GenomicRanges::sort.GenomicRanges(unlist(chr_ranges_list))[1:features]
    
    features_names = paste(as.data.frame(chr_ranges)$seqnames,
                           as.data.frame(chr_ranges)$start,
                           as.data.frame(chr_ranges)$end, sep="_")
  }
  if(featureType[1] == "gene") {
    eval(parse(text = paste0("chr <- ChromSCape::",ref,".GeneTSS"))) #load species chromsizes
    features_names = as.character(sample(chr$gene,features,replace = F))
  }
  vec = rpois(cells*features,0.5) #Add count to values > 0, iteratively
  for(i in 1:10) vec[vec >= i] = vec[vec >= i]  +  i^2*rpois(length(vec[vec >= i]),0.5)
  mat = matrix(vec, nrow = features, ncol = cells, 
               dimnames = list( features_names,cell_names))
  annot = data.frame(cell_id = cell_names,
                     sample_id = sample,
                     batch_id = batches,
                     total_counts = Matrix::colSums(mat))
  if(sparse) return(list("mat" =  as(mat,"dgCMatrix"), "annot" = annot)) else return(list("mat" =  mat, "annot" = annot))
}

out = create_scDataset_raw(featureType = "window")
mat = out$mat
annot = out$annot
