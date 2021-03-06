% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing_filtering_reduction.R
\name{reduce_dims_scExp}
\alias{reduce_dims_scExp}
\title{Reduce dimensions (PCA, TSNE, UMAP)}
\usage{
reduce_dims_scExp(
  scExp,
  dimension_reductions = c("PCA", "TSNE", "UMAP"),
  n = 50,
  batch_correction = FALSE,
  batch_list = NULL,
  verbose = TRUE
)
}
\arguments{
\item{scExp}{A SingleCellExperiment object.}

\item{dimension_reductions}{A character vector of methods to apply.
(c('PCA','TSNE','UMAP'))}

\item{n}{Numbers of dimensions to keep for PCA. (50)}

\item{batch_correction}{Do batch correction ? (FALSE)}

\item{batch_list}{List of characters. Names are batch names, characters are
sample names.}

\item{verbose}{(TRUE)}
}
\value{
A SingleCellExperiment object containing feature spaces. See
?reduceDims().
}
\description{
Reduce dimensions (PCA, TSNE, UMAP)
}
\examples{

scExp = create_scExp(create_scDataset_raw()$mat,create_scDataset_raw()$annot)
scExp = reduce_dims_scExp(scExp,dimension_reductions=c("PCA","UMAP"))
scExp = normalize_scExp(scExp)
scExp = reduce_dims_scExp(scExp,dimension_reductions=c("PCA","UMAP"))
}
