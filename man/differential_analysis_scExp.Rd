% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_analysis_gene_set_enrichment.R
\name{differential_analysis_scExp}
\alias{differential_analysis_scExp}
\title{Runs differential analysis between cell clusters}
\usage{
differential_analysis_scExp(
  scExp,
  de_type = "one_vs_rest",
  method = "wilcox",
  qval.th = 0.01,
  cdiff.th = 1,
  block = NULL
)
}
\arguments{
\item{scExp}{A SingleCellExperiment object containing consclust with selected
number of cluster.}

\item{de_type}{Type of comparisons. Either 'one_vs_rest', to compare each
cluster against all others, or 'pairwise' to make 1 to 1 comparisons.
('one_vs_rest')}

\item{method}{Wilcoxon or edgerGLM}

\item{qval.th}{Adjusted p-value threshold. (0.01)}

\item{cdiff.th}{Fold change threshold. (1)}

\item{block}{Use batches as blocking factors ?}
}
\value{
Returns a SingleCellExperiment object containing a differential list.
}
\description{
Based on clusters of cell defined previously, runs non-parametric Wilcoxon
Rank Sum test to find significantly depleted or enriched features, in
'one_vs_rest' mode or 'pairwise' mode. In pairwise mode, each cluster is
compared to all other cluster individually, and then pairwise comparisons
between clusters are combined to find overall differential features using
combineMarkers function from scran.
}
\details{
This functions takes as input a SingleCellExperiment object with consclust,
the type of comparison, either 'one_vs_rest' or 'pairwise', the adjusted
p-value threshold (qval.th) and the fold-change threshold (cdiff.th). It
outputs a SingleCellExperiment object containing a differential list.
}
\examples{
data("scExp")
scExp_cf = differential_analysis_scExp(scExp)

}
