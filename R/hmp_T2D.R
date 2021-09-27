#' hmp_T2D
#' 
#' This dataset contained a list of data.frames.
#' Raw data is a subset of the data available at: 
#' https://github.com/aametwally/ipop_seasonal
#' The package will be illustrated on longitudinal MO dataset to study the 
#' seasonality of MO expression in patients with diabetes 
#' (see `netOmics` vignette).
#' In this subset we focused on a single individual with 7 timepoints.
#' Briefly 6 different omics were sampled (RNA, proteins, cytokines, 
#' gut microbiome, metabolites and clinical variables).
#' 
#' 
#' @format a list of data.frame
#' \describe{
#'  \item{raw}{data.frame, raw data}
#'  \item{modelled}{data.frame, modelled data}
#'  \item{getCluster.res}{data.frame, clustering results from timeOmics}
#'  \item{getCluster.sparse.res}{data.frame, sparse clustering results from timeOmics}
#'  \item{interaction.biogrid}{data.frame, interactions from BioGRID database}
#'  \item{interaction.TF}{data.frame, TFome interactions from TTrust and TF2DNA}
#'  \item{medlineranker.res.df}{data.frame, medlineRanker enrichment results}
#'  \item{graph.gut}{list of igraph, gut graph obtained with SparCC}
#'  }
#' 
"hmp_T2D"