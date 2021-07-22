#' hmp_T2D
#' 
#' A list of data.frames 
#' 
#' @format a list of data.frame
#' \describe{
#'   \item{raw}{data.frame, raw data}
#'   \item{modelled}{data.frame, modelled data}
#'   \item{getCluster.res}{data.frame, clustering results from timeOmics}
#'   \item{getCluster.sparse.res}{data.frame, sparse clustering results from timeOmics}
#'   \item{interaction.biogrid}{data.frame, interactions from BioGRID database}
#'   \item{interaction.TF}{data.frame, TFome interactions from TTrust and TF2DNA}
#'   \item{medlineranker.res.df}{data.frame, medlineRanker enrichment results}
#'   \item{graph.gut}{list of igraph, gut graph obtained with SparCC}
#'  }
#' 
"hmp_T2D"