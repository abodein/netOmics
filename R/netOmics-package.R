#' netOmics: network-based multi-omics integration and interpretation
#'
#' netOmics is a multi-omics networks builder and explorer.
#' It uses a combination of network inference algorithms and and knowledge-based
#'  graphs to build multi-layered networks. 
#'  
#'  The package can be combined with 
#'  `timeOmics` to incorporate time-course expression data and build 
#'  sub-networks from multi-omics kinetic clusters.
#'  
#' Finally, from the generated multi-omics networks, propagation analyses allow
#'  the identification of missing biological functions (1), 
#' multi-omics mechanisms (2) and molecules between kinetic clusters (3). 
#' This helps to resolve complex regulatory mechanisms.
#' Here are the main functions.
#' 
#' @section Network building:

#' \describe{
#'   \item{`get_grn`}{Based on expression matrix, this function build a gene 
#'   gene regulatory network. Additionally, if clustering information is given,
#'   it builds cluster specific graph.}
#'   \item{`get_interaction_from_database`}{From a database (graph or data.frame 
#'   with interactions between 2 molecules), this function build the induced
#'   graph based on a list of molecules . Alternatively, the function can 
#'   build a graph with the first degree neighbors.}
#'   \item{`get_interaction_from_correlation`}{Compute correlation between two 
#'   dataframe X and Y (or list of data.frame).
#' An incidence graph is returned. A link between two features is produced 
#' if their correlation (absolute value) is above the threshold.}
#'   \item{`combine_layers`}{Combine 2 (or list of) graphs based on given 
#'   intersections.}
#'   }
#'
#' @section Network exploration:
#' \describe{
#'   \item{`random_walk_restart`}{This function performs a propagation analysis 
#'   by random walk with restart
#'  in a multi-layered network from specific seeds.}
#'   \item{`rwr_find_seeds_between_attributes`}{From rwr results, this function 
#'   returns a subgraph if any vertex shares 
#' different attributes value.
#' In biological context, this might be useful to identify vertex shared between
#'  clusters or omics types.}
#'   \item{`rwr_find_closest_type`}{From a rwr results, this function returns 
#'   the closest nodes from a seed with 
#' a given attribute and value.
#' In biological context, it might be useful to get the closest Gene Ontology
#'  annotation nodes from unannotated seeds.}
#'  }
#'  
#' @section Visualisation:
#' \describe{
#'   \item{`summary_plot_rwr_attributes`}{#' Based on the results of 
#' \code{\link[netOmics]{rwr_find_seeds_between_attributes}} which identify the
#'  closest k neighbors from a seed, this function returns a barplot of the node
#'   types (layers) reached for each seed.}
#'   \item{`plot_rwr_subnetwork`}{Display the subgraph from a RWR results. 
#'   This function colors adds a specific
#'  color to each node based on their 'type' attribute.
#' It also adds a legend including the number of vertices/edges and the number 
#' of nodes of specific type.
#' Additionally, the function can display any igraph object.}
#'  }
#' 
#'
#' @docType package
#' @name netOmics
#' 
NULL
#> NULL