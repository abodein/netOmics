#' Gene Regulatory Network
#' 
#' Get Gene Regulatory Network (GRN) from a data.frame.
#' Optionally, if the gene are clustered, sub_network are build for 
#' each cluster.
#' 
#' @param X a \code{data.frame}/\code{matrix} with gene expression 
#' (genes in columns, samples in rows).
#' @param cluster (optional) clustering result from 
#'      \code{\link[timeOmics]{getCluster}}
#' @param method network building method, one of c('aracne')
#' @param type character added to node metadata
#' 
#' 
#' @details 
#' Methods of GRN reconstruction are as follows:
#' 'aracne': use ARACNe algorithm on Mutual Information (MI) adjency matrix 
#' to remove low MI edges in triangles.
#' 
#' @return 
#' An igraph object if no cluster informations are given. 
#' Otherwise, it returns a list of igraph object (\code{list.igraph}) with 
#' a subgraph for each cluster and a global graph with all the genes.
#' 
#' @seealso 
#' \code{\link[minet]{build.mim}}, 
#' \code{\link[minet]{aracne}}, 
#' \code{\link[timeOmics]{getCluster}}
#' 
#' @examples
#' data(hmp_T2D)
#' # grn only on gene
#' cluster.mRNA <- timeOmics::getCluster(hmp_T2D$getCluster.res, 
#'                                       user.block = 'RNA')
#' X <- hmp_T2D$raw$RNA
#' grn.res <- get_grn(X = hmp_T2D$raw$RNA, 
#'                    cluster = cluster.mRNA, 
#'                    method = 'aracne')
#' 
#' 
#' @importFrom minet build.mim aracne
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom purrr map
#' @importFrom igraph set_vertex_attr graph_from_adjacency_matrix as.undirected
#' @export
get_grn <- function(X, 
                    cluster = NULL, 
                    method = c("aracne"),
                    type = "gene"
) {
    
    # check if X
    X <- validate_matrix_X(X, var.name = "'X' ")
    
    # check cluster
    cluster <- check_getCluster(cluster)
    
    # check method, for now only 1
    method <- match.arg(method)
    
    # check type
    type <- check_vector_char(type, X.length = 1, default = "gene")
    
    if (is.null(cluster)) {
        # no clusteing info -> perform grn on all molecules
        mim <- minet::build.mim(X)
        grn.adj <- minet::aracne(mim)
        grn.graph <- igraph::graph_from_adjacency_matrix(grn.adj) %>% 
          igraph::as.undirected()
        
        # add type attribute 'type' <- 'Gene'
        grn.graph <- igraph::set_vertex_attr(graph = grn.graph, 
                                             name = "type", 
                                             value = type)
        grn.graph <- igraph::set_vertex_attr(graph = grn.graph, 
                                             name = "mode", 
                                             value = "core")
        grn.graph <- igraph::set_vertex_attr(graph = grn.graph, 
                                             name = "cluster", 
                                             value = "All")
        
        # res <- list()
        return(grn.graph)
    } else {
        # cluster != NULL we do have cluster info and data are clustered 1. grn
        # for all
        mim <- minet::build.mim(X)
        grn.adj <- minet::aracne(mim)
        grn.graph <- igraph::graph_from_adjacency_matrix(grn.adj) %>% 
          igraph::as.undirected()
                                                         
        grn.graph <- igraph::set_vertex_attr(graph = grn.graph, 
                                             name = "type", 
                                             value = type)
        grn.graph <- igraph::set_vertex_attr(graph = grn.graph, 
                                             name = "mode", 
                                             value = "core")
        
        res <- list()
        res[["All"]] <- grn.graph
        
        # 2. grn by all clusters
        mol_cluster <- cluster %>%
            split(.$cluster) %>%
            purrr::map(~.x$molecule)
        X.by.cluster <- purrr::map(
            mol_cluster, ~{
                dplyr::select(
                    as.data.frame(X, check.names = FALSE),
                    .x
                )
            }
        )
        
        # names_mol_cluster <- check_name_list(mol_cluster)
        for (i in names(mol_cluster)) {
            mim.cluster <- minet::build.mim(X.by.cluster[[i]])
            grn.adj.cluster <- minet::aracne(mim.cluster)
            grn.graph.cluster <- igraph::graph_from_adjacency_matrix(
                grn.adj.cluster) %>% 
              igraph::as.undirected()
            grn.graph.cluster <- igraph::set_vertex_attr(
                graph = grn.graph.cluster, 
                name = "type", 
                value = type)
            grn.graph.cluster <- igraph::set_vertex_attr(
                graph = grn.graph.cluster, 
                name = "mode", 
                value = "core")
            grn.graph.cluster <- igraph::set_vertex_attr(
                graph = grn.graph.cluster, 
                name = "cluster", 
                value = i)
            
            res[[i]] <- grn.graph.cluster
            
            # also add cluster info to 'All' graph
            res[["All"]] <- igraph::set_vertex_attr(
                graph = res[["All"]], name = "cluster", 
                value = i, 
                index = igraph::V(grn.graph.cluster)$name
            )
            
        }
        class(res) <- c("list.igraph")
    }
    return(res)
}


