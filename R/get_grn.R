#' Get Gene Regulatory Network from dataset
#' 
#' A more detailed description
#' 
#' @param X gene expression matrix/data.frame
#' @param cluster optional, clustering result from timeOmics::getCluster
#' @param method network building method, one of c('aracne')
#' 
#' @examples
#' data(HeLa)
#' # grn only on gene
#' cluster.mRNA <- timeOmics::getCluster(HeLa$getCluster, user.block = "mRNA")
#' X <- HeLa$raw$mRNA
#' grn.res <- get_grn(X = HeLa$raw$mRNA, cluster = cluster.mRNA, method = "aracne")
#' 
#' @importFrom minet build.mim aracne
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom purrr map
#' @importFrom igraph set_vertex_attr graph_from_adjacency_matrix
#' @export
get_grn <- function(X, cluster=NULL, method = c("aracne")){
    
    # check if X
    X <- validate_matrix_X(X)
    
    # check cluster 
    cluster <- check_getCluster(cluster)
    
    # check method, for now only 1
    method <- match.arg(method)
    
    if(is.null(cluster)){ # no clusteing info -> perform grn on all molecules
        mim <- minet::build.mim(X)
        grn.adj <- minet::aracne(mim)
        grn.graph <- igraph::graph_from_adjacency_matrix(grn.adj, mode = "undirected")
        
        # add type attribute "type" <- "Gene"
        grn.graph <- igraph::set_vertex_attr(graph = grn.graph, name = "type", value = "Gene")
        #res <- list()
        return(grn.graph)
    } else { # cluster != NULL
        # we do have cluster info and data are clustered
        # 1. grn for all
        mim <- minet::build.mim(X)
        grn.adj <- minet::aracne(mim)
        grn.graph <- igraph::graph_from_adjacency_matrix(grn.adj, mode = "undirected")
        grn.graph <- igraph::set_vertex_attr(graph = grn.graph, name = "type", value = "Gene")
        res <- list()
        res[["All"]] <- grn.graph
        
        #2. grn by all clusters
        mol_cluster <- cluster %>% split(.$cluster) %>% purrr::map(~.x$molecule)
        X.by.cluster <- purrr::map(mol_cluster, ~{dplyr::select(as.data.frame(X, check.names = FALSE), .x)})
        
        for(i in names(mol_cluster)){
            mim.cluster <- minet::build.mim(X.by.cluster[[i]])
            grn.adj.cluster <- minet::aracne(mim.cluster)
            grn.graph.cluster <- igraph::graph_from_adjacency_matrix(grn.adj.cluster, mode = "undirected")
            grn.graph.cluster <- igraph::set_vertex_attr(graph = grn.graph.cluster, name = "type", value = "Gene")
            res[[i]] <- grn.graph.cluster

        }
        class(res) <- c("list.igraph")
    }
    return(res)
}

