#' get graph statistics
#' 
#' 
#' 
#' @examples
#' data(HeLa)
#' # grn only on gene
#' cluster.mRNA <- timeOmics::getCluster(HeLa$getCluster, user.block = "mRNA")
#' X <- HeLa$raw$mRNA
#' grn.res <- get_grn(X = HeLa$raw$mRNA, cluster = cluster.mRNA, method = "aracne")
#' get_graph_stats(grn.res)
#' 
#' 
#' @importFrom minet build.mim aracne
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom purrr imap_dfr
#' @importFrom igraph set_vertex_attr graph_from_adjacency_matrix
#' @export
get_graph_stats <- function(X){
    
    X <- check_graph(X)
    
    if(is(X, "list") | is(X, "list.igraph")){
        as.data.frame(lapply(X, function(x) .get_graph_stats_graph(x)))
        stats <- as.data.frame(purrr::imap_dfr(X, ~.get_graph_stats_graph(.x)))
        rownames(stats) <- names(X)

    } else { # X is not a list
        stats <- .get_graph_stats_graph(X)
    }
    class(stats) <- c("stats", "data.frame")
    return(stats)
}

#' @importFrom igraph degree ecount edge_density
.get_graph_stats_graph <- function(X){
    node.c <- sum(igraph::degree(X) != 0)
    node.i <- sum(igraph::degree(X) == 0)
    edge <- igraph::ecount(X)
    edge.density <- igraph::edge_density(graph = X)
    res <- list(node.c = node.c, node.i = node.i, edge = edge, edge.density = edge.density)
    
    if(any(names(vertex_attr(X)) != "name")){
        item <- names(vertex_attr(X))[names(vertex_attr(X)) != "name"]
        vertex.attr.res <- as.data.frame(vertex_attr(X)) %>% dplyr::select(item)
        for(i in item){
            tmp <- as.list(table(vertex.attr.res[[i]]))
            for(j in names(tmp)){
                name_item <- paste(i, j, sep = ".")
                res[[name_item]] <- tmp[[j]]
            }
        }
    }
    return(as.data.frame(res))
}

