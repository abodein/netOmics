#' Add cluster label to graph nodes
#' 
#' plpl
#' 
#' 
#' @examples 
#' cluster <- timeOmics::getCluster(HeLa$getCluster, user.block = "mRNA")
#' X <- HeLa$raw$mRNA
#' 
#' @export

graph_add_cluster <- function(X, cluster){
    
    # check if X
    if(!is(X, "igraph") & !is(X, "list.igraph")){
        stop("X must be an igraph or list.igraph object")
    }
    
    # check cluster
    cluster <- check_getCluster(cluster)
    
    if(is.null(cluster)){
        
    }else{
        if(is(X, "igrapĥ")){
            
        } else { # X is list.igraph
            for(Xi in X){
                vids <- V(Xi)$name
                vids_intersect <- intersect(vids, cluster$molecule)
                cluster.sub <-  dplyr::filter(cluster, molecule %in% vids_intersect)
                stopifnot(nrow(cluster.sub) > length(vids_intersect))
            }
        }
    }
}
