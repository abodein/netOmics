#' Get interaction from database
#' 
#' Returns an interaction graph from a vector of nodes (or a list of vectors) and an interaction database (data.frame or igraph)
#' 
#' @param X vector of nodes or list of vectors
#' @param db data.frame (with two columns: from, to) or igraph 
#' @param type character added to node metadata
#' @param user.ego logical, if user.ego == TRUE looks for first degree neighbors in db and add 'mode' metadata ('core'/'extended')
#' 
#' @return a subset graph of db from X list of nodes
#' 
#' @examples 
#' data(BIOGRID)
#' db <- BIOGRID
#' X <- V(db)$name[1:10]
#' biogrid.res <- get_interaction_from_database(X, db, type = "db", user.ego = FALSE)

#' @importFrom purrr is_empty map reduce
#' @importFrom igraph induced_subgraph set_vertex_attr adjacent_vertices
#' @export
get_interaction_from_database <- function(X, db, type = "db", user.ego = FALSE) {
    # check db
    db <- check_db(db)
    
    # check X
    X <- check_vector_char(X)
    
    # check type
    type <- check_vector_char(type, X.length = 1, default = "db")
    
    # check user.ego
    user.ego <- return_true_false(user.ego, default = FALSE)
    
    # filter db from X
    if(is(db, "igraph")){
        node.names <- intersect(X, igraph::V(db)$name)
        if(purrr::is_empty(node.names)){
            message("no shared elements between X and db, return empty graph")
        }
        if(isTRUE(ego)){
            ego.neighbors <- igraph::adjacent_vertices(graph = db, v = X, mode = "all")
            ego.neighbors <- unique(purrr::reduce(purrr::map(ego.neighbors, ~names(.x)), union))
            ego.neighbors <- setdiff(ego.neighbors, node.names)
            
            db.subgraph <- igraph::induced_subgraph(graph = db, vids = c(node.names, ego.neighbors))
            db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="mode", index = node.names, value = "core")
            db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="mode", index = ego.neighbors, value = "extended")
        } else { # ego = FALSE
            db.subgraph <- igraph::induced_subgraph(graph = db, vids = c(node.names))
            db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="mode", index = node.names, value = "core")
        }
    }
    
    # return graph
    db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="type", value = type)
    class(db.subgraph) <- c("interaction.igraph", "igraph")
    return(db.subgraph)
}
