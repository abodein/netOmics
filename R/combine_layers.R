#' Combine layers
#' 
#' Return merged graph from two graph layers
#' 
#' @param graph1 
#' @param graph2
#' @param type.graph1
#' @param type.graph2

#' @return a graph
#' 
#' @examples 
#' data(BIOGRID)
#' db <- BIOGRID
#' X <- V(db)$name[1:10]
#' biogrid.res <- get_interaction_from_database(X, db, type = "db", user.ego = FALSE)
#' graph1 <- make_ring(10)
#' graph2 <- make_ring(5); graph2 <- set_vertex_attr(graph2, "name", value = letters[1:5])

#' @importFrom purrr is_empty map reduce
#' @importFrom igraph induced_subgraph set_vertex_attr adjacent_vertices
#' @export
combine_layers <- function(graph1, graph2, type.graph1 = NULL, type.graph2 = NULL) {
    
    # check graph
    if(!is(graph1, "igraph") && !is(graph2, "igraph")){
        stop("graph1 and graph2 must be igraph objects")
    }
    
    # check type
    type.graph1 <- check_vector_char(X = type.graph1, X.length = 1)
    type.graph2 <- check_vector_char(X = type.graph2, X.length = 1)
    
    # check if type.graph1 in graph1 types
    if(!is.null(type.graph1)){
        igraph::vertex_attr(graph1, type)
    }
}


#' graph1 <- make_ring(10)
#' graph2 <- make_ring(5); graph2 <- set_vertex_attr(graph2, "name", value = letters[1:5])
#' 
merge_graphs <- function(graph1, graph2){
    # check if 
}
