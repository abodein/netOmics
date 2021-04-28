#' Combine layers
#' 
#' Return merged graph from two graph layers
#' 
#' @param graph1 
#' @param graph2
#' @param interaction.df

#' @return a graph
#' 
#' @examples 
#' data(BIOGRID)
#' db <- BIOGRID
#' X <- V(db)$name[1:10]
#' biogrid.res <- get_interaction_from_database(X, db, type = "db", user.ego = FALSE)
#' 
#' graph1 <- erdos.renyi.game(10, 0.3)
#' graph1 <- set_vertex_attr(graph1, "name", value = letters[1:10])
#' graph1 <- set_vertex_attr(graph1, "type", value = "shared_g1")
#' graph1 <- set_vertex_attr(graph1, "uniq_g1", value = "uniq_g1")
#' 
#' 
#' graph2 <- erdos.renyi.game(5, 0.3)
#' graph2 <- set_vertex_attr(graph2, "name", value = letters[c(2,4,6,8,10)])
#' graph2 <- set_vertex_attr(graph2, "type", value = "shared_g2")
#' graph2 <- set_vertex_attr(graph2, "uniq_g2", value = "uniq_g2")
#' 
#' 
#' 
#' # not exported 
#' # graph.merged <- merge_graphs(graph1, graph2)
#' # graph.merged <- merge_graphs(graph2, graph1)


#' @importFrom purrr is_empty map reduce
#' @importFrom igraph induced_subgraph set_vertex_attr adjacent_vertices
#' @export
combine_layers <- function(graph1, graph2, interaction.df = NULL) {
    
    # check graph
    if(!is(graph1, "igraph") && !is(graph2, "igraph")){
        stop("graph1 and graph2 must be igraph objects")
    }
    
    # check db
    interaction.df <- check_db(interaction.df)
    
    # check if type.graph1 in graph1 types
    
}



#' @importFrom igraph vertex_attr union delete_vertex_attr set_vertex_attr
merge_graphs <- function(graph1, graph2){
    # shared attr except 'name'
    shared_attr <- intersect(names(igraph::vertex_attr(graph1)), names(igraph::vertex_attr(graph2)))
    shared_attr <- shared_attr[!(shared_attr == "name")]
    # uniq_1 <- setdiff(names(igraph::vertex_attr(graph1)), names(igraph::vertex_attr(graph2)))
    # uniq_2 <- setdiff(names(igraph::vertex_attr(graph2)), names(igraph::vertex_attr(graph1)))
    
    merged_graphs <- igraph::union(graph1, graph2)
    #vertex_attr(merged_graphs) %>% as.data.frame()
    merged_attr <- igraph::vertex_attr(merged_graphs)
    for(sa in shared_attr){
        merged_attr[[sa]] <- vector(length = vcount(merged_graphs))
        for(i in 1:vcount(merged_graphs)){
            # if !is.na _1, return _1 else return _2
            merged_attr[[sa]][i] <- ifelse(!is.na(merged_attr[[paste0(sa, "_1")]][i]), 
                                           merged_attr[[paste0(sa, "_1")]][i],
                                           merged_attr[[paste0(sa, "_2")]][i])
        }
        merged_graphs <- delete_vertex_attr(graph = merged_graphs, name = paste0(sa, "_1"))
        merged_graphs <- delete_vertex_attr(graph = merged_graphs, name = paste0(sa, "_2"))
        merged_graphs <- set_vertex_attr(graph = merged_graphs, name = sa, value = merged_attr[[sa]])
    }
    return(merged_graphs)
}
