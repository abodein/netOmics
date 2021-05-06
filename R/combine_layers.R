#' Combine layers
#'
#' Return merged graph from two graph layers
#'
#' @param graph1 an igraph object or list of igraph (list.igraph)
#' @param graph2 an igraph object or list of igraph (list.igraph) with the same length as graph1.
#' @param interaction.df

#' @return a graph
#'
#' @examples
#' data(BIOGRID)
#' db <- BIOGRID
#' X <- V(db)$name[1:10]
#' biogrid.res <- get_interaction_from_database(X, db, type = "db", user.ego = FALSE)




#' @importFrom purrr is_empty map reduce map2
#' @importFrom igraph induced_subgraph set_vertex_attr adjacent_vertices graph_from_data_frame
#' @export
combine_layers <- function(graph1, graph2 = NULL, interaction.df = NULL) {

    # check graph1
    if(!is(graph1, "igraph") & !is(graph1, "list.igraph")){
        stop("graph1 must be an igraph or list.igraph object")
    }

    if(!is(graph2, "igraph") & !is(graph2, "list.igraph") & !is.null(graph2)){
        stop("graph2 must be an igraph or list.igraph object or NULL")
    }
    if(!is.null(interaction.df)){
        interaction.df <- check_db(interaction.df)
        if(!is(interaction.df, "igraph")){
            interaction.df <- interaction.df %>% dplyr::select(c("from", "to"))
            interaction.graph <- igraph::graph_from_data_frame(interaction.df, directed = FALSE)
        }
    }

    # case1: graph2 = NULL, interaction.df = NULL
    if(is.null(graph2) & is.null(interaction.df)){
        merged.res <- graph1
    }

    # case2: graph1 and graph2 are single graph (+ interaction.df)
    if(is(graph1, "igraph") & is(graph2, "igraph")){
        merged.res <- merge_graphs(graph1, graph2)
        if(!is.null(interaction.df)){ # interaction.graph can be not found, df can be NULL
            merged.res <- merge_graphs(merged.res, interaction.graph)
        }

    # case3: graph1 is a list and graph2 is a single graph (+ interaction.df)
    } else if(is(graph1, "list.igraph") & is(graph2, "igraph")){
        merged.res <- purrr::map(graph1, ~{merge_graphs(.x, graph2)})
        if(!is.null(interaction.df)){ # interaction.graph can be not found, df can be NULL
            merged.res <- purrr::map(merged.res, ~{merge_graphs(.x, interaction.graph)})
        }

    # case4: graph1 and graph2 are list of graph (+ interaction.df)
    } else if(is(graph1, "list.igraph") & is(graph2, "list.igraph")){
        if(length(graph1) != length(graph2)){
            stop("graph1 and graph2 must have the same length")
        }
        if(!is.null(names(graph1)) & !is.null(names(graph2))){
            # graph1 and graph2 have names
            if(!all(names(graph1) %in% names(graph2))){ # same length so reciprocal is TRUE
                # they don't have the same names
                stop("graph1 and graph2 must have the same names")
            } else {
                merged.res <- purrr::map2(graph1, graph2[names(graph1)], ~{merge_graphs(.x, .y)})
            }
        } else {
            # no names, don't care about the order
            merged.res <- purrr::map2(graph1, graph2, ~{merge_graphs(.x, .y)})
            
        }
        if(!is.null(interaction.df)){ # interaction.graph can be not found, df can be NULL
            merged.res <- purrr::map(merged.res, ~{merge_graphs(.x, interaction.graph)})
        }
        
    # case5: inverse of case3 -> error
    } else if(is(graph1, "igraph") & is(graph2, "list.igraph")){
        stop("graph1 and graph2 must have the same length")
        
    # case6: graph1 and interaction.df    
    } else if(is(graph1, "igraph") & is.null(graph2) & !is.null(interaction.df)){
        merged.res <- merge_graphs(graph1, interaction.graph)
        
    # case7: graph1 list and interaction.df
    } else if(is(graph1, "list.igraph") & is.null(graph2) & !is.null(interaction.df)){
        merged.res <- purrr::map(merged.res, ~{merge_graphs(.x, interaction.graph)})
    } 
    
    if(is(merged.res, "list")){
        class(merged.res) <- c("list.igraph", "list.merged.igraph")
    }
    return(merged.res)
}


# graph1 <- erdos.renyi.game(10, 0.3)
# graph1 <- set_vertex_attr(graph1, "name", value = letters[1:10])
# graph1 <- set_vertex_attr(graph1, "type", value = "shared_g1")
# graph1 <- set_vertex_attr(graph1, "uniq_g1", value = "uniq_g1")
#
#
# graph2 <- erdos.renyi.game(5, 0.3)
# graph2 <- set_vertex_attr(graph2, "name", value = letters[c(2,4,6,8,10)])
# graph2 <- set_vertex_attr(graph2, "type", value = "shared_g2")
# graph2 <- set_vertex_attr(graph2, "uniq_g2", value = "uniq_g2")
#
#
#
# # not exported
# # graph.merged <- merge_graphs(graph1, graph2)
# # graph.merged <- merge_graphs(graph2, graph1)

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
    class(merged_graphs) <- c("merged.igraph", "igraph")
    return(merged_graphs)
}
