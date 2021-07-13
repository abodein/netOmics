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
#' 
#' X <- letters[1:5]
#' db <- as.data.frame(list(from = c(letters[1:5], letters[11:15]), to = c(letters[1:10])))

#' @importFrom purrr is_empty map reduce
#' @importFrom igraph induced_subgraph set_vertex_attr adjacent_vertices
#' @export
get_interaction_from_database <- function(X, db = NULL, type = "db", user.ego = FALSE) {
    
    # check X
    if(is(X, "list")){
        X <- lapply(X,function(x)check_vector_char(x))
        if(is.null(names(X))){
            names(X) <- seq_along(X)
        }
    } else {
        X <- check_vector_char(X)
    }
    
    # check db
    db <- check_db(db, var.name = "'db' ")
    
    # check type
    type <- check_vector_char(type, X.length = 1, default = "db")
    
    # check user.ego
    user.ego <- return_true_false(user.ego, default = FALSE)
    

    if(is.null(X)){
        message("X is NULL, returning an empty graph")
        db.subgraph <- igraph::make_empty_graph(directed = FALSE)
        class(db.subgraph) <- c("interaction.igraph", "igraph")
        return(db.subgraph)
    }
    else if(is(X, "list")){
    # filter db from X
        db.subgraph.list <- list()
        if(is(db, "igraph")){
            for(i in names(X)){
                db.subgraph.list[[i]] <-  .interaction_from_igraph(X = X[[i]], db = db, ego = user.ego, type = type)
            }
        } else { # db is a data.frame
            for(i in names(X)){
                db.subgraph.list[[i]] <-  .interaction_from_dataframe(X = X[[i]], db = db, ego = user.ego, type = type)
            }
        }
        class(db.subgraph.list) <- "list.igraph"
        return(db.subgraph.list)
    } else { 
        # X is a single vector
        if(is(db, "igraph")){
            db.subgraph <-  .interaction_from_igraph(X = X, db = db, ego = user.ego, type = type)
        } else { # db is a data.frame
            db.subgraph <-  .interaction_from_dataframe(X = X, db = db, ego = user.ego, type = type)
        }
        return(db.subgraph)
    }
}

.interaction_from_igraph <- function(X, db, ego, type){
    node.names <- intersect(X, igraph::V(db)$name)
    if(purrr::is_empty(node.names)){
        message("no shared elements between X and db, return empty graph")
        db.subgraph <- igraph::make_empty_graph(directed = FALSE)
    }
    else if(isTRUE(ego)){
        ego.neighbors <- igraph::adjacent_vertices(graph = db, v = node.names, mode = "all")
        ego.neighbors <- unique(purrr::reduce(purrr::map(ego.neighbors, ~names(.x)), union))
        ego.neighbors <- setdiff(ego.neighbors, node.names)
        
        db.subgraph <- igraph::induced_subgraph(graph = db, vids = c(node.names, ego.neighbors))
        db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="mode", index = node.names, value = "core")
        db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="mode", index = ego.neighbors, value = "extended")
    } else { # ego = FALSE
        db.subgraph <- igraph::induced_subgraph(graph = db, vids = c(node.names))
        db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="mode", index = node.names, value = "core")
    }
    # return graph
    db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="type", value = type)
    class(db.subgraph) <- c("interaction.igraph", "igraph")
    return(db.subgraph)
}

.interaction_from_dataframe <- function(X, db, ego, type){
    db <- as.data.frame(db) %>% dplyr::select(c("from", "to")) # checked colnames
    db.all.nodes <- unique(c(db$from, db$to))
    node.names <- intersect(X, db.all.nodes)
    if(purrr::is_empty(node.names)){
        message("no shared elements between X and db, return empty graph")
        db.subgraph <- igraph::make_empty_graph(directed = FALSE)
    }
    else if(isTRUE(ego)){
        ego.db <- db %>% dplyr::filter(.$from %in% node.names | .$to %in% node.names)
        #ego.neighbors <- setdiff(db.all.nodes, node.names)
        ego.neighbors <- setdiff(unique(c(ego.db$from, ego.db$to)), node.names)
        
        db.subgraph <- igraph::graph_from_data_frame(ego.db, directed = FALSE)
        db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="mode", index = node.names, value = "core")
        db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="mode", index = ego.neighbors, value = "extended")
    } else {  # ego = FALSE
        ego.db <- db %>% dplyr::filter(.$from %in% node.names & .$to %in% node.names)

        db.subgraph <- igraph::graph_from_data_frame(ego.db, directed = FALSE)
        db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="mode", value = "core")
    }
    
    # return graph
    db.subgraph <- igraph::set_vertex_attr(graph = db.subgraph, name="type", value = type)
    class(db.subgraph) <- c("interaction.igraph", "igraph")
    return(db.subgraph)
}