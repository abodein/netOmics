#' Random Walk with Restart
#' 
#' This function performs a propagation analysis by random walk with restart in a multi-layered network from specific seeds.
#' 
#' @param X an igraph or list.igraph object.
#' @param seed a character vector. Only seeds present in X are considered.
#' @param r a numeric value between 0 and 1. It sets the probability of restarting to a seed node after each step. 
#'  
#' @return 
#' Each element of X returns a list (class = 'rwr') containing the following elements:
#' \item{rwr}{a \code{data.frame}, the RWR results for each valid seed.}
#' \item{seed}{a character vector with the valid seeds}
#' \item{graph}{\code{igraph} object from X}
#' If X is a \code{list.igraph}, the returned object is a \code{list.rwr}.
#' 
#' @seealso 
#' \code{\link[RandomWalkRestartMH]{Random.Walk.Restart.Multiplex}}, \code{\link[netOmics]{rwr_find_seeds_between_attributes}}, \code{\link[netOmics]{rwr_find_closest_type}}
#' 
#' @examples
#' graph1 <- igraph::graph_from_data_frame(list(from = c("A", "B", "A", "D", "C", "A", "C"), 
#' to = c("B", "C", "D", "E", "D", "F", "G")), directed = FALSE)
#' graph1 <- igraph::set_vertex_attr(graph = graph1, name = 'type', index = c("A","B","C"),value = "1")
#' graph1 <- igraph::set_vertex_attr(graph = graph1, name = 'type', index = c("D","E"),value = "2")
#' graph1 <- igraph::set_vertex_attr(graph = graph1, name = 'type', index = c("F", "G"),value = "3")
#' 
#' rwr_res <- random_walk_restart(X = graph1, seed = c("A", "B", "C", "D", "E"))
#' 
#' @import RandomWalkRestartMH
#' @importFrom dplyr mutate left_join
#' @importFrom purrr imap_dfr
#' @importFrom magrittr %>%
#' @export
random_walk_restart <- function(X, seed = NULL, r = 0.7){
    
    # check X is graph or list of graph
    X <- check_graph(X)
    
    # check seed
    seed <- check_vector_char(X = seed, var.name = "'seed' ")
    
    # check r 
    r <- check_single_numeric_value(r, min = 0, max = 1, var.name = "'r' ")
    
    # delta
    delta <- 0.5
    
    res <- list()
    if(is(X, "list.igraph")){
        # apply RWR on each graph
        for(i in seq_along(X)){
            Xi <- X[[i]]
            Xi <- remove_unconnected_nodes(Xi)
            index_name_i <- ifelse(!is.null(names(X)[i]), names(X)[i], i)

            ## possible implementation to benchmark: 
            ## extract graph component and make couples with seeds and matching subgraph
            
            seed_xi <- intersect(seed, igraph::V(Xi)$name) # prevent the error: "Some of the seeds are not nodes of the network"

            # rwr layer names: to change if we include some day multiplex network
            layers_name <- ifelse(!is.null(names(X)[i]), names(X)[i], "graph")
            
            multiplex <- RandomWalkRestartMH::create.multiplex(L1 = Xi,Layers_Name=layers_name)
            adj_matrix <- RandomWalkRestartMH::compute.adjacency.matrix(x = multiplex, delta = delta)
            adj_matrix_norm <- RandomWalkRestartMH::normalize.multiplex.adjacency(x = adj_matrix) # time/RAM consuming
            
            res_tmp <- list()
            for(seed_xi_i in seed_xi){
                rwr_res <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(x = adj_matrix_norm, MultiplexObject = multiplex, Seeds = seed_xi_i, r = r)
                res_tmp[[seed_xi_i]] <- rwr_res
            } 
            if(!is_empty(seed_xi)){
                res[[index_name_i]] <- list()
                res[[index_name_i]][["rwr"]] <- purrr::imap_dfr(res_tmp, ~{.x$RWRM_Results %>% dplyr::mutate(SeedName = .y)}) %>%
                    dplyr::left_join(as.data.frame(vertex_attr(X[[i]])), by = c("NodeNames" = "name"))
                res[[index_name_i]][["graph"]] <- X[[i]]
                res[[index_name_i]][["seed"]] <- seed_xi
                class(res[[index_name_i]]) <- "rwr"
            }
            class(res) <- c("list.rwr")
        }
    } else { # X is a single graph
        Xi <- remove_unconnected_nodes(X)
        
        ## possible implementation to benchmark: 
        ## extract graph component and make couples with seeds and matching subgraph
        
        seed_xi <- intersect(seed, igraph::V(Xi)$name) # prevent the error: Some of the seeds are not nodes of the network    
        
        # rwr layer names: to change if we include some day multiplex network
        #layers_name <- ifelse(!is.null(names(X)[i]), names(X)[i], "graph")
        layers_name <- c("graph")
        
        multiplex <- RandomWalkRestartMH::create.multiplex(L1 = Xi,Layers_Name=layers_name)
        adj_matrix <- RandomWalkRestartMH::compute.adjacency.matrix(x = multiplex, delta = delta)
        adj_matrix_norm <- RandomWalkRestartMH::normalize.multiplex.adjacency(x = adj_matrix) # time/RAM consuming
        
        res_tmp <- list()
        for(seed_xi_i in seed_xi){
            rwr_res <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(x = adj_matrix_norm, MultiplexObject = multiplex, Seeds = seed_xi_i, r = r)
            res_tmp[[seed_xi_i]] <- rwr_res
        }
        # all seeds for a graph X has been computed -> merge result (more efficient than having seperate results + associated graph)
        if(!is_empty(seed_xi)){
            res[["rwr"]] <- purrr::imap_dfr(res_tmp, ~{.x$RWRM_Results %>% dplyr::mutate(SeedName = .y)}) %>%
                dplyr::left_join(as.data.frame(vertex_attr(X)), by = c("NodeNames" = "name"))
            res[["graph"]] <- X
            res[["seed"]] <- seed_xi
        }
        
        class(res) <- c("rwr")
    }
    return(res)
}

#' @importFrom igraph delete.vertices simplify degree
remove_unconnected_nodes <- function(X){
    # remove unconnected nodes but does not simplify
    X.simplified <- igraph::simplify(X)
    isolated_nodes = which(igraph::degree(X.simplified)==0)
    X = igraph::delete.vertices(X, isolated_nodes)
    return(X)
}

#' @importFrom dplyr filter pull top_n
#' @importFrom igraph induced_subgraph set_vertex_attr V
rwr_top_k_graph <- function(X, RWRM_Result_Object, Seed, k = 15){
    Top_Results_Nodes <- RWRM_Result_Object %>% dplyr::filter(SeedName == Seed) %>%
        dplyr::top_n(n = k, wt = Score) %>% dplyr::pull(NodeNames)
        #RWRM_Result_Object$RWRM_Results$NodeNames[seq_len(k)]
    #Query_Nodes <- c(RWRM_Result_Object$Seed_Nodes, Top_Results_Nodes)
    Query_Nodes <- intersect(c(Seed, Top_Results_Nodes), igraph::V(X)$name)
    Target_Nodes <- intersect(Top_Results_Nodes, igraph::V(X)$name)
    
    if(!purrr::is_empty(Query_Nodes)){
        top_k_graph <- igraph::induced_subgraph(graph = X, vids = Query_Nodes)
        top_k_graph <- igraph::set_vertex_attr(graph = top_k_graph, name = "rwr", index = Seed, value = "seed")
        top_k_graph <- igraph::set_vertex_attr(graph = top_k_graph, name = "rwr", index = Target_Nodes, value = "target")
        return(top_k_graph)
    }
    return(NULL)
}



#' RWR Find seeds between attributes
#' 
#' From rwr results, this function returns a subgraph if any vertex shares different attributes value.
#' In biological context, this might be useful to identify vertex shared between clusters or omics types.
#' 
#' @param X a random walk result from \code{random_walk_restart}
#' @param seed a character vector or NULL. If NULL, all the seeds from X are considered.
#' @param attribute a character value or NULL. If NULL, the closest node is returned.
#' @param k a integer, k closest nodes to consider in the search
#' 
#' @return 
#' A list of igraph object for each seed.
#' If X is a list, it returns a list of list of graph.
#' 
#' @examples 
#' graph1 <- igraph::graph_from_data_frame(list(from = c("A", "B", "A", "D", "C", "A", "C"), 
#'                                              to = c("B", "C", "D", "E", "D", "F", "G")), directed = FALSE)
#' graph1 <- igraph::set_vertex_attr(graph = graph1, name = 'type', index = c("A","B","C"),value = "1")
#' graph1 <- igraph::set_vertex_attr(graph = graph1, name = 'type', index = c("D","E"),value = "2")
#' graph1 <- igraph::set_vertex_attr(graph = graph1, name = 'type', index = c("F", "G"),value = "3")
#' 
#' rwr_res <- random_walk_restart(X = graph1, seed = c("A", "B", "C", "D", "E"))
#' rwr_res_type <- rwr_find_seeds_between_attributes(X = rwr_res, attribute = "type", k = 3)
#' 
#' @export
rwr_find_seeds_between_attributes <- function(X, seed = NULL, k = 15, attribute = "type"){
    # check X
    if(!(is(X, "rwr") | is(X, "list.rwr"))){
        stop("X must be a random walk result")
    }
    
    # check k
    if(!is.null(k)){
        k <-  check_single_numeric_value(k, min = 0, max = 200, var.name = "'k' ")

    } else {
        k <- 15
    }
    
    # check seed  # if seed is null, all seeds are considered found in rwr are considered
    if(!is.null(seed)){
        # don't check if all seeds are in vids -> NULL results anyway
        seed <- check_vector_char(X = seed, var.name = "'seed' ", default = NULL)
    } 
    
    # check attribute
    attribute <-  check_vector_char(X = attribute, var.name = "'attribute' ", default = "type", X.length = 1)
    
    if(is(X, "rwr")){
        if(is.null(seed)){ # seed = all seeds 
            seed <- X$seed  # can be NULL
        }
        res <- .rwr_find_seeds_between_attribute(rwr = X, k = k, attribute = attribute, seed = seed)
        class(res) <- "rwr.attributes"
    } else { # X is list.res
        # should not be run on list.res because each item contains a unique cluster
        res <- list()

        for(i in seq_along(X)){
            index_name_i <- ifelse(!is.null(names(X)[i]), names(X)[i], i)
            
            if(is.null(seed)){ # seed = all seeds 
                seed_i <- X[[index_name_i]]$seed  # can be NULL
            } else {
                seed_i <- seed
            }
            
            res[[index_name_i]] <- .rwr_find_seeds_between_attribute(rwr = X[[index_name_i]], k = k, attribute = attribute, seed = seed_i)
            class(res[[index_name_i]]) <- "rwr.attributes"
            
        }
        class(res) <- "list.rwr.attributes"
    }
    return(res)
}

#' @importFrom igraph vertex_attr
.rwr_find_seeds_between_attribute <- function(rwr, k, attribute, seed){
    res <- list()
    for(seed_xi in seed){
        # print(seed_xi)
        top_k_graph <- rwr_top_k_graph(X = rwr$graph, RWRM_Result_Object = rwr$rwr, Seed = seed_xi, k = k)
        
        # find different cluster
        #if(nrow(table(vertex_attr(top_k_graph)$cluster)) >= 2){ # at least 2 != clusters (NA excluded)
        if(!is.null(top_k_graph)){
            if(nrow(table(igraph::vertex_attr(top_k_graph)[[attribute]])) >= 2){  # generic version
                res[[seed_xi]] <- top_k_graph
            }
        }
        # att_seed_value <- igraph::vertex_attr(graph = top_k_graph, index = seed_xi, name = attribute)
        # att_seed_all <- igraph::vertex_attr(graph = top_k_graph, name = attribute) 
        # if(any(att_seed_all != att_seed_all)){
        #     res[[seed_xi]] <- top_k_graph
        # }
    }
    return(res)
}



#' RWR Find closest nodes
#' 
#' From a rwr results, this function returns the closest nodes from a seed with a given attribute and value.
#' In biological context, it might be useful to get the closest Gene Ontology annotation nodes from unannotated seeds.
#' 
#' @param X a random walk result from \code{random_walk_restart}
#' @param seed a character vector or NULL. If NULL, all the seeds from X are considered.
#' @param attribute a character value or NULL. If NULL, the closest node is returned.
#' @param value a character value or NULL. If NULL, the closest node for a given attribute is returned.
#' @param top a numeric value, the top closest nodes to extract
#' 
#' 
#' 
#' @return 
#' A list of \code{data.frame} for each seed containing the closest nodes per seed and their vertex attributes.
#' If X is \code{list.rwr}, the returned value is a list of list. 
#' 
#' 
#' @examples 
#' graph1 <- igraph::graph_from_data_frame(list(from = c("A", "B", "A", "D", "C", "A", "C"), 
#'                                              to = c("B", "C", "D", "E", "D", "F", "G")), directed = FALSE)
#' graph1 <- igraph::set_vertex_attr(graph = graph1, name = 'type', index = c("A","B","C"),value = "1")
#' graph1 <- igraph::set_vertex_attr(graph = graph1, name = 'type', index = c("D","E"),value = "2")
#' graph1 <- igraph::set_vertex_attr(graph = graph1, name = 'type', index = c("F", "G"),value = "3")
#' 
#' rwr_res <- random_walk_restart(X = graph1, seed = c("A", "B", "C", "D", "E"))
#' rwr_find_closest_type(X=rwr_res, attribute = "type", seed = "A")
#' 

#' @export
rwr_find_closest_type <- function(X, seed = NULL, attribute = NULL, value = NULL, top = 1){
    # check X
    if(!(is(X, "rwr") | is(X, "list.rwr"))){
        stop("X must be a random walk result")
    }
    
    # check attribute or replace with default value
    attribute <- check_vector_char(X = attribute, X.length = 1, default = NULL, var.name = "'attribute' ")
    
    # check value or replace with default value
    value <- check_vector_char(X = value, X.length = 1, default = NULL, var.name = "'value' ")
    
    # check top
    top <- check_single_numeric_value(top, var.name = "'top' ")
    
    # check seed  # if seed is null, all seeds are considered found in rwr are considered
    if(!is.null(seed)){
        # don't check if all seeds are in vids -> NULL results anyway
        seed <- check_vector_char(X = seed, var.name = "'seed' ", default = NULL)
    }
    
    if(is(X, "rwr")){
        if(is.null(seed)){ # seed = all seeds 
            seed <- X$seed  # can be NULL
        }
        res <- .rwr_find_closest(rwr = X, user.attribute = attribute, seed = seed, user.value = value, top = top)
        class(res) <- "rwr.closest"
    } else { # X is list.res
        # should not be run on list.res because each item contains a unique cluster
        res <- list()
        
        for(i in seq_along(X)){
            index_name_i <- ifelse(!is.null(names(X)[i]), names(X)[i], i)
            
            if(is.null(seed)){ # seed = all seeds 
                seed_i <- X[[index_name_i]]$seed  # can be NULL
            } else {
                seed_i <- seed
            }
            
            res[[index_name_i]] <- .rwr_find_closest(rwr = X[[index_name_i]], user.attribute = attribute, seed = seed_i, user.value = value, top = top)
            class(res[[index_name_i]]) <- "rwr.closest"
            
        }
        class(res) <- "list.rwr.closest"
    }
    return(res)
}

#' @importFrom dplyr filter top_n left_join select everything
#' @importFrom purrr map_dfr
#' @importFrom tidyr pivot_longer
.rwr_find_closest <- function(rwr, user.attribute, user.value, seed, top){
    res <- list()
    for(seed_xi in seed){
        rwr.res.filtered <- dplyr::filter(rwr$rwr, SeedName == seed_xi) 
        # fix to use pivot_longer with different cast columns (integer/logical, ...)
        rwr.res.filtered <- rwr.res.filtered %>% t %>% t %>% as.data.frame()
        rwr.res.filtered <- tidyr::pivot_longer(rwr.res.filtered, names_to = "attribute", values_to = "value", -c(NodeNames, Score, SeedName),
                                                values_ptypes = list(value=character()))
        if(!is.null(user.attribute)){
            rwr.res.filtered <- dplyr::filter(rwr.res.filtered, attribute == user.attribute)
        }
        if(!is.null(user.value)){
            rwr.res.filtered <- dplyr::filter(rwr.res.filtered, value == user.value)
        }
        rwr.res.filtered <- dplyr::top_n(x = rwr.res.filtered, n = top, wt = Score) %>%
            dplyr::select(c(NodeNames, SeedName)) %>% unique
        if(nrow(rwr.res.filtered) > 0){
            res[[seed_xi]] <- dplyr::left_join(rwr.res.filtered, rwr$rwr, by =  c("NodeNames" = "NodeNames", "SeedName" = "SeedName")) %>% 
                dplyr::select(c(NodeNames, Score, SeedName), dplyr::everything())
        }
    }
    #res <- purrr::map_dfr(res, ~.x)
    return(res)
}



