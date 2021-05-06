#' RWR
#' 
#' @examples
#' 
#' data(HeLa)
#' X <- HeLa$raw$mRNA
#' cluster.mRNA <- timeOmics::getCluster(HeLa$getCluster, user.block = "mRNA")
#' grn.res <- get_grn(X = HeLa$raw$mRNA, cluster = cluster.mRNA, method = "aracne")
#' 
#' X <- grn.res
#' Xi <- grn.res[["Up_Up"]]
#' seed <- igraph::V(Xi)$name[1:5]
#' rwr_res <- random_walk_restart(Xi, seed)
#' 
#' rwr_res <- random_walk_restart(X, seed)
#' 
#' random_walk_restart(X, seed = "A")
#' 
#' @import RandomWalkRestartMH
#' @export
random_walk_restart <- function(X, seed, k = 15, delta = 0.5, r = 0.7){
    
    # check X is graph or list of graph
    X <- check_graph(X)
    
    # check seed
    seed <- check_vector_char(X = seed, var.name = "'seed' ")
    
    # check K
    
    res <- list()
    if(is(X, "list.igraph")){
        # apply RWR on each graph
        for(i in seq_along(X)){
            Xi <- X[[i]]
            Xi <- remove_unconnected_nodes(Xi)
            index_name_i <- ifelse(!is.null(names(X)[i]), names(X)[i], i)
            res[[index_name_i]] <- list()
            
            ## possible implementation to benchmark: 
            ## extract graph component and make couples with seeds and matching subgraph
            
            seed_xi <- intersect(seed, igraph::V(Xi)$name) # prevent the error: Some of the seeds are not nodes of the network    

            # rwr layer names: to change if we include some day multiplex network
            layers_name <- ifelse(!is.null(names(X)[i]), names(X)[i], "graph")
            
            multiplex <- RandomWalkRestartMH::create.multiplex(L1 = Xi,Layers_Name=layers_name)
            adj_matrix <- RandomWalkRestartMH::compute.adjacency.matrix(x = multiplex, delta = delta)
            adj_matrix_norm <- RandomWalkRestartMH::normalize.multiplex.adjacency(x = adj_matrix) # time/RAM consuming
            
            for(seed_xi_i in seed_xi){
                rwr_res <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(x = adj_matrix_norm, MultiplexObject = multiplex, Seeds = seed_xi_i, r = r)
                top_k_graph <- RandomWalkRestartMH::create.multiplexNetwork.topResults(RWRM_Result_Object = rwr_res, MultiplexObject = multiplex, k = k)
                
                res[[index_name_i]][[seed_xi_i]] <- list()
                res[[index_name_i]][[seed_xi_i]][["rwr_res"]] <- rwr_res
                res[[index_name_i]][[seed_xi_i]][["top_k_graph"]] <- top_k_graph
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
        
        for(seed_xi_i in seed_xi){
            rwr_res <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(x = adj_matrix_norm, MultiplexObject = multiplex, Seeds = seed_xi_i, r = r)
            top_k_graph <- RandomWalkRestartMH::create.multiplexNetwork.topResults(RWRM_Result_Object = rwr_res, MultiplexObject = multiplex, k = k)
            
            res[[seed_xi_i]] <- list()
            res[[seed_xi_i]][["rwr_res"]] <- rwr_res
            res[[seed_xi_i]][["top_k_graph"]] <- top_k_graph
            
        }
        class(res) <- c("rwr")
    }
    return(res)
}


remove_unconnected_nodes <- function(X){
    # remove unconnected nodes but does not simplify
    X.simplified <- igraph::simplify(X)
    isolated_nodes = which(igraph::degree(X.simplified)==0)
    X = igraph::delete.vertices(X, isolated_nodes)
    return(X)
}

get_element_composition <- function(graph, seed){
    composition <- igraph::decompose(graph)
    
}



get_component <- function(graph, elt){
    composition <- igraph::decompose(graph)
    index <- lapply(composition, function(x) {elt %in% names(V(x))}) %>% unlist %>% which
    if(!is_null(index)){
        return(composition[[i]])
    }
    return()
}

extract_component <- function(graph, ids){
    log <- igraph::decompose(graph)
    index <- lapply(log, function(x) any(ids %in% V(x)$name)) %>% unlist()
    return(log[index])
}

homo_network_RWR <- function(graph, seed, AdjMatrixNorm = NULL, k = 15){
    seed <- seed[seed %in% V(graph)$name]
    network <- extract_component(graph, seed)[[1]]
    
    MultiplexObject <- RandomWalkRestartMH::create.multiplex(network, Layers_Name = c("component"))
    if(is.null(AdjMatrixNorm)){
        AdjMatrix <- compute.adjacency.matrix(MultiplexObject)
        AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)   
    }
    
    RWR_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm, MultiplexObject,seed)
    TopResults  <- create.multiplexNetwork.topResults(RWR_Results,MultiplexObject,k=k)
    return(TopResults)
}

homo_network_RWR
