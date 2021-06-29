get_correlation_incidence <- function(X, Y, threshold = 0.5){
    res <- list()
    for(i in names(Y)){
        Yi <- Y[[i]]
        res.corr <- cor(x = X, y = Yi, method = "spearman")
        corr.graph <- abs(res.corr) >= threshold
        # res.graph <- graph_from_adjacency_matrix(corr.graph, mode = "undirected") %>% simplify
        res.graph <- graph_from_incidence_matrix(corr.graph, directed = FALSE) %>% simplify
        res[[i]] <-  as.data.frame(get.edgelist(res.graph)) %>% set_names(c("from", "to"))
    }
    return(res)
}

get_correlation_incidence_v2 <- function(X, Y, threshold = 0.5){
    # X and Y are list of the same size with the same names
    res <- list()
    for(i in names(X)){
        Xi <- X[[i]]
        Yi <- Y[[i]]
        res.corr <- cor(x = Xi, y = Yi, method = "spearman")
        corr.graph <- abs(res.corr) >= threshold
        # res.graph <- graph_from_adjacency_matrix(corr.graph, mode = "undirected") %>% simplify
        res.graph <- graph_from_incidence_matrix(corr.graph, directed = FALSE) %>% simplify
        res[[i]] <-  as.data.frame(get.edgelist(res.graph)) %>% set_names(c("from", "to"))
    }
    return(res)
}



