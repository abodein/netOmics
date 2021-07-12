#' Summary Plot RWR attributes
#' 
#' Based on the results of \code{\link[netOmics]{rwr_find_seeds_between_attributes}} which identify the closest k neighbors from a seed, 
#' this function returns a barplot of the node types (layers) reached for each seed.
#' 
#' @param X a 'rwr.attributes' or 'list.rwr.attributes' object from rwr_find_seeds_between_attributes()
#' @param color (optional) a named character vector or list, list of color to apply to each type
#' @param seed.id (optional) a character vector, to filter the results and filter on specific seeds IDs
#' @param seed.type (optional) a character vector, to filter the results and filter on specific seeds types
#' @param plot logical, if TRUE then the plot is produced
#' 
#' @return 
#' a 'ggplot' object
#' 
#' @seealso \code{\link[netOmics]{random_walk_restart}}, \code{\link[netOmics]{rwr_find_seeds_between_attributes}}
#' 
#' @examples 
#' graph1 <- igraph::graph_from_data_frame(list(from = c("A", "B", "A", "D", "C", "A", "C"), 
#'                                              to = c("B", "C", "D", "E", "D", "F", "G")), directed = FALSE)
#' graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("A","B","C"),value = "1")
#' graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("D","E"),value = "2")
#' graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("F", "G"),value = "3")
#' 
#' rwr_res <- random_walk_restart(X = graph1, seed = c("A", "B", "C", "D", "E"))
#' rwr_res_type <- rwr_find_seeds_between_attributes(X = rwr_res, attribute = "type", k = 3)
#' summary_plot_rwr_attributes(rwr_res_type)
#' 
#' 
#' @importFrom tibble rownames_to_column
#' @import ggplot2
#' @importFrom purrr imap_dfr set_names
#' @importFrom igraph vertex_attr
#' @importFrom dplyr filter mutate left_join group_by select summarise n
#' @export
summary_plot_rwr_attributes <- function(X, color = NULL, seed.id = NULL, seed.type = NULL, plot = TRUE){
    stopifnot(is(X, "rwr.attributes") | is(X, "list.rwr.attributes"))
    
    # check seed.id
    seed.id <- check_vector_char(X = seed.id, default = NULL, var.name = "'seed.id' ")
    
    # check seed.type
    seed.type <- check_vector_char(X = seed.type, default = NULL, var.name = "'seed.type' ")
    
    # check color 
    if(!is.null(color)){
        color <- check_named_vector(X = color, var.name = "'color' ")
    }
    # check plot
    plot <- return_true_false(x = plot, default = TRUE)
    
    
    if(is(X, "rwr.attributes")){
        # seed type 
        seed_types <- purrr::imap_dfr(X, ~{vertex_attr(.x) %>% 
                as.data.frame() %>% dplyr::filter(rwr == "seed") %>% 
                dplyr::select(name, type) %>% purrr::set_names(c("name", "seed.type"))})
        # count layer
        va.all <- purrr::imap_dfr(X, ~{igraph::vertex_attr(.x) %>% as.data.frame() %>% 
                dplyr::mutate(seed = .y) %>% 
                dplyr::group_by(seed, type) %>% 
                dplyr::summarise(N = dplyr::n(), .groups = "keep")}) %>% 
            dplyr::left_join(seed_types, by = c("seed"="name"))
    } else { #X is list.rwr.attributes
        seed_types <- lapply(names(X), function(Y){
            purrr::imap_dfr(X[[Y]], ~{igraph::vertex_attr(.x) %>% 
                as.data.frame() %>% dplyr::filter(rwr == "seed") %>% 
                dplyr::select(name, type) %>% purrr::set_names(c("name", "seed.type"))}) %>% 
                dplyr::mutate(sub = Y)}) %>% do.call(what = "rbind")
        
        va.all <- lapply(names(X), function(Y){
            purrr::imap_dfr(X[[Y]], ~{vertex_attr(.x) %>% as.data.frame() %>% 
                    dplyr::mutate(seed = .y) %>% 
                    dplyr::group_by(seed, type) %>% 
                    dplyr::summarise(N = dplyr::n(), .groups = "keep")}) %>% 
                dplyr::mutate(sub = Y)
        }) %>% do.call(what = "rbind") %>%  
            dplyr::left_join(seed_types, by = c("seed"="name", "sub" = "sub"))
    }
    
    # filter seed.id
    if(!is.null(seed.id)){
        va.all <- va.all %>% dplyr::filter(seed %in% seed.id)
    }
    
    # filter seed.type
    if(!is.null(seed.type)){
        user.seed.type <- seed.type
        va.all <- dplyr::filter(va.all, seed.type %in% user.seed.type)
    }
    
    if(!nrow(va.all)){
        return(NULL)
    }
    
    # user color
    if(!is.null(color)){
        user.color <- as.list(color) %>% # named list/vector
            as.data.frame(check.names = FALSE) %>% t %>% as.data.frame(check.names = FALSE) %>% tibble::rownames_to_column("type") %>% 
            purrr::set_names(c("type", "color"))
        
    } else { # color is NULL  -> defined color
        color.tmp <- va.all$type %>% unique %>%  sort()
        user.color <- mixOmics::color.mixo(seq(color.tmp)) %>% purrr::set_names(color.tmp) %>% 
            as.data.frame(check.names = FALSE)  %>% tibble::rownames_to_column("type") %>% 
            purrr::set_names(c("type", "color"))
    }

    # barplot
    # -----------
    gg.tmp <- ggplot2::ggplot(va.all, aes(x = seed, y = N, fill = type)) + geom_bar(stat = "identity") +
        #scale_fill_identity(guide = "legend", labels = user.color$type) 
        scale_fill_manual(values = user.color$color) + ylab("Node Types") + xlab("Seeds") + labs(fill = "Types") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust=1))
    
    if(is(X, "list.rwr.attributes")){
        gg.tmp <- gg.tmp + facet_grid(.~sub)
    }
    if(plot == TRUE){
        print(gg.tmp)
    }
    return(invisible(gg.tmp))
}

#' Plot RWR subnetwork 
#'
#' Display the subgraph from a RWR results. This function colors adds a specific color to each node based on their 'type' attribute.
#' It also adds a legend including the number of vertices/edges and the number of nodes of specific type.
#' Additionally, the function can display any igraph object.
#' 
#' @param X an igraph object
#' @param color (optional) a named character vector or list, list of color to apply to each type
#' @param plot logical, if TRUE then the plot is produced
#' @param legend (optional) logical, if TRUE then the legend is displayed with number of veretices/edges and the number of nodes of specific type.
#' @param ... Arguments to be passed to the plot method
#' 
#' @return 
#' X is returned with additional vertex attributes
#' 
#' @examples
#' graph1 <- igraph::graph_from_data_frame(list(from = c("A", "B", "A", "D", "C", "A", "C"), 
#'                                              to = c("B", "C", "D", "E", "D", "F", "G")), directed = FALSE)
#' graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("A","B","C"),value = "1")
#' graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("D","E"),value = "2")
#' graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("F", "G"),value = "3")
#' 
#' rwr_res <- random_walk_restart(X = graph1, seed = c("A"))
#' rwr_res_type <- rwr_find_seeds_between_attributes(X = rwr_res, attribute = "type")
#' 
#' plot_rwr_subnetwork(rwr_res_type$A)
#' 
#' 
#' @import ggplot2

plot_rwr_subnetwork <- function(X, color = NULL, plot = TRUE, legend = TRUE, ...){
    # check X
    stopifnot(is(X, "igraph"))
    
    # check color 
    if(!is.null(color)){
        color <- check_named_vector(X = color, var.name = "'color' ")
    }

    # check plot
    plot <- return_true_false(x = plot, default = TRUE)
    legend <- return_true_false(x = legend, default = TRUE)
    
    
    va <- igraph::vertex_attr(X) %>% as.data.frame()
    
    # user color
    if(!is.null(color)){
        user.color <- as.list(color) %>% # named list/vector
            as.data.frame(check.names = FALSE) %>% t %>% as.data.frame(check.names = FALSE) %>% tibble::rownames_to_column("type") %>% 
            purrr::set_names(c("type", "color"))
        
    } else { # color is NULL  -> defined color
        color.tmp <- va$type %>% unique %>%  sort()
        user.color <- mixOmics::color.mixo(seq(color.tmp)) %>% purrr::set_names(color.tmp) %>% 
            as.data.frame(check.names = FALSE) %>% tibble::rownames_to_column("type") %>% 
            purrr::set_names(c("type", "color"))
    }
    
    
    va <- va %>% dplyr::left_join(user.color, by = c("type" = "type"))
        #mutate(color = ifelse(rwr == "seed", 'red', color)) %>% 
    if('rwr' %in% names(va)){
        va <- va %>%
            mutate(shape = ifelse(rwr == "seed", 'rectangle', "circle")) %>%
            mutate(frame.color = ifelse(rwr == "seed", 'red', "black"))
    }
    
    igraph::vertex_attr(X) <- va
    
    # graph stats
    legend.graph.stats <- list(leg = c(paste0("V: ",c(igraph::vcount(X))), paste0("E: ",c(igraph::ecount(X)))),
                               pch = c(1, NA), lty = c(NA, 1))

    ## type
    legend.graph.type <- va %>% group_by(type) %>% summarise(N = dplyr::n()) %>% mutate(leg = paste0(type, ": ", N)) %>% 
        mutate(pch = c(19)) %>% left_join(user.color, by = c('type'))
    
    if(plot == TRUE){
        # plot(X, ...)
        plot(X)
        
        if(legend == TRUE){
            # legend.graph.stats
            legend("topleft", 
                   legend = legend.graph.stats$leg, 
                   pch = legend.graph.stats$pch, 
                   lty = legend.graph.stats$lty)
            
            # legend.graph.type
            legend("bottomleft", 
                   legend = legend.graph.type$leg, 
                   pch = legend.graph.type$pch,
                   col = legend.graph.type$color)
         
            if('rwr' %in% names(va)){
                title(main = va %>% filter(rwr == "seed") %>% pull(name))
            }
        }
    }
    return(X)
}
