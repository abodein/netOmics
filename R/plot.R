#' Plot RWR attributes
#' 
summary_plot_rwr_attributes <- function(X, color = NULL, seed.id = NULL, seed.type = NULL){
    stopifnot(is(X, "rwr.attributes") | is(X, "list.rwr.attributes"))
    
    # check seed.id
    seed.id <- check_vector_char(X = seed.id, default = NULL, var.name = "'seed.id' ")
    
    # check seed.type
    seed.type <- check_vector_char(X = seed.type, default = NULL, var.name = "'seed.type' ")
    
    if(is(X, "rwr.attributes")){
        # seed type 
        seed_types <- purrr::imap_dfr(X, ~{vertex_attr(.x) %>% 
                as.data.frame() %>% filter(rwr == "seed") %>% 
                dplyr::select(name, type) %>% purrr::set_names(c("name", "seed.type"))})
        # count layer
        va.all <- purrr::imap_dfr(X, ~{vertex_attr(.x) %>% as.data.frame() %>% 
                mutate(seed = .y) %>% 
                dplyr::group_by(seed, type) %>% 
                dplyr::summarise(N = n(), .groups = "keep")}) %>% 
            dplyr::left_join(seed_types, by = c("seed"="name"))
    } else { #X is list.rwr.attributes
        seed_types <- lapply(names(X), function(Y){
            purrr::imap_dfr(X[[Y]], ~{vertex_attr(.x) %>% 
                as.data.frame() %>% filter(rwr == "seed") %>% 
                dplyr::select(name, type) %>% purrr::set_names(c("name", "seed.type"))}) %>% 
                mutate(sub = Y)}) %>% do.call(what = "rbind")
        
        va.all <- lapply(names(X), function(Y){
            purrr::imap_dfr(X[[Y]], ~{vertex_attr(.x) %>% as.data.frame() %>% 
                    dplyr::mutate(seed = .y) %>% 
                    dplyr::group_by(seed, type) %>% 
                    dplyr::summarise(N = n(), .groups = "keep")}) %>% 
                dplyr::mutate(sub = Y)
        }) %>% do.call(what = "rbind") %>%  
            dplyr::left_join(seed_types, by = c("seed"="name", "sub" = "sub"))
    }
    
    # filter seed.id
    if(!is.null(seed.id)){
        va.all <- va.all %>% dplyr::filter(seed %in% seed.id)
    }
    
    # filter seed.type
    if(!is.null(seed.id)){
        user.seed.type <- seed.type
        va.all <- dplyr::filter(va.all, seed.type %in% user.seed.type)
    }
    
    # check color
    if(!is.null(color)){
        # named list
    } else { # color is NULL
        color <- va.all$type %>% unique %>%  sort()
        user.color <- mixOmics::color.mixo(seq(color)) %>% purrr::set_names(color) %>% 
            as.data.frame()  %>% rownames_to_column("type") %>% 
            set_names(c("type", "color"))
        
        va.all <- va.all %>% 
            dplyr::left_join(user.color, by = c('type' = 'type'))
    }
    
    # barplot
    # -----------
    gg.tmp <- ggplot(va.all, aes(x = seed, y = N, group = type, fill = color)) + geom_bar(stat = "identity") +
        scale_fill_identity()
    
    if(is(X, "list.rwr.attributes")){
        gg.tmp <- gg.tmp + facet_grid(.~sub)
    }
        scale_colour_identity()
}
