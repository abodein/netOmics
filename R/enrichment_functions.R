#' Get interaction from ORA enrichment analysis
#' 
#' Returns results of an ORA analysis as an interaction graph
#' 
#' @param query a vector (or a list) of character with the ID to perform the ORA analysis
#' @param sources (optional) a character in (GO, KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP)
#' @param organism (optional) a character (default = "hsapiens")
#' @param signif.value (optional) a logical, default = ""
#'
#' @return 
#' a graph object (or list of graph) containing the interaction between the query and the target terms.
#'  
#' @seealso \code{\link[gprofiler2]{gost}}  \code{\link[gprofiler2]{gconvert}}
#'
#' @examples
#' query <- c("IL15", "CDHR5", "TGFA", 'C4B')
#' get_interaction_from_ORA(query, sources = "GO")
#' 
#' query <- list("All" = c("IL15", "CDHR5", "TGFA", 'C4B'),
#'               "c1" = c("IL15", "CDHR5", "TGFA"))
#' get_interaction_from_ORA(query, sources = "GO")
#' 
#' @import gprofiler2
#' @importFrom dplyr pull select filter
#' @export
get_interaction_from_ORA <- function(query, sources = "GO", organism = "hsapiens", signif.value = TRUE){
    # validate query (char)
    if(is(query, "list")){
        query <- lapply(query,function(x)check_vector_char(x, var.name = "'query '"))
        if(is.null(names(query))){
            names(query) <- seq_along(query)
        }
    } else {
        query <- check_vector_char(query)
    }
    
    # check organism
    organism = check_vector_char(X = organism, X.length = 1, default = "hsapiens", var.name = "'organism' ")
    
    # check source
    sources <-  match.arg(arg = sources, choices = c("GO", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP"), several.ok = FALSE)
    sources <- check_vector_char(sources, default = "GO")  # default value
    
    # check signif
    signif.value <- return_true_false(signif.value, default = TRUE)
    
    if(is(query, "list")){
        res.ora <- list()
        term_map <- list()
        res.graph <- list()
        for(i in names(query)){
            res.ora[[i]] <- get_ORA(query = query[[i]], sources = sources, organism = organism)
            
            term_map_tmp <- gprofiler2::gconvert(query = query[[i]], organism = organism, target = sources) 
            
            target_id <- (res.ora[[i]] %>% dplyr::filter(significant == signif.value) %>% dplyr::pull(term_id))
            
            term_map[[i]] <- term_map_tmp %>% 
                dplyr::filter(target %in% target_id) %>% 
                dplyr::select(input, target) %>% unique %>% na.omit()
            
            res.graph[[i]] <- igraph::graph_from_data_frame(term_map[[i]], directed = FALSE)
            res.graph[[i]] <- set_vertex_attr(graph = res.graph[[i]], name = "mode", index = term_map[[i]]$input, value = "core")
            res.graph[[i]] <- set_vertex_attr(graph = res.graph[[i]], name = "mode", index = term_map[[i]]$target, value = "extended")
            class(res.graph) <- c("list.interaction.igraph", "list.igraph")
        }
    } else {
        # query is not a list
        res.ora <- get_ORA(query = query, sources = sources, organism = organism)
        
        term_map_tmp <- gprofiler2::gconvert(query = query, organism = organism, target = sources) 
        
        target_id <- (res.ora %>% dplyr::filter(significant == signif.value) %>% dplyr::pull(term_id))
        
        term_map <- term_map_tmp %>% 
            dplyr::filter(target %in% target_id) %>% 
            dplyr::select(input, target) %>% unique %>% na.omit()
        
        res.graph <- igraph::graph_from_data_frame(term_map, directed = FALSE)
        res.graph <- set_vertex_attr(graph = res.graph, name = "mode", index = term_map$input, value = "core")
        res.graph <- set_vertex_attr(graph = res.graph, name = "mode", index = term_map$target, value = "extended")
        
        class(res.graph) <- c("interaction.igraph", "igraph")
        
    }
    return(res.graph)
}


#' ORA enrichment analysis
#' 
#' Returns results of an ORA analysis
#' 
#' @param query a vector of character, a lit of ID
#' @param sources a character or list of character
#' @param organism a character (default = "hsapiens")
#'
#' @return 
#' a data.frame containing the enrichment result
#'  
#' @seealso \code{\link[gprofiler2]{gost}}
#'
#' @import gprofiler2
get_ORA <- function(query, sources = NULL, organism = "hsapiens"){
    
    if (is(query, "list")){
        res <- list()
        for(i in names(query)){
            ORA <- gprofiler2::gost(query = query[[i]], organism = organism, significant = FALSE, sources = sources, multi_query = FALSE)
            ORA.res <- ORA$result
            if(!is.null(ORA.res)){
                ORA.res <- ORA.res %>% mutate(cluster =  i) %>% 
                    dplyr::select("cluster", "term_id","source","term_name","p_value","significant", "term_size", "query_size", "intersection_size", "precision", "recall")
                res[[i]] <- ORA.res
            } 
        }
        RES <- purrr::map_dfr(res, ~.x)
        
    } else {
        ORA <- gprofiler2::gost(query = query, organism = organism, significant = FALSE, sources = sources, multi_query = FALSE)
        ORA.res <- ORA$result
        if(!is.null(ORA.res)){
            ORA.res <- ORA.res %>% mutate(cluster =  "All") %>% 
                dplyr::select("cluster", "term_id","source","term_name","p_value","significant", "term_size", "query_size", "intersection_size", "precision", "recall")
            RES <- ORA.res
        } else {
            RES <- NULL
        }
    }
    return(RES)
}
    
#' Get GO info
#' 
#' From a GO terms (GOID), return definition, ontology and term values from GO.db
#' 
#' @param go a character, GO term
#' 
#' @return 
#' a data.frame with the following columns: "GOID", "DEFINITION", "ONTOLOGY", "TERM"
#' 
#' @import GO.db
#' @importFrom AnnotationDbi keytypes
get_go_info <- function(go){
    res <- AnnotationDbi::select(x = GO.db, keys = go, keytype = 'GOID', columns = keytypes(GO.db))
    return(res)
}

    
