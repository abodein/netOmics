#' ORA enrichment analysis
#' 
#' Returns results of an ORA analysis
#' 
#' @import gprofiler2
get_ORA <- function(query, sources = c("GO"), organism = "hsapiens"){
    # validate query (char)
    
    # check organism
    
    # check source
    ORA <- gprofiler2::gost(query = query, organism = organism, significant = FALSE, sources = sources, multi_query = FALSE)
    ORA.res <- ORA$result
    return(result)
}
    

#' @import GO.db
#' @import AnnotationDbi
get_go_info <- function(go){
    res <- AnnotationDbi::select(x = GO.db, keys = go, keytype = 'GOID', columns = keytypes(GO.db))
    return(res)
}
    
    
    
