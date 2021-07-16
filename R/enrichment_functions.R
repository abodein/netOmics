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
get_ORA <- function(query, sources = c("GO"), organism = "hsapiens"){
    # validate query (char)
    
    # check organism
    
    # check source
    ORA <- gprofiler2::gost(query = query, organism = organism, significant = FALSE, sources = sources, multi_query = FALSE)
    ORA.res <- ORA$result
    return(ORA.res)
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
    
    
    
