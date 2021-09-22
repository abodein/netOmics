# from timeOmics
check_matrix <- function(X){
    # add rownames and colnames if absent, cast into matrix
    if(!(is.matrix(X) || is.data.frame(X))) return(FALSE)
    
    if(is.data.frame(X)){
        X <- as.matrix(X)
    }
    if(is.null(rownames(X))){
        rownames(X) <- seq_len(nrow(X))
    }
    if(is.null(colnames(X))){
        colnames(X) <- paste0("V", seq_len(ncol(X)))
    }
    return(X)
}

validate_matrix_X <- function(X, var.name = "'X' "){
    # X should be a numeric matrix
    X <- check_matrix(X)
    if(!is.numeric(X)){
        stop(var.name,"must be a numeric matrix/data.frame")
    }
    # if(any(!X)) stop("X must be a numeric matrix/data.frame")
    return(X)
}

validate_list_matrix_X <- function(X, var.name = "'X' "){
    if(!is.list(X)){
        stop(var.name, "must be a list of matrix/data.frame")
    }
    X <- lapply(X, validate_matrix_X)
    return(X)
}

# is_almostInteger <- function (X) 
# {
#     if (!is.numeric(X) & !is.vector(X)) 
#         return(FALSE)
#     if (length(X) != 1) 
#         return(FALSE)
#     if (!is.finite(X)) 
#         return(FALSE)
#     X.round <- round(X)
#     if (X == X.round) 
#         return(TRUE)
#     return(FALSE)
# }

#' @importFrom methods is
check_getCluster <- function(X){
    if(!(is(X, "cluster.df") || is.null(X))){
        stop("cluster must be NULL or a result from getCluster()")
    }
    #stopifnot(is(X, "cluster.df") || is.null(X))
    return(X)
}

check_graph <- function(X){
    stopifnot(is(X, "igraph") || is(X, "grn") || is.list(X) || is(X, "list.igraph"))
    
    if(is(X, "list")){
        stopifnot(all(as.logical(lapply(X, function(x) {is(x, "igraph") || is(x, "list.igraph")}))))
        class(X) <- c("list.igraph", class(X))
    }
    return(X)
}

check_db <- function(X, var.name = "'db' "){
    # ADD list of db
    # x is a dataframe with 2 columns (from, to) or igraph
    if(!(is(X, "igraph") || is(X, "data.frame"))){
        stop(var.name, "must be an igraph or data.frame object")
    }
    if(is(X, "data.frame") & !(all(c("from", "to") %in% colnames(X)))){
        stop(var.name, "must contains the columns 'from' and 'to'")
    }
    return(X)
}

#' @importFrom purrr is_empty
#' @importFrom stats na.omit
check_vector_char <- function(X, X.length = NULL, default = NULL, var.name = "'X' "){
    if(is.null(X)){
        return(default)
    }
    
    # remove NA
    X <- na.omit(X)
    if(is_empty(X)){
        return(default)
    } else if(!is.character(X)){
        stop(var.name, "must be a charactor vector")
    } else if(!is.null(X.length)){
        if(length(X) != X.length){
            stop("invalid length")
        } else { # good length
            return(X)
        }
    } else{
        return(X)
    }
#    return(default)
}

return_true_false <- function(x, default){
    if(is.logical(x)){
        if(is.finite(x)){
            return(x)
        } else { #NA
            return(default)
        }
    } else {
        return(default)
    }
}


check_single_numeric_value <- function(x, min = NULL, max = NULL, var.name = "'r' "){
    if(!is.numeric(x) & !is.matrix(x) & length(x) == 1){
        stop(var.name, "must be a numeric value")
    }
    if(!is.null(min) & !is.null(max)){
        if(x < min | x > max){
            # internal, no need to check min and max order
            stop(var.name, "must be a numeric value between ", min, " and ", max)
        }
    }
    return(x)
}

check_named_vector <- function(X, var.name = "'X' "){
    if(!(is(X, 'list') | is(X, "atomic"))){
        stop(var.name, "must be a named verctor or list")
    }
    if(is.null(names(X))){
        stop(var.name, "must be a named verctor or list")
    }
    return(X)
}


