context("utils")

# # from timeOmics
# check_matrix <- function(X){
#     # add rownames and colnames if absent, cast into matrix
#     if(!(is.matrix(X) || is.data.frame(X))) return(FALSE)
#     
#     if(is.data.frame(X)){
#         X <- as.matrix(X)
#     }
#     if(is.null(rownames(X))){
#         rownames(X) <- 1:nrow(X)
#     }
#     if(is.null(colnames(X))){
#         colnames(X) <- paste0("V", 1:ncol(X))
#     }
#     return(X)
# }

test_that("check_matrix", {
    expect_false(check_matrix(c(1,2,3)))
    X = data.frame("a" = c(1,2,3), 'b'= c(2,3,4))
    expect_is(check_matrix(X), 'matrix')
    X = matrix(c(1,2,3,2,3,4), nrow = 2)
    expect_is(check_matrix(X), 'matrix')
    rownames(X) <-  NULL
    expect_is(check_matrix(X), 'matrix')
    colnames(X) <-  NULL
   expect_is(check_matrix(X), 'matrix')
})
    
# validate_matrix_X <- function(X){
#     # X should be a numeric matrix
#     X <- check_matrix(X)
#     if(!is.numeric(X)){
#         stop("'X' must be a numeric matrix/data.frame")
#     }
#     # if(any(!X)) stop("'X' must be a numeric matrix/data.frame")
#     return(X)
# }

test_that("validate_matrix_X", {
    X = data.frame("a" = c("A","B", "C"), 'b'= c(NA))
    expect_error(validate_matrix_X(X), "'X' must be a numeric matrix/data.frame", fixed = TRUE)
})

# validate_list_matrix_X <- function(X){
#     if(!is.list(X)){
#         stop("'X' must be a list of matrix/data.frame")
#     }
#     X <- lapply(X, validate_matrix_X)
#     return(X)
# }

test_that("validate_list_matrix_X", {
    expect_error(validate_list_matrix_X(X = 2), "'X' must be a list of matrix/data.frame", fixed=TRUE)
    X <- list(data.frame("a" = c(1,2,3), 'b'= c(2,3,4)), data.frame("a" = c(1,2,3), 'b'= c(2,3,4)))
    expect_is(validate_list_matrix_X(X), "list")
})
    






