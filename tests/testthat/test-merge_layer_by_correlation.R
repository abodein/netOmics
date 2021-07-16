context('merge_layer_by_correlation')


test_that("interaction_from_correlation fails on invalind input", {
    X <- matrix(rexp(200, rate=.1), ncol=20)
    Y <- matrix(rexp(200, rate=.1), ncol=20)
    expect_error(interaction_from_correlation(X = "",Y), "'X' must be a numeric matrix/data.frame", fixed = TRUE)

    X <- list(matrix(rexp(200, rate=.1), ncol=20), matrix(rexp(200, rate=.1), ncol=2))
    expect_error(interaction_from_correlation(X = X,Y),"'X' must have the same number of rows", fixed = TRUE)
    
    X <- matrix(rexp(200, rate=.1), ncol=20)
    Y <- list(matrix(rexp(200, rate=.1), ncol=20), matrix(rexp(200, rate=.1), ncol=2))
    expect_error(interaction_from_correlation(X = X,Y),"'Y' must have the same number of rows", fixed = TRUE)
    expect_error(interaction_from_correlation(X = X,Y = ""),"'Y' must be a numeric matrix/data.frame", fixed = TRUE)
    
    
})

test_that("interaction_from_correlation works", {
    X <- matrix(rexp(200, rate=.1), ncol=20)
    Y <- matrix(rexp(200, rate=.1), ncol=20)
    expect_is(interaction_from_correlation(X,Y), "igraph")
    
    X <- list(matrix(rexp(200, rate=.1), ncol=20), matrix(rexp(200, rate=.1), ncol=20))
    Y <- matrix(rexp(200, rate=.1), ncol=20)
    expect_is(interaction_from_correlation(X,Y), "igraph")
    
    X <- list(matrix(rexp(200, rate=.1), ncol=20), matrix(rexp(200, rate=.1), ncol=20))
    Y <- list(matrix(rexp(200, rate=.1), ncol=20), matrix(rexp(200, rate=.1), ncol=20))
    expect_is(interaction_from_correlation(X,Y), "igraph")
})