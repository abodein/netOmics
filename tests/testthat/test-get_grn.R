context("get_grn")

#' data(HeLa)
#' #' # grn only on gene
#' cluster.mRNA <- HeLa$getCluster %>% dplyr::filter(block == "mRNA")
#' X <- HeLa$raw$mRNA
#' #grn.res <- get_grn(X = HeLa$raw$mRNA, cluster = cluster.mRNA, method = "aracne")
#' 
data("hmp_T2D")
cluster.mRNA <- hmp_T2D$getCluster.res %>% dplyr::filter(block == "RNA")
X <- hmp_T2D$raw$RNA

test_that("get_grn fails on invalid input - X", {
    # X
    expect_error(get_grn(X = ""), "'X' must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(get_grn(X = 1), "'X' must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(get_grn(X = NA), "'X' must be a numeric matrix/data.frame", fixed = TRUE)
    expect_error(get_grn(X = list()), "'X' must be a numeric matrix/data.frame", fixed = TRUE)
    
})


test_that("get_grn fails on invalid input - X", {
    expect_error(get_grn(X = X, cluster = ""), "cluster must be NULL or a result from getCluster()", fixed = TRUE)
})

test_that("get_grn fails on invalid input - method", {
    expect_error(get_grn(X = X, cluster = NULL, method = ""))
})

test_that("get_grn works", {
    # expect_is(get_grn(X = X, cluster = cluster.mRNA, method = "aracne"), "list.igraph")
    expect_warning(get_grn(X = X, cluster = cluster.mRNA, method = "aracne"))
    
    # expect_is(get_grn(X = X, method = "aracne"), "igraph")
    expect_warning(get_grn(X = X, method = "aracne"))
    
})
