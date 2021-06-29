context("RWR")

data(HeLa)
X <- HeLa$raw$mRNA
cluster.mRNA <- timeOmics::getCluster(HeLa$getCluster, user.block = "mRNA")
grn.res <- get_grn(X = HeLa$raw$mRNA, cluster = cluster.mRNA, method = "aracne")

X <- grn.res
Xi <- grn.res[["Up_Up"]]
seed <- igraph::V(Xi)$name[1:5]
rwr_res <- random_walk_restart(Xi, seed)
rwr_res <- random_walk_restart(X, seed)

test_that("random_walk_with_restart fails with invalid input", {
    # X
    expect_error(random_walk_restart(X=NULL))
    expect_error(random_walk_restart(X=data.frame()))
    expect_error(random_walk_restart(X=list(3)))
    
    # seed
    expect_error(random_walk_restart(X=X, seed = 3))
    
    # r
    expect_error(random_walk_restart(X=X, r = -2), "'r' must be a numeric value between 0 and 1", fixed = TRUE)
    expect_error(random_walk_restart(X=X, r = matrix(2)), "'r' must be a numeric value", fixed = TRUE)
    
})

test_that("random_walk_with_restart works", {
    
    # X
    expect_is(random_walk_restart(X), "list.rwr")
    expect_is(random_walk_restart(Xi), "rwr")
    
    # seed
    expect_is(random_walk_restart(Xi, seed = seed), "rwr")
    expect_is(random_walk_restart(Xi, seed = 'A'), "rwr")
    expect_is(random_walk_restart(X, seed = seed), "list.rwr")
    expect_is(random_walk_restart(Xi, seed = NULL), "rwr")

    # r
    expect_is(random_walk_restart(Xi, seed = NULL, r = 0.4), "rwr")
    
})


X <- grn.res[["All"]]
all_seed <- igraph::V(X)$name
seed <- c("AAK1","ABCC1","ABL2", "AAAAA")

rwr_res <- random_walk_restart(X, seed)
rwr_cluster_res <- rwr_find_seeds_between_attributes(rwr_res, attribute = "cluster")

rwr_res.list <- random_walk_restart(grn.res, seed)
rwr_cluster_res.list <- rwr_find_seeds_between_attributes(rwr_res.list, attribute = "cluster")


test_that("rwr_find_seeds_between_attributes fails with invalid input", {
    # X
    expect_error(rwr_find_seeds_between_attributes(X=NULL))
    expect_error(rwr_find_seeds_between_attributes(X=matrix()))
    
    # seed
    expect_error(rwr_find_seeds_between_attributes(X=rwr_res))
    
})


test_that("rwr_find_seeds_between_attributes works", {
    # X
    expect_error(rwr_find_seeds_between_attributes(X=rwr_res))
    
    
    
})

