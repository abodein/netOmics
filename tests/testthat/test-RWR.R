context("RWR")

graph1 <- igraph::graph_from_data_frame(list(from = c("A", "B", "A", "D", "C", "A", "C"), 
                                             to = c("B", "C", "D", "E", "D", "F", "G")), directed = FALSE)
graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("A","B","C"),value = "1")
graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("D","E"),value = "2")
graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("F", "G"),value = "3")

rwr_res <- random_walk_restart(X = graph1, seed = c("A", "B", "C", "D", "E"))
rwr_res_type <- rwr_find_seeds_between_attributes(X = rwr_res, attribute = "type", k = 3)

graph1.list <- list("X" = graph1, "Y"= graph1)

rwr_res.list <- random_walk_restart(X = graph1.list, seed = c("A", "B", "C", "D", "E"))
rwr_res_type.list <- rwr_find_seeds_between_attributes(X = rwr_res.list, attribute = "type", k = 3)

test_that("random_walk_with_restart fails with invalid input", {
    # X
    expect_error(random_walk_restart(X=NULL))
    expect_error(random_walk_restart(X=data.frame()))
    expect_error(random_walk_restart(X=list(3)))
    
    # seed
    expect_error(random_walk_restart(X=graph1, seed = 3))
    
    # r
    expect_error(random_walk_restart(X=graph1, r = -2), "'r' must be a numeric value between 0 and 1", fixed = TRUE)
    expect_error(random_walk_restart(X=graph1, r = matrix(2)), "'r' must be a numeric value", fixed = TRUE)
    
})

test_that("random_walk_with_restart works", {
    
    # X
    expect_is(random_walk_restart(graph1.list), "list.rwr")
    expect_is(random_walk_restart(graph1), "rwr")
    
    # seed
    expect_is(random_walk_restart(graph1, seed = 'A'), "rwr")
    expect_is(random_walk_restart(graph1.list, seed = "A"), "list.rwr")
    expect_is(random_walk_restart(graph1, seed = NULL), "rwr")

    # r
    expect_is(random_walk_restart(graph1, seed = NULL, r = 0.4), "rwr")
    
})



test_that("rwr_find_seeds_between_attributes fails with invalid input", {
    # X
    expect_error(rwr_find_seeds_between_attributes(X=NULL))
    expect_error(rwr_find_seeds_between_attributes(X=matrix()))
    
    # k
    expect_error(rwr_find_seeds_between_attributes(X=rwr_res, k = -5), "'k' must be a numeric value between 0 and 200", fixed = TRUE)
    expect_error(rwr_find_seeds_between_attributes(X=rwr_res, k = list(3)),"'k' must be a numeric value", fixed = TRUE)
    
    # seed / attribute
    expect_error(rwr_find_seeds_between_attributes(X=rwr_res, seed = matrix(3)), "'seed' must be a charactor vector", fixed = TRUE)
    expect_error(rwr_find_seeds_between_attributes(X=rwr_res, attribute = matrix(3)), "'attribute' must be a charactor vector", fixed = TRUE)
    expect_error(rwr_find_seeds_between_attributes(X=rwr_res, attribute = c("4", "5")), "invalid length", fixed = TRUE)
})


test_that("rwr_find_seeds_between_attributes works", {
    # X
    expect_is(rwr_find_seeds_between_attributes(X=rwr_res), "rwr.attributes")
    expect_is(rwr_find_seeds_between_attributes(X=rwr_res.list), "list.rwr.attributes")
    
    # seed
    expect_is(rwr_find_seeds_between_attributes(X=rwr_res, seed = logical(0)), "rwr.attributes")
    expect_is(rwr_find_seeds_between_attributes(X=rwr_res, attribute = logical(0)), "rwr.attributes")
    expect_is(rwr_find_seeds_between_attributes(X=rwr_res, seed = "Z"), "rwr.attributes")
    
    expect_is(rwr_find_seeds_between_attributes(X=rwr_res.list, seed = "Z"), "list.rwr.attributes")
    
    
    # k 
    expect_is(rwr_find_seeds_between_attributes(X=rwr_res, k = 1) , "rwr.attributes")
    expect_is(rwr_find_seeds_between_attributes(X=rwr_res, k = NULL), "rwr.attributes")
    
})


test_that("rwr_find_closest_type fails with invalid input", {
    # X
    expect_error(rwr_find_closest_type(X=NULL))
    expect_error(rwr_find_closest_type(X=matrix()))
    
    # seed / attribute
    expect_error(rwr_find_closest_type(X=rwr_res, seed = matrix(3)), "'seed' must be a charactor vector", fixed = TRUE)
    expect_error(rwr_find_closest_type(X=rwr_res, attribute = matrix(3)), "'attribute' must be a charactor vector", fixed = TRUE)
    expect_error(rwr_find_closest_type(X=rwr_res, attribute = c("4", "5")), "invalid length", fixed = TRUE)
    
    # value
    expect_error(rwr_find_closest_type(X=rwr_res, attribute = NULL, seed = NULL, value = c("a", "b")), "invalid length", fixed = TRUE)
    expect_error(rwr_find_closest_type(X=rwr_res, attribute = NULL, seed = NULL, value = matrix(1,2)), "'value' must be a charactor vector", fixed = TRUE)
    
})


test_that("rwr_find_closest_type works", {
    # X
    expect_is(rwr_find_closest_type(X=rwr_res), "rwr.closest")
    expect_is(rwr_find_closest_type(X=rwr_res.list), "list.rwr.closest")
    
    # seeds / attributes 
    expect_is(rwr_find_closest_type(X=rwr_res, seed = logical(0)), "rwr.closest")
    expect_is(rwr_find_closest_type(X=rwr_res, attribute = logical(0)), "rwr.closest")
    expect_is(rwr_find_closest_type(X=rwr_res, attribute = NULL, seed = NULL), "rwr.closest")
    expect_is(rwr_find_closest_type(X=rwr_res, attribute = "NULL", seed = NULL), "rwr.closest")
    
    #value 
    expect_is(rwr_find_closest_type(X=rwr_res, attribute = "NULL", seed = NULL, value = NULL), "rwr.closest")
    expect_is(rwr_find_closest_type(X=rwr_res, attribute = NULL, seed = NULL, value = "test"), "rwr.closest")
    expect_is(rwr_find_closest_type(X=rwr_res, value = "test"), "rwr.closest")
    
    expect_is(rwr_find_closest_type(X=rwr_res, seed = "Z"), "rwr.closest")
    
    expect_is(rwr_find_closest_type(X=rwr_res.list, seed = "Z"), "list.rwr.closest")
    
})

# test_that("top_k_graph returns NULL", {
#     rwr_top_k_graph
#     rwr_top_k_graph(X = rwr_res$graph, RWRM_Result_Object = rwr_res$rwr, Seed = "A", k = 1)
# })
    