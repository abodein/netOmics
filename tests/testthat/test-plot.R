context("plot")

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

summary_plot_rwr_attributes(rwr_res_type)

test_that("summary_plot_rwr_attributes fails on invalid input", {
    #X
    expect_error(summary_plot_rwr_attributes(X = c()))
    expect_error(summary_plot_rwr_attributes(X = NULL))
    # color must be a named vector/list
    expect_error(summary_plot_rwr_attributes(X = rwr_res_type, color = 3), "'color' must be a named verctor or list", fixed = TRUE)
    expect_error(summary_plot_rwr_attributes(X = rwr_res_type, color = data.frame()), "'color' must be a named verctor or list", fixed = TRUE)
    
    # seed.id 
    expect_error(summary_plot_rwr_attributes(X = rwr_res_type, seed.id = c(1)), "'seed.id' must be a charactor vector", fixed = TRUE)
    # seed.type
    expect_error(summary_plot_rwr_attributes(X = rwr_res_type, seed.type = c(1)), "'seed.type' must be a charactor vector", fixed = TRUE)
    
    # plot TRUE/FALSE -> default
})

test_that("summary_plot_rwr_attributes works", {
    expect_is(summary_plot_rwr_attributes(rwr_res_type), "ggplot")
    
    expect_is(summary_plot_rwr_attributes(rwr_res_type,color = list("1" = "red", "2"="blue", "3"="pink")),"ggplot")
    expect_is(summary_plot_rwr_attributes(rwr_res_type, seed.id = c("A","B")), "ggplot")
    expect_is(summary_plot_rwr_attributes(rwr_res_type, seed.type = c("1")), "ggplot")
    
    expect_is(summary_plot_rwr_attributes(rwr_res_type, seed.type = c("1"), plot = NA), "ggplot")
    expect_is(summary_plot_rwr_attributes(rwr_res_type, seed.type = c("1"), plot = FALSE), "ggplot")
    
    expect_is(summary_plot_rwr_attributes(rwr_res_type.list), "ggplot")
    
    expect_null(summary_plot_rwr_attributes(rwr_res_type, seed.type = c("f"), plot = FALSE))
    
})

test_that("plot_rwr_subnetwork fails on invalid input", {
    # X
    expect_error(plot_rwr_subnetwork(c()))
    
    # color
    expect_error(plot_rwr_subnetwork(X = rwr_res_type$A, color = 3), "'color' must be a named verctor or list", fixed = TRUE)
    expect_error(plot_rwr_subnetwork(X = rwr_res_type$A, color = data.frame()), "'color' must be a named verctor or list", fixed = TRUE)
    
    # plot, legend -> TRUE/FALSE 
    # ...
})

test_that("plot_rwr_subnetwork works", {
    # X
    expect_is(plot_rwr_subnetwork(X = rwr_res_type$A), "igraph")
    
    #color
    expect_is(plot_rwr_subnetwork(X = rwr_res_type$A, color = list("1" = "red", "2"="blue", "3"="green")), "igraph")
    expect_is(plot_rwr_subnetwork(X = rwr_res_type$A, color = list("1" = "red", "4"="blue", "3"="green")), "igraph")
    
    #plot
    expect_is(plot_rwr_subnetwork(X = rwr_res_type$A, plot = FALSE), "igraph")
    expect_is(plot_rwr_subnetwork(X = rwr_res_type$A, plot = NA), "igraph")
    
    expect_is(plot_rwr_subnetwork(X = rwr_res_type$A, legend = FALSE), "igraph")
    expect_is(plot_rwr_subnetwork(X = rwr_res_type$A, legend = NA), "igraph")
    expect_is(plot_rwr_subnetwork(X = rwr_res_type$A, legend = NULL), "igraph")
    
})

    