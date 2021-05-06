context("combine_layers")

graph1 <- igraph::graph_from_data_frame(list(from = c("A", "B"), to = c("B", "C")), directed = FALSE)
graph2 <- igraph::graph_from_data_frame(list(from = c(1), to = c(2)), directed = FALSE)
graph3 <- igraph::make_empty_graph(directed = FALSE)

interaction.df1 <- as.data.frame(list(from = c("C", "B"), to = c(1, 2)))
interaction.df2 <- as.data.frame(list(from = c("D", "D"), to = c(3, 4)))


graph1.list <- list(graph1, graph1)
graph2.list <- list(graph2, graph2)
graph3.list <- list(graph1, graph1, graph1)
class(graph1.list) <- class(graph2.list) <- class(graph3.list) <- "list.igraph"

graph1.list.named <- list(graph1 = graph1, graph2 = graph1)
graph2.list.named <- list(graph1 = graph2, graph2 = graph2)
graph3.list.named <- list(graph2 = graph2, graph1 = graph2)
graph4.list.named <- list(graph3 = graph2, graph2 = graph2)
class(graph1.list.named) <- class(graph2.list.named) <- class(graph3.list.named) <- class(graph4.list.named)  <- "list.igraph"




test_that("get_grn fails on invalid input", {
    expect_error(combine_layers(graph1 = ""), "graph1 must be an igraph or list.igraph object", fixed = TRUE)
    expect_error(combine_layers(graph1, graph2 = ""), "graph2 must be an igraph or list.igraph object or NULL", fixed = TRUE)
    expect_error(combine_layers(graph1, graph2, interaction.df = ""))
    expect_error(combine_layers(graph1, graph2 = graph2.list), "graph1 and graph2 must have the same length", fixed = TRUE)
    expect_error(combine_layers(graph1.list, graph2 = graph3.list), "graph1 and graph2 must have the same length", fixed = TRUE)
    
    expect_error(combine_layers(graph1.list, graph2 = graph3.list), "graph1 and graph2 must have the same length", fixed = TRUE)
    expect_error(combine_layers(graph1.list.named, graph2 = graph4.list.named), "graph1 and graph2 must have the same names", fixed = TRUE)
    
})


test_that("combine_layers works", {
    # graph1 and graph2
    expect_is(combine_layers(graph1 = graph1, graph2 = graph2), "merged.igraph")
    
    # graph1 and interaction.df
    expect_is(combine_layers(graph1 = graph1, interaction.df = interaction.df1), "merged.igraph")
    expect_is(combine_layers(graph1 = graph1, interaction.df = interaction.df2), "merged.igraph")

    expect_is(combine_layers(graph1 = graph1, graph2 = graph2, interaction.df = interaction.df1), "merged.igraph")
    expect_is(combine_layers(graph1 = graph1, graph2 = graph2, interaction.df = interaction.df2), "merged.igraph")
    
    expect_is(combine_layers(graph1 = graph1.list, graph2 = graph2), "list.merged.igraph")
    expect_is(combine_layers(graph1 = graph1.list, graph2 = graph2, interaction.df = interaction.df1), "list.merged.igraph")
    
    expect_is(combine_layers(graph1 = graph1.list, graph2 = graph2.list), "list.merged.igraph")
    expect_is(combine_layers(graph1 = graph1.list, graph2 = graph2.list, interaction.df = interaction.df1), "list.merged.igraph")
    
    #named list
    expect_is(combine_layers(graph1.list.named, graph2 = graph2.list.named), "list.merged.igraph")
    expect_is(combine_layers(graph1.list.named, graph2 = graph3.list.named), "list.merged.igraph")
    
    # with empty graph
    expect_warning(combine_layers(graph1, graph2 = graph3), "Some, but not all graphs are named, not using vertex names")
})
    