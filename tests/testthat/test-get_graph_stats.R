context("get_graph_stats")

graph1 <- igraph::graph_from_data_frame(list(from = c("A", "B", "A", "D", "C", "A", "C"), 
                                             to = c("B", "C", "D", "E", "D", "F", "G")), directed = FALSE)
graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("A","B","C"),value = "1")
graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("D","E"),value = "2")
graph1 <- set_vertex_attr(graph = graph1, name = 'type', index = c("F", "G"),value = "3")


graph1.list <- list("X" = graph1, "Y"= graph1)

test_that("get_graph_stats fails on invalid input", {
    expect_error(get_graph_stats(X = ""))
})

test_that("get_graph_stats works", {
    expect_is(get_graph_stats(X = graph1), "data.frame")
    expect_is(get_graph_stats(X = graph1.list), "data.frame")
})