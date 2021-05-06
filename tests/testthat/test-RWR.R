context("RWR")

graph1 <- igraph::graph_from_data_frame(as.data.frame(list(from = c("A", "B", "E", "A"), 
                                    to = c("B", "C", "F", "C"))), directed = FALSE)

graph2 <- igraph::induced_subgraph(graph1, vids = c("A", "B", "C"))

graph3 <- igraph::graph_from_data_frame(as.data.frame(list(from = c("A", "B", "E", "A", "Z"), 
                                                           to = c("B", "C", "F", "C", "Z"))), directed = FALSE) %>% 
    igraph::simplify()


MultiplexObject <- RandomWalkRestartMH::create.multiplex(graph1,Layers_Name=c("graph"))
AdjMatrix <- RandomWalkRestartMH::compute.adjacency.matrix(MultiplexObject)
AdjMatrixNorm <- RandomWalkRestartMH::normalize.multiplex.adjacency(AdjMatrix)
Seed <- c("A")
RWR_Results <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(AdjMatrixNorm, MultiplexObject,Seed)
TopResults  <- RandomWalkRestartMH::create.multiplexNetwork.topResults(RWR_Results,MultiplexObject,k=15)

MultiplexObject <- RandomWalkRestartMH::create.multiplex(graph3,Layers_Name=c("graph"))
AdjMatrix <- RandomWalkRestartMH::compute.adjacency.matrix(MultiplexObject)
AdjMatrixNorm <- RandomWalkRestartMH::normalize.multiplex.adjacency(AdjMatrix)
Seed <- c("A")
RWR_Results <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(AdjMatrix, MultiplexObject,Seed)
TopResults  <- RandomWalkRestartMH::create.multiplexNetwork.topResults(RWR_Results,MultiplexObject,k=15)

data(PPI_Network)
PPI_MultiplexObject <- create.multiplex(PPI_Network,Layers_Name=c("PPI"))
AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
SeedGene <- c("PIK3R1")
RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_MultiplexObject,SeedGene)
TopResults_PPI  <- create.multiplexNetwork.topResults(RWR_PPI_Results,PPI_MultiplexObject,k=15)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

