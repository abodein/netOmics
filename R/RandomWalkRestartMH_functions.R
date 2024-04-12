# Credit to Valdeolivas et al.
# All the functions below come from the package RandomWalkRestartMH. 
# Due to build errors in Bioconducteur, and to avoid depreciation of the netOmics package, all functions imported by the package have been repatriated here.

# need ti keep:
# - RandomWalkRestartMH::create.multiplex
# - RandomWalkRestartMH::compute.adjacency.matrix
# - RandomWalkRestartMH::normalize.multiplex.adjacency
# - RandomWalkRestartMH::Random.Walk.Restart.Multiplex

#' @importFrom igraph set_vertex_attr
create.multiplex <- function(LayersList,...){
  
  if (!class(LayersList) == "list"){
    stop("The input object should be a list of graphs.")
  }
  
  
  Number_of_Layers <- length(LayersList)
  SeqLayers <- seq(Number_of_Layers)
  Layers_Name <- names(LayersList)
  
  # if (!all(sapply(SeqLayers, function(x) is.igraph(LayersList[[x]])))){
  #   stop("Not igraph objects")
  # }
  
  Layer_List <- lapply(SeqLayers, function (x) {
    if (is.null(V(LayersList[[x]])$name)){
      LayersList[[x]] <- 
        igraph::set_vertex_attr(LayersList[[x]],"name", 
                        value=seq(1,vcount(LayersList[[x]]),by=1))
    } else {
      LayersList[[x]]
    }
  })
  
  ## We simplify the layers 
  Layer_List <- 
    lapply(SeqLayers, function(x) simplify.layers(Layer_List[[x]]))
  
  ## We set the names of the layers. 
  
  if (is.null(Layers_Name)){
    names(Layer_List) <- paste0("Layer_", SeqLayers)
  } else {
    names(Layer_List) <- Layers_Name
  }
  
  ## We get a pool of nodes (Nodes in any of the layers.)
  Pool_of_Nodes <- 
    sort(unique(unlist(lapply(SeqLayers, 
                              function(x) V(Layer_List[[x]])$name))))
  
  Number_of_Nodes <- length(Pool_of_Nodes)
  
  Layer_List <-
    lapply(Layer_List, add.missing.nodes,Number_of_Layers,Pool_of_Nodes)
  
  # We set the attributes of the layer
  counter <- 0 
  Layer_List <- lapply(Layer_List, function(x) { 
    counter <<- counter + 1; 
    igraph::set_edge_attr(x,"type",igraph::E(x), value = names(Layer_List)[counter])
  })
  
  
  MultiplexObject <- c(Layer_List,list(Pool_of_Nodes=Pool_of_Nodes,
                                       Number_of_Nodes_Multiplex=Number_of_Nodes, 
                                       Number_of_Layers=Number_of_Layers))
  
  class(MultiplexObject) <- "Multiplex"
  
  return(MultiplexObject)
}



# internal 
#' @importFrom igraph as.undirected is_weighted E simplify
simplify.layers <- function(Input_Layer){
  
  ## Undirected Graphs
  Layer <- igraph::as.undirected(Input_Layer, mode = c("collapse"),
                         edge.attr.comb = igraph::igraph_opt("edge.attr.comb"))
  
  ## Unweighted or Weigthed Graphs
  if (igraph::is_weighted(Layer)){
    b <- 1
    weigths_layer <- igraph::E(Layer)$weight
    if (min(weigths_layer) != max(weigths_layer)){
      a <- min(weigths_layer)/max(weigths_layer)
      range01 <- (b-a)*(weigths_layer-min(weigths_layer))/
        (max(weigths_layer)-min(weigths_layer)) + a
      igraph::E(Layer)$weight <- range01
    } else {
      igraph::E(Layer)$weight <- rep(1, length(weigths_layer))
    }
  } else {
    igraph::E(Layer)$weight <- rep(1, ecount(Layer))
  }
  
  ## Simple Graphs
  Layer <- 
    igraph::simplify(Layer,remove.multiple = TRUE,remove.loops = TRUE, 
                     edge.attr.comb=mean)
  
  return(Layer)
}

#' @importFrom igraph add_vertices
add.missing.nodes <- function (Layers,Nr_Layers,NodeNames) {
  
  igraph::add_vertices(Layers,
               length(NodeNames[which(!NodeNames %in% igraph::V(Layers)$name)]),
               name=NodeNames[which(!NodeNames %in%  igraph::V(Layers)$name)])
}


#' @importFrom Matrix Diagonal bdiag
#' @importFrom igraph as_adjacency_matrix is_weighted

compute.adjacency.matrix <- function(x,delta = 0.5)
{
  if (!isMultiplex(x) & !isMultiplexHet(x)) {
    stop("Not a Multiplex or Multiplex Heterogeneous object")
  }
  if (delta > 1 || delta <= 0) {
    stop("Delta should be between 0 and 1")
  }
  
  N <- x$Number_of_Nodes_Multiplex
  L <- x$Number_of_Layers
  
  ## We impose delta=0 in the monoplex case.
  if (L==1){
    delta = 0
  }
  
  Layers_Names <- names(x)[seq(L)]
  
  ## IDEM_MATRIX.
  Idem_Matrix <- Matrix::Diagonal(N, x = 1)
  
  counter <- 0 
  Layers_List <- lapply(x[Layers_Names],function(x){
    
    counter <<- counter + 1;    
    if (igraph::is_weighted(x)){ 
      Adjacency_Layer <-  igraph::as_adjacency_matrix(x,sparse = TRUE, 
                                              attr = "weight")
    } else {
      Adjacency_Layer <-  igraph::as_adjacency_matrix(x,sparse = TRUE)
    }
    
    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),
                                       order(colnames(Adjacency_Layer))]
    colnames(Adjacency_Layer) <- 
      paste0(colnames(Adjacency_Layer),"_",counter)
    rownames(Adjacency_Layer) <- 
      paste0(rownames(Adjacency_Layer),"_",counter)
    Adjacency_Layer
  })
  
  MyColNames <- unlist(lapply(Layers_List, function (x) unlist(colnames(x))))
  MyRowNames <- unlist(lapply(Layers_List, function (x) unlist(rownames(x))))
  names(MyColNames) <- c()
  names(MyRowNames) <- c()
  SupraAdjacencyMatrix <- (1-delta)*(Matrix::bdiag(unlist(Layers_List)))
  colnames(SupraAdjacencyMatrix) <-MyColNames
  rownames(SupraAdjacencyMatrix) <-MyRowNames
  
  offdiag <- (delta/(L-1))*Idem_Matrix
  
  i <- seq_len(L)
  Position_ini_row <- 1 + (i-1)*N
  Position_end_row <- N + (i-1)*N
  j <- seq_len(L)
  Position_ini_col <- 1 + (j-1)*N
  Position_end_col <- N + (j-1)*N
  
  for (i in seq_len(L)){
    for (j in seq_len(L)){
      if (j != i){
        SupraAdjacencyMatrix[(Position_ini_row[i]:Position_end_row[i]),
                             (Position_ini_col[j]:Position_end_col[j])] <- offdiag
      }    
    }
  }
  
  SupraAdjacencyMatrix <- as(SupraAdjacencyMatrix, "dgCMatrix")
  return(SupraAdjacencyMatrix)
}
#' @importFrom Matrix t colSums
normalize.multiplex.adjacency <- function(x)
{
  if (!is(x,"dgCMatrix")){
    stop("Not a dgCMatrix object of Matrix package")
  }
  
  Adj_Matrix_Norm <- Matrix::t(Matrix::t(x)/(Matrix::colSums(x, na.rm = FALSE, dims = 1,
                                             sparseResult = FALSE)))
  
  return(Adj_Matrix_Norm)
}

Random.Walk.Restart.Multiplex <- function(x, MultiplexObject, Seeds, 
                                                  r=0.7,tau,MeanType="Geometric", DispResults="TopScores",...){

  
  L <- MultiplexObject$Number_of_Layers
  N <- MultiplexObject$Number_of_Nodes
  
  Seeds <- as.character(Seeds)
  if (length(Seeds) < 1 | length(Seeds) >= N){
    stop("The length of the vector containing the seed nodes is not 
         correct")
  } else {
    if (!all(Seeds %in% MultiplexObject$Pool_of_Nodes)){
      stop("Some of the seeds are not nodes of the network")
      
    }
  }
  
  if (r >= 1 || r <= 0) {
    stop("Restart partameter should be between 0 and 1")
  }
  
  if(missing(tau)){
    tau <- rep(1,L)/L
  } else {
    tau <- as.numeric(tau)
    if (sum(tau)/L != 1) {
      stop("The sum of the components of tau divided by the number of 
           layers should be 1")
    }
  }
  
  if(!(MeanType %in% c("Geometric","Arithmetic","Sum"))){
    stop("The type mean should be Geometric, Arithmetic or Sum")
  }
  
  if(!(DispResults %in% c("TopScores","Alphabetic"))){
    stop("The way to display RWRM results should be TopScores or
         Alphabetic")
  }
  
  ## We define the threshold and the number maximum of iterations for
  ## the random walker.
  Threeshold <- 1e-10
  NetworkSize <- ncol(x)
  
  ## We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1
  
  ## We compute the scores for the different seeds.
  Seeds_Score <- get.seed.scoresMultiplex(Seeds,L,tau)
  
  ## We define the prox_vector(The vector we will move after the first RWR
  ## iteration. We start from The seed. We have to take in account
  ## that the walker with restart in some of the Seed nodes, depending on
  ## the score we gave in that file).
  prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)
  
  prox_vector[which(colnames(x) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])
  
  prox_vector  <- prox_vector/sum(prox_vector)
  restart_vector <-  prox_vector
  
  while(residue >= Threeshold){
    
    old_prox_vector <- prox_vector
    prox_vector <- (1-r)*(x %*% prox_vector) + r*restart_vector
    residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
    iter <- iter + 1;
  }
  
  NodeNames <- character(length = N)
  Score = numeric(length = N)
  
  rank_global <- data.frame(NodeNames = NodeNames, Score = Score)
  rank_global$NodeNames <- gsub("_1", "", row.names(prox_vector)[seq_len(N)])
  
  if (MeanType=="Geometric"){
    rank_global$Score <- geometric.mean(as.vector(prox_vector[,1]),L,N)    
  } else {
    if (MeanType=="Arithmetic") {
      rank_global$Score <- regular.mean(as.vector(prox_vector[,1]),L,N)    
    } else {
      rank_global$Score <- sumValues(as.vector(prox_vector[,1]),L,N)    
    }
  }
  
  if (DispResults=="TopScores"){
    ## We sort the nodes according to their score.
    Global_results <- 
      rank_global[with(rank_global, order(-Score, NodeNames)), ]
    
    ### We remove the seed nodes from the Ranking and we write the results.
    Global_results <- 
      Global_results[which(!Global_results$NodeNames %in% Seeds),]
  } else {
    Global_results <- rank_global    
  }
  
  rownames(Global_results) <- c()
  
  RWRM_ranking <- list(RWRM_Results = Global_results,Seed_Nodes = Seeds)
  
  class(RWRM_ranking) <- "RWRM_Results"
  return(RWRM_ranking)
}

#' @method print RWRM_Results
#' @export
print.RWRM_Results <- function(x,...)
{
  cat("Top 10 ranked Nodes:\n")
  print(head(x$RWRM_Results,10))
  cat("\nSeed Nodes used:\n")
  print(x$Seed_Nodes)
}

get.seed.scoresMultiplex <- function(Seeds,Number_Layers,tau) {
  
  Nr_Seeds <- length(Seeds)
  
  Seeds_Seeds_Scores <- rep(tau/Nr_Seeds,Nr_Seeds)
  Seed_Seeds_Layer_Labeled <- 
    paste0(rep(Seeds,Number_Layers),sep="_",rep(seq(Number_Layers), 
                                                length.out = Nr_Seeds*Number_Layers,each=Nr_Seeds))
  
  Seeds_Score <- data.frame(Seeds_ID = Seed_Seeds_Layer_Labeled,
                            Score = Seeds_Seeds_Scores, stringsAsFactors = FALSE)
  
  return(Seeds_Score)
}


geometric.mean <- function(Scores, L, N) {
  
  FinalScore <- numeric(length = N)
  
  for (i in seq_len(N)){
    FinalScore[i] <- prod(Scores[seq(from = i, to = N*L, by=N)])^(1/L)
  }
  
  return(FinalScore)
}

regular.mean <- function(Scores, L, N) {
  
  FinalScore <- numeric(length = N)
  
  for (i in seq_len(N)){
    FinalScore[i] <- mean(Scores[seq(from = i, to = N*L, by=N)])
  }
  
  return(FinalScore)
}

sumValues <- function(Scores, L, N) {
  
  FinalScore <- numeric(length = N)
  
  for (i in seq_len(N)){
    FinalScore[i] <- sum(Scores[seq(from = i, to = N*L, by=N)])
  }
  
  return(FinalScore)
}

isMultiplex <- function (x) 
{
  is(x, "Multiplex")
}

isMultiplexHet <- function (x) 
{
  is(x, "MultiplexHet")
}
