## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(fig.align = "center")


## ----eval=FALSE---------------------------------------------------------------
## # install the package via BioConductor
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## 
## BiocManager::install("netOmics")


## ----eval=FALSE---------------------------------------------------------------
## # install the package via github
## library(devtools)
## install_github("abodein/netOmics")


## ----eval=TRUE, message=FALSE-------------------------------------------------
# load the package
library(netOmics)


## ----eval=TRUE, message=FALSE-------------------------------------------------
# usefull packages to build this vignette
library(timeOmics)
library(tidyverse)
library(igraph)


## ----load_data----------------------------------------------------------------
# load data
data("hmp_T2D")


## ----timeOmics_1, eval=FALSE--------------------------------------------------
## # not evaluated in this vignette
## 
## #1 filter fold-change
## remove.low.cv <- function(X, cutoff = 0.5){
##     # var.coef
##     cv <- unlist(lapply(as.data.frame(X),
##                     function(x) abs(sd(x, na.rm = TRUE)/mean(x, na.rm= TRUE))))
##     return(X[,cv > cutoff])
## }
## fc.threshold <- list("RNA"= 1.5, "CLINICAL"=0.2, "GUT"=1.5, "METAB"=1.5,
##                      "PROT" = 1.5, "CYTO" = 1)
## 
## # --> hmp_T2D$raw
## data.filter <- imap(raw, ~{remove.low.cv(.x, cutoff = fc.threshold[[.y]])})
## 
## #2 scale
## data <- lapply(data.filter, function(x) log(x+1))
## # --> hmp_T2D$data
## 
## 
## #3 modelling
## lmms.func <- function(X){
##     time <- rownames(X) %>% str_split("_") %>%
##       map_chr(~.x[[2]]) %>% as.numeric()
##     lmms.output <- lmms::lmmSpline(data = X, time = time,
##                                    sampleID = rownames(X), deri = FALSE,
##                                    basis = "p-spline", numCores = 4,
##                                    keepModels = TRUE)
##     return(lmms.output)
## }
## data.modelled <- lapply(data, function(x) lmms.func(x))
## 
## # 4 clustering
## block.res <- block.pls(data.modelled, indY = 1, ncomp = 1)
## getCluster.res <- getCluster(block.res)
## # --> hmp_T2D$getCluster.res
## 
## 
## # 5 signature
## list.keepX <- list("CLINICAL" = 4, "CYTO" = 3, "GUT" = 10, "METAB" = 3,
##                    "PROT" = 2,"RNA" = 34)
## sparse.block.res  <- block.spls(data.modelled, indY = 1, ncomp = 1, scale =TRUE,
##                                 keepX =list.keepX)
## getCluster.sparse.res <- getCluster(sparse.block.res)
## # --> hmp_T2D$getCluster.sparse.res


## ----timeOmics_2--------------------------------------------------------------
# clustering results
cluster.info <- hmp_T2D$getCluster.res


## ----graph.rna, warning=FALSE-------------------------------------------------
cluster.info.RNA <- timeOmics::getCluster(cluster.info, user.block = "RNA")
graph.rna <- get_grn(X = hmp_T2D$data$RNA, cluster = cluster.info.RNA)

# to get info about the network
get_graph_stats(graph.rna)


## ----PROT_graph, warning=FALSE------------------------------------------------
# Utility function to get the molecules by cluster
get_list_mol_cluster <- function(cluster.info, user.block){
  require(timeOmics)
    tmp <- timeOmics::getCluster(cluster.info, user.block) 
    res <- tmp %>% split(.$cluster) %>% 
        lapply(function(x) x$molecule)
    res[["All"]] <- tmp$molecule
    return(res)
}

cluster.info.prot <- get_list_mol_cluster(cluster.info, user.block = 'PROT')
graph.prot <-  get_interaction_from_database(X = cluster.info.prot, 
                                             db = hmp_T2D$interaction.biogrid, 
                                             type = "PROT", user.ego = TRUE)
# get_graph_stats(graph.prot)


## ----GUT_graph, eval = FALSE--------------------------------------------------
## # not evaluated in this vignette
## library(SpiecEasi)
## 
## get_sparcc_graph <- function(X, threshold = 0.3){
##     res.sparcc <- sparcc(data = X)
##     sparcc.graph <- abs(res.sparcc$Cor) >= threshold
##     colnames(sparcc.graph) <-  colnames(X)
##     rownames(sparcc.graph) <-  colnames(X)
##     res.graph <- graph_from_adjacency_matrix(sparcc.graph,
##                                              mode = "undirected") %>% simplify
##     return(res.graph)
## }
## 
## gut_list <- get_list_mol_cluster(cluster.info, user.block = 'GUT')
## 
## graph.gut <- list()
## graph.gut[["All"]] <- get_sparcc_graph(hmp_T2D$raw$GUT, threshold = 0.3)
## graph.gut[["1"]] <- get_sparcc_graph(hmp_T2D$raw$GUT %>%
##                                        dplyr::select(gut_list[["1"]]),
##                                      threshold = 0.3)
## graph.gut[["-1"]] <- get_sparcc_graph(hmp_T2D$raw$GUT %>%
##                                         dplyr::select(gut_list[["-1"]]),
##                                       threshold = 0.3)
## class(graph.gut) <- "list.igraph"


## ----GUT----------------------------------------------------------------------
graph.gut <- hmp_T2D$graph.gut
# get_graph_stats(graph.gut)


## ----CYTO_graph, warning=FALSE------------------------------------------------
# CYTO -> from database (biogrid)
cyto_list = get_list_mol_cluster(cluster.info = cluster.info, 
                                 user.block = "CYTO")
graph.cyto <-  get_interaction_from_database(X = cyto_list,
                                             db = hmp_T2D$interaction.biogrid, 
                                             type = "CYTO", user.ego = TRUE)
# get_graph_stats(graph.cyto)

# METAB -> inference
cluster.info.metab <-  timeOmics::getCluster(X = cluster.info, 
                                             user.block = "METAB")
graph.metab <-  get_grn(X = hmp_T2D$data$METAB, 
                        cluster = cluster.info.metab)
# get_graph_stats(graph.metab)

# CLINICAL -> inference
cluster.info.clinical <- timeOmics::getCluster(X = cluster.info, 
                                               user.block = 'CLINICAL')
graph.clinical <- get_grn(X = hmp_T2D$data$CLINICAL,
                          cluster = cluster.info.clinical)
# get_graph_stats(graph.clinical)


## ----merged_0-----------------------------------------------------------------
full.graph <- combine_layers(graph1 = graph.rna, graph2 = graph.prot)
full.graph <- combine_layers(graph1 = full.graph, graph2 = graph.cyto)

full.graph <- combine_layers(graph1 = full.graph,
                             graph2 = hmp_T2D$interaction.TF)
# get_graph_stats(full.graph)


## ----merged_1_gut, warning=FALSE----------------------------------------------
all_data <- reduce(hmp_T2D$data, cbind)

# omic = gut
gut_list <- get_list_mol_cluster(cluster.info, user.block = "GUT")
omic_data <- lapply(gut_list, function(x)dplyr::select(hmp_T2D$data$GUT, x))

# other data = "RNA", "PROT", "CYTO"
other_data_list <- get_list_mol_cluster(cluster.info,
                                        user.block = c("RNA", "PROT", "CYTO"))
other_data <- lapply(other_data_list, function(x)dplyr::select(all_data, x))

# get interaction between gut data and other data
interaction_df_gut <- get_interaction_from_correlation(X = omic_data,
                                                       Y = other_data,
                                                       threshold = 0.99)

# and merge with full graph
full.graph <- combine_layers(graph1 = full.graph,
                             graph2 = graph.gut,
                             interaction.df = interaction_df_gut$All)


## ----merged_2_clinical, warning=FALSE-----------------------------------------
# omic =  Clinical
clinical_list <- get_list_mol_cluster(cluster.info, user.block = "CLINICAL")
omic_data <- lapply(clinical_list, 
                    function(x)dplyr::select(hmp_T2D$data$CLINICAL, x))

# other data = "RNA", "PROT", "CYTO", "GUT"
other_data_list <- get_list_mol_cluster(cluster.info,
                                        user.block = c("RNA", "PROT", 
                                                       "CYTO", "GUT"))
other_data <- lapply(other_data_list, function(x)dplyr::select(all_data, x))


# get interaction between gut data and other data
interaction_df_clinical <- get_interaction_from_correlation(X = omic_data
                                                            , Y = other_data,
                                                            threshold = 0.99)

# and merge with full graph
full.graph <- combine_layers(graph1 = full.graph,
                             graph2 = graph.clinical, 
                             interaction.df = interaction_df_clinical$All)


## ----merged_3_metab, warning=FALSE--------------------------------------------
# omic =  Metab
metab_list <- get_list_mol_cluster(cluster.info, user.block = "METAB")
omic_data <- lapply(metab_list, function(x)dplyr::select(hmp_T2D$data$METAB, x))

# other data = "RNA", "PROT", "CYTO", "GUT", "CLINICAL"
other_data_list <- get_list_mol_cluster(cluster.info,
                                        user.block = c("RNA", "PROT", "CYTO", 
                                                       "GUT", "CLINICAL"))
other_data <- lapply(other_data_list, function(x)dplyr::select(all_data, x))

# get interaction between gut data and other data
interaction_df_metab <- get_interaction_from_correlation(X = omic_data,
                                                         Y = other_data, 
                                                         threshold = 0.99)

# and merge with full graph
full.graph <- combine_layers(graph1 = full.graph, 
                             graph2 = graph.metab, 
                             interaction.df = interaction_df_metab$All)


## -----------------------------------------------------------------------------
# ORA by cluster/All
mol_ora <- get_list_mol_cluster(cluster.info, 
                                user.block = c("RNA", "PROT", "CYTO"))

# get ORA interaction graph by cluster
graph.go <- get_interaction_from_ORA(query = mol_ora,
                                     sources = "GO",
                                     organism = "hsapiens",
                                     signif.value = TRUE)

# merge
full.graph <- combine_layers(graph1 = full.graph, graph2 = graph.go)


## -----------------------------------------------------------------------------
# medlineRanker -> database
medlineranker.res.df <- hmp_T2D$medlineranker.res.df %>% 
  dplyr::select(Disease, symbol) %>% 
  set_names(c("from", "to"))
  
mol_list <-  get_list_mol_cluster(cluster.info = cluster.info,
                                  user.block = c("RNA", "PROT", "CYTO"))
graph.medlineranker <-  get_interaction_from_database(X = mol_list,
                                                      db = medlineranker.res.df, 
                                                      type = "Disease",
                                                      user.ego = TRUE)
# get_graph_stats(graph.medlineranker)

# merging
full.graph <- combine_layers(graph1 = full.graph, graph2 = graph.medlineranker)


## -----------------------------------------------------------------------------
# graph cleaning
graph_cleaning <- function(X, cluster.info){
    # no reusability
    X <- igraph::simplify(X)
    va <- vertex_attr(X)
    viewed_mol <- c()
    for(omic in unique(cluster.info$block)){
        mol <- intersect(cluster.info %>% dplyr::filter(.$block == omic) %>%
                           pull(molecule), V(X)$name)
        viewed_mol <- c(viewed_mol, mol)
        X <- set_vertex_attr(graph = X, 
                             name = "type", 
                             index = mol, 
                             value = omic)
        X <- set_vertex_attr(graph = X, 
                             name = "mode",
                             index = mol,
                             value = "core")
    }
    # add medline ranker and go
    mol <- intersect(map(graph.go, ~ as_data_frame(.x)$to) %>%
                       unlist %>% unique(), V(X)$name) # only GO terms
    viewed_mol <- c(viewed_mol, mol)
    X <- set_vertex_attr(graph = X, name = "type", index = mol, value = "GO")
    X <- set_vertex_attr(graph = X, name = "mode", 
                         index = mol, value = "extended")
    
    mol <- intersect(as.character(medlineranker.res.df$from), V(X)$name)
    viewed_mol <- c(viewed_mol, mol)
    X <- set_vertex_attr(graph = X, name = "type",
                         index = mol, value = "Disease")
    X <- set_vertex_attr(graph = X, name = "mode",
                         index = mol, value = "extended")
    
    other_mol <- setdiff(V(X), viewed_mol)
    if(!is_empty(other_mol)){
        X <- set_vertex_attr(graph = X, name = "mode",
                             index = other_mol, value = "extended")
    }
    X <- set_vertex_attr(graph = X, name = "mode", 
                         index = intersect(cluster.info$molecule, V(X)$name), 
                         value = "core")
    
    # signature
    mol <-  intersect(V(X)$name, hmp_T2D$getCluster.sparse.res$molecule)
    X <- set_vertex_attr(graph = X, name = "sparse", index = mol, value = TRUE)
    mol <-  setdiff(V(X)$name, hmp_T2D$getCluster.sparse.res$molecule)
    X <- set_vertex_attr(graph = X, name = "sparse", index = mol, value = FALSE)
    
    return(X)
}


## -----------------------------------------------------------------------------
FULL <- lapply(full.graph, function(x) graph_cleaning(x, cluster.info))
get_graph_stats(FULL)


## ----eval = FALSE-------------------------------------------------------------
## # degree analysis
## d <- degree(FULL$All)
## hist(d)
## d[max(d)]
## 
## # modularity # Warnings: can take several minutes
## res.mod <- walktrap.community(FULL$All)
## # ...
## 
## # modularity
## sp <- shortest.paths(FULL$All)


## -----------------------------------------------------------------------------
# seeds = all vertices -> takes 5 minutes to run on regular computer
# seeds <- V(FULL$All)$name
# rwr_res <- random_walk_restart(FULL, seeds)

# seed = some GO terms
seeds <- head(V(FULL$All)$name[V(FULL$All)$type == "GO"])
rwr_res <- random_walk_restart(FULL, seeds)


## -----------------------------------------------------------------------------
rwr_type_k15 <- rwr_find_seeds_between_attributes(X = rwr_res, 
                                                  attribute = "type", k = 15)

# a summary plot function
summary_plot_rwr_attributes(rwr_type_k15)
summary_plot_rwr_attributes(rwr_type_k15$All)


## -----------------------------------------------------------------------------
rwr_type_k15 <- rwr_find_seeds_between_attributes(X = rwr_res$All, 
                                                  attribute = "cluster", k = 15)
summary_plot_rwr_attributes(rwr_type_k15)


## -----------------------------------------------------------------------------
sub_res <- rwr_type_k15$`GO:0005737`
sub <- plot_rwr_subnetwork(sub_res, legend = TRUE, plot = TRUE)


## -----------------------------------------------------------------------------
rwr_res <- random_walk_restart(FULL$All, seed = "ZNF263")

# closest GO term
rwr_find_closest_type(rwr_res, seed = "ZNF263", attribute = "type", 
                      value = "GO", top = 5)

# closest Disease
rwr_find_closest_type(rwr_res, seed = "ZNF263", attribute = "type", 
                      value = "Disease", top = 5)

# closest nodes with an attribute "cluster" and the value "-1"
rwr_find_closest_type(rwr_res, seed = "ZNF263", attribute = "cluster",
                      value = "-1", top = 5)


## ----eval = FALSE-------------------------------------------------------------
## seeds <- V(FULL$All)$name[V(FULL$All)$type %in% c("GO", "Disease")]


## -----------------------------------------------------------------------------
sessionInfo()

