load("~/Documents/TO2/CS1/analysis/data_filtered.rda")

# from CS1
# data.lmms <- list("mRNA" = mRNA.filter, "trans" = trans.filter, "prot"=prot.filter)
# # #cluster.info.updown.mean <- imap_dfr(data.lmms, ~{up_down_cluster(.x) %>% mutate(block = .y)})
# # 
# cluster.info.updown.lmms <- timeOmics::getUpDownCluster(X = data.lmms)
# cluster.info.updown.lmms <- timeOmics::getCluster(cluster.info.updown.lmms)

#usethis::use_data(cluster.info.updown.lmms, internal = FALSE)



colnames(mRNA.filter)  <- colnames(mRNA.filter) %>% stringr::str_remove("mRNA_")
colnames(trans.filter)  <- colnames(trans.filter) %>% stringr::str_remove("trans_")
colnames(prot.filter)  <- colnames(prot.filter) %>% stringr::str_remove("prot_")

data.lmms <- list("mRNA" = mRNA.filter, "trans" = trans.filter, "prot"=prot.filter)
cluster.info.updown.lmms <- timeOmics::getUpDownCluster(X = data.lmms)
cluster.info.updown.lmms <- timeOmics::getCluster(cluster.info.updown.lmms)

colnames(mRNA.fc)  <- colnames(mRNA.fc) %>% stringr::str_remove("mRNA_")
colnames(trans.fc)  <- colnames(trans.fc) %>% stringr::str_remove("trans_")
colnames(prot.fc)  <- colnames(prot.fc) %>% stringr::str_remove("prot_")


HeLa <- list("raw" = list("mRNA" = mRNA.fc, "prot" = prot.fc, "tran" = trans.fc), 
              "modelled" = list("mRNA" = mRNA.filter, "trans" = trans.filter, "prot"=prot.filter),
              "getCluster" = cluster.info.updown.lmms)

# grn for mRNA only
cluster.mRNA <- timeOmics::getCluster(HeLa$getCluster, user.block = "mRNA")
grn.res <- get_grn(X = HeLa$raw$mRNA, cluster = cluster.mRNA, method = "aracne")
HeLa[["GRN"]] <- grn.res

# prot
data("human.BIOGRID") # biogrid v3.5.187

# only proteins
cluster.proteins <- timeOmics::getCluster(HeLa$getCluster, user.block = "prot") 

# proteins name by cluster
proteins <- lapply(split(cluster.proteins, f = cluster.proteins$cluster), function(x) x$molecule)
biogrid.res <- get_interaction_from_database(X = proteins, db = human.BIOGRID, type = "db", user.ego = TRUE)

biogrid.res[["Up_Down"]] <- igraph::make_empty_graph(directed = FALSE)
biogrid.res[["All"]] <- get_interaction_from_database(X = unlist(proteins), db = human.BIOGRID, type = "prot", user.ego = TRUE)
get_graph_stats(biogrid.res)
merged.graph <- combine_layers(graph1 = grn.res, graph2 = biogrid.res)

HeLa[["prot.graph"]] <- biogrid.res
HeLa[["grn.prot.graph"]] <- merged.graph

#vertex_attr(HeLa[["grn.prot.graph"]]$All)

# miss tfome

usethis::use_data(HeLa, internal = FALSE, overwrite = TRUE)

data("HeLa")

