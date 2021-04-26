load("~/Documents/TO2/CS1/analysis/data_filtered.rda")

# from CS1
data.lmms <- list("mRNA" = mRNA.filter, "trans" = trans.filter, "prot"=prot.filter)
#cluster.info.updown.mean <- imap_dfr(data.lmms, ~{up_down_cluster(.x) %>% mutate(block = .y)})

cluster.info.updown.lmms <- timeOmics::getUpDownCluster(X = data.lmms)
cluster.info.updown.lmms <- timeOmics::getCluster(cluster.info.updown.lmms)

data.lmms
#usethis::use_data(cluster.info.updown.lmms, internal = FALSE)

HeLa <- list("raw" = list("mRNA" = mRNA.fc, "prot" = prot.fc, "tran" = trans.fc), 
             "modelled" = list("mRNA" = mRNA.filter, "trans" = trans.filter, "prot"=prot.filter),
             "getCluster" = cluster.info.updown.lmms)
usethis::use_data(HeLa, internal = FALSE, overwrite = TRUE)

pca.res <- pca(data.lmms$mRNA)
getCluster(pca.res)$molecule[46]

debug(timeOmics::getUpDownCluster)
