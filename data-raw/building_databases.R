# General database setup

library(tidyverse) 
library(igraph)

########################################
##               BIOGRID              ##
########################################
 
# Release Version: 4.3.196
# File: BIOGRID-ALL-4.3.196.tab3.zip 
biogrid.tab <- read_tsv("./data-raw/databases/BIOGRID-ALL-4.3.196.tab3.zip")
biogrid.tab.filtered <- biogrid.tab %>% 
    dplyr::select("Official Symbol Interactor A", "Official Symbol Interactor B") %>% 
    unique %>%
    set_names(c("Symbol A", "Symbol B"))

BIOGRID <- graph_from_data_frame(biogrid.tab.filtered, directed = FALSE) %>% simplify()

## TTRUST
