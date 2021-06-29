# General database setup

library(tidyverse) 
library(igraph)

########################################
##               BIOGRID              ##
########################################
 
# Release Version: 4.3.196
# File: BIOGRID-ALL-4.3.196.tab3.zip 
biogrid.tab <- read_tsv("./data-raw/databases/BIOGRID-ALL-4.3.196.tab3.zip")

biogrid.tab <- read_tsv("/home/antoine/Documents/TO2/netOmics-case-studies/HeLa_Cell_Cycling/data/BIOGRID-ALL-3.5.187.tab3.txt")
biogrid.tab.filtered <- biogrid.tab %>% filter(.$`Organism Interactor A` == 9606 & .$`Organism Interactor B` == 9606) %>%
#    dplyr::select("Official Symbol Interactor A", "Official Symbol Interactor B") %>%
    dplyr::select("SWISS-PROT Accessions Interactor A", "SWISS-PROT Accessions Interactor B") %>%
    unique %>%
    set_names(c("from", "to")) 

human.BIOGRID <- graph_from_data_frame(biogrid.tab.filtered, directed = FALSE) %>% simplify()
usethis::use_data(human.BIOGRID)

########################################
##               TF2DNA              ##
########################################

### bash 
# http://fiserlab.org/pscan_files.tar.gz
# for i in `ls **/*.pscan`; do cut -f1,2 $i | tail -n+2; done > TF2DNA_exp.txt
# sort TF2DNA_exp.txt | uniq > TF2DNA_exp_uniq.txt
###

### Gene SYMBOL
TF2DNA <- read_tsv("~/work/WORK:analysis/human_tf2DNA/TF2DNA_exp_uniq.txt", col_names = F)
names(TF2DNA) <- c("TF","Target")

## TTRUST
TTRUST <- read_tsv("~/Documents/TO2/LOREAL/Data/Proteomic/trrust_rawdata.human.tsv", col_names = F)
colnames(TTRUST) <- c("TF", "Target", "Regulation", 'PMID')
TFome <- rbind(TF2DNA, dplyr::select(TTRUST, c("TF", "Target")))

saveRDS(TFome, file = "~/Documents/TO2/TFome.Rds")


library(org.Hs.eg.db)

AnnotationDbi::select(org.Hs.eg.db, keytype = "UNIPROT", keys = V()$name, 
                                  columns = c("SYMBOL")) %>%
    unique %>% na.omit


# use TF(uniprot) -> Target(Symbol)
TFome.filtered <- TFome %>% left_join(conv.ego, by = c("TF"="SYMBOL")) %>% na.omit %>%
    dplyr::rename("TF.uniprot" = "UNIPROT") %>%
    dplyr::select(TF.uniprot, Target) %>% # UNIPROT -> SYMBOL
    unique %>%
    set_names(c("Prot", "Gene")) %>%
    mutate(from = Prot, to = Gene)
# mutate(interaction_type = "TF")
coding.filtered <- conv.ego %>% set_names(c("Prot", "Gene")) %>%
    mutate(from = Gene, to = Prot)
#%>%
#mutate(interaction_type = "coding")

