
library(tidyverse)
library(timeOmics)
library(lubridate)
library(lmms)

rm(list=ls())

# RData from:
# Sailani, M. R., Metwally, A. A., Zhou, W., Rose, S. M. S. F., Ahadi, S., Contrepois, K., ... & Snyder, M. P. (2020). 
# Deep longitudinal multiomics profiling reveals two biological seasonal patterns in California. Nature communications, 11(1), 1-12.
load("/home/antoine/Documents/timeomics_analysis/HMP_seasoning/Multi_Omics_Seasonal.RData")

# 0. DATA CLEANING

Gut_annotation_colData <- Gut_annotation_colData %>%
    mutate(YMD = lubridate::dmy(IRIS)) %>%
    mutate(Date = IRIS) %>%
    mutate(Time = yday(YMD)) %>%
    mutate(omics = "Gut") %>% dplyr::select(-Date, -BMI, -IRIS)

list_lab <- list("RNA" = RNA_annotation_colData,
                 "Metabo" = Metabolomics_annotation_colData,
                 #"Gut" =Gut_annotation_colData,
                 "Clinical" = Clinical_labs_annotation_colData)

list_lab_df <- imap_dfr(list_lab, ~{.x %>% 
        mutate("omics" = .y) %>% 
        mutate(YMD = lubridate::ymd(as.Date(Date))) %>% 
        dplyr::select(-Date)
})

IRIS_BMI <- list_lab_df[c(1:3)] %>% unique() 
IRIS_only <- IRIS_BMI %>% dplyr::select(-BMI) %>% unique %>% filter(!is.na(IRIS), !is.na(SubjectID)) %>%
    unique
IRIS_1 <- IRIS_only %>% group_by(SubjectID) %>%
    dplyr::summarise(N = n()) %>%
    filter(N == 1) %>% pull(SubjectID) %>% as.character()
IRIS_only <- IRIS_only %>% filter(SubjectID %in% IRIS_1)

# 1. DATA PREPARATION
# GUT
###########################
GUT_sample <- Gut_annotation_colData %>% 
    mutate(Year = ifelse(year(YMD) < 2000, year(YMD) +2000, year(YMD))) %>%
    left_join(IRIS_only) %>%
    mutate(SampleID = paste0(SubjectID, "_", Time, "_",  Year, "_", IRIS, "_", rownames(.)))

GUT <- gut_df_Data
rownames(GUT) <- GUT_sample$SampleID

# CLINICAL
###########################
CLINICAL_sample <- Clinical_labs_annotation_colData %>% 
    mutate(Year = ifelse(year(Date) < 2000, year(Date) +2000, year(Date))) %>%
    left_join(IRIS_only) %>%
    mutate(SampleID = paste0(SubjectID, "_", Time, "_", Year, "_", IRIS, "_", rownames(.))) 

CLINICAL <- clinical_labs_Data
rownames(CLINICAL) <- CLINICAL_sample$SampleID
index.na <- CLINICAL %>% lapply(function(x) is.na(x) %>% sum) %>% unlist
CLINICAL <- CLINICAL[index.na<=11] %>% na.omit

# RNA
###########################
RNA_sample <- RNA_annotation_colData %>% 
    mutate(Year = ifelse(year(Date) < 2000, year(Date) +2000, year(Date))) %>%
    left_join(IRIS_only) %>%
    mutate(SampleID = paste0(SubjectID, "_", Time, "_", Year, "_", IRIS, "_", rownames(.))) 

index.NA <- (!is.na(RNA_annotation_colData$Time) & !is.na(RNA_annotation_colData$Time))
RNA_sample <- RNA_sample[index.NA,]
RNA <- RNA_df_Data[index.NA,]
rownames(RNA) <- RNA_sample$SampleID

# NOSE
###########################
Nasal_sample <- Nasal_annotation_colData %>% 
    mutate(Year = ifelse(year(Date) < 2000, year(Date) +2000, year(Date))) %>%
    left_join(IRIS_only) %>%
    mutate(SampleID = paste0(SubjectID, "_", Time, "_",  Year, "_", IRIS, "_", rownames(.)))

NASAL <- Nasal_df_Data
rownames(NASAL) <- Nasal_sample$SampleID

# PROTEIN
###########################
PROT_sample <- Proteomics_annotation_colData %>% 
    mutate(Year = ifelse(year(Date) < 2000, year(Date) +2000, year(Date))) %>%
    left_join(IRIS_only) %>%
    mutate(SampleID = paste0(SubjectID, "_", Time, "_",  Year, "_", IRIS, "_", rownames(.)))

PROT <- Proteomics_df_Data
rownames(PROT) <- PROT_sample$SampleID

# METABOLITE
###########################
METAB_sample <- Metabolomics_annotation_colData %>% 
    mutate(Year = ifelse(year(Date) < 2000, year(Date) +2000, year(Date))) %>%
    left_join(IRIS_only) %>%
    mutate(SampleID = paste0(SubjectID, "_", Time, "_",  Year, "_", IRIS, "_", rownames(.)))

METAB <- Metabolomics_df_Data
rownames(METAB) <- METAB_sample$SampleID

# CYTOKINE
###########################
CYTO_sample <- Cytokines_annotation_colData %>% 
    mutate(Year = ifelse(year(Date) < 2000, year(Date) +2000, year(Date))) %>%
    left_join(IRIS_only) %>%
    mutate(SampleID = paste0(SubjectID, "_", Time, "_",  Year, "_", IRIS, "_", rownames(.)))

CYTO <- Cytokines_df_Data
rownames(CYTO) <- CYTO_sample$SampleID

############################

# DATA: only RNA/CLINICAL/GUT/METAB
# split by IR/IS
DATA <- list("RNA.IR" = RNA[str_split(rownames(RNA),"_") %>% map_chr(~.x[[4]]) == "IR",],
             "GUT.IR" = GUT[str_split(rownames(GUT),"_") %>% map_chr(~.x[[4]]) == "IR",],
             
             "CLINICAL.IR" = CLINICAL[str_split(rownames(CLINICAL),"_") %>% map_chr(~.x[[4]]) == "IR",],
             "RNA.IS" = RNA[str_split(rownames(RNA),"_") %>% map_chr(~.x[[4]]) == "IS",],
             
             "GUT.IS" = GUT[str_split(rownames(GUT),"_") %>% map_chr(~.x[[4]]) == "IS",],
             "CLINICAL.IS" = CLINICAL[str_split(rownames(CLINICAL),"_") %>% map_chr(~.x[[4]]) == "IS",],
             
             "METAB.IR" = METAB[str_split(rownames(METAB),"_") %>% map_chr(~.x[[4]]) == "IR",],
             "METAB.IS" = METAB[str_split(rownames(METAB),"_") %>% map_chr(~.x[[4]]) == "IS",],
             
             "PROT.IS" = PROT[str_split(rownames(PROT),"_") %>% map_chr(~.x[[4]]) == "IS",],
             "PROT.IR" = PROT[str_split(rownames(PROT),"_") %>% map_chr(~.x[[4]]) == "IR",],
             
             "CYTO.IS" = CYTO[str_split(rownames(CYTO),"_") %>% map_chr(~.x[[4]]) == "IS",],
             "CYTO.IR" = CYTO[str_split(rownames(CYTO),"_") %>% map_chr(~.x[[4]]) == "IR",]
) 

COMBINED <- list("RNA" = RNA, CLINICAL = CLINICAL, GUT = GUT, METAB = METAB, PROT = PROT, CYTO = CYTO)
save(DATA, COMBINED, file = "/home/antoine/Documents/timeomics_analysis/HMP_seasoning/netomics/RAW_DATA.RDA")
############################################################

stat_raw_data <- lapply(list(RNA=RNA, GUT=GUT, METAB=METAB, CLINICAL=CLINICAL, PROT = PROT, CYTO = CYTO), dim) %>%
    as.data.frame() %>% t %>% as.data.frame() %>%
    setNames(c("sample", "feature"))
lapply(list(RNA=RNA, GUT=GUT, METAB=METAB, CLINICAL=CLINICAL, PROT = PROT, CYTO = CYTO), function(x){
    rownames(x) %>% str_remove("_.*") %>% unique %>% length()}) %>% 
    as.data.frame() %>%  t %>% as.data.frame() %>% setNames("uniqueID") %>%
    rownames_to_column("omic") %>%
    left_join(stat_raw_data %>% rownames_to_column("omic")) %>% column_to_rownames("omic") %>% t %>%
    as.data.frame() %>% knitr::kable()

# 2. DATA FILTERING

# 1. coef. of var
cv.data <- lapply(DATA, function(X){
    unlist(lapply(as.data.frame(X), 
                  function(x) abs(sd(x, na.rm = TRUE)/mean(x, na.rm= TRUE))))
})

fc.data <- list("RNA.IR"= 1.5, "RNA.IS"=1.5,
                "CLINICAL.IR"=0, "CLINICAL.IS"=0,
                "GUT.IR"=1.5, "GUT.IS"=1.5,
                "METAB.IR"=1.5 , "METAB.IS"=1.5, 
                "PROT.IR"=1.5 , "PROT.IS"=1.5, 
                "CYTO.IR" = 1.5, "CYTO.IS"=1.5)

par(mfrow = c(6,2))
#for(i in c("RNA.IR","CLINICAL.IR", "GUT.IR", "METAB.IR", "RNA.IS","CLINICAL.IS", "GUT.IS", "METAB.IS", "PROT.IR", "PROT.IS", "CYTO.IR", "CYTO.IS")){
for(i in c("RNA.IR","RNA.IS", "CLINICAL.IR", "CLINICAL.IS", "GUT.IR","GUT.IS", "METAB.IS", "METAB.IR", "PROT.IR", "PROT.IS", "CYTO.IR", "CYTO.IS")){
    
    hist(cv.data[[i]], breaks = 20, main =i)
    abline(v = fc.data[[i]], col = "red")
    legend("topright", legend = paste0("FC = ",fc.data[[i]]), col = "red", lty = 1) 
}
par(mfrow = c(1,1))


# 2. Remove low cv features
remove.low.cv <- function(X, cutoff = 0.5){
    # var.coef
    cv <- unlist(lapply(as.data.frame(X), 
                        function(x) abs(sd(x, na.rm = TRUE)/mean(x, na.rm= TRUE))))
    return(X[,cv > cutoff])
}

DATA.filtered <- list("RNA.IR" = remove.low.cv(DATA$RNA.IR, 1.5),
                      "RNA.IS" = remove.low.cv(DATA$RNA.IS, 1.5),
                      "GUT.IR" = remove.low.cv(DATA$GUT.IR, 1.5),
                      "GUT.IS" = remove.low.cv(DATA$GUT.IS, 1.5),
                      # "CLINICAL.IR" = remove.low.cv(DATA$CLINICAL.IR, 0),
                      # "CLINICAL.IS" = remove.low.cv(DATA$CLINICAL.IS, 0),
                      "CLINICAL.IR" = DATA$CLINICAL.IR,
                      "CLINICAL.IS" = DATA$CLINICAL.IS,
                      
                      "METAB.IR" = remove.low.cv(DATA$METAB.IR, 1.5),
                      "METAB.IS" = remove.low.cv(DATA$METAB.IS, 1.5), 
                      "PROT.IS" = remove.low.cv(DATA$PROT.IS, 1.5),
                      "PROT.IR" = remove.low.cv(DATA$PROT.IR, 1.5),
                      "CYTO.IS" = remove.low.cv(DATA$CYTO.IS, 1),
                      "CYTO.IR" = remove.low.cv(DATA$CYTO.IR, 1))
lapply(DATA.filtered, dim)

# 3. scale filtered value (log, scale, CLR)

# scale for OTU
norm_OTU <- function(DF, AR = F){
    DF <- DF + 0.0001
    
    data.TSS.clr = mixOmics::logratio.transfo(DF, logratio = 'CLR')
    
    # reconstrcuct dataframe
    data.good <- as.data.frame(matrix(ncol = ncol(data.TSS.clr), 
                                      nrow = nrow( data.TSS.clr)))
    rownames(data.good) <- rownames(data.TSS.clr)
    colnames(data.good) <- colnames(data.TSS.clr)
    for( i in c(1:nrow(data.TSS.clr))){
        for( j in c(1:ncol(data.TSS.clr))){
            data.good[i,j] <- data.TSS.clr[i,j]
        }
    }
    return(data.good)
}


DATA.filtered.scale <- list(
    "RNA.IR" = log(DATA.filtered$RNA.IR + 1) %>% scale,
    "RNA.IS" = log(DATA.filtered$RNA.IS + 1) %>% scale,
    
    "CLINICAL.IR" = log(DATA.filtered$CLINICAL.IR +1)%>% scale,
    "CLINICAL.IS" = log(DATA.filtered$CLINICAL.IS +1)%>% scale,
    
    "GUT.IR" = norm_OTU(DATA.filtered$GUT.IR),
    "GUT.IS" = norm_OTU(DATA.filtered$GUT.IS),
    
    "METAB.IR" = log(DATA.filtered$METAB.IR +1)%>% scale,
    "METAB.IS" = log(DATA.filtered$METAB.IS +1)%>% scale,
    
    # "PROT.IR" = log(DATA.filtered$PROT.IR +1)%>% scale,
    # "PROT.IS" = log(DATA.filtered$PROT.IS +1)%>% scale
    
    "PROT.IR" = DATA.filtered$PROT.IR,  # already scale
    "PROT.IS" = DATA.filtered$PROT.IS,
    
    "CYTO.IR" = log(DATA.filtered$CYTO.IR +1),  
    "CYTO.IS" = log(DATA.filtered$CYTO.IS +1)
    
)

lapply(DATA.filtered, dim) %>%
    as.data.frame() %>% t %>% as.data.frame() %>%
    setNames(c("sample", "feature")) %>%
    rownames_to_column("OMIC") %>%
    mutate(IRIS = str_extract(OMIC,"..$"), OMIC = str_remove(OMIC, "...$"))  %>%
    gather(meta, value, -c(OMIC, IRIS)) %>%
    spread(OMIC, value) %>% arrange(IRIS) %>%
    dplyr::select(IRIS, meta, RNA, GUT, METAB, CLINICAL, PROT, CYTO)

save(DATA.filtered.scale, DATA.filtered, file = "/home/antoine/Documents/timeomics_analysis/HMP_seasoning/netomics/DATA_FILTERED.RDA")
############################################################

fc.data.combined <- list("RNA"= 1.5, 
                         "CLINICAL"=0.2,
                         "GUT"=1.5,
                         "METAB"=1.5,
                         "PROT" = 1.5,
                         "CYTO" = 1)
cv.data.combined <- lapply(COMBINED, function(X){
    unlist(lapply(as.data.frame(X), 
                  function(x) abs(sd(x, na.rm = TRUE)/mean(x, na.rm= TRUE))))
})
fc.color <- list("RNA"= color.mixo(4), 
                 "CLINICAL"=color.mixo(1),
                 "GUT"=color.mixo(2),
                 "METAB"=color.mixo(3),
                 "PROT"=color.mixo(5),
                 "CYTO" = color.mixo(6))

par(mfrow = c(3,2))
for(i in c("RNA","CLINICAL", "GUT", "METAB", "PROT", "CYTO")){
    hist(cv.data.combined[[i]], breaks = 20, main =i, xlab = paste0("Var. Coef. (", i, ")"), 
         col = fc.color[[i]])
    abline(v = fc.data.combined[[i]], col = "red")
    legend("topright", legend = paste0("CV = ",fc.data.combined[[i]]), col = "red", lty = 1) 
}
par(mfrow = c(1,1))

# 3. MODELLING

lmms.func <- function(X, mode = "p-spline"){
    time <- rownames(X) %>% str_split("_") %>% map_chr(~.x[[2]]) %>% as.numeric()
    lmms.output <- lmms::lmmSpline(data = X, time = time,
                                   sampleID = rownames(X), deri = FALSE,
                                   basis = mode, numCores = 4, 
                                   keepModels = TRUE)
    return(lmms.output)
}

# only one ID/Year
ID_u <- "ZLZNCLZ"
Year_u <- 2015

# ID_u <- "ZOZOW1T"
# Year_u <- 2015


tmp <-   imap_dfr(DATA.filtered.scale, ~{
    .x %>% as.data.frame() %>% rownames_to_column("sample") %>%
        gather(feature, value, -sample) %>%
        mutate(ID = str_split(sample, "_") %>% map_chr(~.x[[1]])) %>%
        mutate(year = str_split(sample, "_") %>% map_chr(~.x[[3]])) %>%
        mutate(OMIC = str_remove(.y, "...$")) %>% 
        mutate(IRIS = str_extract(.y, "..$"))
}) 
tmp %>%  dplyr::select(-feature, -value) %>% unique() %>% na.omit %>% 
    group_by(ID, year, OMIC, IRIS) %>% summarise(N = n()) %>% spread(OMIC, N) %>% 
    na.omit() %>% split(.$year)
tmp %>%  dplyr::select(-feature, -value) %>% unique() %>% na.omit %>% 
    group_by(ID, year, OMIC, IRIS) %>% summarise(N = n()) %>% spread(OMIC, N) %>% 
    filter(ID == ID_u) %>% dplyr::select(-IRIS) %>% na.omit

# just a filter to get only the selected ID/Year
DATA.GOOD <- imap_dfr(DATA.filtered.scale, ~{
    .x %>% as.data.frame() %>% rownames_to_column("sample") %>%
        gather(feature, value, -sample) %>%
        mutate(ID = str_split(sample, "_") %>% map_chr(~.x[[1]])) %>%
        mutate(year = str_split(sample, "_") %>% map_chr(~.x[[3]])) %>%
        mutate(OMIC = str_remove(.y, "...$")) %>%
        dplyr::filter(ID == ID_u) %>% 
        dplyr::filter(year == Year_u)
}) %>% split(.$OMIC) %>%
    purrr::map(~{
        .x %>% 
            dplyr::select(sample, feature, value) %>%
            spread(feature, value) %>% 
            column_to_rownames("sample")
    })


# only 1 year
MODELLED <- lapply(DATA.GOOD, function(x) lmms.func(x))

MODELLED %>% lapply(function(x)x@predSpline %>% dim) %>%
    as.data.frame() %>% t %>% as.data.frame() %>%
    setNames(c("sample", "feature")) %>% t %>%as.data.frame() %>%
    dplyr::select(RNA, GUT, METAB, CLINICAL, PROT)

MODELLED %>% imap_dfr(~.x@modelsUsed %>% table %>% as.data.frame  %>%
                          column_to_rownames(".") %>% t %>% as.data.frame %>% 
                          mutate(omic = .y)) %>% remove_rownames() %>%  
    column_to_rownames('omic') %>% t %>%
    as.data.frame() %>% dplyr::select(RNA, GUT, METAB, CLINICAL, PROT)

# 4. STRAIGHT LINE FILTERING

filterlmms.func <- function(modelled.data, lmms.output){
    time = modelled.data %>% rownames() %>% str_split("_") %>% map_chr(~.x[[2]]) %>% as.numeric()
    #time = rownames(modelled.data) %>% as.numeric()
    filter.res <- lmms.filter.lines(data = modelled.data,
                                    lmms.obj = lmms.output, time = time,
                                    homoskedasticity.cutoff=0.05)$filtered
}

FILTER <- lapply(names(DATA.GOOD), function(x) filterlmms.func(modelled.data = DATA.GOOD[[x]], lmms.output = MODELLED[[x]]))
names(FILTER) <- names(MODELLED)

FILTER %>% lapply(dim) %>%
    as.data.frame() %>% t %>% as.data.frame() %>%
    setNames(c("sample", "feature")) %>%
    t %>%as.data.frame() %>%
    dplyr::select(RNA, GUT, METAB, CLINICAL)

FINAL.FILTER <- FILTER[c("CLINICAL", "GUT", "METAB", "RNA", "PROT")]
rownames(FINAL.FILTER[["GUT"]]) <- rownames(FINAL.FILTER[["RNA"]]) # change 86 par 85

DATA.LMMS <- lapply(MODELLED, function(x)x@predSpline %>% t %>% as.data.frame)
rownames(DATA.LMMS[["GUT"]]) <- rownames(DATA.LMMS[["RNA"]]) # change 86 par 85

save(FINAL.FILTER, MODELLED, DATA.GOOD, DATA.LMMS, file = "/home/antoine/Documents/timeomics_analysis/HMP_seasoning/netomics/LMMS.RDA")

lapply(DATA.LMMS, dim)
# 5. MULTI-OMICS CLUSTERING

block.res <- block.pls(DATA.LMMS, indY = 1, ncomp = 5)
getNcomp.res <- getNcomp(block.res, X = DATA.LMMS, indY = 1)

# block.res <- block.pls(FINAL.FILTER, indY = 1, ncomp = 3)
# getNcomp.res <- getNcomp(block.res, X = FINAL.FILTER, indY = 1)

plot(getNcomp.res)

# ncomp = 2
block.res <- block.pls(DATA.LMMS, indY = 1, ncomp = 1, scale =FALSE) 

# block.res <- block.pls(FINAL.FILTER, indY = 1, ncomp = 1, scale =FALSE) 

plotLong(object = block.res, title = "Block-PLS Clusters, scale = TRUE", legend = TRUE)

getCluster(block.res) %>% group_by(block, cluster) %>% summarise(N = n()) %>%
    spread(block, N) %>%
    dplyr::select(cluster, RNA, GUT, METAB, CLINICAL)

save(block.res, file = "/home/antoine/Documents/timeomics_analysis/HMP_seasoning/netomics/timeomics_res_block.rda")


# elagage
test.list.keepX <- list(
    "CLINICAL" = seq(2,8,by=1),
    "GUT" = seq(2,10,by=1),
    "METAB" = seq(2,9,by=1),
    "RNA" = seq(10,50,by=2)
)

tune.block.res <- tuneCluster.block.spls(X= FINAL.FILTER, indY = 1,
                                         test.list.keepX=test.list.keepX, 
                                         scale=FALSE, 
                                         mode = "canonical", ncomp = 1)
tune.block.res$choice.keepX 
final.block <- block.spls(FINAL.FILTER, indY = 1, ncomp = 1, scale =FALSE, 
                          keepX = tune.block.res$choice.keepX) 
plotLong(final.block, legend = TRUE)

getCluster(final.block) %>% group_by(block, cluster) %>% summarise(N = n()) %>%
    spread(block, N) %>%
    dplyr::select(cluster, RNA, GUT, METAB, CLINICAL)

library("openxlsx")
cluster_comp <- getCluster(final.block) %>% dplyr::select(molecule, block, cluster, comp, contribution) %>% 
    mutate(cluster = ifelse(cluster == -1, "Cluster 1", "Cluster 2")) %>%
    split(.$cluster)
write.xlsx(cluster_comp, file = "cluster_composition.xlsx")


# final data for netOmics package // shrink
#   DATA.GOOD -> ! individual, 7 timepoints but no modelisation
#   DATA RAW -> filter based on DATA.GOOD molecules (for OTU -> sparcc needs RAW)
# DATA$RNA.IS %>% rownames() %>% str_remove("_.*") %>% str_detect("ZLZNCLZ") %>% any()
# [1] TRUE

hmp_T2D <- list()
hmp_T2D$raw <- list() 
hmp_T2D$data <- list()
for(i in names(DATA.GOOD)){
    # hmp_diabetes$raw[[i]] <- DATA[[paste0(i, ".IS")]][str_detect(rownames(DATA[[paste0(i, ".IS")]]), "ZLZNCLZ"), colnames(DATA.GOOD[[i]])]
    hmp_T2D$raw[[i]] <- DATA[[paste0(i, ".IS")]][rownames(DATA.GOOD[[i]]), colnames(DATA.GOOD[[i]])]
    rownames(hmp_T2D$raw[[i]]) <- rownames(DATA.GOOD$RNA)
    hmp_T2D$raw[[i]] <- hmp_T2D$raw[[i]] %>% rownames_to_column("Rownames") %>% 
        mutate(Rownames = Rownames %>% str_split("_") %>% map_chr(~.x[[2]]) %>% as.numeric()) %>% 
        arrange(Rownames) %>% column_to_rownames(var = "Rownames")
    
    #data
    hmp_T2D$data[[i]] <- DATA.GOOD[[i]]
    rownames(hmp_T2D$data[[i]]) <- rownames(hmp_T2D$raw[[i]])
}

hmp_T2D$data <- DATA.GOOD
lmms.func <- function(X){
    # time <- rownames(X) %>% str_split("_") %>% map_chr(~.x[[2]]) %>% as.numeric()
    time <- rownames(X) %>% as.numeric()
    lmms.output <- lmms::lmmSpline(data = X, time = time,
                                   sampleID = rownames(X), deri = FALSE,
                                   basis = "p-spline", numCores = 4, 
                                   keepModels = TRUE)
    return(lmms.output)
}

MODELLED <- lapply(hmp_T2D$data, function(x) lmms.func(x))

MODELLED %>% imap_dfr(~.x@modelsUsed %>% table %>% as.data.frame  %>%
                          column_to_rownames(".") %>% t %>% as.data.frame %>% 
                          mutate(omic = .y)) %>% remove_rownames() %>%  
    column_to_rownames('omic') %>% t %>%
    as.data.frame()

data.MODELLED <-  lapply(MODELLED, function(x) x@predSpline %>% t)
rownames(data.MODELLED$GUT) <- rownames(data.MODELLED$RNA)


# .$timeOmics
# block.res.no.model  <- block.pls(hmp_T2D$data, indY = 1, ncomp = 5, scale =FALSE) 
# getNcomp.res <- getNcomp(block.res.no.model, X = hmp_T2D$data, indY = 1)
# plot(getNcomp.res)
block.res.no.model  <- block.pls(hmp_T2D$data, indY = 1, ncomp = 1, scale =TRUE) 
hmp_T2D$getCluster.res <- getCluster(block.res.no.model)

block.res.w.model  <- block.pls(data.MODELLED, indY = 1, ncomp = 5, scale =TRUE) 
getNcomp.res <- getNcomp(block.res.w.model, X = data.MODELLED, indY = 1)
block.res.w.model  <- block.pls(data.MODELLED, indY = 1, ncomp = 1, scale =TRUE) 

hmp_T2D$getCluster.res <- getCluster(block.res.w.model)

# .$sparse
test.list.keepX <- list(
    "CLINICAL" = seq(2,39,by=2),
    "CYTO" = seq(2,10,b=2),
    "GUT" = seq(2,50,by=2),
    "METAB" = seq(2,50,by=2),
    "PROT" = seq(2,30,b=2),
    "RNA" = seq(10,100,by=3))

# tune.block.res <- tuneCluster.block.spls(X= hmp_T2D$data, indY = 1,
#                                          test.list.keepX=test.list.keepX, 
#                                          scale=FALSE, 
#                                          mode = "canonical", ncomp = 1)
# too long

list.keepX <- list(
    "CLINICAL" = 4,
    "CYTO" = 3,
    "GUT" = 10,
    "METAB" = 3,
    "PROT" = 2,
    "RNA" = 34)

sparse.block.res.w.model  <- block.spls(data.MODELLED, indY = 1, ncomp = 1, scale =TRUE, keepX = list.keepX) 
plotLong(sparse.block.res.w.model, scale = TRUE, legend = TRUE)


hmp_T2D$getCluster.sparse.res <-  getCluster(sparse.block.res.w.model)
timeOmics::getSilhouette(sparse.block.res.w.model)
timeOmics::getSilhouette(block.res.w.model)

## Biogrid ## DATABASES
biogrid <- read_tsv("/home/antoine/Documents/TO2/netOmics-case-studies/HeLa_Cell_Cycling/data/BIOGRID-ALL-3.5.187.tab3.txt")
biogrid.filtered <- biogrid %>% dplyr::select("Official Symbol Interactor A", "Official Symbol Interactor B") %>% unique %>%
    set_names(c("from", "to"))

biogrid.filtered.tmp <- biogrid.filtered %>% filter(from %in% cluster.info$molecule | to %in% cluster.info$molecule)

hmp_T2D$interaction.biogrid <- biogrid.filtered.tmp

TFome <- readRDS( "~/Documents/TO2/TFome.Rds")
tf.1 <- TFome %>% filter(TF %in% hmp_T2D$getCluster.res$molecule | Target %in% hmp_T2D$getCluster.res$molecule) %>% 
    unlist() %>% unique
hmp_T2D$interaction.TF <- TFome %>% filter(TF %in% tf.1 | Target %in% tf.1)

hmp_T2D$interaction.TF <- TFome.igraph

hmp_T2D$interaction.TF <- TF.interact
usethis::use_data(hmp_T2D, overwrite = TRUE)


medlineranker.res <- read_tsv("/home/antoine/Documents/timeomics_analysis/HMP_seasoning/netomics/res_medlineranker_all.csv")
tmp <- dplyr::select(medlineranker.res, c("Disease", 'Gene symbols'))  %>% set_names(c("Disease", "symbol")) 
# separate_rows(tmp, Disease, symbol, sep = "|", convert = TRUE)

library(splitstackshape)
tmp[c(1,2),] %>% separate_rows(Disease, symbol)
tmp[c(2,3),] %>% 
    rownames_to_column("id") %>% 
    cSplit(., sep = "|", splitCols = 2:ncol(.), direction = "long", makeEqual = T) %>% 
    as_tibble() %>% 
    group_by(id) %>% 
    fill(2:ncol(.)) %>% 
    unique() %>% 
    ungroup(id) %>% 
    select(-id)

medlineranker.res.df <- tmp %>% rownames_to_column("id") %>% 
    cSplit(., sep = "|", splitCols = 3:ncol(.), direction = "long", makeEqual = T) %>% 
    as_tibble() %>% group_by(id) %>% fill(2:ncol(.)) %>% unique() %>% ungroup(id) %>% select(-id)

medlineranker.res.df <- left_join(medlineranker.res.df, medlineranker.res) %>% mutate(symbol = as.character(symbol))
hmp_T2D$medlineranker.res.df <- medlineranker.res.df 
usethis::use_data(hmp_T2D, overwrite = TRUE)