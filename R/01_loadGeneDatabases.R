##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
# Updated on January 09, 2025
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_settings.R") # load baseDir <- "/WORKING_DIR"

dataDir <- file.path(baseDir, "GeneLists")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(stringr)
library(CellChat)
library(msigdbr)
library(GSEABase)

# Gene-sets from the mSigDB (assessed on Jan 09, 2024)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark <- data.frame(
        Database = "Hallmark",
        Category = gsub("HALLMARK_", "", hallmark$gs_name),
        Gene = hallmark$gene_symbol
)
table(hallmark$Category)

genesetFilename <- list.files(dataDir, pattern = "kegg_legacy", full.names = T)
geneset <- getGmt(con = genesetFilename)
genesetNames <- names(geneset)
kegg <- c()
for (genesetName in genesetNames) {
        df <- data.frame(
                Database = "KEGG",
                Category = gsub("KEGG_", "", genesetName), 
                Gene = unlist(geneIds(geneset[genesetName]))
        )
        kegg <- rbind(kegg, df)
}
table(kegg$Category)

genesetFilename <- list.files(dataDir, pattern = "3ca", full.names = T)
geneset <- getGmt(con = genesetFilename)
genesetNames <- names(geneset)
cancer <- c()
for (genesetName in genesetNames) {        
        df <- data.frame(
                Database = "CancerCellAtlas",
                Category = gsub("GAVISH_3CA_", "", genesetName), 
                Gene = unlist(geneIds(geneset[genesetName]))
        )
        cancer <- rbind(cancer, df)
}
table(cancer$Category)

oncogenic <- msigdbr(species = "Homo sapiens", category = "C6")
oncogenic <- data.frame(
        Database = "OncogenicSignature",
        Category = oncogenic$gs_name,
        Gene = oncogenic$gene_symbol
)
table(oncogenic$Category)

celltype <- msigdbr(species = "Homo sapiens", category = "C8")
celltype <- data.frame(
        Database = "CellType",
        Category = celltype$gs_name,
        Gene = celltype$gene_symbol
)
celltypeStats <- table(celltype$Category) # exclude too small/large gene-sets
include <- names(celltypeStats)[intersect(which(celltypeStats >= 10), which(celltypeStats < 1000))]
celltype <- celltype[which(celltype$Category %in% include),]

# Ligand-Receptor pairs from CellChat
cellchat <- CellChatDB.human$interaction

symbolsL <- strsplit(as.character( cellchat[, "ligand.symbol"]), ", ", fixed=TRUE)
dupCellchat <- data.frame(pathway_name = rep(cellchat[,"pathway_name"], sapply(symbolsL, length)), ligand.symbol = unlist(symbolsL))
ligands <- dupCellchat[!duplicated(dupCellchat), ]

cellchat1 <- data.frame(
        Database = "CellChat",
        Category = "Ligand",
        Gene = unique(ligands$ligand.symbol)
)

symbolsL <- strsplit(as.character( cellchat[, "receptor.symbol"]), ", ", fixed=TRUE)
dupCellchat <- data.frame(pathway_name = rep(cellchat[,"pathway_name"], sapply(symbolsL, length)), receptor.symbol = unlist(symbolsL))
receptors <- dupCellchat[!duplicated(dupCellchat), ]

cellchat2 <- data.frame(
        Database = "CellChat",
        Category = "Receptor",
        Gene = unique(receptors$receptor.symbol)
)

cellchat <- rbind(cellchat1, cellchat2)
table(cellchat$Category)
# Ligand Receptor 
#    786      720 

# GAD (Genetic Association Database)
fileNames <- list.files(dataDir, pattern = "GAD_", full.names = T)
gad <- c()
for (filename in fileNames) {
        category <- gsub("GAD_", "", gsub(".txt", "", basename(filename)))
        buff <- read.table(filename, header=F, stringsAsFactor=F)
        df <- data.frame(
                Database = "GAD",
                Category = category,
                Gene = buff$V1
        )
        gad <- rbind(gad, df)
}
table(gad$Ca)
#           AGING          CANCER  CARDIOVASCULAR  CHEMDEPENDENCY   DEVELOPMENTAL 
#              31             275             286              54              54 
#   HEMATOLOGICAL          IMMUNE       INFECTION       METABOLIC    NEUROLOGICAL 
#              40             313              99             370             167 
# NORMALVARIATION PHARMACOGENOMIC           PSYCH           RENAL    REPRODUCTION 
#              22              19             184              72              79 
#          VISION 
#              34 

# OMIM (Online Mendelian Inheritance in Man)
fileNames <- list.files(dataDir, pattern = "OMIM_", full.names = T)
omim <- c()
for (filename in fileNames) {
        category <- gsub("OMIM_", "", gsub(".txt", "", basename(filename)))
        buff <- read.table(filename, header=F, stringsAsFactor=F)
        df <- data.frame(
                Database = "OMIM",
                Category = category,
                Gene = buff$V1
        )
        omim <- rbind(omim, df)
}
table(omim$Category)
# DENOVO  DOMINANT_NEGATIVE HAPLOINSUFFICIENCY           RCESSIVE 
#   467                364                175                817 

# Drug-related, ADME and PharmGKB
fileNames <- list.files(dataDir, pattern = "DRUG_", full.names = T)
drug <- c()
for (filename in fileNames) {
        category <- gsub("DRUG_", "", gsub(".txt", "", basename(filename)))
        buff <- read.table(filename, header=F, stringsAsFactor=F)
        df <- data.frame(
                Database = "DRUG",
                Category = category,
                Gene = buff$V1
        )
        drug <- rbind(drug, df)
}
table(drug$Category)
# ADME PharmGKB 
#  299      367

# Protein Family: enzymes, carriers, transporters
fileNames <- list.files(dataDir, pattern = "FAMILY_", full.names = T)
family <- c()
for (filename in fileNames) {
        category <- gsub("FAMILY_", "", gsub(".txt", "", basename(filename)))
        buff <- read.table(filename, header=F, stringsAsFactor=F)
        df <- data.frame(
                Database = "FAMILY",
                Category = category,
                Gene = buff$V1
        )
        family <- rbind(family, df)
}
table(family$Category)
# CYP SLC UGT 
# 113 431  32

## Altogether
out <- rbind(
        hallmark,
        kegg,
        cancer,
        oncogenic,
        celltype,
        cellchat,
        gad,
        omim,
        drug,
        family
)

table(out$Database)
length(unique(out$Gene))
saveRDS(out, file.path(dataDir, "GeneList.RDS")) # Available at https://ShinyApps.UCalgary.ca/Coverage

q("no")
