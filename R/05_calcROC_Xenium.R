##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
# Updated on January 13, 2025
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_settings.R") # load baseDir <- "/WORKING_DIR"

dataDir <- file.path(baseDir, "Xenium", "data")
outDir <- file.path(baseDir, "Xenium", "results")
resourceDir <- file.path(baseDir, "GeneLists")

plotting <- TRUE

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(stringr)
library(tidyverse)
library(pROC)
library(dplyr)

panel <- read.csv(file.path(dataDir, "XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv"))
geneSymbols <- panel$gene_name
names(geneSymbols) <- panel$gene_id

stats <- readRDS(file.path(dataDir, "Stats.Xenium.RDS"))
stats <- stats[which(stats$Gene %in% panel$gene_id), c(1:8)] # human
stats$Gene <- as.character(geneSymbols[stats$Gene])

geneList <- readRDS(file.path(resourceDir, "GeneList.RDS"))

AUCs <- c()
nohit <- c()
for (database in unique(geneList$Database)) {
        message(database)
        ggList <- geneList[which(geneList$Database == database),]

        for (categ in unique(ggList$Category)) {
                message(categ)
                gList <- ggList[which(ggList$Category == categ),]

                for (idx in c(2:ncol(stats))) {
                        dName <- colnames(stats)[idx]
                        message(dName)

                        df <- as.data.frame(stats[,c(1,idx)])
                        colnames(df) <- c("Gene", "Count")
                        if (length(which(is.na(df$Count))) > 0) {
                                df <- df[-which(is.na(df$Count)),]
                        }
                        df <- df[order(df$Count),]

                        df$Hit <- "No"
                        df$Hit[which(df$Gene %in% gList$Gene)] <- "Yes"
                        df$Hit <- factor(df$Hit, levels=c("Yes", "No"))
                        table(df$Hit)

                        if (all(table(df$Hit) > 2)) {
                                p_roc <- roc(Hit ~ Count, data = df, plot=F, direction = ">", levels=c("Yes", "No"))
                                auc <- unlist(strsplit(toString(p_roc$auc),":"))
                                AUCs <- rbind(AUCs, c(database, categ, dName, auc))

                                if (plotting) {
                                        pdf(file.path(outDir, "PDFs", paste0(database, "__", categ, "__", dName, ".pdf")), width=14, height=7)
                                        par(mfrow=c(1,2))
                                        b <- barplot(log10(df$Count + 1), horiz = T, col = "grey80", border = "grey80", xlab="Count, log10", ylab="Gene", main=dName)
                                        text(log10(df$Count[which(df$Hit == "Yes")] + 1), b[which(df$Hit == "Yes")], labels = df$Gene[which(df$Hit == "Yes")], col="red", pos=2)
                                        legend("bottomright", legend=paste0(length(which(df$Hit == "Yes")), " genes"), title = paste0(database, ", ", categ), bty="n")

                                        roc(df$Hit ~ df$Count, plot=T, col="red", smooth=F, print.auc=F, lwd=2, direction = ">", levels=c("Yes", "No"))
                                        auc_f <- format(round(as.numeric(auc),4), nsmall = 4)
                                        legend("bottomright", c(paste("AUC=", auc_f, sep="")), lty=1, bty="n", col="red", lwd=2)
                                        dev.off()
                                }
                        } else {
                                nohit <- rbind(nohit, c(database, categ, dName))
                        }
                }
        }
}

AUCs
summary(as.numeric(AUCs[,4]))
write.table(AUCs, file.path(outDir, "AUCs.txt"), row.names=F, col.names=T, quote=F, sep="\t")

outDf <- data.frame(
        DB = paste0(AUCs[,1], "$", AUCs[,2]),
        Dataset = AUCs[,3],
        AUC = as.numeric(AUCs[,4])
)
outDf <- outDf %>% group_by(DB) %>%  summarise(AvgAUC = mean(AUC, na.rm=T))
write.table(outDf, file.path(outDir, "AUCs_summary.txt"), row.names=F, col.names=T, quote=F, sep="\t")

write.table(nohit, file.path(outDir, "No_hit.txt"), row.names=F, col.names=F, quote=F, sep="\t")

q("no")
