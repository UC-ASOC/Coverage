##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
# Updated on January 09, 2025
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_settings.R") # load baseDir <- "/WORKING_DIR"

dataDir <- file.path(baseDir, "GeoMx", "Human")
resourceDir <- file.path(baseDir, "GeneLists")
outDir <- file.path(dataDir, "results")

plotting <- TRUE

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(GeomxTools)
library(stringr)
library(tidyverse)
library(pROC)
library(dplyr)
library(corrplot)

if (file.exists(file.path(outDir, "ObjL.GeoMx.RDS"))) {
        objL <- readRDS(file.path(outDir, "ObjL.GeoMx.RDS"))
} else {
        stop("Run loadGeoMx.R script")
}

sumL <- lapply(seq_along(objL), function(idx) {
        obj <- objL[[idx]]

        gSet <- GeomxTools::aggregateCounts(obj)
        colSum <- apply(as.matrix(assayData(gSet)$exprs), 1, sum)
        df <- data.frame(
                Gene = names(colSum),
                Count = colSum
        )
        return(df)
})
names(sumL) <- names(objL)

stats <- sumL %>% reduce(full_join, by = "Gene")
colnames(stats) <- c("Gene", names(objL))

counts <- stats[, c(2:ncol(stats))]
rownames(counts) <- stats[, 1]
counts <- counts[rowSums(is.na(counts))==0,]

coef <- cor(counts, method = "spearman")
pdf(file.path(outDir, "Correlation.pdf"), width=15, height=15)
corrplot(coef, type='lower', order = 'hclust', addCoef.col = 'black', tl.pos = 'd', tl.col = 'black', cl.pos = 'n')
dev.off()

geneList <- readRDS(file.path(resourceDir, "GeneList.RDS"))

AUCs <- c()
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
                                        legend("bottomright", c(paste("AUC=", auc_f, sep="")), lty=1, bty="n", col="blue", lwd=2)
                                        dev.off()
                                }
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

q("no")
