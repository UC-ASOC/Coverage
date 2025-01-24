##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
# Updated on October 29, 2024
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_settings.R") # load baseDir <- "/WORKING_DIR"

dataDir <- file.path(baseDir, "Xenium", "data")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(HDF5Array)
library(tidyverse)

h5Filenames <- list.files(file.path(dataDir), pattern = "h5", full.names = T)
h5Names <- gsub(".h5", "", basename(h5Filenames))

sumL <- lapply(h5Filenames, function(h5Filename) {
        expr <- TENxMatrix(h5Filename, "matrix")
        colSum <- apply(as.matrix(expr), 1, sum)
        df <- data.frame(
                Gene = names(colSum),
                Count = colSum
        )
        return(df)
})
names(sumL) <- h5Names

stats <- sumL %>% reduce(full_join, by = "Gene")
colnames(stats) <- c("Gene", h5Names)

saveRDS(stats, file.path(dataDir, "Stats.Xenium.RDS")) # Available at https://ShinyApps.UCalgary.ca/Coverage

q("no")
