##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
# Updated on January 06, 2025
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
setwd(file.path(getwd(), "R"))
source("00_settings.R") # load baseDir <- "/WORKING_DIR"

dataDir <- file.path(baseDir, "GeoMx", "Human")
outDir <- file.path(dataDir, "results")

dir.create(outDir, showWarnings = FALSE)

folderNames <- c(
        "SOA_Brain", "SOA_Colon", "SOA_Kidney", "SOA_Liver", "SOA_LymphNode", "SOA_Pancreas", "NS_Kidney", 
        "GSE208747", "GSE244117", "GSE254145", "GSE263897", "GSE272995", "GSE274938", "GSE275677", "GSE278670", "GSE281193"
)

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(GeomxTools)

objL <- lapply(seq_along(folderNames), function(idx) {
        dName <- folderNames[idx]
        message(dName)

        dccFiles <- dir(file.path(dataDir, dName, "DCC"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
        pkcFile <- dir(file.path(dataDir, dName, "PKC"), pattern = ".pkc$", full.names = TRUE, recursive = TRUE)
        annotFile <- dir(file.path(dataDir, dName, "Annot"), pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)
        initSet <- readNanoStringGeoMxSet(
                dccFiles = dccFiles,
                pkcFiles = pkcFile,
                phenoDataFile = annotFile,
                phenoDataSheet = "Annotation",
                phenoDataDccColName = "FileName",
                protocolDataColNames = c("ROI", "AOI"),
                experimentDataColNames = c("Panel")
        )
        initSet$Tag <- dName
        initSet

        saveRDS(initSet, file.path(outDir, paste0(dName, ".GeoMx.RDS")))
        
        return(initSet)
})
names(objL) <- folderNames
saveRDS(objL, file.path(outDir, "ObjL.GeoMx.RDS")) # Available at https://ShinyApps.UCalgary.ca/Coverage

q("no")

## SOA_Kidney, 19 Mar 2021 - *078 (052858A2, hu_kidney_001)
# Not all probes are found within PKC probe metadata. The following probes are ignored from analysis and were most likely removed from metadata while resolving multiple module PKC version conflicts.
# RTS0021749RTS0021916RTS0021976RTS0021979RTS0022732RTS0027943RTS0031708RTS0036612RTS0036631RTS0050279RTS0050433RTS0050547RTS0050580RTS0051214RTS0051219RTS0051232RTS0051449RTS0052009RTS0029660

## GSE274938
# A custom panel for the COVID19
