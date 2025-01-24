set.seed(1118)
baseDir <- "~/Documents/ASOC/Analysis/Coverage"

# User defined functions
overlapGroups <- function(listInput, sort = TRUE) {
        listInputmat <- fromList(listInput) == 1
        listInputunique <- unique(listInputmat)
        grouplist <- list()
        for (i in 1:nrow(listInputunique)) {
                currentRow <- listInputunique[i, ]
                myelements <- which(apply(listInputmat, 1, function(x) all(x == currentRow)))
                attr(myelements, "groups") <- currentRow
                grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
                myelements
        }
        if (sort) {
                grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
        }
        attr(grouplist, "elements") <- unique(unlist(listInput))
        return(grouplist)
}

