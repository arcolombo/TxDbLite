load("repeatsByClass.rda")
rptsByClass <- unlist(rptsByClass) 
repClass <- sort(unique(paste(rptsByClass$name, 
                              rptsByClass$family, 
                              rptsByClass$class, 
                              sep=":")))
repClass <- as.data.frame(do.call(rbind, strsplit(repClass, ":")))
names(repClass) <- c("name", "family", "class")
repClass <- repClass[ repClass$family != repClass$class, ] 
rownames(repClass) <- make.unique(repClass$name)

names(repClass) <- c("name", "tx_biotype", "gene_biotype")
repClass[ repClass$gene_biotype == "SINE?", "gene_biotype" ] <- "SINE"
repClass[ repClass$gene_biotype == "RC?", "gene_biotype" ] <- "RC"
repClass[ !(repClass$gene_biotype %in% c("LINE","SINE","DNA","LTR")),
         "gene_biotype" ] <- "other_repeat"
repClass[ repClass$gene_biotype == "DNA", "gene_biotype" ] <- "DNA_element"
repClass[ repClass$gene_biotype == "LTR", "gene_biotype" ] <- "LTR_element"
repClass$class <- "repeat"
saveRDS(repClass, file="repClass.rds")

repeattype <- repClass
save(repeattype, 
     file="~/Dropbox/artemis.ttriche/data/repeattype.rda", 
     compress="xz")

bestMatch <- function(query, corpus) {
  if (length(query) > 1) return(sapply(query, bestMatch, corpus=corpus))
  hits <- grep(query, corpus, ignore.case=TRUE, value=TRUE)
  if (length(hits) == 0) {
    hits <- agrep(query, 
                  corpus, 
                  max.distance=100,
                  costs=c(ins=1, del=2, sub=3),
                  ignore.case=TRUE, 
                  value=TRUE)
  }
  mindist <- which.min(adist(query, hits, ignore.case=TRUE, 
                             costs=c(ins=1, del=2, sub=3)))
  hits[mindist] 
}

repBaseToRmsk <- function(repeats, repClass, ...) { 
  perfectMatches <- names(repeats)[ names(repeats) %in% rownames(repClass) ]
  names(perfectMatches) <- perfectMatches
  imperfectMatches <- sapply(setdiff(names(repeats), perfectMatches),
                             bestMatch, corpus=repClass$name)
  allMatches <- c(perfectMatches, imperfectMatches)[names(repeats)] 
  res <- repClass[allMatches, ]
  rownames(res) <- names(repeats)
  res$family <- gsub("\\?$", "", res$family)
  res$class <- gsub("\\?$", "", res$class)
  res$class[ res$family == "SVA" ] <- "SINE"
  res$class[ grep("tRNA", res$family) ] <- "Small_RNA"
  res$class[ res$class == "RC" ] <- "DNA"
  res$class[ res$class == "DNA" ] <- "Transposon"
  names(res) <- c("name", "tx_biotype", "gene_biotype")
  res <- res[, c("tx_biotype", "gene_biotype")]
  res$class <- "repeat"
  mcols(repeats) <- DataFrame(res)
  repeats$tx_id <- NA
  repeats$gene_id <- NA
  repeats$gene_name <- NA
  repeats$entrezid <- rep(NA, length(repeats))
  cols <- c("tx_id", "gene_id", "gene_name", "entrezid", 
            "tx_biotype", "gene_biotype", "class")
  return(repeats[, cols])
}
