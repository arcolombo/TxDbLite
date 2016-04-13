# needs cleanup
if (FALSE) {
  
  library(biomaRt) # only to get current ENSG/ENST-to-REACTOME
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") 

  columns <- list(gene=c(ENSG="ensembl_gene_id",
                         ENST="ensembl_transcript_id",
                         HUGO="hgnc_symbol",
                         CCDS="ccds",
                         REACTOME="reactome"),
                  transcript=c(ENSG="ensembl_gene_id",
                               ENST="ensembl_transcript_id",
                               HUGO="hgnc_transcript_name",
                               CCDS="ccds",
                               REACTOME="reactome"))
  filters <- list(gene="ensembl_gene_id",
                  transcript="ensembl_transcript_id")

  # by gene, for example, but we can do this for all txs in an ENSEMBL release
  genes <- c(CBL="ENSG00000110395", 
             WT1="ENSG00000184937", 
             DNMT3A="ENSG00000119772")

  geneLevel <- getBM(attributes=columns$gene,
                     filters=filters$gene, 
                     values=genes, 
                     mart=ensembl)

  # at transcript level, which is probably preferable for TxDbLite mappings:
  txs <- unique(geneLevel$ensembl_transcript_id)
  txLevel <- getBM(attributes=columns$transcript,
                   filters=filters$transcript,
                   values=txs,
                   mart=ensembl)

  # last step: pull the Reactome pathway names and FI IDs out of the text files:
  reactomeNames <- read.table("ReactomePathways.txt", sep="\t") 
  names(reactomeNames) <- c("reactome", "name", "species") 
  reactomeNames <- reactomeNames[!duplicated(reactomeNames$reactome), ] # weird
  rownames(reactomeNames) <- reactomeNames$reactome

  # fill in the pathway names for those that have them
  hasReactome <- which(txLevel$reactome != "")
  txLevel$reactomePathway <- ""
  txLevel$reactomePathway[hasReactome] <- 
    reactomeNames[txLevel$reactome[hasReactome], "name"]

  # add in Reactome FI support (gene level)
  reactomeFI <- read.table("FIsInGene_031516_with_annotations.txt", 
                           head=TRUE, sep="\t") # need to map by ENSG, so...

}
