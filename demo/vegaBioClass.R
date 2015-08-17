library(ensembldb)


vegaBioClass <- function(txmap) {
  gene_biotypes <- names(sort(table(txmap$gene_biotype)))
  tx_biotypes <- names(sort(table(txmap$tx_biotype)))
  biotypes <- DataFrame(biotype=sort(union(gene_biotypes, tx_biotypes)))
  biotypes$class <- rep("Unknown", nrow(biotypes))
  rownames(biotypes) <- biotypes$biotype

  biotypes$class[grep("lnc", biotypes$biotype)] <- "long_noncoding"
  biotypes$class[grep("ncrna", biotypes$biotype)] <- "long_noncoding"
  biotypes$class[grep("lincRNA", biotypes$biotype)] <- "long_noncoding"
  biotypes$class[grep("(IG|TR|CDS|decay|coding)", 
                      biotypes$biotype)] <- "protein_coding"
  biotypes$class[grep("seudogene", biotypes$biotype)] <- "pseudogene"
  biotypes[grep("Mt_", biotypes$biotype), "class"] <- "mitochondrial"
  biotypes[grep("sense_", biotypes$biotype), "class"] <- "long_noncoding"
  biotypes["LRG_gene", "class"] <- "exclude_from_analysis"
  biotypes["TEC", "class"] <- "exclude_from_analysis"
  biotypes["vaultRNA", "class"] <- "short_noncoding"
  biotypes["ribozyme", "class"] <- "long_noncoding"
  biotypes["snoRNA", "class"] <- "short_noncoding"
  biotypes["scaRNA", "class"] <- "short_noncoding"
  biotypes["miRNA", "class"] <- "short_noncoding"
  biotypes["snRNA", "class"] <- "short_noncoding"
  biotypes["sRNA", "class"] <- "short_noncoding"
  biotypes["rRNA", "class"] <- "short_noncoding"
  biotypes["misc_RNA", "class"] <- "short_noncoding"
  biotypes["antisense", "class"] <- "long_noncoding"
  biotypes["polymorphic_pseudogene", "class"] <- "protein_coding"
  biotypes["processed_transcript", "class"] <- "long_noncoding"
  biotypes["retained_intron", "class"] <- "long_noncoding"
  
  return(biotypes[, "class", drop=FALSE])
}

if (FALSE) {
  transcriptome <- "EnsDb.Hsapiens.v80"
  library(transcriptome, character.only=TRUE)
  txcolumns <- c("gene_id","gene_name","entrezid","tx_biotype","gene_biotype")
  txmap <- transcripts(get(transcriptome), columns=txcolumns)
  seqlevelsStyle(txmap) <- "UCSC"
  biotypes <- vegaBioClass(txmap)
}
