#' create an ERCC annotation DB from the FASTA file
#' NOTE: we probably shouldn't even bother exporting this.
#' If ERCC spike-ins change radically, we'll need to build an entirely new one.
#'
#' @import Rsamtools
#'
#' @export
erccDbLiteFromFasta <- function(fastaFile, verbose=TRUE) { 

  if (verbose) cat("Extracting spike-in associations...")
  faFile <- FaFile(fastaFile)
  if (!file.exists(index(faFile))) indexFa(path(faFile))
  txs <- scanFaIndex(faFile)
  names(txs) <- seqnames(txs)
  mcols(txs) <- DataFrame(tx_length=width(txs),
                          gc_content=GCcontent(scanFa(fastaFile)),
                          tx_id=names(txs),
                          gene_id= names(txs),
                          gene_name=rep(NA, length(txs)),
                          entrezid=rep(NA, length(txs)))
  txs <- .addSpikeInSubgroup(txs)
  if (verbose) cat("done.\n")

  if (verbose) cat("Creating the database...") # {{{
 
  #for purpsoses of building the package, the name can not have an underscore

   


 outstub <- getTxDbLiteName(fastaFile)
  dbname <- paste(outstub, "sqlite", sep=".") 
  con <- dbConnect(dbDriver("SQLite"), dbname=dbname)
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Writing the spike-in tables...") # {{{
  txcols <- c("seqnames", "start", "end", "strand",
              "tx_length", "gc_content", 
              "tx_id", "gene_id", "gene_name", "entrezid", 
              "tx_biotype", "gene_biotype", "biotype_class")
  tx <- as(txs, "data.frame")[, txcols] 
  dbWriteTable(con, name="tx", tx, overwrite=T, row.names=F)
  rm(tx)
  if (verbose) cat("done.\n") # }}}

  Metadata <- .erccDbLiteMetadata(packageName=outstub, sourceFile=fastaFile)
  dbWriteTable(con, name="metadata", Metadata, overwrite=TRUE, row.names=FALSE)
  dbGetQuery(con, "create index tx_id_idx on tx (tx_id);")

  ## finish 
  dbDisconnect(con)
  return(dbname)

}

.addSpikeInSubgroup <- function(txs) { # {{{

  data(ERCC_annotated, package="TxDbLite")
  mcols(txs)$tx_biotype <- rep("SpikeIn_unannotated", length(txs))
  mcols(txs)$gene_biotype <- rep("SpikeIn", length(txs))

  named <- intersect(as.character(seqnames(txs)), rownames(ERCC_annotated))
  mcols(txs[named])$tx_biotype <- paste0("SpikeIn_", 
                                         ERCC_annotated[named, "subgroup"])
  mcols(txs)$biotype_class <- rep("SpikeIn", length(txs))
  
  return(txs)

} # }}}

.erccDbLiteMetadata <- function(packageName, sourceFile) { # {{{

  MetaData <- data.frame(matrix(ncol=2, nrow=8))
  colnames(MetaData) <- c("name", "value")
  MetaData[1,] <- c("package_name", packageName)
  MetaData[2,] <- c("db_type", "ErccDbLite")
  MetaData[3,] <- c("type_of_gene_id", "N/A")
  MetaData[4,] <- c("created_by", paste("TxDbLite", packageVersion("TxDbLite")))
  MetaData[5,] <- c("creation_time", date())
  MetaData[6,] <- c("organism", "N/A")
  MetaData[7,] <- c("genome_build", "N/A")
  MetaData[8,] <- c("source_file", sourceFile)
  return(MetaData)

} # }}}
