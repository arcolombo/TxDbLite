#' create an ERCC annotation DB from the FASTA file
#' NOTE: we probably shouldn't even bother exporting this.
#' If ERCC spike-ins change radically, we'll need to build an entirely new one.
#'
#' @importFrom Rsamtools indexFa index
#' @importFrom Rsamtools scanFaIndex FaFile scanFa path
#' @importFrom GenomeInfoDb seqnames
#' @importFrom DBI dbConnect dbDriver dbWriteTable dbGetQuery dbDisconnect
#' @importFrom S4Vectors DataFrame
#' @param fastaFile a fasta file containing ERCC sequence information
#' @param verbose boolean if true will print process messages 
#' @param dryRun  boolean if false a sql-lite data base is saved, if true, no database is saved.
#' @export
erccDbLiteFromFasta <- function(fastaFile, verbose=TRUE, dryRun=FALSE) { 

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
  txs<-as.data.frame(txs)
  txs$gene<-txs$seqnames
  txs$tx_biotype_id<-as.numeric(as.factor(txs$tx_biotype))
  txs$gene_biotype_id<-as.numeric(as.factor(txs$gene_biotype))
  txs$median_length<-txs$tx_length
  if (verbose) cat("done.\n")

  if (verbose) cat("Creating the database...") # {{{
 
  #for purpsoses of building the package, the name can not have an underscore

   


 outstub <- getTxDbLiteName(basename(fastaFile))
  dbname <- paste(outstub, "sqlite", sep=".") 
  if(dryRun==FALSE){
  con <- dbConnect(dbDriver("SQLite"), dbname=dbname)
  }
  if (verbose) cat("done.\n") # }}}

 if(verbose) cat("Writing the gene table...")
   gxcols <- c("seqnames", "start", "end", "strand",
              "gene", "gene_id", "gene_biotype_id",
              "entrezid", "gene_name", "median_length","copyNumber")
   if(dryRun==FALSE){
   gxs<-as.data.frame(txs[,gxcols])
   dbWriteTable(con,name="gene",gxs,overwrite=T,row.names=F)
  #}}}
  }
  if(verbose) cat("Tabulating gene biotypes...")
  gene_biotypes<-data.frame(id=seq_along(levels(as.factor(txs$gene_biotype))),
                            gene_biotype=levels(as.factor(txs$gene_biotype)))

  if(dryRun==FALSE){
  dbWriteTable(con,name="gene_biotype",gene_biotypes,overwrite=T,row.names=F)
   }
  if(verbose) cat ("done. \n")

  if (verbose) cat("Writing the spike-in tables...") # {{{
  txcols <- c("start", "end", "tx_id",
              "tx_length", "gc_content", 
               "tx_biotype_id","gene","copyNumber")
  tx <- txs[, txcols] 
  if(dryRun==FALSE){
  dbWriteTable(con, name="tx", tx, overwrite=T, row.names=F)
  }
  rm(tx)
  if (verbose) cat("done.\n") # }}}
  
  if(verbose) cat("Tabulating transcript biotypes...")
  if(dryRun==FALSE){
  tx_biotypes<-data.frame(id=seq_along(levels(as.factor(txs$tx_biotype))),
                         tx_biotype=levels(as.factor(txs$tx_biotype)))
  dbWriteTable(con,name="tx_biotype",tx_biotypes,overwrite=T,row.names=F)
  }

  if(verbose) cat("Writing the biotype_class table...")
  if(dryRun==FALSE){
   biotype_class<-data.frame(biotype=levels(as.factor(txs$gene_biotype)),
                             class="SpikeIn")
  dbWriteTable(con,name="biotype_class",biotype_class,overwrite=T,row.names=F)
   }

  Metadata <- .erccDbLiteMetadata(packageName=outstub, sourceFile=fastaFile)
  if(dryRun==FALSE){
  dbWriteTable(con, name="metadata", Metadata, overwrite=TRUE, row.names=FALSE)
  dbGetQuery(con, "create index tx_id_idx on tx (tx_id);")
  
  ## finish 
  dbDisconnect(con)
  }
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
  mcols(txs)$copyNumber<-rep(1,length(txs))
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
