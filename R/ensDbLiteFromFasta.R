#' Function to create an EnsDbLite object (a stripped-down SQLite db) from
#' specified Ensembl FASTA files for a given version, species, and type(s).
#' Note that this is shockingly easier than the same process for the GTF...
#' 
#' @param organism      the organism involved (Homo_sapiens, Mus_musculus, etc.)
#' @param fastaFiles    the FASTA files to collate into a EnsDbLite instance
#' @param version       the Ensembl build (e.g. 81)
#' @param verbose       make a lot of noise? (TRUE) 
#' 
#' @import GenomicRanges
#' @import OrganismDbi
#' @import Biostrings
#' 
#' @export
#'
ensDbLiteFromFasta <- function(organism, fastaFiles, version, verbose=TRUE){#{{{

  require(Biostrings) 
  require(GenomicRanges)
  options(useFancyQuotes=FALSE)
  options(stringsAsFactors=FALSE)
  organism <- sub("\\.", "_", ## try & be robust
                  sub("Mmusculus", "Mus_musculus", 
                      sub("Hsapiens", "Homo_sapiens", organism)))

  ## utility functions for parsing ENSEMBL FASTAs
  splt <- function(x, y=" ") strsplit(x, y)[[1]]
  popx <- function(x) x[length(x)] 
  pop <- function(x, y=" ") popx(splt(x, y))
  grab <- function(x, y=" ", i=1) splt(x, y)[i]
  shift <- function(x, y=" ") grab(x, y, i=1)

  if (organism == "Homo_sapiens") { ## {{{ prefixes differ by organism 
    txpre <- "ENST"
    gxpre <- "ENSG"
  } else if (organism == "Mus_musculus") {
    txpre <- "ENSMUST"
    gxpre <- "ENSMUSG"
  } else {
    stop("Currently supporting Homo_sapiens & Mus_musculus... patches welcome!")
  } # }}}

  if (verbose) cat("Extracting transcript lengths...") # {{{
  txLen <- do.call(c, lapply(fastaFiles, fasta.seqlengths))
  names(txLen) <- sapply(sapply(names(txLen), shift), pop, y="\\.")
  if (verbose) cat("done.\n") # }}}
  
  if (verbose) cat("Extracting transcript descriptions...") # {{{
  txDesc <- do.call(c, lapply(lapply(fastaFiles, fasta.seqlengths), names))
  names(txDesc) <- sapply(txDesc, shift)
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Extracting genomic coordinates...") # {{{
  txCoords <- as.data.frame(strsplit(sapply(txDesc, grab, i=3), ":"))[-1,]
  rownames(txCoords) <- c("genomeVersion","seqnames","start","end","strand")
  genomeVersion <- unique(as.character(txCoords["genomeVersion",]))
  if (length(genomeVersion) > 1) {
    message("These fasta files are from different genome assemblies!")
    message(paste(genomeVersion, collapse=", "))
    stop("Aborting construction of annotation database.")
  }
  txs <- GRanges(seqnames=as.character(txCoords["seqnames",]),
                 ranges=IRanges(start=as.integer(txCoords["start",]),
                                end=as.integer(txCoords["end",])),
                 strand=ifelse(txCoords["strand",] == "1", "+", "-"))
  names(txs) <- names(txCoords)
  genome(txs) <- genomeVersion
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Extracting gene and biotype associations...") # {{{
  mcols(txs) <- DataFrame(tx_id=names(txs),
                          tx_length=txLen[names(txs)],
                          gene_id=sapply(sapply(txDesc,grab,i=4), pop, ":"),
                          tx_biotype=sapply(sapply(txDesc,grab,i=6), pop,":"),
                          gene_biotype=sapply(sapply(txDesc,grab,i=5), pop,":"))
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Tabulating transcript biotypes...") # {{{
  tx_biotypes <- data.frame(id=seq_along(levels(as.factor(txs$tx_biotype))),
                            tx_biotype=levels(as.factor(txs$tx_biotype)))
  txs$tx_biotype_id <- as.numeric(as.factor(txs$tx_biotype))
  ## retrieved as factors, with levels set from db
  if (verbose) cat("done.\n") # }}}

  message("FIXME: add pre-stored biotype classes!")

  squash <- function(x, by, FUN) { # {{{
    xx <- aggregate(x=x, by=by, FUN=FUN)
    x <- xx[,2]
    names(x) <- xx[,1]
    x
  } # }}}

  if (verbose) cat("Tabulating genes...") # {{{
  gxChr <- squash(as.character(seqnames(txs)), by=list(txs$gene_id), FUN=unique)
  gxStart <- squash(start(txs), by=list(txs$gene_id), FUN=min)
  stopifnot(identical(names(gxChr), names(gxStart)))
  gxEnd <- squash(end(txs), by=list(txs$gene_id), FUN=max)
  stopifnot(identical(names(gxChr), names(gxEnd)))
  gxStrand <- squash(as.character(strand(txs)), by=list(txs$gene_id),FUN=unique)
  stopifnot(identical(names(gxChr), names(gxStrand)))
  gxs <- GRanges(seqnames=gxChr,
                 ranges=IRanges(start=gxStart, end=gxEnd),
                 strand=gxStrand)
  genome(gxs) <- genomeVersion
  names(gxs) <- names(gxChr)
  gxs$gene_id <- names(gxs)
  gxs$gene <- as.integer(sub(gxpre, "", names(gxs)))
  gxBiotype <- squash(txs$gene_biotype, by=list(txs$gene_id), FUN=unique)
  gxs$gene_biotype <- gxBiotype[names(gxs)] 
  gxs$gene_biotype_id <- as.numeric(as.factor(gxs$gene_biotype))
  gxMedLen <- squash(txs$tx_length, by=list(txs$gene_id), FUN=median)
  gxs$median_length <- gxMedLen[names(gxs)]
  gxs$entrezid <- getEntrezIDs(gxs)[names(gxs)]
  gxs$gene_name <- getSymbols(gxs)[names(gxs)]
  txs$gene <- as.integer(sub(gxpre, "", txs$gene_id))
  if (verbose) cat("...done.\n") # }}}

  if (verbose) cat("Creating the database...") # {{{
  txVersion <- paste0("v", version)
  abbr <- c(Homo_sapiens="Hsapiens", Mus_musculus="Mmusculus")
  outstub <- paste("EnsDbLite", abbr[organism], txVersion, sep=".")
  dbname <- paste(outstub, "sqlite", sep=".") 
  con <- dbConnect(dbDriver("SQLite"), dbname=dbname)
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Writing the gene table...") # {{{
  gxcols <- c("seqnames", "start", "end", "strand",
              "gene", "gene_id", "gene_biotype_id",
              "entrezid", "gene_name", "median_length")
  gx <- as(gxs, "data.frame")[, gxcols] 
  dbWriteTable(con, name="gene", gx, overwrite=T, row.names=F)
  rm(gx)
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Tabulating gene biotypes...") # {{{
  gene_biotypes <- data.frame(id=seq_along(levels(as.factor(txs$gene_biotype))),
                              gene_biotype=levels(as.factor(txs$gene_biotype)))
  gxs$gene_biotype_id <- as.numeric(as.factor(gxs$gene_biotype))
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Writing the gene_biotype table...") # {{{
  dbWriteTable(con, name="gene_biotype", gene_biotypes, 
               overwrite=T, row.names=F)
  rm(gene_biotypes) 
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Writing the tx table...") # {{{
  txcols <- c("start", "end", "tx_id", "tx_length", "tx_biotype_id", "gene")
  tx <- as(txs, "data.frame")[, txcols]
  dbWriteTable(con, name="tx", tx, overwrite=TRUE, row.names=FALSE)
  rm(tx)
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Tabulating transcript biotypes...") # {{{
  gene_biotypes <- data.frame(id=seq_along(levels(as.factor(txs$gene_biotype))),
                              gene_biotype=levels(as.factor(txs$gene_biotype)))
  gxs$gene_biotype_id <- as.numeric(as.factor(gxs$gene_biotype))
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Writing the tx_biotype table...") # {{{
  dbWriteTable(con, name="tx_biotype", tx_biotypes, overwrite=T, row.names=F)
  rm(tx_biotypes)
  if (verbose) cat("done.\n") # }}}

  ## write metadata table # {{{ 
  Metadata <- ensDbLiteMetadata(packageName=outstub, 
                                genomeVersion=txVersion,
                                sourceFile=paste(sub(".gz$", "", fastaFiles),
                                                     collapse="+"))
  dbWriteTable(con, name="metadata", Metadata, overwrite=TRUE, row.names=FALSE)
  # }}}

  ## create indices  # {{{
  dbGetQuery(con, "create index tx_id_idx on tx (tx_id);")
  dbGetQuery(con, "create index gene_idx on gene (gene);")
  dbGetQuery(con, "create index txb_id_idx on tx_biotype (id);")
  dbGetQuery(con, "create index gxb_id_idx on gene_biotype (id);")
  # }}}

  ## finish 
  dbDisconnect(con)
  return(dbname)

} # }}}


#' @describeIn ensDbLiteFromFasta
#' 
#' add EntrezGene IDs for Ensembl genes
#'
#' @param gxs       a GRanges of genes
#' 
#' @return  entrez_id values for the genes, where found
#'
getEntrezIDs <- function(gxs) { # {{{
  if (substr(unique(genome(gxs)), 1, 4) == "GRCh") {
    orgDb <- "Homo.sapiens"
  } else if (substr(unique(genome(gxs)), 1, 4) == "GRCm") {
    orgDb <- "Mus.musculus"
  } else {
    stop("Could not determine organism from genome of gxs.  Patches welcome!")
  } 
  library(orgDb, character.only=TRUE) 
  mapIds(get(orgDb), keys=gxs$gene_id, column="ENTREZID", keytype="ENSEMBL") 
} # }}}


#' @describeIn ensDbLiteFromFasta
#' 
#' add symbols for Ensembl genes
#'
#' @param gxs       a GRanges of genes
#' 
#' @return  symbols for the genes, where found 
#'
getSymbols <- function(gxs) { # {{{
  if (substr(unique(genome(gxs)), 1, 4) == "GRCh") {
    orgDb <- "Homo.sapiens"
  } else if (substr(unique(genome(gxs)), 1, 4) == "GRCm") {
    orgDb <- "Mus.musculus"
  } else {
    stop("Could not determine organism from genome of gxs.  Patches welcome!")
  } 
  library(orgDb, character.only=TRUE) 
  mapIds(get(orgDb), keys=gxs$gene_id, column="SYMBOL", keytype="ENSEMBL") 
} # }}}


#' @describeIn ensDbLiteFromFasta
#' 
#' create metadata for an EnsDbLite instance
#'
#' @param packageName   the name of the annotation package to be built 
#' @param genomeVersion name of genome assembly for coordinates, e.g. "GRCh38"
#' @param sourceFile    name of FASTA file(s) whence it was built, as a string
#' 
#' @return a data.frame of metadata suitable for cramming into the database
#'
ensDbLiteMetadata <- function(packageName, genomeVersion, sourceFile) { # {{{

  tokens <- strsplit(packageName, "\\.")[[1]]
  organism <- c(Hsapiens="homo_sapiens", Mmusculus="mus_musculus")[tokens[2]]  
  ensemblVersion <- sub("v", "", tokens[3])

  MetaData <- data.frame(matrix(ncol=2, nrow=11))
  colnames(MetaData) <- c("name", "value")
  MetaData[1,] <- c("package_name", packageName)
  MetaData[2,] <- c("db_type", "EnsDbLite")
  MetaData[3,] <- c("type_of_gene_id", "Ensembl Gene ID")
  MetaData[4,] <- c("created_by", paste("TxDbLite", packageVersion("TxDbLite")))
  MetaData[5,] <- c("creation_time", date())
  MetaData[6,] <- c("organism", organism )
  MetaData[7,] <- c("genome_build", genomeVersion)
  MetaData[8,] <- c("ensembl_version", ensemblVersion)
  MetaData[9,] <- c("source_file", sourceFile)
  return(MetaData)

} # }}}
