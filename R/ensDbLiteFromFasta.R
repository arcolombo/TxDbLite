#' Function to create an EnsDbLite object (a stripped-down SQLite db) from
#' a specified Ensembl FASTA file for a given version, species, and type(s).
#' Note that this is shockingly easier than the same process for the GTF...
#' 
#' @param fastaFile     the FASTA file to collate into a EnsDbLite instance
#' @param verbose       make a lot of noise? (TRUE) 
#' 
#' @import GenomicRanges
#' @import OrganismDbi
#' @importFrom Rsamtools indexFa scanFaIndex scanFa
#' 
#' @export
ensDbLiteFromFasta <- function(fastaFile, verbose=TRUE){#{{{

  options(useFancyQuotes=FALSE)
  options(stringsAsFactors=FALSE)

  ## utility functions for parsing ENSEMBL FASTAs
  splt <- function(x, y=" ") strsplit(x, y)[[1]]
  popx <- function(x) x[length(x)] 
  pop <- function(x, y=" ") popx(splt(x, y))
  grab <- function(x, y=" ", i=1) splt(x, y)[i]
  shift <- function(x, y=" ") grab(x, y, i=1)

  txDbLiteName <- getTxDbLiteName(fastaFile)
  genomeVersion <- strsplit(fastaFile, "\\.")[[1]][1]
  tokens <- strsplit(txDbLiteName, "\\.")[[1]]
  organism <- tokens[2] 
  version <- tokens[3]
  org <- getOrgDetails(organism) 
  if (!require(org$package, character.only=TRUE)) {
    stop("Please install the", org$package, "package, then try again. Thanks!")
  } 

  if (verbose) cat("Extracting transcript lengths...") # {{{
  txLen <- fasta.seqlengths(fastaFile)
  names(txLen) <- sapply(names(txLen), idStub)
  if (verbose) cat("done.\n") # }}}
  
  if (verbose) cat("Extracting transcript descriptions...") # {{{
  txDesc <- names(fasta.seqlengths(fastaFile))
  names(txDesc) <- sapply(txDesc, idStub)
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

  matrixStrand <- ifelse(txCoords["strand",] == "1", "+", "-") # 1xN matrix 
  #converting to named character class for proper strand conversion. 
  characterStrand <- as.character(matrixStrand)
  names(characterStrand) <- sapply(colnames(matrixStrand), idStub)
  stopifnot(identical(length(matrixStrand),length(characterStrand)))

  # a sanity check enforcing strand conservation:
  stopifnot(identical(Rle(strand(characterStrand)),
                      strand(Rle(characterStrand)))) 

  txs <- GRanges(seqnames=as.character(txCoords["seqnames",]),
                 ranges=IRanges(start=as.integer(txCoords["start",]),
                                end=as.integer(txCoords["end",])),
                 strand=characterStrand)
  names(txs) <- names(txCoords)
  genome(txs) <- genomeVersion
  if (verbose) cat("done.\n") # }}}

  # from here on out, all of the ENSG/ENST fuckery is fixed by idStub()

  if (verbose) cat("Extracting gene and biotype associations...") # {{{
  mcols(txs) <- DataFrame(tx_id=names(txs),
                          tx_length=txLen[names(txs)],
                          gene_id=sapply(sapply(txDesc,grab,i=4), pop, ":"),
                          tx_biotype=sapply(sapply(txDesc,grab,i=6), pop,":"),
                          gene_biotype=sapply(sapply(txDesc,grab,i=5), pop,":"))
  mcols(txs)$gene_id <- sapply(mcols(txs)$gene_id, idStub)
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Tabulating GC content...") # {{{
  txs$gc_content <- GCcontent(scanFa(fastaFile))
  if (verbose) cat("done.\n") # }}}

  if (verbose) cat("Tabulating transcript biotypes...") # {{{
  tx_biotypes <- data.frame(id=seq_along(levels(as.factor(txs$tx_biotype))),
                            tx_biotype=levels(as.factor(txs$tx_biotype)))
  txs$tx_biotype_id <- as.numeric(as.factor(txs$tx_biotype))
  ## retrieved as factors, with levels set from db
  if (verbose) cat("done.\n") # }}}

  squash <- function(x, by, FUN) { # {{{
    xx <- aggregate(x=x, by=by, FUN=FUN)
    if(all(sapply(xx[,2], length) == 1) == "TRUE") { 
      x <- xx[,2]
      names(x) <- xx[,1]
    }

    if(!all(sapply(xx[,2],length) == 1) == "TRUE"){
      idx <- which(sapply(xx[,2],length) > 1)
      xx[idx,2] <- sapply(xx[idx, 2], 
                          function(x) x <- "*") # does this do anything?
      stopifnot(all(sapply(xx[,2], length) == 1))
      x <- xx[,2]
      names(x) <- xx[,1] 
      }
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
  gxCharacterStrand <- as.character(gxStrand)
  names(gxCharacterStrand) <- names(gxStrand)
  stopifnot(identical(Rle(strand(gxCharacterStrand)),
                      strand(Rle(gxCharacterStrand))))  
  stopifnot(identical(names(gxCharacterStrand),names(gxChr)))

  gxs <- GRanges(seqnames=gxChr,
                 ranges=IRanges(start=gxStart, end=gxEnd),
                 strand=gxCharacterStrand)
  genome(gxs) <- genomeVersion
  names(gxs) <- names(gxChr)
  gxs$gene_id <- names(gxs)
  gxs$gene <- as.integer(sub(org$gxpre, "", names(gxs)))
  gxBiotype <- squash(txs$gene_biotype, by=list(txs$gene_id), FUN=unique)
  gxs$gene_biotype <- gxBiotype[names(gxs)] 
  gxs$gene_biotype_id <- as.numeric(as.factor(gxs$gene_biotype))
  gxMedLen <- squash(txs$tx_length, by=list(txs$gene_id), FUN=median)
  gxs$median_length <- gxMedLen[names(gxs)]
  gxs$entrezid <- getEntrezIDs(gxs, organism)[names(gxs)]
  gxs$gene_name <- getSymbols(gxs, organism)[names(gxs)]
  txs$gene <- as.integer(sub(org$gxpre, "", txs$gene_id))
  gxs$gene <- as.integer(sub(org$gxpre,"",gxs$gene_id))
  if (verbose) cat("...done.\n") # }}}

  if (verbose) cat("Creating the database...") # {{{
  txVersion <- gsub("^v", "", version)
  outstub <- getTxDbLiteName(fastaFile)
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
  txcols <- c("start", "end", 
              "tx_id", "tx_length", "gc_content", "tx_biotype_id", "gene")
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

  if (verbose) cat("Writing the biotype_class table...") # {{{
  data(ensembl_biotypes, package="TxDbLite")
  if (class(ensembl_biotypes) != "data.frame") browser()
  dbWriteTable(con, name="biotype_class", ensembl_biotypes,
               overwrite=T, row.names=F)
  if (verbose) cat("done.\n") # }}}

  # write metadata table # {{{ 
  Metadata <- ensDbLiteMetadata(packageName=outstub, 
                                genomeVersion=txVersion,
                                sourceFile=fastaFile)
  dbWriteTable(con, name="metadata", Metadata, overwrite=TRUE, row.names=FALSE)
  # }}}

  # create indices  # {{{
  dbGetQuery(con, "create index tx_id_idx on tx (tx_id);")
  dbGetQuery(con, "create index gene_idx on gene (gene);")
  dbGetQuery(con, "create index txb_id_idx on tx_biotype (id);")
  dbGetQuery(con, "create index gxb_id_idx on gene_biotype (id);")
  dbGetQuery(con, "create index class_id_idx on biotype_class (class);")
  dbGetQuery(con, "create index biotype_id_idx on biotype_class (biotype);")
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
#' @param organism  what kind of organism these genes are from
#' 
#' @return  entrez_id values for the genes, where found
#'
#' @export
getEntrezIDs <- function(gxs, organism) { # {{{
  org <- getOrgDetails(organism)
  library(org$package, character.only=TRUE) 
  res <- try(mapIds(get(org$package), keys=gxs$gene_id, 
                column="ENTREZID", keytype=org$keytype), silent=TRUE)
  if (inherits(res, "try-error")) {
    warning("No ENTREZID mappings were found for these genes...")
    return(rep(NA, length(gxs)))
  } else { 
    return(res)
  }
} # }}}


#' @describeIn ensDbLiteFromFasta
#' 
#' add symbols for Ensembl genes
#'
#' 
#' 
#' @return  symbols for the genes, where found 
#'
#' @export
getSymbols <- function(gxs, organism) { # {{{
  org <- getOrgDetails(organism)
  library(org$package, character.only=TRUE) 

  ## needed since Ens83...
  if (any(grepl("\\.", gxs$gene_id))){
   gxs$gene_id<-gsub("\\.","",gxs$gene_id)
  }#no period in gene names.... 
     
  res <- try(mapIds(get(org$package), keys=gxs$gene_id, 
                    column=org$symbol, keytype=org$keytype), silent=TRUE)
  if (inherits(res, "try-error")) {
    warning("No SYMBOLS were found for these genes...")
    return(rep(NA, length(gxs)))
  } else { 
    return(res)
  }
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
#' @export
ensDbLiteMetadata <- function(packageName, genomeVersion, sourceFile) { # {{{

  tokens <- strsplit(getFastaStub(sourceFile), "\\.")[[1]]
  names(tokens)[1:3] <- c("organism", "genome", "version")
  organism <- getOrgDetails(tokens["organism"])
  build <- tokens["genome"]

  MetaData <- data.frame(matrix(ncol=2, nrow=8))
  colnames(MetaData) <- c("name", "value")
  MetaData[1,] <- c("package_name", packageName)
  MetaData[2,] <- c("db_type", "EnsDbLite")
  MetaData[3,] <- c("type_of_gene_id", "Ensembl Gene ID")
  MetaData[4,] <- c("created_by", paste("TxDbLite", packageVersion("TxDbLite")))
  MetaData[5,] <- c("creation_time", date())
  MetaData[6,] <- c("organism", organism$name)
  MetaData[7,] <- c("genome_build", build)
  MetaData[8,] <- c("source_file", sourceFile)
  return(MetaData)

} # }}}

