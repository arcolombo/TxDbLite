#' functions to create a RepDbLite object (a stripped-down SQLite db) from
#' a specified RepBase FASTA file for a given version, species, and type(s)
#' 
#' @param organism      the organism involved (Homo_sapiens, Mus_musculus, etc.)
#' @param fastaFile     the FASTA file to collate into a RepDbLite instance
#' @param version       the RepBase build (e.g. 20_07)
#' 
#' @export
#'
repDbLiteFromFasta <- function(organism, fastaFile, version) {
 
  verbose <- TRUE 
  options(useFancyQuotes=FALSE)
  organism <- sub("\\.", "_", ## try & be robust
                  sub("Mmusculus", "Mus_musculus", 
                      sub("Hsapiens", "Homo_sapiens", organism)))
  sapply(fastaFiles, indexFa)
  names(fastaFiles) <- fastaFiles
  faFiles <- FaFileList(lapply(fastaFiles, FaFile))

  ## this is only really relevant for RepBase repeat annots
  combinedSeqinfo <- Reduce(merge, lapply(faFiles, seqinfo)) 

  ## create the SQLite database... 
  repVersion <- gsub("Repbase","", ignore.case=TRUE, version)
  outstub <- getTxDbLiteName(fastaFile)
  dbname <- paste(outstub, "sqlite", sep=".") 
  con <- dbConnect(dbDriver("SQLite"), dbname=dbname)

  ## write repeat table 
  rpt <- data.frame(seqnames=seqlevels(combinedSeqinfo),
                    start=rep(0, length(combinedSeqinfo)), 
                    end=as.data.frame(combinedSeqinfo)$seqlengths,
                    repeat_id=seqlevels(combinedSeqinfo))
  dbWriteTable(con, name="repeat", as.data.frame(rpt), overwrite=T, row.names=F)

  ## write metadata table 
  Metadata <- repDbLiteMetadata(outstub, genomeVersion=repVersion,
                                sourceFile=fastaFile)

  dbWriteTable(con, name="metadata", Metadata, overwrite=TRUE, row.names=FALSE)

  ## create indices 
  dbGetQuery(con, "create index repeat_id_idx on repeat (repeat_id);")

  ## write repeat_family table? (how?)
  message("repeat_family table not written")
  # 

  if (verbose) cat("Writing the biotype_class table...") # {{{
  data(repeat_biotypes, package="TxDbLite")
  dbWriteTable(con, name="biotype_class", as(repeat_biotypes, "data.frame"), 
               overwrite=T, row.names=F)
  rm(repeat_biotypes)
  if (verbose) cat("done.\n") # }}}

  ## more indices
  # 
  # dbGetQuery(con, paste("create index repeat_family_id_idx",
  #                       "on repeat_family (family_id);")
  # dbGetQuery(con, paste("create index repeat_family_repeat_id_idx",
  #                       "on repeat_family (repeat_id);")

  ## finish 
  dbDisconnect(con)
  return(dbname)

} # }}}


#' @describeIn repDbLiteFromFasta
#' 
#' create metadata for a RepDbLite instance
#'
#' @param packageName   the name of the annotation package to be built 
#' @param genomeVersion name of genome assembly for coordinates, e.g. "20_07"
#' @param sourceFile    name of FASTA file(s) whence it was built, as a string
#' 
#' @return a data.frame of metadata suitable for cramming into the database
#'
repDbLiteMetadata <- function(packageName, genomeVersion, sourceFile) { # {{{

  tokens <- strsplit(packageName, "\\.")[[1]]
  organism <- getOrganismAbbreviation(tokens[2])

  MetaData <- data.frame(matrix(ncol=2, nrow=8))
  colnames(MetaData) <- c("name", "value")
  MetaData[1,] <- c("package_name", packageName)
  MetaData[2,] <- c("db_type", "RepDbLite")
  MetaData[3,] <- c("type_of_gene_id", "RepBase identifiers")
  MetaData[4,] <- c("created_by", paste("TxDbLite", packageVersion("TxDbLite")))
  MetaData[5,] <- c("creation_time", date())
  MetaData[6,] <- c("organism", organism )
  MetaData[7,] <- c("genome_build", genomeVersion)
  MetaData[8,] <- c("source_file", sourceFile)
  return(MetaData)

} # }}}
