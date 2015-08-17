#' functions to create a RepDbLite object (a stripped-down SQLite db) from
#' specified RepBase FASTA files for a given version, species, and type(s)
#' 
#' @param organism      the organism involved (Homo_sapiens, Mus_musculus, etc.)
#' @param fastaFiles    the FASTA files to collate into a RepDbLite instance
#' @param version       the RepBase build (e.g. 20_07)
#' 
#' @export
#'
repDbLiteFromFasta <- function(organism, fastaFiles, version) {
 
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
  repVersion <- paste0("RepBase", gsub("Repbase","", ignore.case=TRUE, version))
  abbr <- c(Homo_sapiens="Hsapiens", Mus_musculus="Mmusculus")
  outstub <- paste("RepDbLite", abbr[organism], repVersion, sep=".")
  dbname <- paste(outstub, "sqlite", sep=".") 
  con <- dbConnect(dbDriver("SQLite"), dbname=dbname)

  ## write repeat table 
  rpt <- data.frame(seqnames=seqlevels(combinedSeqinfo),
                    start=rep(0, length(combinedSeqinfo)), 
                    end=as.data.frame(combinedSeqinfo)$seqlengths,
                    repeat_id=seqlevels(combinedSeqinfo))
  dbWriteTable(con, name="repeat", as.data.frame(rpt), overwrite=T, row.names=F)

  ## write metadata table 
  Metadata <- buildMetadataLite(organism, version, host="unknown",
                                sourceFile=paste(fastaFiles, collapse="+"),
                                genomeVersion=repVersion,
                                packageName=outstub)
  dbWriteTable(con, name="metadata", Metadata, overwrite=TRUE, row.names=FALSE)

  ## create indices 
  dbGetQuery(con, "create index repeat_id_idx on repeat (repeat_id);")

  ## write repeat_family table? (how?)
  message("repeat_family table not written")
  # 
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
