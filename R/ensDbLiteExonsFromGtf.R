#' as the name suggests, this fetches transcript-to-exon mappings from a GTF.
#' 
#' @param ensdb         an (existing) EnsDbLite object (with con, etc.)
#' @param gtf           the file name for the GTF file (perhaps compressed)
#'
ensDbLiteExonsFromGtf <- function(ensdb, gtf) { # {{{
 
  require(rtracklayer)
  options(useFancyQuotes=FALSE)
  if (!file.exists(gtf)) stop(paste("Could not find GTF file", gtf))
  header <- fetchHeader(gtf)
  if (header[header[,1] == "genome-version", "value"] != genomeVersion) {
    stop("The header from", gtf, "does not match the specified genome version!")
  }

  ## get tx2exon and exon tables
  exs <- fetchExons(gtf, verbose=verbose)
  tx2exon <- as(mcols(exs), "data.frame")
  swap <- c(transcript_id="tx_id", exon_number="exon_idx", exon_id="exon_id")
  names(tx2exon) <- swap[names(tx2exon)] 
  keepCols <- c("exon_id","start","end")
  exs <- as(exs[!duplicated(exs$exon_id)], "data.frame")[, keepCols]
  tbls <- list(exs=exs, tx2exon=tx2exon)
  return(tbls)

} # }}}

fetchHeader <- function(gtf, verbose=FALSE) { # {{{
  if(verbose) cat("Reading header from", gtf, "...")
  tmp <- readLines(gtf, n=10)
  tmp <- tmp[grep(tmp, pattern="^#")]
  if(length(tmp) > 0){
    tmp <- gsub(tmp, pattern="^#", replacement="")
    tmp <- gsub(tmp, pattern="^!", replacement="")
    Header <- do.call(rbind, strsplit(tmp, split=" ", fixed=TRUE))
    colnames(Header) <- c("name", "value")
    if(verbose) cat("done.\n")
    return(Header)
  } else {
    if(verbose) cat("failed!\n")
    stop("Unable to read header from", gtf)
  }
} # }}}

fetchExons <- function(gtf, verbose=FALSE) { # {{{
  if(verbose) cat("Importing exons from", gtf, "...")
  cols <- c("transcript_id", "exon_number", "exon_id")
  exs <- import(con=gtf, format="gtf", feature.type="exon")
  mcols(exs) <- mcols(exs)[, cols]
  if(verbose) cat("done.\n")
  return(exs)
} # }}}
