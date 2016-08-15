#' HUGO transcript and gene names for ENSEMBL IDs of the corresponding type
#' 
#' @param   useCache  boolean (TRUE, and should stay that way, save for updates)
#'
#' @return  a list of HUGO names associated with ENSEMBL transcript and gene IDs
#' 
#' @importFrom  biomaRt useEnsembl
#' @importFrom biomaRt getBM
#'
#' @export
getHugoCache <- function(useCache=TRUE) {
  if (useCache == TRUE) { 
    if (!exists("hugoCache")) data(hugoCache, package="TxDbLite")
    return(hugoCache)
  } else { 
    # this should only be used to repopulate the cache!
    columns <- list(transcript=c("ensembl_transcript_id",
                                 "hgnc_transcript_name"),
                    gene=c("ensembl_gene_id",
                           "hgnc_symbol"))
    filters <- list(transcript=list(with_ens_hs_transcript=TRUE, 
                                    with_hgnc_transcript_name=TRUE),
                    gene=list(with_hgnc=TRUE))
    dset <- "hsapiens_gene_ensembl"
    caches <- lapply(c(transcript="transcript", gene="gene"), function(x) 
                     fetchFromEnsembl(columns[[x]], filters[[x]], dataset=dset))
    cache <- do.call(c, 
                     lapply(caches, function(cache) split(cache[,2],cache[,1])))
    names(cache) <- gsub("transcript\\.", "", names(cache))
    names(cache) <- gsub("gene\\.", "", names(cache))
    return(cache)
  }
}
