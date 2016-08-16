#' @rdname TxDbLite-class
#' @export
setMethod("dbconn", "TxDbLite", function(x) return(x@con))

#' @rdname TxDbLite-class
setMethod("show", "TxDbLite", function(object) { # {{{
  if(is.null(object@con)) stop(paste("Invalid", class(object), "instance!"))
  info <- metadata(object)
  cat(class(object), ":\n")
  catrow <- function(x) cat(paste0("|", x["name"], ": ", x["value"], "\n"))
  for(i in 1:nrow(info)) catrow(info[i,])
}) # }}}

#' @rdname TxDbLite-class
setMethod("metadata", "TxDbLite", function(x, ...) { # {{{
  md <- dbGetQuery(dbconn(x), "select * from metadata")
  rownames(md) <- md$name
  return(md)
}) # }}}

setMethod("transcripts", "TxDbLite", function(x) { # {{{
  res <- makeGRangesFromDataFrame(dbGetQuery(dbconn(x), "select * from tx"),
                                  keep.extra.columns=TRUE)
  genome(res) <- metadata(x)["genome_build","value"]
  names(res) <- res$tx_id
  return(res)
}) # }}}

setMethod("promoters", "TxDbLite", function(x,upstream=2000,downstream=200,...){ # {{{
  trim(suppressWarnings(promoters(transcripts(x, ...),
                                  upstream=upstream,
                                  downstream=downstream)))
}) # }}}

#' 
#' 
#' get transcripts by gene, promoter, tx_biotype, gene_biotype, or biotype_class
#'
#' @param x   the TxDbLite instance
#' @param by  how to split the transcripts 
#'
#' @return a GRangesList
#' @rdname TxDbLite-class
#' @export
setMethod("transcriptsBy", "TxDbLite", function(x, # {{{
                                                by=c("gene",
                                                     "promoter",
                                                     "tx_biotype",
                                                     "gene_biotype",
                                                     "biotype_class")) {
  by <- match.arg(by)
  txs <- transcripts(x)
  name <- function(x) paste0(seqnames(x),":",start(x),"-",end(x),":",strand(x)) 
  proms <- promoters(x)
  seqlevelsStyle(proms) <- "UCSC" # else names can potentially break downstream
  promoterNames <- name(proms)
  switch(by,
         gene=split(txs, txs$gene_id),
         promoter=split(txs, promoterNames),
         tx_biotype=split(txs, txs$tx_biotype),
         gene_biotype=split(txs, txs$gene_biotype),
         biotype_class=split(txs, txs$biotype_class))
}) # }}} 

setMethod("genes", "TxDbLite", function(x) { # {{{
  sql <- paste("select seqnames, start, end, strand, ",
               "       tx_length, 'NA' as gc_content, 'NA' as tx_id,",
               "       tx_id as gene_id, 'NA' as gene_name, 'NA' as entrezid, ",
               "       'NA' as tx_biotype, gene_biotype, biotype_class", 
               "  from tx")
  res <- makeGRangesFromDataFrame(dbGetQuery(dbconn(x), sql), 
                                  keep.extra.columns=TRUE)
  genome(res) <- metadata(x)["genome_build","value"]
  names(res) <- res$gene_id
  return(res)
}) # }}}

#' Generic for querying genes
#' @name genesBy
#' @rdname TxDbLite-class
#' @export
setGeneric("genesBy", function(x,by=c("gene_biotype","biotype_class"), ...){#{{{
             standardGeneric("genesBy")
           }) # }}}

#' 
#' 
#' get genes by gene_biotype or biotype_class
#'
#' @param x   the TxDbLite instance
#' @param by  how to split the genes 
#' @aliases genesBy TxDbLite-method
#' @return a GRangesList
#' @rdname TxDbLite-class
#' @export
setMethod("genesBy", "TxDbLite", function(x, by=c("gene_biotype","biotype_class")) { # {{{
  by <- match.arg(by)
  gxs <- genes(x)
  switch(by,
         gene_biotype=split(gxs, gxs$gene_biotype),
         biotype_class=split(gxs, gxs$biotype_class))
}) # }}} 

## EnsDbLite methods

setMethod("genes", "EnsDbLite", function(x) { # {{{
  sql <- paste("select seqnames, start, end, strand, ",
               "       median_length as tx_length, 'NA' as gc_content,",
               "       'NA' as tx_id, gene_id, gene_name, entrezid,",
               "       'NA' as tx_biotype, gene_biotype,",
               "       class as biotype_class", 
               "  from gene, gene_biotype, biotype_class",
               " where gene.gene_biotype_id = gene_biotype.id",
               "   and gene_biotype.gene_biotype = biotype_class.biotype",
               " order by gene_id asc")
  res <- makeGRangesFromDataFrame(dbGetQuery(dbconn(x), sql),
                                  keep.extra.columns=TRUE)
  genome(res) <- metadata(x)["genome_build","value"]
  names(res) <- res$gene_id
  return(res)
}) # }}}

setMethod("transcripts", "EnsDbLite", function(x) { # {{{
  sql <- paste("select gene.seqnames, tx.start, tx.end, gene.strand,",
               "       tx_length, gc_content, tx_id, gene_id, gene_name,",
               "       entrezid, tx_biotype, gene_biotype,",
               "       class as biotype_class",
               "  from gene, tx, gene_biotype, tx_biotype, biotype_class",
               " where gene.gene = tx.gene",
               "   and tx.tx_biotype_id = tx_biotype.id",
               "   and gene.gene_biotype_id = gene_biotype.id",
               "   and gene_biotype.gene_biotype = biotype_class.biotype",
               " order by tx_id asc")
  res <- makeGRangesFromDataFrame(dbGetQuery(dbconn(x), sql),
                                  keep.extra.columns=TRUE)
  genome(res) <- metadata(x)["genome_build","value"]
  names(res) <- res$tx_id
  return(res)
}) # }}}

setMethod("listGenebiotypes", "EnsDbLite", function(x, ...){ # {{{
  return(dbGetQuery(dbconn(x), "select * from gene_biotype")[,"gene_biotype"])
}) # }}}

setMethod("listTxbiotypes", "EnsDbLite", function(x, ...){ # {{{
  return(dbGetQuery(dbconn(x), "select * from tx_biotype")[,"tx_biotype"])
}) # }}}

setMethod("show", "EnsDbLite", function(object) { # {{{
  callNextMethod() # TxDbLite show method -- basic information on the db  
  genesql <- "select count(distinct gene) from gene"
  g <- dbGetQuery(dbconn(object), genesql)[1,1]
  txsql <- "select count(distinct tx_id) from tx"
  tx <- dbGetQuery(dbconn(object), txsql)[1,1]
  cat(paste0("| ", tx, " transcripts from ", g, " bundles (genes).\n"))
}) # }}}


## RepDbLite show method
setMethod("show", "RepDbLite", function(object) { # {{{
  callNextMethod() # TxDbLite show method -- basic information on the db  
  repsql <- "select count(distinct tx_id) from tx"
  famsql <- "select count(distinct tx_biotype) from tx"
  rpts <- dbGetQuery(dbconn(object), repsql)[1,1]
  fam <- dbGetQuery(dbconn(object), famsql)[1,1]
  cat(paste0("| ", rpts, " repeat exemplars from ", 
                   fam, " repeat families (no known genes).\n"))
}) # }}}

## RepDbLite objects have no genes in them
setMethod("genes", "RepDbLite", function(x) callNextMethod()[0] ) ## no genes
setMethod("promoters", "RepDbLite", function(x) callNextMethod()[0] ) ## none


## ErccDbLite show method
setMethod("show", "ErccDbLite", function(object) { # {{{
  callNextMethod() # TxDbLite show method -- basic information on the db  
  ctlsql <- "select count(distinct tx_id) from tx"
  grpsql <- "select count(distinct tx_biotype) from tx"
  ctl <- dbGetQuery(dbconn(object), ctlsql)[1,1]
  grp <- dbGetQuery(dbconn(object), grpsql)[1,1]
  ## subtract 1 from the number of subgroups as "unannotated" is in there
  cat(paste0("| ", ctl, " spike-in controls from ", 
                   grp - 1, " subgroups (no known genes).\n"))
}) # }}}

## ErccDbLite objects have no genes in them
setMethod("genes", "ErccDbLite", function(x) callNextMethod()[0] ) ## no genes
setMethod("promoters", "ErccDbLite", function(x) callNextMethod()[0] ) ## none
