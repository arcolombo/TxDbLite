## TxDbLite methods (inherited by all subclasses)

setMethod("dbconn", "TxDbLite", function(x) return(x@con))

setMethod("show", "TxDbLite", function(object) { # {{{
  if(is.null(object@con)) stop(paste("Invalid", class(object), "instance!"))
  info <- metadata(object)
  cat(class(object), ":\n")
  catrow <- function(x) cat(paste0("|", x["name"], ": ", x["value"], "\n"))
  for(i in 1:nrow(info)) catrow(info[i,])
}) # }}}

setMethod("metadata", "TxDbLite", function(x, ...) { # {{{
  md <- dbGetQuery(dbconn(x), "select * from metadata")
  rownames(md) <- md$name
  return(md)
}) # }}}

setMethod("transcripts", "TxDbLite", function(x) { # {{{
  res <- makeGRangesFromDataFrame(dbGetQuery(dbconn(x), "select * from tx"),
                                  keep.extra.columns=TRUE)
  names(res) <- res$tx_id
  return(res)
}) # }}}


## EnsDbLite methods

setMethod("genes", "EnsDbLite", function(x) { # {{{
  sql <- paste("select seqnames, start, end, strand, median_length,",
               "       gene_id, gene_name, gene_biotype, ",
               "       class as biotype_class, entrezid", 
               "  from gene, gene_biotype, biotype_class",
               " where gene.gene_biotype_id = gene_biotype.id",
               "   and gene_biotype.gene_biotype = biotype_class.biotype",
               " order by gene_id asc")
  res <- makeGRangesFromDataFrame(dbGetQuery(dbconn(x), sql),
                                  keep.extra.columns=TRUE)
  names(res) <- res$gene_id
  return(res)
}) # }}}

setMethod("transcripts", "EnsDbLite", function(x) { # {{{
  sql <- paste("select gene.seqnames, tx.start, tx.end, gene.strand,",
               "       tx_length, tx_id, gene_id, gene_name, entrezid,",
               "       tx_biotype, gene_biotype, class as biotype_class",
               "  from gene, tx, gene_biotype, tx_biotype, biotype_class",
               " where gene.gene = tx.gene",
               "   and tx.tx_biotype_id = tx_biotype.id",
               "   and gene.gene_biotype_id = gene_biotype.id",
               "   and gene_biotype.gene_biotype = biotype_class.biotype",
               " order by tx_id asc")
  res <- makeGRangesFromDataFrame(dbGetQuery(dbconn(x), sql),
                                  keep.extra.columns=TRUE)
  names(res) <- res$tx_id
  return(res)
}) # }}}

setMethod("transcriptsBy", "EnsDbLite", function(x, by=c("gene")) { # {{{
  txs <- transcripts(x)
  split(txs, txs$gene_id)
}) # }}} 

setMethod("listGenebiotypes", "EnsDbLite", function(x, ...){ # {{{
  return(dbGetQuery(dbconn(x), "select * from gene_biotype")[,"gene_biotype"])
}) # }}}

setMethod("listTxbiotypes", "EnsDbLite", function(x, ...){ # {{{
  return(dbGetQuery(dbconn(x), "select * from tx_biotype")[,"tx_biotype"])
}) # }}}

setMethod("promoters", "EnsDbLite", function(x, upstream=2000, downstream=200, ...) { # {{{
  trim(suppressWarnings(promoters(transcripts(x, ...),
                                  upstream=upstream,
                                  downstream=downstream)))
}) # }}}

setMethod("show", "EnsDbLite", function(object) { # {{{
  callNextMethod() ## TxDbLite show method -- basic information on the db  
  genesql <- "select count(distinct gene) from gene"
  g <- dbGetQuery(dbconn(object), genesql)[1,1]
  txsql <- "select count(distinct tx_id) from tx"
  tx <- dbGetQuery(dbconn(object), txsql)[1,1]
  cat(paste0("| ", tx, " transcripts from ", g, " bundles (genes).\n"))
}) # }}}


## RepDbLite show method

setMethod("show", "RepDbLite", function(object) { # {{{
  callNextMethod() ## TxDbLite show method -- basic information on the db  
  famsql <- "select count(distinct tx_biotype) from tx"
  repsql <- "select count(distinct tx_id) from tx"
  rpts <- dbGetQuery(dbconn(object), repsql)[1,1]
  fam <- dbGetQuery(dbconn(object), famsql)[1,1]
  cat(paste0("| ", rpts, " repeat exemplars from ", fam, " repeat families.\n"))
}) # }}}


## ErccDbLite show method

setMethod("show", "ErccDbLite", function(object) { # {{{
  callNextMethod() ## TxDbLite show method -- basic information on the db  
  ctlsql <- "select count(distinct tx_id) from tx"
  grpsql <- "select count(distinct tx_biotype) from tx"
  ctl <- dbGetQuery(dbconn(object), ctlsql)[1,1]
  grp <- dbGetQuery(dbconn(object), grpsql)[1,1]
  ## subtract 1 from the number of subgroups as "unannotated" is in there
  cat(paste0("| ", ctl, " spike-in controls from ", grp - 1, " subgroups.\n"))
}) # }}}
