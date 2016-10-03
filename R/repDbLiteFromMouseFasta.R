#' functions to create a RepDbLite object (a stripped-down SQLite db) from
#' a specified RepBase FASTA file for a given version, species, and type(s). this is the mouse rep biotype to be used.
#' @param fastaFile     the FASTA file to collate into a RepDbLite instance
#' @param verbose       make a lot of noise? (TRUE) 
#' @param dryRun        boolean if false, a sql database is created
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Rsamtools indexFa scanFaIndex scanFa FaFile
#' @importFrom Biostrings fasta.seqlengths
#' @importFrom DBI dbConnect dbDriver dbWriteTable dbGetQuery dbDisconnect
#' @export
repDbLiteFromMouseFasta <- function(fastaFile, verbose=TRUE, dryRun=FALSE) {

#FIX ME support mouse uncataloged
# FIX ME add three columns for reactome Ensembl support
 
  verbose <- TRUE 
  options(useFancyQuotes=FALSE)
  fastaStub <- getFastaStub(fastaFile)
  tokens <- strsplit(fastaStub, "\\.")[[1]]
  version <- tokens[3]
  organism <- tokens[1] 
  organism <- sub("\\.", "_", ## try & be robust
                  sub("Mmusculus", "Mus_musculus", 
                      sub("Hsapiens", "Homo_sapiens", organism)))
 



  ## utility functions for parsing FASTAs
  splt <- function(x, y=" ") strsplit(x, y)[[1]]
  popx <- function(x) x[length(x)] 
  pop <- function(x, y=" ") popx(splt(x, y))
  grab <- function(x, y=" ", i=1) splt(x, y)[i]
  shift <- function(x, y=" ") grab(x, y, i=1)
  shiftall <- function(x, y=" ") sapply(x, shift, y=y)

  indexFa(FaFile(fastaFile))
  txs <- scanFaIndex(FaFile(fastaFile))
  names(txs) <- seqnames(txs)

  if (verbose) cat("Extracting repeat lengths...") # {{{
  txLen <- fasta.seqlengths(fastaFile)
  names(txLen) <- sapply(sapply(names(txLen), shift), shift, y="\\t")
  stopifnot(identical(names(txs), names(txLen)))
  txs$tx_length <- txLen
  if (verbose) cat("done.\n") # }}}
 
  if (verbose) cat("Extracting repeat descriptions...") # {{{
  data(mouse_repeat_biotypes, package="TxDbLite")
  uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
  if (length(uncataloged) > 1) { 

      message(length(uncataloged), " uncataloged repeat biotypes, fix case...")
      hit <- !is.na(match(tolower(uncataloged), 
                          tolower(rownames(mouse_repeat_biotypes))))
      rows <- rownames(mouse_repeat_biotypes)
      hits <- uncataloged[hit] 
      cased <- mouse_repeat_biotypes[match(tolower(hits), tolower(rows)),]
      rownames(cased) <- hits
      cased$name <- hits
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, cased)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))

      message(length(uncataloged), 
              " uncataloged repeat biotypes, fix Tiggers...")
      Tigger <- mouse_repeat_biotypes["Tigger1", , drop=FALSE]
      Tiggers <- uncataloged[which(toupper(substr(uncataloged, 1, 6)) == 
                                   "TIGGER")]
      makeTigger <- function(name) return(data.frame(name=name, Tigger[,-1]))

      if (length(Tiggers)>0){
      TiggerRows <- do.call(rbind, lapply(Tiggers, makeTigger))
      rownames(TiggerRows) <- Tiggers
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, TiggerRows)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
     }
      message(length(uncataloged), " uncataloged mouse repeat biotypes, fix Alus...")
      Alu <- mouse_repeat_biotypes["Alu", , drop=FALSE]
      Alus <- uncataloged[which(substr(uncataloged, 1, 3) == "Alu")]
      makeAlu <- function(name) return(data.frame(name=name, Alu[,-1]))
      if(length(Alus)==0){
     message("Alus were not found in uncataloged mouse repeats, skipping ...")
     } 
     if(length(Alus)>0) {
      AluRows <- do.call(rbind, lapply(Alus, makeAlu))
      rownames(AluRows) <- Alus
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, AluRows)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
         }#if Alus are found

  #    message(length(uncataloged), " uncataloged repeat biotypes, fix HERVs...")
   #   HERVs <- uncataloged[which(substr(uncataloged, 1, 4) == "HERV")]
   #   HERVi <- paste0(HERVs[paste0(HERVs, "-int") %in% 
  #                    rownames(repeat_biotypes)], "-int")
   #   HERVis <- repeat_biotypes[HERVi, ]
    #  rownames(HERVis) <- sub("-int", "", rownames(HERVis))
  #    HERVis$name <- sub("-int", "", HERVis$name)
  #    repeat_biotypes <- rbind(repeat_biotypes, HERVis)
  #    repeat_biotypes <- repeat_biotypes[order(repeat_biotypes$name),] 
  #    uncataloged <- setdiff(names(txs), rownames(repeat_biotypes))
#
  #    message(length(uncataloged)," uncataloged repeat biotypes, more HERVs...")
  #    HERVs <- HERVs[!paste0(HERVs, "-int") %in% rownames(repeat_biotypes)]
  #    HERV <- repeat_biotypes["HERV1_I-int", , drop=FALSE]
  #    makeHERV <- function(name) return(data.frame(name=name, HERV[,-1]))
  #    HERVRows <- do.call(rbind, lapply(HERVs, makeHERV))
  #    rownames(HERVRows) <- HERVs
  #    repeat_biotypes <- rbind(repeat_biotypes, HERVRows)
  #    repeat_biotypes <- repeat_biotypes[order(repeat_biotypes$name),] 
  #    uncataloged <- setdiff(names(txs), rownames(repeat_biotypes))
  # there ain't no HERVs for mouse
    
      message(length(uncataloged), " uncataloged mouse repeat biotypes, fix LINE1...")
      L1s <- uncataloged[which(substr(uncataloged, 1, 2) == "L1")]
      LINE1 <- mouse_repeat_biotypes["L1M",]
      makeLINE1 <- function(name) return(data.frame(name=name, LINE1[,-1]))
      if(length(L1s)==0){
     message("uncataloged mouse repeats did not contain L1s...skipping")
     } 
     if(length(L1s)>0) {
      LINE1Rows <- do.call(rbind, lapply(L1s, makeLINE1))
      rownames(LINE1Rows) <- L1s
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, LINE1Rows)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
      }
 
      message(length(uncataloged), " uncataloged mouse repeat biotypes, fix MERs...")
      MERs <- uncataloged[which(substr(uncataloged, 1, 3) == "MER")]
      MER <- mouse_repeat_biotypes["MER101",]
      makeMER <- function(name) return(data.frame(name=name, MER[,-1]))
      if(length(MERs)>0) {
      MERRows <- do.call(rbind, lapply(MERs, makeMER))
      rownames(MERRows) <- MERs
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, MERRows)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
   }


      message(length(uncataloged), " uncataloged mouse repeat biotypes, fix LTRs...")
      LTRs <- uncataloged[which(substr(uncataloged, 1, 3) == "LTR")]
      LTR <- mouse_repeat_biotypes["LTR1",]
      makeLTR <- function(name) return(data.frame(name=name, LTR[,-1]))
     if(length(LTRs)==0) {
      message("No LTRs were found in the uncataloged repeats...skipping")
       }
      if(length(LTRs)>0) {
      LTRRows <- do.call(rbind, lapply(LTRs, makeLTR))
      rownames(LTRRows) <- LTRs
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, LTRRows)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
      }

      message(length(uncataloged), " uncataloged mouse repeat biotypes, fix SVAs...")
      SVAs <- uncataloged[which(substr(uncataloged, 1, 3) == "SVA")]
      SVA <- mouse_repeat_biotypes["SVA_A",]
      makeSVA <- function(name) return(data.frame(name=name, SVA[,-1]))
     if(length(SVAs)==0){ 
     message("SVAs were not found in the uncataloged mouse ...")
     } 
     if(length(SVAs)>0){
      SVARows <- do.call(rbind, lapply(SVAs, makeSVA))
      rownames(SVARows) <- SVAs
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, SVARows)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
     }

      message(length(uncataloged), " uncataloged mouse repeat biotypes, fix SINEs...")
      SINEs <- uncataloged[grep("(FLA|FAM|FLAM|FRAM)", uncataloged)]
      SINE <- mouse_repeat_biotypes["FRAM",]
      makeSINE <- function(name) return(data.frame(name=name, SINE[,-1]))
      if(length(SINEs)==0){
      message("SINEs were not found in the repeat fasta for mouse ...")
     }
     if(length(SINEs)>0){
      SINERows <- do.call(rbind, lapply(SINEs, makeSINE))
      rownames(SINERows) <- SINEs
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, SINERows)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
      }


      message(length(uncataloged),
              " uncataloged mouse repeat biotypes, fix Mariners...")
      MARs <- uncataloged[which(substr(uncataloged, 1, 3) == "MAR")]
      MAR <- mouse_repeat_biotypes["MARNA",] 
      makeMAR <- function(name) return(data.frame(name=name, MAR[,-1]))
     if(length(MARs)>0){
       MARRows <- do.call(rbind, lapply(MARs, makeMAR))
      rownames(MARRows) <- MARs
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, MARRows)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
   }
      ## fall back on hinting...?
      message(length(uncataloged), " uncataloged mouse repeat biotypes... hinting...")
      txDesc <- names(fasta.seqlengths(fastaFile))
      names(txDesc) <- shiftall(names(txLen), y="\\t")
      txUncat <- txDesc[uncataloged] 
      txHints <- sapply(sapply(txUncat, strsplit, "\\t"), `[`, 2)
      hinted <- data.frame(name=names(txHints),
                           tx_biotype=txHints,
                           gene_biotype=NA,
                           class="repeat")
      hinted$gene_biotype[grep("DNA", hinted$name)] <- "DNA_element"
      hinted$gene_biotype[grep("DNA", hinted$tx_biotype)] <- "DNA_element"
      hinted$gene_biotype[grep("hAT", hinted$tx_biotype)] <- "DNA_element"
      hinted$gene_biotype[grep("Mariner", hinted$tx_biotype)] <- "DNA_element"
      hinted$gene_biotype[grep("Transposable", 
                               hinted$tx_biotype)] <- "DNA_element"
      hinted$gene_biotype[grep("EUT",ignore.case=T, 
                               hinted$tx_biotype)] <- "DNA_element"
      hinted$gene_biotype[grep("SAT",ignore.case=T, 
                               hinted$tx_biotype)] <- "other_repeat"
      hinted$gene_biotype[grep("L1", hinted$name)] <- "LINE"
      hinted$gene_biotype[grep("L2", hinted$name)] <- "LINE"
      hinted$gene_biotype[grep("LINE", hinted$name)] <- "LINE"
      hinted$gene_biotype[grep("L1", hinted$tx_biotype)] <- "LINE"
      hinted$gene_biotype[grep("L2", hinted$tx_biotype)] <- "LINE"
      hinted$gene_biotype[grep("LINE", hinted$tx_biotype)] <- "LINE"
      hinted$gene_biotype[grep("SINE", hinted$name)] <- "SINE"
      hinted$gene_biotype[grep("SINE", hinted$tx_biotype)] <- "SINE"
      hinted$gene_biotype[grep("Repetitive", 
                               hinted$tx_biotype)] <- "other_repeat"
      hinted$gene_biotype[grep("Endogenous", 
                               hinted$tx_biotype)] <- "LTR_element"
      hinted$gene_biotype[grep("ERV", hinted$tx_biotype)] <- "LTR_element"
      hinted$gene_biotype[grep("LTR", hinted$tx_biotype)] <- "LTR_element"
      hinted$gene_biotype[is.na(hinted$gene_biotype)] <- "other_repeat"
      mouse_repeat_biotypes <- rbind(mouse_repeat_biotypes, hinted)
      mouse_repeat_biotypes <- mouse_repeat_biotypes[order(mouse_repeat_biotypes$name),] 
      uncataloged <- setdiff(names(txs), rownames(mouse_repeat_biotypes))
      message(length(uncataloged), 
              " uncataloged repeat biotypes after hinting.")

      ## collapse tx_biotypes with "?" into their parent & tidy up the rest 
      mouse_repeat_biotypes$tx_biotype <- shiftall(mouse_repeat_biotypes$tx_biotype, y="-")
      mouse_repeat_biotypes$tx_biotype <- sub("\\?$", "", mouse_repeat_biotypes$tx_biotype)
      mouse_repeat_biotypes$tx_biotype <- sub("acro$", "acromeric", 
                                        mouse_repeat_biotypes$tx_biotype)
      mouse_repeat_biotypes$tx_biotype <- sub("centr$", "centromeric", 
                                        mouse_repeat_biotypes$tx_biotype)
  }
  if (verbose) cat("done.\n") # }}}

  if(verbose) cat("Extracting copy number...")
   if(grepl("musculus",organism)==TRUE){
   data("mouseCopyNumber.mm10",package="TxDbLite")
      for(i in 1:nrow(repeat_biotypes)){
        id<-which(rownames(repeat_biotypes)[i]==rownames(mouseCopyNumber.mm10))
         if(length(id)>0){
         repeat_biotypes$copyNumber[i]<-mouseCopyNumber.mm10$count[id]
         } else if(length(id)<=0){
           repeat_biotypes$copyNumber[i]<-1
            }
        } #for each biotype
    }
  ## create the SQLite database... 
  repVersion <- gsub("Repbase","", ignore.case=TRUE, version)
  outstub <- getTxDbLiteName(basename(fastaFile))
  dbname <- paste(outstub, "sqlite", sep=".") 
  if(dryRun==FALSE){
  con <- dbConnect(dbDriver("SQLite"), dbname=dbname)
  }
  if (verbose) cat("Creating the database...") # {{{
  rpt <- as.data.frame(txs)
  rpt$tx_length <- txLen
  rpt$gc_content <- GCcontent(scanFa(fastaFile))
  rpt$tx_id <- names(txLen)
  rpt$gene_id <- names(txLen)
  rpt$gene_name <- NA 
  rpt$entrezid <- NA 
  rpt$tx_biotype <- mouse_repeat_biotypes[rownames(rpt), "tx_biotype"]
  rpt$tx_biotype_id<-as.numeric(as.factor(rpt$tx_biotype))
  rpt$gene_biotype <- mouse_repeat_biotypes[rownames(rpt), "gene_biotype"]
  rpt$gene_biotype_id<-as.numeric(as.factor(rpt$gene_biotype))
  rpt$biotype_class <- "repeat"
  rpt$copyNumber<-repeat_biotypes[rownames(rpt),"copyNumber"]
  rpt$gene<-rpt$seqnames
  rpt$median_length<-rpt$tx_length

  if(verbose) cat("Writing the gene table...\n")
  gxcols <- c("seqnames", "start", "end", "strand",
              "gene","gene_id","gene_biotype_id",
               "entrezid","gene_name","median_length","copyNumber")
  if(dryRun==FALSE){
  dbWriteTable(con, name="gene", rpt[,gxcols], overwrite=T, row.names=F)
  } 
  if (verbose) cat("done.\n") # }}}

  if(verbose) cat("Tabulating gene biotypes...\n")
   gene_biotypes<-data.frame(id=seq_along(levels(as.factor(rpt$gene_biotype))),
                            gene_biotype=levels(as.factor(rpt$gene_biotype)))
  if(verbose) cat("Writing the gene_biotype table...")
  if(dryRun==FALSE){
  dbWriteTable(con,name="gene_biotype",gene_biotypes,overwrite=T,row.names=F)
  }
  if(verbose) cat("done.\n")

 if(verbose) cat("Writing the tx table...\n") #{{{
  txcols <- c("start", "end", "tx_id","tx_length","gc_content",
              "tx_biotype_id", "gene","copyNumber")
  if(dryRun==FALSE){
  dbWriteTable(con, name="tx", rpt[,txcols], overwrite=T, row.names=F)
   }
  if (verbose) cat("done.\n") # }}}

 if(verbose) cat("Tabulating transcript biotypes...\n")
  if(dryRun==FALSE){
    if(verbose) cat("Writing the transcript bio types...") #{{{
   tx_biotypes<-data.frame(id=seq_along(levels(as.factor(rpt$tx_biotype))),
                           tx_biotype=levels(as.factor(rpt$tx_biotype)))
    dbWriteTable(con,name="tx_biotype",tx_biotypes,overwrite=T,row.names=F)
   #}}}
   }

 if(verbose) cat("Writing the biotype_class table...\n")
  if(dryRun==FALSE){
   biotype_class=data.frame(biotype=levels(as.factor(rpt$gene_biotype)),
                           class="repeat")
   dbWriteTable(con,name="biotype_class",biotype_class,overwrite=T,row.names=F)
 }


  ## write metadata table 
  Metadata <- repDbLiteMetaMousedata(outstub, sourceFile=fastaFile)
  if(dryRun==FALSE){
  dbWriteTable(con, name="metadata", Metadata, overwrite=TRUE, row.names=FALSE)

  ## create indices 
  dbGetQuery(con, "create index tx_id_idx on tx (tx_id);")
  dbGetQuery(con, "create index gene_idx on gene (gene);")
  dbGetQuery(con, "create index txb_id_idx on tx_biotype (id);")
  dbGetQuery(con, "create index tx_biotype_idx on tx (tx_biotype_id);")
  dbGetQuery(con, "create index gxb_id_idx on gene_biotype (id);")
  dbGetQuery(con, "create index class_id_idx on biotype_class (class);")
  dbGetQuery(con, "create index biotype_id_idx on biotype_class (biotype);")
  dbGetQuery(con, "create index copyNumber_idx on tx (copyNumber);")

  ## finish 
  dbDisconnect(con)
  }
  return(dbname)

}


#' @describeIn repDbLiteFromMouseFasta
#' 
#' create metadata for a RepDbLite instance
#'
#' @param packageName   the name of the annotation package to be built 
#' @param sourceFile    name of FASTA file(s) whence it was built, as a string
#' 
#' @return a data.frame of metadata suitable for cramming into the database
#'
repDbLiteMetaMousedata <- function(packageName, sourceFile) { # {{{

  tokens <- strsplit(getFastaStub(sourceFile), "\\.")[[1]]
  names(tokens)[1:3] <- c("organism", "genome", "version")
  organism <- getOrgDetails(tokens["organism"])
  build <- paste0(tokens["genome"], tokens["version"])
  version <- tokens["version"]

  MetaData <- data.frame(matrix(ncol=2, nrow=8))
  colnames(MetaData) <- c("name", "value")
  MetaData[1,] <- c("package_name", packageName)
  MetaData[2,] <- c("db_type", "RepDbLite")
  MetaData[3,] <- c("type_of_gene_id", "RepBase identifiers")
  MetaData[4,] <- c("created_by", paste("TxDbLite", packageVersion("TxDbLite")))
  MetaData[5,] <- c("creation_time", date())
  MetaData[6,] <- c("organism", organism$name)
  MetaData[7,] <- c("genome_build", build)
  MetaData[8,] <- c("source_file", sourceFile)
  return(MetaData)

} # }}}
