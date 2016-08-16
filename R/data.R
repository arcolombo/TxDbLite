#' Structure for annotating transcript biotypes from ENSEMBL
#'
#' holds 48 transcript biotypes and 6 levels of classes for Homo Sapiens.
#'
#' 
#' @source \url{http://www.ensembl.org/}
"ensembl_biotypes"


#' Data frame of ENSEMBL transcript IDs and corresponding reactome pathway IDs, Reactome URL pathway-browser, reactome pathway name, and organism for numerous species including Mus musculus, Homo Sapiens and many more.  This is used for enrichment activation analysis.
#' @format A data frame full of reactome data
#'   \itemize{
#'    \item code: names of ENSEMBL IDs of transcripts and genes
#'    \item term: Reactome ID pathway
#'   \item URL: Reactome Pathway URL
#'  \item name: Pathway name
#' \item from: Pathway Origin
#' \item organism: organism of pathway
#' }
#' @source \url{www.reactome.org}
"ensemblToReactome"


#' Data frame of External RNA-Control-Consortium (ERCC) Ambion Spike-Ins Concentration Mixes
#' @source \url{https://www.thermofisher.com/order/catalog/product/4456740}
"ERCC_annotated"

#' A list of HUGO Gene Symbols and associated ENSEMBL Transcripts and Gene ID
#' @format A list of two, 
#' \itemize{
#'  \item transcript: a list of ENSEMBL Transcript IDs and associated HUGO Gene Symbols
#'  \item gene: a list of ENSEMBL Gene IDs and associated HUGO Gene Symbols.
#' }
#' @source \url{www.genenames.org}
"hugoCache"

#' A Data frame of repetitive elements relating to Mus Musculus
#' @format  A Data frame filled with mouse repeats
#' \itemize{
#'   \item name: names of repeat elements
#'   \item tx_biotype: transcript biotype of given repeat
#'   \item gene_biotype: gene biotype
#'   \item class: class of transcript
#'   }
#' @source \url{www.girinst.org}
"mouse_repeat_biotypes"


#' A list of ENSEMBL Transcripts and associated Reactome Pathway IDs
#' @source \url{www.reactome.org}
"reactomeCache"

#' A list of Reactome IDs and associated names for each reactome ID
#' @source \url{www.reactome.org}
"reactomePathways"

#' A Data frame of repetitive elements relating to Homo Sapiens
#' @format  A Data frame filled with human repeats
#' \itemize{
#'   \item name: names of repeat elements
#'   \item tx_biotype: transcript biotype of given repeat
#'   \item gene_biotype: gene biotype
#'   \item class: class of transcript
#'   }
#' @source \url{www.girinst.org}
"repeat_biotypes"

#' A Data Frame of all supported organisms and supported annotations available by TxDbLite
#' @format A Data Frame of information of supported species
#' \itemize{
#' \item name: names of species
#' \item abbreviation: supported species abbreviation
#' \item package: associated package for given species
#' \item keytype: species key
#' \item symbol: species symbol
#' \item txpre: transcrpit prefix for species
#' \item gxpre: gene prefix
#' \item reactome: reactome abbreviation
#' }
#' @source \url{www.ensembl.org}
"supportedOrganismsForTxDbLite"

#' GRanges object for entries without annotations
#'
#' @format A GRanges object for unannotated transcripts with 9 metadata columns
#' \itemize{
#'   \item seqnames: names of transcripts
#'   \item ranges: coordinates if known
#'   \item strand:  strand info 
#'   \item tx_length: transcript length not effective length
#'   \item gc_content: content of GC bases 
#'   \item tx_id: transcript id, or repeat name
#'   \item gene_id: ENSEMBL gene id or repeat name
#'   \item gene_name: HUGO gene name from entrez conversion
#'   \item entrez_id: entrez id 
#'   \item tx_biotype: transcript biotype
#'   \item gene_biotype:  gene biotype
#'   \item biotype_class: biotype class
#' }
#' @source \url{http://www.zeroskateboards.com/}
"unannotatedTranscript"

