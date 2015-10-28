library(TxDbLite)
library(artemisData)

fasta_dir <- system.file("extdata", "fasta", package="artemisData")
repbase_files <- list.files(fasta_dir, pattern="^Homo_sapiens.RepBase.*.gz$")
findDupes(repbase_files)
