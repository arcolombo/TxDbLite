repPrep <- function(repeats) { 
  sapply(repeats, indexFa)
  merged <- do.call(c, lapply(humanRepeats, scanFaIndex))
  names(merged) <- make.unique(as.character(seqnames(merged)))
  return(merged)
}
