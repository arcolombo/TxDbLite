require(rtracklayer)
require(GenomicRanges)

#setwd("~/BAMs/")
rpts <- read.table("rmsk.hg38.txt.gz")
#colnames(rpts) <- 
 # c('chrom','chromStart','chromEnd','strand','name','class','family','NA1','NA2','NA3','NA4','NA5','NA6','NA7','NA8','NA9')

colnames(rpts)<-c('bin','swScore','milliDiv','milliDel','millIns','chrom','chromStart','chromEnd','genoLeft','strand','name','class','family','repStart','repEnd','repLeft','id')
#the hg19, UCSC is a 0-based genome, where we need the GrCh38 1-based genome 
rptsByClass <- GRangesList(lapply(split(rpts, rpts$class), 
                                  makeGRangesFromDataFrame,
                                  keep.extra.columns=TRUE,
                                  starts.in.df.are.0based=TRUE))

require(Homo.sapiens)
seqinfo(rptsByClass) <- seqinfo(Homo.sapiens)[seqlevels(rptsByClass)]
genome(seqinfo(rptsByClass))<-"hg38"
#saveRDS(rptsByClass, file="~/Dropbox/Giris_Data/stemness/repeatsByClass.rds")
