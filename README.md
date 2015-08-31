TxDbLite is a (very) minimalist annotation package generator, designed to extract annotations from FASTA files.  

The underlying assumption is that users want to quantify transcript-level RNA expression with Kallisto, may or may not have a GTF file for each transcriptome, and would like to extract as much information as possible about the transcripts from what they do have.  This in turn allows the artemis package to automatically determine what transcriptomes a dataset came from, whether those are already known to the package, and how to generate certain types of annotation for specific transcriptomes from Ensembl, RepBase, and ERCC (spike-in controls from Ambion). 

The user presumably already has the fasta file in question; typically the TxDbLite package will be automatically invoked during the creation of a new Kallisto index.  If a fasta file "smells like" a recognizable class of transcripts, and the user does not already appear to have installed or created annotations for the transcriptome, artemis will then attempt to do so.  This can also be performed manually:

The below example takes a recent Ensembl build for Drosophila melanogaster (the fruit fly) and processes it.   
The FASTA file is at 
ftp://ftp.ensembl.org/pub/release-81/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.ncrna.fa.gz and should be downloaded into your working directory. 

```r 
fasta <- "Drosophila_melanogaster.BDGP6.ncrna.fa.gz"
ensdb <- ensDbLiteFromFasta(fasta)
show(EnsDbLite(ensdb))
```

```
EnsDbLite :
|package_name: EnsDbLite.Drosophila_melanogaster.BDGP6.ncrna
|db_type: EnsDbLite
|type_of_gene_id: Ensembl Gene ID
|created_by: TxDbLite 1.9.16
|creation_time: Sun Aug 30 15:48:41 2015
|organism: Dmelanogaster
|genome_build: ncrna
|source_file: Drosophila_melanogaster.BDGP6.ncrna.fa.gz
| 4098 transcripts from 3384 bundles (genes).
```

Installation: 

```r
biocLite("devtools")
biocLite("RamsinghLab/TxDbLite")
```
