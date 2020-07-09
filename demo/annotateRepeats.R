library(arkas) # EISAI release or later
library(arkasData) # any release should do 
library(RepDbLite.Hsapiens.2005) # or generate
library(ComplexHeatmap)

data(NS) # normal vs senescent HSPCs from Capone, Connor, et al 
DE_repeats <- read.csv(system.file("extdata", "DE_repeats.csv", 
                                   package="TxDbLite", mustWork=TRUE),
                       row.names=1)
annotations <- transcripts(RepDbLite.Hsapiens.2005)[rownames(DE_repeats)]

DE_repeat_tpm <- tpm(NS)[rownames(DE_repeats), ]
repeat_colors <- list(repeat_class=c(DNA_element="purple",
                                     LINE="red", SINE="orange",
                                     LTR_element="green",
                                     other_repeat="gray50"),
                      repeat_family=c(Alu="orange",
                                      "DNA transposon"="violet",
                                      "Endogenous Retrovirus"="green",
                                      ERVK="darkgreen",
                                      ERV1="forestgreen",
                                      ERVL="green",
                                      hAT="purple",
                                      L1="red",
                                      "LTR Retrotransposon"="seagreen",
                                      TcMar="violet",
                                      SAT="gray40",
                                      SVA="gray60"))

rpt_anno <- 
  rowAnnotation(data.frame(repeat_class=annotations$gene_biotype,
                           repeat_family=annotations$tx_biotype),
                col=repeat_colors)
anno <- HeatmapAnnotation(data.frame(senescent=ifelse(NS$treatment, "S", "N")),
                          col=list(senescent=c(S="gray50", N="white")))

# plot 
rpt_anno + 
Heatmap(log1p(DE_repeat_tpm), 
        row_title_side="left",
        show_row_dend=FALSE, 
        top_annotation=anno,
        row_names_gp=gpar(fontsize=9),
        name="log(1+TPM)")

bySd <- function(x, k=25) { 
  sds <- apply(x, 1, sd, na.rm=TRUE)
  x[rev(order(sds))[seq_len(k)], ]
}

top10sd <- bySd(log1p(DE_repeat_tpm), k=10)
annoTop10 <- annotations[rownames(top10sd)]
rpt_anno_10 <- 
  rowAnnotation(data.frame(repeat_class=annoTop10$gene_biotype,
                           repeat_family=annoTop10$tx_biotype),
                col=repeat_colors)

# plot again
rpt_anno_10 + 
Heatmap(top10sd, 
        row_title_side="left",
        show_row_dend=FALSE, 
        top_annotation=anno, 
        row_names_gp=gpar(fontsize=9),
        name="log(1+TPM)")
