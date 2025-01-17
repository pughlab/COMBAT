library(dplyr)
library(ggplot2)
library(readr)
library(BoutrosLab.plotting.general);
library(annotatr);
library(GenomicRanges);
load("/Users/arnavaz/Library/CloudStorage/OneDrive-UHN/Projects/COMBAT/hg38_gene_annotations.RData")
load("/Users/arnavaz/Library/CloudStorage/OneDrive-UHN/Projects/COMBAT/hg38_cpg_annotations.RData")
COMBAT_meth = read_tsv("/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/COMBAT_rpkm.txt")
COMBAT_meth = COMBAT_meth %>% mutate(interval = paste0(chr, "_", start, "_", end))
combat_meth = COMBAT_meth[,4:22]


intervals = read_tsv()
hg38.annotations <- annotate_regions(
  regions = hg38_gene_annotations,
  annotations = hg38_cpg_annotations
)


genes = c("TLR5","C1orf61","NR2F1-AS1" ,
          "ZDHHC1" ,"OTP" ,"FEZF2" ,"METAP1D" ,"TSHZ3",
           "THEG5" ,"FLJ45513","DLX4" ,"POU3F1" ,"NKX2-1-AS1",
            "NKX2-8" ,"PNPLA1" ,"TIFAB","NEUROG1" ,"RADIL",
             "EFCAB10" , "MAST1" ,"ZFHX4-AS1" ,"POLR1A",
           "CAMKMT", "SIX3-AS1" ,"HOXB13","TTLL6" ,
           "LHX5-AS1","LECT1" ,"ALX3" ,"BMP7" ,"PRDM13" ,
           "ULBP1","RAET1K")

#hg38.annotations = as.data.frame(hg38.annotations)

hg38.annotations = hg38.annotations[!is.na(hg38.annotations$symbol),] 

gene_intervals = hg38.annotations[hg38.annotations$symbol %in% genes,]

gene_df = as.data.frame(gene_intervals)

gene_df$intervals = paste(gene_df$seqnames, gene_df$annot.start, gene_df$annot.end, sep = "_")


combat_gr = GRanges(seqnames = COMBAT_meth$chr, ranges = IRanges(start = COMBAT_meth$start, end = COMBAT_meth$end), mcols = combat_meth)

overlaps <- findOverlaps(gene_intervals, combat_gr)
overlap_results <- data.frame(
  query_interval = gene_intervals[queryHits(overlaps)],     # Intervals from table
  subject_interval = combat_gr[subjectHits(overlaps)], # Intervals from matrix
  mcols(combat_gr[subjectHits(overlaps)])              # Overlapping values
)



pdf("test.pdf")
ggplot(overlap_results, aes(x = query_interval.annot.id, y = subject_interval.mcols.COMBAT_0001)) +
   geom_boxplot() +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")
dev.off()
