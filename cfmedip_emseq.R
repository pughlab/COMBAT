



### Test ###

library(dplyr)
library(ggplot2)
library(readr)
library(BoutrosLab.plotting.general);
library(annotatr);
library(GenomicRanges);
load("/Users/arnavaz/Library/CloudStorage/OneDrive-UHN/Projects/COMBAT/hg38_gene_annotations.RData")

# load("/Users/arnavaz/Library/CloudStorage/OneDrive-UHN/Projects/COMBAT/cfMeDIP_src_Kui/black_bin_v2.RData")
# Define directories
directory1 = "/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/MethylDackel/cfmedip"

directory2 = "/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/MethylDackel/em-seq"

load("/Users/arnavaz/Library/CloudStorage/OneDrive-UHN/Projects/COMBAT/hg38_cpg_annotations.RData")

# Get list of files in each directory
files1 <- list.files(directory1, pattern = "count.txt", full.names = TRUE)
files2 <- list.files(directory2, pattern = "summed_intervals", full.names = TRUE)

# Ensure both directories have the same number of files
if (length(files1) != length(files2)) {
  stop("The number of files in the directories does not match!")
}

# Initialize an empty list for plots
island <- subset(hg38_cpg_annotations, type == "hg38_cpg_islands")
shores = subset(hg38_cpg_annotations, type == "hg38_cpg_shores")
shelf = subset(hg38_cpg_annotations, type == "hg38_cpg_shelves")

# Loop through the files in both directories
for (i in seq_along(files1)) {
  # Read the tab-separated files
  data1 <- read_tsv(files1[i], col_names = F)
  data2 <- read_tsv(files2[i], col_names = FALSE)
  
  data2 =  data2 %>% mutate(interval_id = paste(X1, X2, X3, sep = "_"))
  interval = data2 %>% select(interval_id)
  data1$interval_id = interval$interval_id
  # data2 = data2 %>% filter(!interval_id %in% black_bin$black_bin)
  # data1 = data1 %>% filter(!interval_id %in% black_bin$black_bin)
  combat_gr = GRanges(seqnames = data1$X1, ranges = IRanges(start = data1$X2, data1$X3), mcols = data1$X4)
  combat_emseq = GRanges(seqnames = data2$X1, ranges = IRanges(start = data2$X2, data2$X3), mcols = data2$X4)
  #Island Shores Shelf for Cfmedip
  
  overlaps <- findOverlaps(island, combat_gr)
  overlap_results_island <- data.frame(
    query_interval = island[queryHits(overlaps)],     # Intervals from table
    subject_interval = combat_gr[subjectHits(overlaps)], # Intervals from matrix
    mcols(combat_gr[subjectHits(overlaps)])              # Overlapping values
  )
  
  
  overlaps <- findOverlaps(shores, combat_gr)
  overlap_results_shores <- data.frame(
    query_interval = shores[queryHits(overlaps)],     # Intervals from table
    subject_interval = combat_gr[subjectHits(overlaps)], # Intervals from matrix
    mcols(combat_gr[subjectHits(overlaps)])              # Overlapping values
  )
  overlaps <- findOverlaps(shelf, combat_gr)
  overlap_results_shelf <- data.frame(
    query_interval = shelf[queryHits(overlaps)],     # Intervals from table
    subject_interval = combat_gr[subjectHits(overlaps)], # Intervals from matrix
    mcols(combat_gr[subjectHits(overlaps)])              # Overlapping values
  )
  
  
 
  ##### Island Shores Shelf for Emseq
  overlaps <- findOverlaps(island, combat_emseq)
  overlap_results_island2 <- data.frame(
    query_interval = island[queryHits(overlaps)],     # Intervals from table
    subject_interval = combat_emseq[subjectHits(overlaps)], # Intervals from matrix
    mcols(combat_emseq[subjectHits(overlaps)])              # Overlapping values
  )
  
  
  overlaps <- findOverlaps(shores, combat_emseq)
  overlap_results_shores2 <- data.frame(
    query_interval = shores[queryHits(overlaps)],     # Intervals from table
    subject_interval = combat_emseq[subjectHits(overlaps)], # Intervals from matrix
    mcols(combat_emseq[subjectHits(overlaps)])              # Overlapping values
  )
  overlaps <- findOverlaps(shelf, combat_emseq)
  overlap_results_shelf2 <- data.frame(
    query_interval = shelf[queryHits(overlaps)],     # Intervals from table
    subject_interval = combat_emseq[subjectHits(overlaps)], # Intervals from matrix
    mcols(combat_emseq[subjectHits(overlaps)])              # Overlapping values
  )
  
  
  
  values1 <- overlap_results_island$mcols
  values2 <- overlap_results_island2$mcols
  
  
  values1_shore <- overlap_results_shores$mcols
  values2_shore <- overlap_results_shores2$mcols
  
  
  values1_shelf <- overlap_results_shelf$mcols
  values2_shelf <- overlap_results_shelf2$mcols
  
  # Create a data frame for plotting
  comparison_df <- data.frame(
    Value1 = values1/sum(values1),
    Value2 = values2/sum(values2)
  )

  
  comparison_df_shore <- data.frame(
    Value1 = values1_shore/sum(values1_shore),
    Value2 = values2_shore/sum(values2_shore)
  )
  
  
  comparison_df_shelf <- data.frame(
    Value1 = values1_shelf/sum(values1_shelf),
    Value2 = values2_shelf/sum(values2_shelf)
  )
  
  
fname = gsub("_count.txt","", basename(files1[i]))


comparison_df <- comparison_df

corval = cor(comparison_df$Value1, comparison_df$Value2, method = "pearson")
pdf(paste0(fname, "_Island.pdf"))

  # Generate a scatter plot
  ggplot(comparison_df, aes(x = log10(Value1), y = log10(Value2))) +
    geom_point(alpha = 0.7, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = paste0("Cfmedip and Em-seq methylation correlation on CpG Islands:", fname ),
      x = "cfMedip methylation normalized counts (log10)" ,
      y = "Em-seq methylation normalized counts (log10)"
    ) +
     # xlim(-9,-5)+
     # ylim(-9,-5)+
    annotate(
      "text", 
      x = log10(1e-06) , y = max(log10(comparison_df$Value2)),  # Position (top-left)
      label = paste0("Correlation: ", round(corval, 2)), 
      hjust = 0, size = 5, color = "red"
    )+
    theme_minimal()
 
dev.off()



}

file1 = "/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/MethylDackel/cfmedip/COMBAT_0006_CpG_Island_cfmedip.tsv"
file2 = "/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/MethylDackel/em-seq/COMBAT_0006_CpG_Island_emseq.tsv"

data1 <- read_tsv(file1, col_names = F)
data2 <- read_tsv(file2, col_names = F)

values1 <- data1[,4]
values2 <- data2[,4]

comparison_df <- data.frame(
  Value1 = values1/sum(values1),
  Value2 = values2/sum(values2)
)


colnames(comparison_df) <- c("Value1", "Value2")

ggplot(comparison_df, aes(x = log10(Value1), y = log10(Value2))) +
  geom_point(alpha = 0.7, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = paste0("Comparison cfmedip and Em-seq methylation at CpG Islands:  ", "COMBAT 0006" ),
    x = "cfMedip methylation normalized counts (log10)" ,
    y = "Em-seq methylation normalized counts (log10)"
  ) +
  
  theme_minimal()
