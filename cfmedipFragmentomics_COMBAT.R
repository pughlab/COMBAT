# Load necessary library
library(dplyr)
library(ggplot2)
library(readr)
library(ggplot2)
library(purrr)
library(gridExtra)
COMBAT = read_tsv ("/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/MEDIPIPE/fragment_profile_GC_corrected_1mb.tsv")
HBC = read_tsv ("/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/MEDIPIPE/aggregated_HBC/fragment_profile_GC_corrected_1mb.tsv")


armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

# Read the data

data <- HBC
mat = data[,6:33]
mat = rowMeans(mat)

data = cbind(data[1:5], mat)
# Create a combined column for chromosome arm
data$chrom_arm_bin <- paste(data$chr , data$arm,data$bin_id, sep = "_")
data$chrom_arm  <- paste(data$chr, data$arm, sep = "_")



# Plot short.corrected2 across chromosome arms
unique_labels <- data$chrom_arm
unique_labels[duplicated(data$chrom_arm)] <- "" 
data2 <- COMBAT

data2$chrom_arm_bin <- paste(data2$chr , data2$arm,data2$bin_id, sep = "_")
data2$chrom_arm  <- paste(data2$chr, data2$arm, sep = "_")
# Plot short.corrected2 across chromosome arms
unique_labels <- data2$chrom_arm
unique_labels[duplicated(data2$chrom_arm)] <- ""

data = data %>% filter(chr == "chr17")
data2 = data2 %>% filter(chr == "chr17")

unique_labels <- data2$chrom_arm
unique_labels[duplicated(data2$chrom_arm)] <- ""
unique_labels <- data$chrom_arm
unique_labels[duplicated(data$chrom_arm)] <- "" 




# Read the TSV file
data <- read.delim("COMBAT_fragRation.tsv", header = TRUE, sep = "\t")

# Function to create a single plot
create_plot <- function(column, d) {
  ggplot() +
    
    geom_line(data = data, aes(x = chrom_arm_bin, y = mat ), color = "blue", linetype = "solid", group = 1) +
    geom_line(data = d, aes(x = chrom_arm_bin, y = !!sym(column)), color = "red", linetype = "dashed", group = 1)+
    scale_x_discrete(labels = unique_labels) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = column, x = "Chromosome Arm Bin", y = "Value")
}

# Create a list of plots
plots <- map(names(data2)[6:24], ~create_plot(.x, data2))

# Save plots to PDF
pdf("chromosome_plots.pdf", width = 40, height = 8)
walk(plots, print)
dev.off()

