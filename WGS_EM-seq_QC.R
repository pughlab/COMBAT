
### QC for EM-seq## Basic plots 


bd.g.path =  "/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/MethylDackel/em-seq"
wgs.met.path = "/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/WgsMetrics"
bd.g.flist = list.files(bd.g.path, full.names = T, pattern ="CpG.bedGraph")
# wgs.metrics = list.files(wgs.met.path, full.names = T, pattern ="alignment_metrics.txt")

for (file in bd.g.flist){

  tab = read.table(file,  sep = "\t", stringsAsFactors = F, header = T, skip=1)

}




# Load necessary libraries
library(ggplot2)

# Read the GC Bias Detail Metrics output from Picard
# Replace "GcBiasDetailMetrics.txt" with your actual file name
f = "/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/WgsMetrics/COMBAT_0006_09_T1_bwamem_aligned_sorted_merged_markdup.bam_gc_bias_metrics.txt"
gc_bias_data <- read.table(f, comment.char = "#", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Explore the data (Optional)
head(gc_bias_data)

# Key columns typically include:
# - "WINDOWS": Number of windows at each GC content
# - "READS_USED": Number of reads contributing to GC content
# - "NORMALIZED_COVERAGE": Coverage normalized to average

# Create a plot of GC content vs. normalized coverage
gc_bias_plot <- ggplot(gc_bias_data, aes(x = GC, y = NORMALIZED_COVERAGE)) +
  geom_line(color = "blue") +                          # Line plot for GC bias trend
  geom_point(size = 1, color = "red", alpha = 0.6) +   # Points for individual GC windows
  theme_minimal() +                                    # Clean plot theme
  labs(
    title = "GC Bias Visualization",
    x = "GC Content (%)",
    y = "Normalized Coverage",
    caption = "Source: Picard GcBiasDetailMetrics"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10))         # Adjust X-axis ticks

# Display the plot
print(gc_bias_plot)

# Save the plot as an image (Optional)
ggsave("GC_Bias_Plot.png", gc_bias_plot, width = 8, height = 6)







# Load required libraries
library(ggplot2)
library(dplyr)

# List of file paths for CollectAlignmentSummaryMetrics output
file_paths <- list.files(path = "/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/WgsMetrics", pattern = ".*alignment_metrics.txt", full.names = TRUE)

# Function to read and extract the relevant columns
read_picard_metrics <- function(file_path) {
  data <- read.delim(file_path, skip = 6, header = TRUE)  # Skip the header rows
  sample_name <- gsub(".*\\/|_alignment_metrics.txt", "", file_path)  # Extract the sample name
  data$Sample <- sample_name
  return(data)
}

# Apply this function to all the files and bind them together
all_samples_data <- do.call(rbind, lapply(file_paths, read_picard_metrics))




