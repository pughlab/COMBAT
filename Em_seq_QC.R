# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# File paths (adjust these if necessary)
setwd("/Users/arnavaz/Desktop")
markdup_file <- "Sample1_markdup.txt"
wgs_file <- "Sample1_WGS_metrics.txt"

# Read the files, skipping the comment lines
markdup_data <- read.delim(markdup_file, comment.char = "#", header = TRUE, sep = "\t")
wgs_data <- read.delim(wgs_file, comment.char = "#", header = TRUE, sep = "\t")

# Extract metrics sections and histograms
# MARKDUP metrics
markdup_metrics <- markdup_data[1, ]
tmp <- markdup_data[-c(1, 2),]
colnames(tmp) <- NULL
# First row contains the metrics
markdup_histogram <- tmp[, 1:2]# Remaining rows are the histogram
colnames(markdup_histogram) <- c("BIN", "VALUE")

# WGS metrics
wgs_metrics <- wgs_data[1, ]  # First row contains the metrics
wgs_histogram <- wgs_data[-1, ]  # Remaining rows are the histogram

# Convert histogram data to numeric
markdup_histogram <- markdup_histogram %>%
  mutate(BIN = as.numeric(BIN), VALUE = as.numeric(VALUE))


tmp <- wgs_data[-c(1, 2),]
colnames(tmp) <-  NULL
tmp <- tmp[, 1:2]
colnames(tmp) <- c("coverage","high_quality_coverage_count")
wgs_histogram <- tmp %>%
  mutate(coverage = as.numeric(coverage),
         high_quality_coverage_count = as.numeric(high_quality_coverage_count))

# 1. Bar Plot: Percent Duplication from MarkDuplicates
ggplot(markdup_metrics, aes(x = LIBRARY, y = PERCENT_DUPLICATION)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_minimal() +
  labs(
    title = "Percent Duplication by Library",
    x = "Library",
    y = "Percent Duplication"
  )

# 2. Line Plot: Histogram from MarkDuplicates
ggplot(markdup_histogram, aes(x = BIN, y = VALUE)) +
  geom_line(color = "darkgreen") +
  theme_minimal() +
  labs(
    title = "MarkDuplicates Histogram",
    x = "BIN",
    y = "VALUE"
  )

# 3. Coverage Histogram: WGS Metrics
ggplot(wgs_histogram, aes(x = coverage, y = high_quality_coverage_count)) +
  geom_bar(stat = "identity", fill = "darkred") +
  theme_minimal() +
  labs(
    title = "Coverage Distribution (WGS)",
    x = "Coverage",
    y = "High Quality Coverage Count"
  )

# 4. Metrics Summary Table: WGS Metrics
wgs_summary <- wgs_metrics %>%
  select(GENOME_TERRITORY, MEAN_COVERAGE, SD_COVERAGE, MEDIAN_COVERAGE,
         PCT_EXC_MAPQ, PCT_EXC_DUPE, PCT_EXC_UNPAIRED, PCT_EXC_BASEQ, PCT_EXC_TOTAL,
         PCT_1X, PCT_5X, PCT_10X, PCT_20X, PCT_30X) %>%
  gather(key = "Metric", value = "Value")

# 5. Bar Plot: WGS Metrics Summary
ggplot(wgs_summary, aes(x = reorder(Metric, Value), y = Value, fill = Metric)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "WGS Metrics Summary",
    x = "Metrics",
    y = "Value"
  )

# 6. Combined Histogram and Line Plot for Coverage vs. Duplication
combined_histogram <- wgs_histogram %>%
  inner_join(markdup_histogram, by = c("coverage" = "BIN")) %>%
  rename(Duplication_Value = VALUE, Coverage_Count = high_quality_coverage_count)



###### Visualizaing Flagstat output from Samtools #####


# Load necessary libraries
library(ggplot2)

# Function to read the flagstat output and extract relevant values
extract_flagstat_values <- function(file_path) {
  # Read the file
  flagstat_output <- readLines(file_path)
  
  # Initialize a named vector for the rates
  rates <- c("primary_mapped_rate" = NA, "secondary_mapped_rate" = NA, "duplicate_rate" = NA)
  
  # Extract the values from the output
  for (line in flagstat_output) {
    if (grepl("primary", line) && grepl("mapped", line)) {
      # Extract primary mapped values
      primary_mapped <- as.numeric(strsplit(strsplit(line, "\\+")[1], " ")[[1]][1])
      total_reads <- as.numeric(strsplit(flagstat_output[1], "\\+")[1])
      rates["primary_mapped_rate"] <- primary_mapped / total_reads
    }
    
    if (grepl("secondary", line) && grepl("alignments", line)) {
      # Extract secondary mapped values
      secondary_mapped <- as.numeric(strsplit(strsplit(line, "\\+")[1], " ")[[1]][1])
      rates["secondary_mapped_rate"] <- secondary_mapped / total_reads
    }
    
    if (grepl("duplicates", line)) {
      # Extract duplicates
      duplicates <- as.numeric(strsplit(strsplit(line, "\\+")[1], " ")[[1]][1])
      rates["duplicate_rate"] <- duplicates / total_reads
    }
  }
  
  return(rates)
}

# File path to your samtools flagstat output
file_path <- "flagstat_output.txt"

# Extract rates from the file
rates <- extract_flagstat_values(file_path)

# Prepare data for plotting
plot_data <- data.frame(
  Rate = names(rates),
  Value = rates
)

# Plot the rates as a bar plot
ggplot(plot_data, aes(x = Rate, y = Value, fill = Rate)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Mapping Rates", y = "Rate", x = "Metric") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text(aes(label = scales::percent(Value)), vjust = -0.5)

ggplot(combined_histogram, aes(x = coverage)) +
  geom_bar(aes(y = Coverage_Count), stat = "identity", fill = "blue", alpha = 0.5) +
  scale_y_continuous(
    name = "Coverage Count",
    sec.axis = sec_axis(~ . / 1000, name = "Duplication Value")
  ) +
  theme_minimal() +
  labs(
    title = "Coverage Count vs. Duplication Value",
    x = "Coverage"
  )

# Save plots
ggsave("Percent_Duplication.png")
ggsave("MarkDuplicates_Histogram.png")
ggsave("Coverage_Distribution.png")
ggsave("WGS_Metrics_Summary.png")
ggsave("Combined_Coverage_Duplication.png")

# 7. Save Merged Metrics Table (Optional)
merged_metrics <- cbind(markdup_metrics, wgs_metrics)
write.table(merged_metrics, "Merged_Metrics.txt", sep = "\t", row.names = FALSE, quote = FALSE)

