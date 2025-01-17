rm(list = ls())
setwd("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/3_fragmentomics/1_fragment_ratio")

library(dplyr)
library(ggplot2)
library(ggpubr)

##############
# ggplot theme
{
  mytheme <- theme_classic(base_size=12) + theme(
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    strip.text.x = element_text(size=11),
    strip.text.y = element_text(size=12),
    axis.title.x = element_text(face="bold", size=17),
    axis.title.y = element_text(size=15),
    axis.text.y = element_text(size=15),
    plot.title = element_text(size=15),
    legend.position = "none",
    legend.title = element_text(size=10),
    legend.text = element_text(size=10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_rect(fill="white", color="white"),
    panel.spacing.x=unit(0.1, "lines"))
}

## healthy mean 
fr <- read.table("TCGE-CFMe-HBC/fragment_profile_GC_corrected_5mb.tsv",
                 comment.char = "!", header = T, check.names = F)
fr_m = data.matrix(fr[, -(1:5)]) 
#fr_m_z = t(scale(t(fr_m), center = TRUE, scale = TRUE))   ## z-score transform
ratio_mean_hbc = rowMeans(fr_m, na.rm = T)

fr_combo <- fr[, 1:5]      ## aggregated ratio matix without scaling for all projects


###############
## per project 
##############
project_id <- c("TCGE-CFMe-AML", "TCGE-CFMe-CRC", "TCGE-CFMe-LFS", "TCGE-CFMe-HBC",
            "TCGE-CFMe-LTX", "TCGE-CFMe-PDAC", "TCGE-CFMe-BCA", "TCGE-CFMe-HNSC",
            "TCGE-CFMe-PRAD", "TCGE-CFMe-BRCA", "TCGE-CFMe-HPCC", "TCGE-CFMe-SCLC")
s_info <- read.csv("/Users/yong/OneDrive - UHN/Projects/TCGE/cfEpigenomics/Resource/0_sample_meta_info/TCGE-CFMe-Only-Samples-light.csv")

fr_aggr <- data.frame()   ## aggregated data.frame after scaling per project


for (cohort in project_id)
{
print(cohort)
fr <- read.table("/Users/arnavaz/Desktop/HPC/pughlab/projects/COMBAT/MEDIPIPE/aggregated/fragment_profile_GC_corrected_5mb.tsv",
                 comment.char = "!", header = T, check.names = F)

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

fr$arm <- factor(fr$arm, levels=armlevels)

## 
arm <- fr %>% 
       group_by(arm) %>%
       summarize(n=n()) %>%
       mutate(arm = as.character(arm))

## adjust labels for visualization 
## show long arms (q) for some of chromosomes
small.arms <- setNames(c("", "12q", "", "16q",
                         "", "17q", "", "18q",
                         "", "19", "", "20",
                         "21", "22"),
                       c("12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

###########################################
fr_m = data.matrix(fr[, -(1:5)])     ## ratio martix

## keep cfMe only 
#idx_s <- match(colnames(fr_m), s_info$sequencing_id)
fr_m <- na.omit(fr_m)
fr_m_z = t(scale(t(fr_m), center = TRUE, scale = TRUE))   ## z-score transform

## merge fr matirx 



N = nrow(fr_m)
M <- ncol(fr_m)

bin = rep(fr$bin_id, M)
arm = rep(fr$arm, M)

project = rep("COMBAT", N*M)
group = rep(group_s, each = N)
sample = rep(colnames(fr_m), each = N)
ratio = as.vector(fr_m)
ratio_mean_project = rep(rowMeans(fr_m, na.rm = T), M)
ratio_mean_healthy = rep(ratio_mean_hbc, M)

## z scores
ratio_z = as.vector(fr_m_z)

## sum z scores for whole genome
ratio_z_genome = rep(colSums(fr_m_z, na.rm = TRUE), each = N)

fr_dat = data.frame(bin, arm, sample, group, project, 
                    ratio, ratio_mean_project, ratio_mean_healthy,
                    ratio_z, ratio_z_genome)

fr_aggr = rbind(fr_aggr, fr_dat)

########################
## using original ratios
{
# Generate Fragmentation and Coverage plots
g <- ggplot(fr_dat, aes(x = bin, y = ratio, group = sample))
g <- g + geom_line(linewidth = 0.5, alpha = 0.5, color = "gray40")
## add mean line per project
g <- g + geom_line(aes(x = bin, y = ratio_mean_project),
                   size=0.75, alpha=0.5, color="black")
## add mean line baed the HBC samples
g <- g + geom_line(aes(x = bin, y = ratio_mean_healthy),
                   size=0.75, alpha=0.5, color="green")
g <- g + facet_grid(~arm, switch = "x",space = "free_x", scales = "free_x", 
                    labeller = labeller(arm = arm.labels))
g <- g + labs(x = cohort, y = "Short/Long Fragment Ratio\n", color="")
g <- g + mytheme
ggsave(paste0(cohort, "_fragment_ratio_in_5Mb.png"), width = 16, height = 4)

## subgroups per each individual projects
g <- ggplot(fr_dat, aes(x = bin, y = ratio, group = sample))
g <- g + geom_line(linewidth = 0.5, alpha = 0.5, color = "gray40")
## add mean line per project
g <- g + geom_line(aes(x = bin, y = ratio_mean_project),
                   size=0.75, alpha=0.5, color="black")
## add mean line baed the HBC samples
g <- g + geom_line(aes(x = bin, y = ratio_mean_healthy),
                   size=0.75, alpha=0.5, color="green")
g <- g + facet_grid(group ~ arm, switch = "x",space = "free_x", scales = "free_x", 
                    labeller = labeller(arm = arm.labels))
g <- g + labs(x = cohort, y = "Short/Long Fragment Ratio\n", color="")
g <- g + mytheme
ggsave(paste0(cohort, "_fragment_ratio_in_5Mb_per_group.png"), width = 15, height = 8)
}

#########################################
## using z transformed ratios per project
{
  # Generate Fragmentation and Coverage plots
  g <- ggplot(fr_dat, aes(x = bin, y = ratio_z, group = sample))
  g <- g + geom_line(linewidth = 0.5, alpha = 0.5, color = "gray40")
  g <- g + geom_hline(yintercept = 0, linetype="dashed", color = "red", linewidth = 0.5)
  g <- g + facet_grid(~arm, switch = "x",space = "free_x", scales = "free_x", 
                      labeller = labeller(arm = arm.labels))
  g <- g + labs(x = cohort, y = "Short/Long Fragment Ratio Z score\n", color="")
  g <- g + mytheme
  ggsave(paste0(cohort, "_fragment_ratio_Z_score_in_5Mb.png"), width = 16, height = 4)
  
  ## subgroups per each individual projects
  g <- ggplot(fr_dat, aes(x = bin, y = ratio_z, group = sample))
  g <- g + geom_line(linewidth = 0.5, alpha = 0.5, color = "gray40")
  g <- g + geom_hline(yintercept = 0, linetype="dashed", color = "red", linewidth = 0.5)
  g <- g + facet_grid(group ~ arm, switch = "x",space = "free_x", scales = "free_x", 
                      labeller = labeller(arm = arm.labels))
  g <- g + labs(x = cohort, y = "Short/Long Fragment Ratio Z score\n", color="")
  g <- g + mytheme
  ggsave(paste0(cohort, "_fragment_ratio_Z_score_in_5Mb_per_group.png"), width = 15, height = 8)
}


}

saveRDS(fr_aggr, "Aggreated_fr_after_scaling_per_project.RDS")
saveRDS(fr_combo, "Aggreated_fr_without_scaling.RDS")


#########################################
## projects aggregated and z transformed
########################################
{
  
  ###########################################
  ## aggregated projects with original ratios
  {
    # Generate Fragmentation and Coverage plots
    g <- ggplot(fr_aggr, aes(x = bin, y = ratio, group = sample))
    g <- g + geom_line(linewidth = 0.5, alpha = 0.5, color = "gray40")
    ## add mean line per project
    g <- g + geom_line(aes(x = bin, y = ratio_mean_project),
                       size=0.75, alpha=0.5, color="black")
    ## add mean line baed the HBC samples
    g <- g + geom_line(aes(x = bin, y = ratio_mean_healthy),
                       size=0.75, alpha=0.5, color="green")
    g <- g + facet_grid(project ~ arm, switch = "x",space = "free_x", scales = "free_x", 
                        labeller = labeller(arm = arm.labels))
    g <- g + labs(x = "All datasets", y = "Short/Long Fragment Ratio\n", color="")
    g <- g + mytheme
    ggsave("Aggregated_fragment_ratio_in_5Mb.png", width = 16, height = 20)
    
    ## healthy samples only 
    fr_aggr_h <- filter(fr_aggr, group == "BRCA_control" | 
                                 group == "Normal" |
                                 group == "SCLC_Normal" |
                                 group == "HBC")
    
    # Generate Fragmentation and Coverage plots
    g <- ggplot(fr_aggr_h, aes(x = bin, y = ratio, group = sample))
    g <- g + geom_line(linewidth = 0.5, alpha = 0.5, color = "gray40")
    ## add mean line per project
    g <- g + geom_line(aes(x = bin, y = ratio_mean_project),
                       size=0.75, alpha=0.5, color="black")
    ## add mean line baed the HBC samples
    g <- g + geom_line(aes(x = bin, y = ratio_mean_healthy),
                       size=0.75, alpha=0.5, color="green")
    g <- g + facet_grid(project ~ arm, switch = "x",space = "free_x", scales = "free_x", 
                        labeller = labeller(arm = arm.labels))
    g <- g + labs(x = "All Healthy samples", y = "Short/Long Fragment Ratio\n", color="")
    g <- g + mytheme
    ggsave("Aggregated_fragment_ratio_in_5Mb_Healthy.png", width = 16, height = 10)
    
  }
  
  #################################################
  ## aggregated projects ratio_z_genome per project
  {
    fr_aggr_p <- fr_aggr %>% 
                 select(sample, project, group, ratio_z_genome) %>% 
                 distinct(sample, .keep_all = TRUE)
    
    ## per cancer type
    g <- ggboxplot(fr_aggr_p, x = "group", y = "ratio_z_genome", fill = "project")
    g <- g + geom_hline(yintercept = 0, linetype="dashed", color = "red", linewidth = 0.5)
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0("Aggreated_fragment_ratio_Z_Scores_per_project_group.png"), width = 16, height = 6)
  }
  
  ############################################################
  ## using the combo fr matrix and do the Z transform globally
  {
    ## merge fr matirx 
    fr_combo_m <- data.matrix(fr_combo[, -(1:5)]) 
    fr_combo_m_z = t(scale(t(fr_combo_m), center = TRUE, scale = TRUE))   ## z-score transform
    ratio_z_score = colSums(fr_combo_m_z, na.rm = TRUE)
    
    idx_ss <- match(colnames(fr_combo_m_z), s_info$sequencing_id)
    s_info_ss <- s_info[idx_ss, ]
    
    fr_combo_z_socre <- cbind(s_info_ss, ratio_z_score)
    
    ## genome-wide score boxplot 
    g <- ggboxplot(fr_combo_z_socre, x = "cancer_type", y = "ratio_z_score", fill = "project_id")
    g <- g + geom_hline(yintercept = 0, linetype="dashed", color = "red", linewidth = 0.5)
    g <- g + labs(y = "Sum of Fragment Ratio Z socres")
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0("Aggreated_fragment_ratio_Z_Scores_all_projects_cancer_type.png"), width = 16, height = 6)
    
    g <- ggboxplot(fr_combo_z_socre, x = "group", y = "ratio_z_score", fill = "project_id")
    g <- g + geom_hline(yintercept = 0, linetype="dashed", color = "red", linewidth = 0.5)
    g <- g + labs(y = "Sum of Fragment Ratio Z socres")
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0("Aggreated_fragment_ratio_Z_Scores_all_projects_group.png"), width = 16, height = 6)
    
    ## sex vs cancer_type
    g <- ggboxplot(fr_combo_z_socre, x = "sex", y = "ratio_z_score", fill = "cancer_type")
    g <- g + geom_hline(yintercept = 0, linetype="dashed", color = "red", linewidth = 0.5)
    g <- g + facet_wrap(~ project_id, nrow = 2)
    g <- g + labs(y = "Sum of Fragment Ratio Z socres")
    g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0("Aggreated_fragment_ratio_Z_Scores_all_projects_sex_cancer_type.png"), width = 16, height = 10)
    
    ## association with ages
    g <- ggscatter(fr_combo_z_socre, y = "ratio_z_score", x = "age", color = "project_id",
                   add = "reg.line") 
    g <- g + facet_wrap(. ~ project_id, nrow =2) + theme(legend.position = "none")
    g <- g + stat_cor(label.y = 4000)
    g <- g + stat_regline_equation(label.y = 3600)
    g <- g + labs(y = "Sum of Fragment Ratio Z socres")
    ggsave(paste0("Aggreated_fragment_ratio_Z_Scores_all_projects_vs_age.png"), width = 16, height = 8)
    
  }
  
  
  
}
