#####################################
###     mesothelioma SV ratio     ###
#####################################

#####################################
###      libraries & options      ###
#####################################
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)
library(gridExtra)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(GenomicRanges)
library(cowplot)
library(gtable)
library(ggbeeswarm)
library(ggpubr)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

#####################################
###           file paths          ###
#####################################
sample_table_path  <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
muttable_path      <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/patient_mutburden.RDS"
SV_table_path      <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/SV/input_files/input_32matrix_patient.rds"

output_path        <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/GenomicOverview/"

#####################################
###              MAIN             ###
#####################################
sample_table <- read.table(sample_table_path, head = TRUE, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

SV_table     <- readRDS(SV_table_path)
SV_table <- data.frame(SV_table)
SV_table <- data.frame(rowSums(SV_table))
SV_table$sample <- rownames(SV_table)
SV_table$sample <- gsub("X", "", SV_table$sample)
colnames(SV_table) <- c("SVcount", "sample")

# check whether mesotheliomas have more SVs than the other subtypes
SV_table$sample <- as.character(SV_table$sample)
sample_table$participant_id <- as.character(sample_table$participant_id)
SV_table <- left_join(SV_table, sample_table[, c("participant_id", "histology")], by = c("sample" = "participant_id"))

SV_table <- SV_table[complete.cases(SV_table), ]

SV_table[SV_table$histology == "MESOTHELIOMA", "group"] <- "Mesothelioma"
SV_table[is.na(SV_table$group), "group"] <- "all_other_subtypes"

p <- ggboxplot(SV_table, x = "group", y = "SVcount",
               color = "group", outlier.shape = NA) +
  geom_quasirandom(data = SV_table, aes(x = group, y = SVcount, color = group), alpha = 0.3) +
  scale_color_manual(values = c("Mesothelioma" = "#e3309e",
                                "all_other_subtypes" = "black"),
                     breaks = c("Mesothelioma", "all_other_subtypes"),
                     labels = c("Mesothelioma" = "mesothelioma",
                                "all_other_subtypes" = "all other subtypes")) +
  scale_y_log10() +
  scale_x_discrete(labels = c("all other subtypes", "Mesothelioma")) +
  xlab("") +
  ylab("# SVs") +
  stat_compare_means() +
  theme_bw() +
  theme(legend.position = "none")


pdf(paste0(output_path, "SV_Meso_hist_comparison.pdf"), width = 4, height = 4)
p
dev.off()

# do mesotheliomas have a higher SV/SBS ratio
muttable <- readRDS(muttable_path)

# get all the SBSs
patient_hist_df <- muttable
patient_hist_df$SV_SBS_ratio <- patient_hist_df$SVcount / patient_hist_df$SBScount

patient_hist_df[patient_hist_df$histology == "MESOTHELIOMA", "group"] <- "Mesothelioma"
patient_hist_df[is.na(patient_hist_df$group), "group"] <- "all_other_subtypes"

p_ratio <- ggboxplot(patient_hist_df, x = "group", y = "SV_SBS_ratio", color = "group", outlier.shape = NA, palette = c("black", "#e3309e") ) +
                geom_quasirandom(data = patient_hist_df, aes(x = group, y = SV_SBS_ratio, color = group), alpha = 0.3) +
                stat_compare_means() +
                scale_x_discrete(labels = c("all other subtypes", "Mesothelioma")) +
                xlab("") +
                ylab("SV:SBS ratio") +
                scale_y_log10() +
                theme_bw() +
                theme(legend.position = "none")

pdf(paste0(output_path, "SV_SBS_ratio_Meso_hist_comparison.pdf"), width = 4, height = 4)
p_ratio
dev.off()
