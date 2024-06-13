# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #           analyse tumours without drivers assigned                    # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#################################################################################
####################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggpubr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(survival)
library(survminer)
library(ggbeeswarm)
library(cowplot)
library(stringr)

options(stringsAsFactors = F)
options(bitmapType = "cairo")


#################################################################################
#######################################################################file paths

driver_count_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/patient_driver_binary.txt"
cluster_assignment_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/new_signature_cluster_patient_assignment.txt"
sample_table_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
mut_burden_path         <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/patient_mutburden.RDS"
smoking_path            <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/smoking_data.txt"

output_path             <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/export/"

#################################################################################
#######################################################################     MAIN

# get the tumours driver counts 

driver_count <- read.table(driver_count_path, head = T, sep = "\t")
colnames(driver_count)[1] <- "patient"
driver_count$patient <- as.character(driver_count$patient)
driver_count <- left_join(driver_count, cluster_assignment, by = c("patient" = "variable"))
driver_count$cluster_assignment <- as.character(driver_count$cluster_assignment)

# add additional data of these patients
sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

data <- left_join(driver_count, sample_table[, c("participant_id", "histology")], by = c("patient" = "participant_id"))

mut_burden <- readRDS(mut_burden_path)

data <- left_join(data, mut_burden,  by = c("patient" = "sample"))

# make a classifier of having a driver or not
data[data$coding == FALSE & data$non_coding == FALSE & data$SVdriver == FALSE & data$CNdriver == FALSE, "driver_present"] <- FALSE
data[is.na(data$driver_present), "driver_present"] <- TRUE

driver_mutburden <- data[data$histology.x %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), c("patient", "histology", "driver_present", "SBScount", "IDcount", "DBScount", "SVcount"), ]
driver_mutburden <- melt(driver_mutburden, id.vars = c("patient", "histology", "driver_present"))

p_mutburden <- ggplot(driver_mutburden, aes(driver_present, value + 1, fill = driver_present)) +
  geom_boxplot() +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), label.y = log10(50000), size = 3) +
  scale_y_log10() +
  scale_fill_manual(values = c("#542788", "#f1a340")) +
  facet_grid(histology.x~variable) +
  xlab("cancer gene aberration present") +
  ylab("mutation burden") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill="white"))

pdf(paste0(output_path, "tumours_no_driver_mutburden.pdf"), width = 5, height = 3.5)
p_mutburden
dev.off()

# how about smoking
smoking <- read.table(smoking_path, head = T, sep = "\t")
smoking$participant_id <- as.character(smoking$participant_id)
data <- left_join(data, smoking[, c("participant_id", "classifier")], by = c("patient" = "participant_id"))

df_smoking <- data.frame(table(data[data$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), c("driver_present", "classifier", "histology")]))
df_smoking[df_smoking$driver_present == TRUE & df_smoking$histology == "ADENOCARCINOMA", "driver_count"] <- driver_prop[driver_prop$driver_present == TRUE & driver_prop$histology == "ADENOCARCINOMA", "Freq"]
df_smoking[df_smoking$driver_present == TRUE & df_smoking$histology == "SQUAMOUS_CELL", "driver_count"] <- driver_prop[driver_prop$driver_present == TRUE & driver_prop$histology == "SQUAMOUS_CELL", "Freq"]
df_smoking[df_smoking$driver_present == FALSE & df_smoking$histology == "ADENOCARCINOMA", "driver_count"] <- driver_prop[driver_prop$driver_present == FALSE & driver_prop$histology == "ADENOCARCINOMA", "Freq"]
df_smoking[df_smoking$driver_present == FALSE & df_smoking$histology == "SQUAMOUS_CELL", "driver_count"] <- driver_prop[driver_prop$driver_present == FALSE & driver_prop$histology == "SQUAMOUS_CELL", "Freq"]

df_smoking$prop <- df_smoking$Freq / df_smoking$driver_count
df_smoking$prop <- df_smoking$prop*100

p_smoking <- ggplot(df_smoking, aes(driver_present, prop, fill = classifier)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#7ad8fa","#46788a")) +
  facet_wrap(~histology.x, nrow = 1) +
  xlab("cancer gene aberration present") +
  ylab("% tumours") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 7))

pdf(paste0(output_path, "tumours_no_driver_smoking.pdf"), width = 4, height = 3)
p_smoking
dev.off()

# test for significance

# lusc
dat <- df_smoking[df_smoking$histology.x == "SQUAMOUS_CELL", ]
dat <- dat[, c("driver_present", "classifier", "Freq")]
dat <- dcast(driver_present ~ classifier, data = dat)
rownames(dat) <- dat$driver_present
dat$driver_present <- NULL

fisher.test(dat)

# luad
dat <- df_smoking[df_smoking$histology.x == "ADENOCARCINOMA", ]
dat <- dat[, c("driver_present", "classifier", "Freq")]
dat <- dcast(driver_present ~ classifier, data = dat)
rownames(dat) <- dat$driver_present
dat$driver_present <- NULL

fisher.test(dat)
