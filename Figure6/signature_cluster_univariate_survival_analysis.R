# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #            univariate signature clustering survival analysis        # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
###################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(dplyr)
library(reshape2)
library(survival)
library(survminer)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

################################################################################
######################################################################file paths

sample_table_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
survival_path         <- "/re_gecip/cancer_lung/shared/pancan_survival_release_v16.RDS"
smoking_path          <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/smoking_data.txt"
clin_data_path        <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/sampleTables/v14_lung_patient_clinical_data_2022-07-05_useful.tsv"
mut_burden_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/patient_mutburden.RDS"
CPI_response_path     <- "/re_gecip/cancer_lung/bsimpson/CPI_cohort/Scripts/Analysis_version_lock_v13/Mastersheet_V31.csv"
cluster_path          <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/new_signature_cluster_patient_assignment.txt"
output_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/"


TCRA_path             <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/5.Immune/summary_files/VDJ_GEL_segmodel_summary_adjusted_20220201.RDS"


################################################################################
######################################################################      MAIN

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

cluster_assignment <- read.table(cluster_path, head = T, sep = "\t")
colnames(cluster_assignment) <- c("cluster_assignment", "participant_id")
cluster_assignment$cluster_assignment <- as.character(cluster_assignment$cluster_assignment)
cluster_assignment$participant_id <- as.character(cluster_assignment$participant_id)

survival <- readRDS(survival_path)
survival$participant_id <- as.character(survival$participant_id)
survival$survival <- sub(" days", "", survival$survival)
survival$survival <- as.numeric(survival$survival)
survival <- survival[survival$participant_id %in% sample_table$participant_id, ]
survival <- survival[survival$disease_type == "LUNG", ]

data <- left_join(cluster_assignment, survival[, c("participant_id", "status", "survival")])
data <- left_join(data, sample_table[, c("participant_id", "histology")])


# add stage
clinical <- read.table(clin_data_path, head = T, sep = "\t")
clinical$participant_id <- as.character(clinical$participant_id)
stage <- clinical[, c("participant_id", "stage_stage_best")]
stage$stage <- stage$stage_stage_best
stage[which(stage$stage == ""), "stage"] <- NA
stage[which(stage$stage == "?"), "stage"] <- NA
stage[which(stage$stage == "U"), "stage"] <- NA
stage[which(stage$stage %in% c("1A", "1A_1", "1A3", "1A_3B_", "1A1", "1A2", "1A2", "1B")), "stage"] <- 1
stage[which(stage$stage %in% c("2", "2A", "2B")), "stage"] <- 2
stage[which(stage$stage %in% c("3", "3A", "3B", "3C")), "stage"] <- 3
stage[which(stage$stage %in% c("4", "4A", "4B")), "stage"] <- 4
stage$stage <- as.character(stage$stage)
stage[stage$stage %in% c("1", "2", "3"), "stage_category"] <- "early"
stage[stage$stage %in% c("4"), "stage_category"] <- "late"

data <- left_join(data, stage)

# test the groups individually vs all other groups
data[data$cluster_assignment == 1, "in1"] <- 1
data[data$cluster_assignment != 1, "in1"] <- 0

data[data$cluster_assignment == 2, "in2"] <- 1
data[data$cluster_assignment != 2, "in2"] <- 0

data[data$cluster_assignment == 3, "in3"] <- 1
data[data$cluster_assignment != 3, "in3"] <- 0

data[data$cluster_assignment == 4, "in4"] <- 1
data[data$cluster_assignment != 4, "in4"] <- 0

data[data$cluster_assignment == 5, "in5"] <- 1
data[data$cluster_assignment != 5, "in5"] <- 0

data[data$cluster_assignment == 6, "in6"] <- 1
data[data$cluster_assignment != 6, "in6"] <- 0

data[data$cluster_assignment == 7, "in7"] <- 1
data[data$cluster_assignment != 7, "in7"] <- 0

data[data$cluster_assignment == 8, "in8"] <- 1
data[data$cluster_assignment != 8, "in8"] <- 0

data[data$cluster_assignment == 9, "in9"] <- 1
data[data$cluster_assignment != 9, "in9"] <- 0

data[data$cluster_assignment == 10, "in10"] <- 1
data[data$cluster_assignment != 10, "in10"] <- 0

data[data$cluster_assignment == 11, "in11"] <- 1
data[data$cluster_assignment != 11, "in11"] <- 0

data[data$cluster_assignment == 12, "in12"] <- 1
data[data$cluster_assignment != 12, "in12"] <- 0

# kick out tumours with a survival that's longer than 10 years
data <- data[data$survival <= (10*365), ]

# thise tumours with a survival of more than 5 years censor them at 5 years
data[which(data$survival >= (5*365)), "status"] <- 1
data[which(data$survival >= (5*365)), "survival"] <- (5*365)


fit1 <- survfit(Surv(survival, status) ~ in1, data = data)
ggsurvplot(fit1, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)


s <- coxph(Surv(survival, status) ~ in2  + stage + histology, data = data[data$survival >= 90, ])
s <- coxph(Surv(survival, status) ~ in2  + stage, data = data[data$survival >= 90, ])
s <- coxph(Surv(survival, status) ~ in2  + stage, data = data[data$survival >= 90 & data$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), ])

fit2 <- survfit(Surv(survival, status) ~ in2, data = data[data$survival >= 90, ])
p <- ggsurvplot(fit2, data = data[data$survival >= 90, ],
           conf.int = TRUE,
           pval = TRUE,
           risk.table = FALSE,
           palette = c("black", "#1f78b4"))

pdf(paste0(output_path, "cluster2_survivial.pdf"), width = 3, height = 3)
p
dev.off()

fit3 <- survfit(Surv(survival, status) ~ in3, data = data)
ggsurvplot(fit3, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

fit4 <- survfit(Surv(survival, status) ~ in4, data = data)
ggsurvplot(fit4, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

fit5 <- survfit(Surv(survival, status) ~ in5, data = data)
ggsurvplot(fit5, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

fit6 <- survfit(Surv(survival, status) ~ in6, data = data)
ggsurvplot(fit6, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

fit7 <- survfit(Surv(survival, status) ~ in7, data = data)
ggsurvplot(fit7, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

fit8 <- survfit(Surv(survival, status) ~ in8, data = data)
ggsurvplot(fit8, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

fit9 <- survfit(Surv(survival, status) ~ in9, data = data)
ggsurvplot(fit9, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

fit10 <- survfit(Surv(survival, status) ~ in10, data = data)
ggsurvplot(fit10, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

fit11 <- survfit(Surv(survival, status) ~ in11, data = data)
ggsurvplot(fit11, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

fit12 <- survfit(Surv(survival, status) ~ in12, data = data)
ggsurvplot(fit12, data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)


# do individual pairings
data[data$cluster_assignment == 2, "in2in12"] <- "in2"
data[data$cluster_assignment == 12, "in2in12"] <- 'in12'

fit2_12 <- survfit(Surv(survival, status) ~ in2in12, data = data[!is.na(data$in2in12), ])
ggsurvplot(fit2_12, data = data[!is.na(data$in2in12), ],
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)


data[data$cluster_assignment == 2, "in2in3"] <- "in2"
data[data$cluster_assignment == 3, "in2in3"] <- 'in3'

fit2_3 <- survfit(Surv(survival, status) ~ in2in3, data = data[!is.na(data$in2in3), ])
ggsurvplot(fit2_3, data = data[!is.na(data$in2in3), ],
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

data[data$cluster_assignment == 1, "in1in3"] <- "in1"
data[data$cluster_assignment == 3, "in1in3"] <- 'in3'

fit1_3 <- survfit(Surv(survival, status) ~ in1in3, data = data[!is.na(data$in1in3), ])
ggsurvplot(fit1_3, data = data[!is.na(data$in1in3), ],
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

data[data$cluster_assignment == 10, "in10in12"] <- "in10"
data[data$cluster_assignment == 12, "in10in12"] <- 'in12'

fit10_12 <- survfit(Surv(survival, status) ~ in10in12, data = data[!is.na(data$in10in12), ])
ggsurvplot(fit10_12, data = data[!is.na(data$in10in12), ],
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

data[data$cluster_assignment == 1, "in1in12"] <- "in1"
data[data$cluster_assignment == 12, "in1in12"] <- 'in12'

fit1_12 <- survfit(Surv(survival, status) ~ in1in12, data = data[!is.na(data$in1in12), ])
ggsurvplot(fit1_12, data = data[!is.na(data$in1in12), ],
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

data[data$cluster_assignment == 11, "in11in12"] <- "in11"
data[data$cluster_assignment == 12, "in11in12"] <- 'in12'

fit11_12 <- survfit(Surv(survival, status) ~ in11in12, data = data[!is.na(data$in11in12), ])
ggsurvplot(fit11_12, data = data[!is.na(data$in11in12), ],
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

data[data$cluster_assignment == 11, "in11in9"] <- "in11"
data[data$cluster_assignment == 9, "in11in9"] <- 'in9'

fit11_9 <- survfit(Surv(survival, status) ~ in11in9, data = data[!is.na(data$in11in9), ])
ggsurvplot(fit11_9, data = data[!is.na(data$in11in9), ],
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

data[data$cluster_assignment == 12, "in12in9"] <- "in12"
data[data$cluster_assignment == 9, "in12in9"] <- 'in9'

fit12_9 <- survfit(Surv(survival, status) ~ in12in9, data = data[!is.na(data$in12in9), ])
ggsurvplot(fit12_9, data = data[!is.na(data$in12in9), ],
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE)

