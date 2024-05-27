# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # #             signature clustering survival analysis                  # # #
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

# assemble data

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

# add histology, timepoint, age, sex
data <- left_join(data, sample_table[, c("participant_id", "histology", "time_point", "patient_age", "sex")])

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

# add smoking status
smoking <- read.table(smoking_path, head = T, sep = "\t")
smoking$participant_id <- as.character(smoking$participant_id)

data <- left_join(data, smoking[, c("participant_id", "classifier")])

# add whether patients had treatment
treatment <- clinical[, c("participant_id", "chemo_start_date_of_regimen", "RT_treatmentstartdate")]
treatment[is.na(treatment$chemo_start_date_of_regimen), "hadCHEMO"] <- FALSE
treatment[!is.na(treatment$chemo_start_date_of_regimen), "hadCHEMO"] <- TRUE
treatment[is.na(treatment$RT_treatmentstartdate), "hadRT"] <- FALSE
treatment[!is.na(treatment$RT_treatmentstartdate), "hadRT"] <- TRUE

# add CPI treatmemt
CPI <- read.csv(CPI_response_path, head = T)
CPI$CPI <- TRUE
CPI$Participant.Id <- as.character(CPI$Participant.Id)

treatment <- left_join(treatment, CPI[, c("Participant.Id", "CPI")], by = c("participant_id" = "Participant.Id"))
treatment[is.na(treatment$CPI), "CPI"] <- FALSE
treatment[treatment$hadCHEMO & treatment$hadRT & treatment$CPI, "treatment_group"] <- "chemo_RT_CPI"
treatment[treatment$hadCHEMO & treatment$hadRT & treatment$CPI == FALSE, "treatment_group"] <- "chemo_RT"
treatment[treatment$hadCHEMO & treatment$hadRT == FALSE & treatment$CPI == TRUE, "treatment_group"] <- "chemo_CPI"
treatment[treatment$hadCHEMO & treatment$hadRT == FALSE & treatment$CPI == FALSE, "treatment_group"] <- "chemo"
treatment[treatment$hadCHEMO == FALSE & treatment$hadRT == FALSE & treatment$CPI == FALSE, "treatment_group"] <- "none"
treatment[treatment$hadCHEMO == FALSE & treatment$hadRT == TRUE & treatment$CPI == FALSE, "treatment_group"] <- "RT"
treatment[treatment$hadCHEMO == FALSE & treatment$hadRT == TRUE & treatment$CPI == TRUE, "treatment_group"] <- "RT_CPI"
treatment[treatment$hadCHEMO == FALSE & treatment$hadRT == FALSE & treatment$CPI == TRUE, "treatment_group"] <- "CPI"
treatment[treatment$treatment_group == "none", "treatment_binary"] <- "none"
treatment[treatment$treatment_group != "none", "treatment_binary"] <- "had_treatment"

data <- left_join(data, treatment[, c("participant_id", "treatment_group", "treatment_binary")])

# remove all tumours that died within 90 days, all mets and all stage 4
data <- data[data$survival >= 90, ]
data <- data[-which(data$stage == 4), ]
data <- data[data$time_point == "PRIMARY", ]

data <- data[!is.na(data$status), ]
data[data$status == 1, "new_status"] <- 0
data[data$status == 2, "new_status"] <- 1
data_df <- data

# perform survival analysis
s <- coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age + stage, data = data_df)
p <- ggforest(coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age, data = data_df[data_df$histology != "CARCINOID", ]))

# # make a survival plot of samples in cluster 9 and not in cluster 9
# 
# data_df[data_df$cluster_assignment == 9, "in_cluster9"] <- TRUE
# data_df[data_df$cluster_assignment != 9, "in_cluster9"] <- FALSE
# 
# fit <- survfit(Surv(survival, new_status) ~ in_cluster9 + histology + patient_age + stage + treatment_group + sex, data = data_df)
# ggsurvplot(fit,  
#            size = 1,
#            pval = TRUE,
#            risk.table = TRUE)




################################################################################

# only perform survival analysis on the three biggest clusters
data_bigcluster <- data_df[data_df$cluster_assignment %in% c(2, 10, 12), ]

s_bigcluster <- coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age + stage + treatment_group + sex, data = data_bigcluster)


################################################################################
# check survival within cluster 1

# add TCRA scores to data
TCRA <- readRDS(TCRA_path)
TCRA <- TCRA[TCRA$type == "tumour", ]
data_df <- left_join(data_df, TCRA[, c("patient", "TCRA.fraction.gc.adj")], by = c("participant_id" = "patient"))

cluster1_df <- data_df[data_df$cluster_assignment == 1, ]

# classify TCRA by high and low
cluster1_df[which(cluster1_df$TCRA.fraction.gc.adj >= median(cluster1_df$TCRA.fraction.gc.ad, na.rm = T)), "TCRA_group"] <- "high"
cluster1_df[which(cluster1_df$TCRA.fraction.gc.adj < median(cluster1_df$TCRA.fraction.gc.ad, na.rm = T)), "TCRA_group"] <- "low"

s_cluster1 <- coxph(Surv(survival, new_status) ~ TCRA_group + histology + patient_age + stage + treatment_group + sex, data = cluster1_df)


fit <- survfit(Surv(survival, new_status) ~ TCRA_group, data = cluster1_df)
p <- ggsurvplot(fit, data = cluster1_df,
                conf.int = TRUE,
                pval = TRUE,
                risk.table = FALSE,
                palette = c("#0215de", "#a7cee3"))

pdf(paste0(output_path, "cluster1_survival_tcra.pdf"), width = 3, height = 3)
p
dev.off()


# is this also the case for any of the other groups
data_df <- data_df %>% group_by(cluster_assignment) %>% mutate(cluster_TCRA_group = ifelse(TCRA.fraction.gc.adj >= median(TCRA.fraction.gc.adj, na.rm = T), "high", "low"))
data_df <- data.frame(data_df)

summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "1", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "2", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "3", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "4", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "5", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "6", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "7", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "8", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "9", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "10", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "11", ]))
summary(coxph(Surv(survival, new_status) ~ cluster_TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df[data_df$cluster_assignment == "12", ]))

# check this for all samples
data_df[which(data_df$TCRA.fraction.gc.adj >= median(data_df$TCRA.fraction.gc.ad, na.rm = T)), "TCRA_group"] <- "high"
data_df[which(data_df$TCRA.fraction.gc.adj < median(data_df$TCRA.fraction.gc.ad, na.rm = T)), "TCRA_group"] <- "low"

summary(coxph(Surv(survival, new_status) ~ TCRA_group + histology + patient_age + stage + treatment_group + sex, data = data_df))



################################################################################
# do survival analysis on all clusters, only LUAD and LUSC, only clusters with more than 10 patients, no mets, more than 80 days survival adjust for stage

data_df_ll <- data_df[data_df$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL"), ]

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll)
p <- ggsurvplot(fit, data = data_df_ll,
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

summary(coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age + stage, data = data_df_ll))




fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("5", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("5", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("6", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("6", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("1", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("1", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("2", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("2", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("4", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("4", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("7", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("7", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("9", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("9", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("10", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("10", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("11", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("11", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("12", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("12", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("12", "10"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("12", "10"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("6", "4"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("6", "4"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("11", "8"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("11", "8"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("1", "9"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("1", "9"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("4", "9"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("4", "9"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)


fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("2", "6"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("2", "6"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("2", "5"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("2", "5"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "7"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "7"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "7", "3"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "7", "3"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "7", "3", "5", "11"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "7", "3", "5", "11"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "5"),])
p <- ggsurvplot(fit, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "5"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

summary(coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age + stage, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "5"),]))

summary(coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age + stage, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "7", "3", "5", "11"),]))
summary(coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age + stage, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "7"),]))
summary(coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age + stage, data = data_df_ll[data_df_ll$cluster_assignment %in% c("8", "3"),]))

# check for survival differences inn different groups of stage 1 tumours

summary(coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age, data = data_df_ll[which(data_df_ll$stage == "1"),]))

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[which(data_df_ll$stage == "1"),])
p <- ggsurvplot(fit, data = data_df_ll[which(data_df_ll$stage == "1"),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)

summary(coxph(Surv(survival, new_status) ~ cluster_assignment + histology + patient_age, data = data_df_ll[which(data_df_ll$stage %in% c("2", "3")),]))

fit <- survfit(Surv(survival, new_status) ~ cluster_assignment, data = data_df_ll[which(data_df_ll$stage %in% c("2", "3")),])
p <- ggsurvplot(fit, data = data_df_ll[which(data_df_ll$stage %in% c("2", "3")),],
                conf.int = FALSE,
                pval = TRUE,
                risk.table = FALSE)


