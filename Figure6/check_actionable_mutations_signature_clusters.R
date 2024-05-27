# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #         check KRAS muttations in signature cluster 1                # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
##################################################################### libraries

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
##################################################################### file paths

sample_table_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
cluster_assignment_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/new_signature_cluster_patient_assignment.txt"
muttable_path           <- "/re_gecip/cancer_lung/kthol/landscape/input/muttable.RData"
survival_path           <- "/re_gecip/cancer_lung/shared/pancan_survival_release_v16.RDS"
clin_data_path          <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/sampleTables/v14_lung_patient_clinical_data_2022-07-05_useful.tsv"

output_path             <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/survival/"

################################################################################
#####################################################################       Main

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)
sample_table[sample_table$histology == "MET_OTHER", "histology"] <- "OTHER"

cluster_assignment <- read.table(cluster_assignment_path, head = T, sep = "\t")
samples <- cluster_assignment[cluster_assignment$cluster_assignment == 12, "variable"]

muttable <- get(load(muttable_path))
muttable <- muttable[muttable$Use.For.Plots == TRUE, ]
muttable <- muttable[muttable$sample %in% samples, ]

KRAS_muttable <- muttable[grep("KRAS", muttable$Gene.refGene), ]

# how many are KRAS g12c?
KRAS_g12c_muttable <- KRAS_muttable[grep("g12c", KRAS_muttable$AAChange.refGene, ignore.case = TRUE), ]

# how many tumours in the whole cohort have KRAS g12c mutations?
all_muttable <- get(load(muttable_path))
all_muttable <- all_muttable[muttable$Use.For.Plots == TRUE, ]
all_KRAS_g12c_samples <- all_muttable[grep("KRAS", all_muttable$Gene.refGene), ]
all_KRAS_g12c_muttable <- all_KRAS_g12c_samples[grep("g12c", all_KRAS_g12c_samples$AAChange.refGene, ignore.case = TRUE), ]
cluster_assignment$variable <- as.character(cluster_assignment$variable)
all_KRAS_g12c_muttable <- left_join(all_KRAS_g12c_muttable, cluster_assignment, by = c("sample" = "variable"))

KRAS_samples <- unique(all_KRAS_g12c_muttable$sample)

cluster_assignment[cluster_assignment$variable %in% KRAS_samples, "KRAS"] <- TRUE
cluster_assignment[!cluster_assignment$variable %in% KRAS_samples, "KRAS"] <- FALSE
 
# do a fishers test on this

dat <- data.frame(
  "KRAS_G12C_yes" = c(length(unique(cluster_assignment[cluster_assignment$cluster_assignment == "12" & cluster_assignment$KRAS == TRUE, "variable"])), length(unique(cluster_assignment[cluster_assignment$cluster_assignment != "12" & cluster_assignment$KRAS == TRUE, "variable"]))),
  "KRAS_G12C_no" = c(length(unique(cluster_assignment[cluster_assignment$cluster_assignment == "12" & cluster_assignment$KRAS == FALSE, "variable"])), c(length(unique(cluster_assignment[cluster_assignment$cluster_assignment != "12" & cluster_assignment$KRAS == FALSE, "variable"])))),
  row.names = c("in_cluster12", "other_cluster"),
  stringsAsFactors = FALSE
)

chisq.test(dat)

# do we also have info on met exon 14 skipping?
met_muttable <- all_muttable[grep("MET", all_muttable$Gene.refGene), ]


# do survival analysis of tumours in signature group 12 with KRASG12C vs those in other groups with KRASG12C

KRAS_samples <- unique(all_KRAS_g12c_muttable$sample)

cluster_assignment[cluster_assignment$variable %in% KRAS_samples, "KRAS"] <- TRUE
cluster_assignment[!cluster_assignment$variable %in% KRAS_samples, "KRAS"] <- FALSE
cluster_assignment$variable <- as.character(cluster_assignment$variable)

#add histology
cluster_assignment <- left_join(cluster_assignment, sample_table[, c("participant_id", "histology")], by = c("variable" = "participant_id"))

survival <- readRDS(survival_path)
survival$participant_id <- as.character(survival$participant_id)
survival$survival <- sub(" days", "", survival$survival)
survival$survival <- as.numeric(survival$survival)
survival <- survival[survival$participant_id %in% sample_table$participant_id, ]
survival <- survival[survival$disease_type == "LUNG", ]

cluster_assignment <- left_join(cluster_assignment, survival[, c("participant_id", "status", "survival")], by = c("variable" = "participant_id"))

# make an identifier waysing wheter the KRAS tumour is in cluster 12 or not in cluster12
KRAS_tumours <- cluster_assignment[cluster_assignment$KRAS == TRUE, ]
KRAS_tumours[KRAS_tumours$cluster_assignment == "12", "cluster12"] <- TRUE
KRAS_tumours[KRAS_tumours$cluster_assignment != "12", "cluster12"] <- FALSE

# kick out tumours with a survival that's longer than 10 years
KRAS_tumours <- KRAS_tumours[KRAS_tumours$survival <= (10*365), ]

# thise tumours with a survival of more than 5 years censor them at 5 years
KRAS_tumours[KRAS_tumours$survival >= (5*365), "status"] <- 1
KRAS_tumours[KRAS_tumours$survival >= (5*365), "survival"] <- (5*365)


# univariate analysis

fit <- survfit(Surv(survival, status) ~ cluster12, data = KRAS_tumours)
p <- ggsurvplot(fit, data = KRAS_tumours,
                conf.int = TRUE,
                pval = TRUE,
                risk.table = TRUE,
                palette = c("#0215de", "#a7cee3"))

pdf(paste0(output_path, "cluster12_KRASG12C_survival.pdf"), width = 7, height = 7)
p
dev.off()

# only do for luad
fit <- survfit(Surv(survival, status) ~ cluster12, data = KRAS_tumours[KRAS_tumours$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL") & KRAS_tumours$survival >= 90, ])
p <- ggsurvplot(fit, data = KRAS_tumours[KRAS_tumours$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL") & KRAS_tumours$survival >= 90, ],
                conf.int = TRUE,
                pval = TRUE,
                risk.table = TRUE,
                palette = c("#0215de", "#a7cee3"))

pdf(paste0(output_path, "cluster12_KRASG12C_survival_LUAD_LUSC.pdf"), width = 7, height = 7)
p
dev.off()

# only do for luad and survival >90days
fit <- survfit(Surv(survival, status) ~ cluster12, data = KRAS_tumours[KRAS_tumours$histology == "ADENOCARCINOMA" & KRAS_tumours$survival >= 90, ])
p <- ggsurvplot(fit, data = KRAS_tumours[KRAS_tumours$histology == "ADENOCARCINOMA", ],
                conf.int = TRUE,
                pval = TRUE,
                risk.table = TRUE,
                palette = c("#0215de", "#a7cee3"))

pdf(paste0(output_path, "cluster12_KRASG12C_survival_LUAD_90days.pdf"), width = 7, height = 7)
p
dev.off()

# adjust for stage and histology
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



KRAS_tumours <- left_join(KRAS_tumours, stage, by = c("variable" = "participant_id"))

s <- coxph(Surv(survival, status) ~ cluster12  + stage, data = KRAS_tumours[KRAS_tumours$histology %in% c("ADENOCARCINOMA") & KRAS_tumours$survival >= 90, ])
s <- coxph(Surv(survival, status) ~ cluster12  + stage + histology, data = KRAS_tumours[KRAS_tumours$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL") & KRAS_tumours$survival >= 90, ])
s <- coxph(Surv(survival, status) ~ cluster12  + stage + histology, data = KRAS_tumours[KRAS_tumours$survival >= 90, ])

############################################################################
# also check for other KRASg mutations

# how many tumours in the whole cohort have KRAS g12 mutations?

all_KRAS_g_samples <- all_muttable[grep("KRAS", all_muttable$Gene.refGene), ]
all_KRAS_g_samples <- all_KRAS_g_samples[grep("p\\.g", all_KRAS_g_samples$AAChange.refGene, ignore.case = TRUE), ]
cluster_assignment$variable <- as.character(cluster_assignment$variable)
all_KRAS_g_samples <- left_join(all_KRAS_g_samples, cluster_assignment, by = c("sample" = "variable"))

KRAS_samples <- unique(all_KRAS_g_samples$sample)

cluster_assignment[cluster_assignment$variable %in% KRAS_samples, "KRAS"] <- TRUE
cluster_assignment[!cluster_assignment$variable %in% KRAS_samples, "KRAS"] <- FALSE


# make an identifier waysing wheter the KRAS tumour is in cluster 12 or not in cluster12
KRAS_tumours <- cluster_assignment[cluster_assignment$KRAS == TRUE, ]
KRAS_tumours[KRAS_tumours$cluster_assignment == "12", "cluster12"] <- TRUE
KRAS_tumours[KRAS_tumours$cluster_assignment != "12", "cluster12"] <- FALSE

# kick out tumours with a survival that's longer than 10 years
KRAS_tumours <- KRAS_tumours[KRAS_tumours$survival <= (10*365), ]

# thise tumours with a survival of more than 5 years censor them at 5 years
KRAS_tumours[KRAS_tumours$survival >= (5*365), "status"] <- 1
KRAS_tumours[KRAS_tumours$survival >= (5*365), "survival"] <- (5*365)


# univariate analysis

fit <- survfit(Surv(survival, status) ~ cluster12, data = KRAS_tumours)
p <- ggsurvplot(fit, data = KRAS_tumours,
                conf.int = TRUE,
                pval = TRUE,
                risk.table = TRUE,
                palette = c("#0215de", "#a7cee3"))

pdf(paste0(output_path, "cluster12_KRASG_survival.pdf"), width = 7, height = 7)
p
dev.off()


col <- c("black", "#fff199")

# only do for luad
fit <- survfit(Surv(survival, status) ~ cluster12, data = KRAS_tumours[KRAS_tumours$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL") & KRAS_tumours$survival >= 90, ])
p <- ggsurvplot(fit, data = KRAS_tumours[KRAS_tumours$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL") & KRAS_tumours$survival >= 90, ],
                conf.int = TRUE,
                pval = TRUE,
                risk.table = TRUE,
                palette = col)

pdf(paste0(output_path, "cluster12_KRASG_survival_LUAD_LUSC.pdf"), width = 7, height = 7)
p
dev.off()

# only do for luad and survival >90days
fit <- survfit(Surv(survival, status) ~ cluster12, data = KRAS_tumours[KRAS_tumours$histology == "ADENOCARCINOMA" & KRAS_tumours$survival >= 90, ])
p <- ggsurvplot(fit, data = KRAS_tumours[KRAS_tumours$histology == "ADENOCARCINOMA", ],
                conf.int = TRUE,
                pval = TRUE,
                risk.table = TRUE,
                palette = c("#0215de", "#a7cee3"))

pdf(paste0(output_path, "cluster12_KRASG_survival_LUAD_90days.pdf"), width = 7, height = 7)
p
dev.off()

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



KRAS_tumours <- left_join(KRAS_tumours, stage, by = c("variable" = "participant_id"))

s <- coxph(Surv(survival, status) ~ cluster12  + stage, data = KRAS_tumours[KRAS_tumours$histology %in% c("ADENOCARCINOMA") & KRAS_tumours$survival >= 90, ])
s <- coxph(Surv(survival, status) ~ cluster12  + stage + histology, data = KRAS_tumours[KRAS_tumours$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL") & KRAS_tumours$survival >= 90, ])
s <- coxph(Surv(survival, status) ~ cluster12  + stage + histology, data = KRAS_tumours[KRAS_tumours$survival >= 90, ])
