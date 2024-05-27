# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #              smoking prediction based on density distribution         # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
###################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggplot2)
library(dplyr)
library(reshape2)
library(grid)
library(viridis)
library(gridExtra)
library(ggbeeswarm)
library(ggpubr)
library(signal)
library(pracma)
library(mixtools)
library(grid)
library(cowplot)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

################################################################################
##################################################################### file paths
signatures_path        <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
sample_table_path      <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
output_path            <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/"
smoking_path           <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/sampleTables/smoking_data.txt"
cancer_gene_mut_path   <- "/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Panlung--hg19.csv"
patient_mutburden_path <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/patient_mutburden.RDS"
signature_cluster_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/new_signature_cluster_patient_assignment.txt"
survival_path          <- "/re_gecip/cancer_lung/shared/pancan_survival_release_v16.RDS"
clin_data_path         <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/sampleTables/v14_lung_patient_clinical_data_2022-07-05_useful.tsv"

output_file            <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/"

################################################################################
##########################################################################  MAIN

# read in smoking data
sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)
sample_table[sample_table$histology == "MET_OTHER", "histology"] <- "OTHER"

smoking_table <- read.table(smoking_path, head = T, sep = "\t")
smoking_table <- smoking_table[, c("participant_id", "smoking_status")]
smoking_table <- unique(smoking_table)
smoking_table$participant_id <- as.character(smoking_table$participant_id)


# read in signature data
load(signatures_path)
SBS_exposure$sample <- as.character(SBS_exposure$sample)
DBS_exposure$sample <- as.character(DBS_exposure$sample)
ID_exposure$sample  <- as.character(ID_exposure$sample)

# we don't have DBS and ID signatures for all samples, make those 0
all_sample <- data.frame(sample_table$participant_id)
colnames(all_sample) <- "sample"

# keep smoking signatures only
SBS_smoking <- SBS_exposure[SBS_exposure$label %in% c("SBS4:smoking", "SBS92:smoking"),]
SBS_smoking_mat <- dcast(SBS_smoking, sample ~ label, value.var = "exposure")

DBS_smoking <- DBS_exposure[DBS_exposure$label == "DBS2:smoking", ]
DBS_smoking_mat <- dcast(DBS_smoking, sample ~ label, value.var = "exposure")

ID_smoking <- ID_exposure[ID_exposure$label == "ID3:smoking", ]
ID_smoking_mat <- dcast(ID_smoking, sample ~ label, value.var = "exposure")

sig_df <- left_join(SBS_smoking_mat, DBS_smoking_mat)
sig_df <- left_join(sig_df, ID_smoking_mat)
sig_df[is.na(sig_df)] <- 0

sig_df$total_smoking_exposure <- sig_df$`SBS4:smoking` + sig_df$`SBS92:smoking` + sig_df$`DBS2:smoking` + sig_df$`ID3:smoking`

# add histology
sig_df <- left_join(sig_df, sample_table[, c("participant_id", "histology")], by = c("sample" = "participant_id"))

# let's have a look at the distributions summed up smoking signatures exposures
sig_df[sig_df$total_smoking_exposure == 0, "total_smoking_exposure"] <- 1
sig_df_ex_mdl <- normalmixEM(log(sig_df$total_smoking_exposure))

pdf(paste0(output_path, "supp_fig2_summed_smoking_exposure_distribution.pdf"), width = 5, height = 5)
plot(sig_df_ex_mdl, which = 2, xlab2 = "# mutations (log)", main2 = "")
dev.off()

# get the posterior probabilities of each sample for each of the distributions
probabilities <- sig_df_ex_mdl$posterior
probabilities <- data.frame(probabilities)

# add sample names to this
probabilities$sample <- sig_df$sample

# give this more intuitive column names
distribution_medians <- sig_df_ex_mdl$mu

colnames(probabilities)[match(max(distribution_medians), distribution_medians)] <- "probability_smoker"
colnames(probabilities)[match(min(distribution_medians), distribution_medians)] <- "probability_non_smoker"

# binary classify if smoker or not
# but under the condition that th absolute differences between the probability smoker and non smoker needs to be bigger than 0.5
# if that is not the case i will classify them as a smoker, because these samples are difficult to classify, usually because they are from ex smokers
# but for our purposes and ex smoker is a smoker
probabilities <- probabilities %>% mutate(classifier = ifelse(probability_non_smoker > probability_smoker, "non_smoker", "smoker"))

# calculate abosulte difference between probablilites
probabilities$abs_difference_prop <- abs(probabilities$probability_smoker - probabilities$probability_non_smoker)
probabilities[probabilities$abs_difference_prop <= 0.95, "classifier"] <- "smoker"

# add the smoking signature values
probabilities <- left_join(probabilities, sig_df)

# now lets add the ground truth smoking data to see how well we're doing
probabilities <- left_join(probabilities, smoking_table[, c("participant_id", "smoking_status")], by = c("sample" = "participant_id"))

probabilities_gt <- probabilities[probabilities$smoking_status != "unknown", ]
probabilities_gt_count <- data.frame(table(probabilities_gt[, c("classifier", "smoking_status")]))

p <- ggplot(probabilities_gt_count, aes(smoking_status, Freq, fill = classifier, label = Freq)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#7ad8fa","#46788a")) +
  ylab("number of samples") +
  xlab("ground truth smoking status") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12))

pdf(paste0(output_path, "supp_fig2_smoking_prediction_ground_truth_all_samples.pdf"), width = 5, height = 5)
p
dev.off()


all <- probabilities

all$corrected_smoking_status <- all$classifier

# add sex
all <- left_join(all, sample_table[, c("participant_id", "sex")], by = c("sample" = "participant_id"))

# calculate some proportions
counts <- all[, c("sample", "corrected_smoking_status", "histology", "sex")]

count_hist <- data.frame(table(counts[, c("corrected_smoking_status", "histology")]))
hist_count <- data.frame(table(sample_table$histology))
colnames(hist_count) <- c("histology", "total_hist_count")
count_hist <- left_join(count_hist, hist_count)
count_hist$prop <- (count_hist$Freq / count_hist$total_hist_count)*100

# plot
hist_order <- c("CARCINOID", "MESOTHELIOMA", "NEUROENDOCRINE_CARCINOMA", "ADENOCARCINOMA", "MET_ADENOCARCINOMA",
                "OTHER", "ADENOSQUAMOUS", "SQUAMOUS_CELL", "MET_SQUAMOUS_CELL",
                "SMALL_CELL", "MET_SMALL_CELL", "LARGE_CELL")

count_hist$histology <- factor(count_hist$histology, levels = hist_order)

histology_names <- c("Carcinoid",
                     "Mesothelioma",
                     "Neuroendocrine\ncarcinoma",
                     "Adenocarcinoma",
                     "Adenocarcinoma\nmetastasis",
                     "Other",
                     "Adenosquamous",
                     "Squamous cell",
                     "Squamonus cell\nmetastasis",
                     "Small cell",
                     "Small cell\nmetastasis",
                     "Large cell")

names(histology_names) <- hist_order

# p_hist <- ggplot(count_hist, aes(corrected_smoking_status, prop, fill = corrected_smoking_status)) +
#   geom_bar(stat = "identity") +
#   ylab("% samples") +
#   xlab(" ") +
#   scale_fill_manual(values = c("#7ad8fa","#46788a")) +
#   scale_x_discrete(labels = c("low\nsmoking\nsignature", "high\nsmoking\nsignatures")) +
#   facet_wrap(~histology, ncol = 1, strip.position = "right", labeller = labeller(histology = histology_names)) +
#   theme_bw() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 14),
#         axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
#         panel.spacing = unit(0.1, 'lines'),
#         plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
#         strip.text.y = element_text(angle = 0),
#         legend.position = "none",
#         strip.background =element_rect(fill="white"))

p_hist <- ggplot(count_hist, aes(histology, prop, fill = corrected_smoking_status)) +
  geom_bar(stat = "identity") +
  ylab("% samples") +
  xlab(" ") +
  scale_fill_manual(name = "inferred smoking status", values = c("#7ad8fa","#46788a"), labels = c("inferred non-smoker", "inferred smoker")) +
  scale_x_discrete(labels = c("CARCINOID" = "Carcinoid",
                             "MESOTHELIOMA" = "Mesothelioma",
                             "NEUROENDOCRINE_CARCINOMA" = "Neuroendocrine\ncarcinoma",
                             "ADENOCARCINOMA" = "Adenocarcinoma",
                             "MET_ADENOCARCINOMA" = "Adenocarcinoma\nmetastasis",
                             "OTHER" = "Other",
                             "ADENOSQUAMOUS" = "Adenosquamous",
                             "SQUAMOUS_CELL" = "Squamous cell",
                             "MET_SQUAMOUS_CELL" = "Squamonus cell\nmetastasis",
                             "SMALL_CELL" = "Small cell",
                             "MET_SMALL_CELL" = "Small cell\nmetastasis",
                             "LARGE_CELL" = "Large cell")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        panel.spacing = unit(0.1, 'lines'),
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.text.y = element_text(angle = 0))

pdf(paste0(output_path, "smoking_proportion_histology.pdf"), width = 9, height = 3)
p_hist
dev.off()

count_sex <- data.frame(table(counts[, c("corrected_smoking_status", "sex")]))
smoke_count <- data.frame(table(counts$corrected_smoking_status))
colnames(smoke_count) <- c("corrected_smoking_status", "total_smoke_count")
count_sex <- left_join(count_sex, smoke_count)
count_sex$prop <- (count_sex$Freq / count_sex$total_smoke_count)*100
count_sex$biological_sex <- "sex"

p_sex <- ggplot(count_sex, aes(corrected_smoking_status, prop, fill = sex)) +
  geom_bar(stat = "identity") +
  ylab("% samples") +
  xlab(" ") +
  scale_fill_manual(values = c("pink", "lightblue")) +
  facet_wrap(~biological_sex, nrow = 1, strip.position = "right") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, angle = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, 'lines'),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
        strip.text.y = element_text(angle = 0),
        strip.background =element_rect(fill="white"),
        legend.position = "top")

# make into one plot
p_all <- plot_grid(p_sex, p_hist, align = "v", axis = "lr", rel_heights = c(1, 7), ncol = 1)

pdf(paste0(output_path, "supp_fig2_smoking_characterisation.pdf"), width = 4, height = 10)
p_all
dev.off()

# get driver gene mutations
cancer_gene_muts <- read.csv(cancer_gene_mut_path, head = T, sep = "\t")
cancer_gene_muts <- cancer_gene_muts[cancer_gene_muts$var_class != "amp", ]
cancer_gene_muts <- cancer_gene_muts[cancer_gene_muts$var_class != "hom.del.", ]
cancer_gene_muts <- unique(cancer_gene_muts[, c("participant_id", "gene_name", "gr_id", "key")])
cancer_gene_muts$gene_cat <- paste0(cancer_gene_muts$gene_name, "_", cancer_gene_muts$gr_id)

# only take those genes that have mutations in more than 10 patients
cancer_gene_count <- cancer_gene_muts[, c("gene_cat", "participant_id")]
cancer_gene_count <- unique(cancer_gene_count)
cancer_gene_count <- data.frame(table(cancer_gene_count$gene_cat))
cancer_gene_count <- cancer_gene_count[cancer_gene_count$Freq >= 10, ]

cancer_gene_muts <- cancer_gene_muts[cancer_gene_muts$gene_cat %in% cancer_gene_count$Var1, ]

EGFR_muts <- cancer_gene_muts[, c("participant_id", "gene_cat") ]
EGFR_muts <- unique(EGFR_muts)
EGFR_muts$participant_id <- as.character(EGFR_muts$participant_id)

# add smoking status
EGFR_muts <- left_join(EGFR_muts, all[, c("sample", "classifier")], by = c("participant_id" = "sample"))

# add mutation burden

gene_mat <- EGFR_muts
gene_mat$val <- TRUE
gene_mat <- dcast(gene_mat, participant_id ~ gene_cat, value.var = "val")
gene_mat[is.na(gene_mat)] <- FALSE

# add the patients that don't have mutations in any of these genes
missing_samples <- sample_table[-which(sample_table$participant_id %in% gene_mat$participant_id), "participant_id"]
add <- matrix(data = FALSE, nrow = length(missing_samples), ncol = ncol(gene_mat), )
colnames(add) <- colnames(gene_mat)
rownames(add) <- missing_samples
add <- data.frame(add)
add$participant_id <- missing_samples
colnames(add) <- sub("\\.", "-", colnames(add))

gene_mat <- rbind(gene_mat, add)

mutburden <- readRDS(patient_mutburden_path)
mutburden$SBS_DBS_ID_mutburden <- mutburden$SBScount + mutburden$DBScount + mutburden$IDcount

gene_mat <- left_join(gene_mat, mutburden[, c("sample", "SBS_DBS_ID_mutburden")], by = c("participant_id" = "sample"))

# add smoking classifier

gene_mat <- left_join(gene_mat, unique(all[, c("sample", "classifier")]), by = c("participant_id" = "sample"))

gene_test <- list()
# for each of these genes do a fishers test of proportion of smoker vs nons-mokers
for(gene in unique(EGFR_muts$gene_cat)){
  
  print(gene)
  df <- gene_mat[, c("participant_id", "SBS_DBS_ID_mutburden", "classifier", gene)]
  colnames(df) <- c("participant_id", "SBS_DBS_ID_mutburden", "classifier", "mutation")
  
  count <- data.frame(table(df[, c("classifier", "mutation")]))
  
  mutation_smoker <-  count[count$classifier == "smoker" & count$mutation == TRUE, "Freq"]
  if(length(mutation_smoker) == 0){mutation_smoker <- 0}
  
  # how many of the samples with this gene mutation are non_smoker
  mutation_non_smoker <- count[count$classifier == "non_smoker" & count$mutation == TRUE, "Freq"]
  if(length(mutation_non_smoker) == 0){mutation_non_smoker <- 0}
  
  # how many of all samples that don't have this mutation are smoker
  no_mutation_smoker <- count[count$classifier == "smoker" & count$mutation == FALSE, "Freq"]
  if(length(no_mutation_smoker) == 0){no_mutation_smoker <- 0}
  
  # how many of all samples that don't have this mutation are non smoker
  no_mutation_non_smoker <- count[count$classifier == "non_smoker" & count$mutation == FALSE, "Freq"]
  if(length(no_mutation_non_smoker) == 0){no_mutation_non_smoker <- 0}

  dat <- data.frame(smoker = c(no_mutation_smoker, mutation_smoker),
                    non_smoker = c(no_mutation_non_smoker, mutation_non_smoker))
  
  rownames(dat) <- c("no_mutation", "mutation")
  
  test <- fisher.test(dat)
  
  # also do a glm on this
  glmtest <- glm(data = df, as.factor(classifier) ~ as.factor(mutation) + SBS_DBS_ID_mutburden,
                 family = binomial)
  
  data_out <- data.frame(gene_element = gene,
                         p_value = test$p.value,
                         glm_p_value = coef(summary(glmtest))[,'Pr(>|z|)'][2],
                         odds_ratio = test$estimate,
                         loCI = test$conf.int[1],
                         hiCI = test$conf.int[2],
                         smoker = mutation_smoker,
                         non_smoker = mutation_non_smoker)
  
  gene_test[[gene]] <- data_out
}

gene_test <- do.call(rbind, gene_test)
gene_test$p_adjust <- p.adjust(gene_test$glm_p_value, method = "fdr")
gene_test$p_fisher_adjust <- p.adjust(gene_test$p_value, method = "fdr")

gene_test$significance <- NULL
gene_test[which(gene_test$p_adjust < 0.0001), "significance"] <- "****"
gene_test[which(gene_test$p_adjust < 0.001 & gene_test$p_adjust >= 0.0001), "significance"] <- "***"
gene_test[which(gene_test$p_adjust < 0.01 & gene_test$p_adjust >= 0.001), "significance"] <- "**"
gene_test[which(gene_test$p_adjust < 0.05 & gene_test$p_adjust >= 0.01), "significance"] <- "*"
gene_test[which(gene_test$p_adjust > 0.05), "significance"] <- ""

# plot this
gene_test <- gene_test[order(gene_test$significance,
                             gene_test$odds_ratio,
                             decreasing = T), ]
# gene_test_p <- gene_test[, c("gene_element", "smoker", "non_smoker", "significance")]
# gene_test_p <- melt(gene_test_p, id.vars = c("gene_element", "significance"))
gene_test$gene_element <- factor(gene_test$gene_element, levels = unique(gene_test$gene_element))
# gene_test_p$value_log <- log10(gene_test_p$value)

# only plot significant ones
# gene_test_p <- gene_test_p[gene_test_p$significance != "", ]

# p <- ggplot(data = gene_test_p, aes(gene_element, variable, fill = variable, alpha = value_log)) +
#             geom_tile(color = "grey") +
#             geom_text(aes(label = significance), color = "black", size = 3) +
#             scale_fill_manual(values=c(non_smoker="#7ad8fa", smoker="#46788a"), labels = c("low\nsmoking\nsignature", "high\nsmoking\nsignatures")) +
#             scale_y_discrete(labels=c("smoker" = "high\nsmoking\nsignatures", "non_smoker" = "low\nsmoking\nsignature")) +
#             theme_bw() +
#             theme(axis.text.x = element_text(size = 12),
#                   axis.text.y = element_text(size = 12),
#                   panel.grid.minor = element_blank(),
#                   panel.grid.major = element_blank(),) +
#             xlab("") +
#             ylab("") +
#             coord_flip()

gene_test[gene_test$odds_ratio < 1, "enriched_in"] <- "smoker"
gene_test[gene_test$odds_ratio >= 1, "enriched_in"] <- "non_smoker"
gene_test[gene_test$significance != "", "is_significant"] <- "significant"
gene_test[gene_test$significance == "", "is_significant"] <- "not_significant"

p <- ggplot(gene_test, aes(x = gene_element, y = odds_ratio, color = enriched_in, alpha = is_significant)) + 
        geom_pointrange(aes(ymin = loCI, ymax = hiCI)) +
        scale_color_manual(values = c("#7ad8fa","#46788a"), labels = c("inferred non-smoker", "inferred smoker"), name = "enriched in") +
        scale_y_log10() +
        scale_alpha_manual(values = c(0.3, 1), labels = c("not_significant", "significant"), name = "") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "blue", linewidth = 0.5) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
              axis.text.y = element_text(size = 12),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),) +
        xlab("") +
        ylab("odds ratio") 


pdf(paste0(output_path, "supp_fig2_smoking_driver.pdf"), width = 14, height = 4)
p
dev.off()

# is this still significant when adjusting for mutation burden
# gene_mat <- EGFR_muts
# gene_mat$val <- TRUE
# gene_mat <- dcast(gene_mat, participant_id ~ gene_cat, value.var = "val")
# gene_mat[is.na(gene_mat)] <- FALSE
# 
# # add the patients that don't have mutations in any of these genes
# missing_samples <- sample_table[-which(sample_table$participant_id %in% gene_mat$participant_id), "participant_id"]
# add <- matrix(data = FALSE, nrow = length(missing_samples), ncol = ncol(gene_mat), )
# colnames(add) <- colnames(gene_mat)
# rownames(add) <- missing_samples
# add <- data.frame(add)
# add$participant_id <- missing_samples
# colnames(add) <- sub("\\.", "-", colnames(add))
# 
# gene_mat <- rbind(gene_mat, add)
# 
# mutburden <- readRDS(patient_mutburden_path)
# mutburden$SBS_DBS_ID_mutburden <- mutburden$SBScount + mutburden$DBScount + mutburden$IDcount
# 
# gene_mat <- left_join(gene_mat, mutburden[, c("sample", "SBS_DBS_ID_mutburden")], by = c("participant_id" = "sample"))
# 
# # add smoking classifier
# 
# gene_mat <- left_join(gene_mat, unique(all[, c("sample", "classifier")]), by = c("participant_id" = "sample"))
#                       
# gene_mat[gene_mat$EGFR_CDS == FALSE, "EGFR_CDS"] <- "no_mutation"
# gene_mat[gene_mat$EGFR_CDS == TRUE, "EGFR_CDS"] <- "mutation"             
# 
# gene_mat[gene_mat$KRAS_CDS == FALSE, "KRAS_CDS"] <- "no_mutation"
# gene_mat[gene_mat$KRAS_CDS == TRUE, "KRAS_CDS"] <- "mutation" 
# 
# gene_mat[gene_mat$TP53_CDS == FALSE, "TP53_CDS"] <- "no_mutation"
# gene_mat[gene_mat$TP53_CDS == TRUE, "TP53_CDS"] <- "mutation" 
# 
# gene_mat[gene_mat$KEAP1_CDS == FALSE, "KEAP1_CDS"] <- "no_mutation"
# gene_mat[gene_mat$KEAP1_CDS == TRUE, "KEAP1_CDS"] <- "mutation" 
# 
# gene_mat[gene_mat$RB1_CDS == FALSE, "RB1_CDS"] <- "no_mutation"
# gene_mat[gene_mat$RB1_CDS == TRUE, "RB1_CDS"] <- "mutation" 
# 
# gene_mat[gene_mat$NOTCH1_CDS == FALSE, "NOTCH1_CDS"] <- "no_mutation"
# gene_mat[gene_mat$NOTCH1_CDS == TRUE, "NOTCH1_CDS"] <- "mutation" 
# 
# gene_mat[gene_mat$STK11_CDS == FALSE, "STK11_CDS"] <- "no_mutation"
# gene_mat[gene_mat$STK11_CDS == TRUE, "STK11_CDS"] <- "mutation" 
# 
# gene_mat[gene_mat$NFE2L2_CDS == FALSE, "NFE2L2_CDS"] <- "no_mutation"
# gene_mat[gene_mat$NFE2L2_CDS == TRUE, "NFE2L2_CDS"] <- "mutation" 
# 
# gene_mat[gene_mat$ARID1A_CDS == FALSE, "ARID1A_CDS"] <- "no_mutation"
# gene_mat[gene_mat$ARID1A_CDS == TRUE, "ARID1A_CDS"] <- "mutation" 
# 
# gene_mat[gene_mat$ZNF578_enhancer == FALSE, "ZNF578_enhancer"] <- "no_mutation"
# gene_mat[gene_mat$ZNF578_enhancer == TRUE, "ZNF578_enhancer"] <- "mutation" 
# 
#                       
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(EGFR_CDS) + SBS_DBS_ID_mutburden,
#              family = "binomial")
# 
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(KRAS_CDS) + SBS_DBS_ID_mutburden,
#              family = "binomial")
# 
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(TP53_CDS) + SBS_DBS_ID_mutburden,
#              family = "binomial")
# 
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(KEAP1_CDS) + SBS_DBS_ID_mutburden,
#              family = "binomial")
# 
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(RB1_CDS) + SBS_DBS_ID_mutburden,
#              family = "binomial")
# 
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(NOTCH1_CDS) + SBS_DBS_ID_mutburden,
#              family = "binomial")
# 
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(STK11_CDS) + SBS_DBS_ID_mutburden,
#              family = "binomial")
# 
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(NFE2L2_CDS) + SBS_DBS_ID_mutburden,
#              family = "binomial")
# 
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(ARID1A_CDS) + SBS_DBS_ID_mutburden,
#              family = "binomial")
# 
# model <- glm(data = gene_mat, as.factor(classifier) ~ as.factor(ZNF578_enhancer) + SBS_DBS_ID_mutburden,
#              family = "binomial")


# add the corrected smoking status to overall smoking table
# smoking_table <- left_join(smoking_table, all, by = c("participant_id" = "sample"))
# smoking_table$smoking_status.y <- NULL
# colnames(smoking_table)[2] <- "reported_smoking_status"
# 
# write.table(smoking_table, paste0(output_file, "smoking_data.txt"), sep = "\t", row.names = F, quote = F)

# how many smoking associated mutations are there in the 9 non-smoking LUSCs?
non_smoker_LUSCs <- all[all$classifier == "non_smoker" & all$histology == "SQUAMOUS_CELL", ]
median_non_smoker_LUSC <- median(non_smoker_LUSCs$total_smoking_exposure)
  
smoker_LUSCs <- all[all$classifier == "smoker" & all$histology == "SQUAMOUS_CELL", ]
median_smoker_LUSCs <- median(smoker_LUSCs$total_smoking_exposure)

# how many mutations are there in total in the non-smoker LUSCs
patient_mutburden <- readRDS(patient_mutburden_path)
patient_mutburden$total_mutburden <- patient_mutburden$SBScount + patient_mutburden$IDcount + patient_mutburden$DBScount

med_non_smoker_LUSC_mutburden <-   median(patient_mutburden[patient_mutburden$sample %in% non_smoker_LUSCs$sample, "total_mutburden"])
med_smoker_LUSC_mutburden     <-   median(patient_mutburden[patient_mutburden$sample %in% smoker_LUSCs$sample, "total_mutburden"])

#################################################################################
# test for survival differencesf of patients with EGFR mutations in group 2 and not i group2 
cancer_gene_muts <- read.csv(cancer_gene_mut_path, head = T, sep = "\t")
cancer_gene_muts <- cancer_gene_muts[cancer_gene_muts$var_class != "amp", ]
cancer_gene_muts <- cancer_gene_muts[cancer_gene_muts$var_class != "hom.del.", ]
cancer_gene_muts <- unique(cancer_gene_muts[, c("participant_id", "gene_name", "gr_id", "key")])
cancer_gene_muts$gene_cat <- paste0(cancer_gene_muts$gene_name, "_", cancer_gene_muts$gr_id)
EGFR_gene_muts <- cancer_gene_muts[cancer_gene_muts$gene_name == "EGFR", ]
EGFR_gene_muts <- unique(EGFR_gene_muts[, c("participant_id", "gene_name"), ])

# add sigantuer clusters
cluster_assignment <- read.table(signature_cluster_path, head = T, sep = "\t")
EGFR_gene_muts <- left_join(EGFR_gene_muts, cluster_assignment, by = c("participant_id" = "variable"))
EGFR_gene_muts$participant_id <- as.character(EGFR_gene_muts$participant_id)

# add survival data
survival <- readRDS(survival_path)
survival$participant_id <- as.character(survival$participant_id)
survival$survival <- sub(" days", "", survival$survival)
survival$survival <- as.numeric(survival$survival)
survival <- survival[survival$participant_id %in% sample_table$participant_id, ]
survival <- survival[survival$disease_type == "LUNG", ]

EGFR_gene_muts <- left_join(EGFR_gene_muts, survival)

# make a flag of being in lcuster 2 or not
EGFR_gene_muts[EGFR_gene_muts$cluster_assignment == 2, "cluster2"] <- TRUE
EGFR_gene_muts[EGFR_gene_muts$cluster_assignment != 2, "cluster2"] <- FALSE

# add histology
EGFR_gene_muts <- left_join(EGFR_gene_muts, sample_table[, c("participant_id", "histology")])

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

EGFR_gene_muts <- left_join(EGFR_gene_muts, stage)


fit <- survfit(Surv(survival, status) ~ cluster2, data = EGFR_gene_muts[EGFR_gene_muts$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")
                                                                        & EGFR_gene_muts$stage_stage_best %in% c("1B", "2A", "3A", "2B") & EGFR_gene_muts$survival >= 90, ])
p <- ggsurvplot(fit, data = EGFR_gene_muts[EGFR_gene_muts$histology %in% c("ADENOCARCINOMA", "SQUAMOUS_CELL")
                                           & EGFR_gene_muts$stage_stage_best %in% c("1B", "2A", "3A", "2B") & EGFR_gene_muts$survival >= 90, ],
                conf.int = TRUE,
                pval = TRUE,
                risk.table = TRUE,
                palette = c("#0215de", "#a7cee3"))



