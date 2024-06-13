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

gene_test$gene_element <- factor(gene_test$gene_element, levels = unique(gene_test$gene_element))


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
