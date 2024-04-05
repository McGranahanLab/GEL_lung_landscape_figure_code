# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #             chromothripsis signature cluster exlploration           # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
##################################################################### libraries
library(dplyr)
library(reshape2)
library(ggplot2)

################################################################################
##################################################################### file paths

signatures_path         <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
sample_table_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
cluster_assignment_path <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/new_signature_cluster_patient_assignment.txt"
purity_path             <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/CopyNumber/combined_ascatTable_segments_joined.rds"
ecdna_path              <- "/re_gecip/cancer_lung/CBailey1/pancan_ecdna_analysis/data/summary_objects/ALL_CANCER_TYPES_BED_LIST_2023-01-18.RDS"
output_path             <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/clustering/ALL_new/" 

num_cluster <- 12

################################################################################
#####################################################################       MAIN

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)
cluster_assignment <- read.table(cluster_assignment_path, head = T, sep = "\t")
cluster_assignment$variable <- as.character(cluster_assignment$variable)

amp_data <- readRDS(ecdna_path)
amp_data <- amp_data[amp_data$participant_id %in% sample_table$participant_id, ]
amp_data <- left_join(amp_data, sample_table[, c("participant_id", "histology")])

# classify patients by having ecdna or not
amp_patient <- amp_data[amp_data$feature_type == "ecDNA", ]
amp_patient <- unique(amp_patient[, c("participant_id", "feature_type")])
amp_patient <- left_join(amp_patient, cluster_assignment, by = c("participant_id" = "variable")) 

data <- cluster_assignment
colnames(data) <- c("cluster_assignment", "participant_id")
data <- left_join(data, amp_patient)
data[is.na(data$feature_type), "feature_type"] <- "no_ecdna"

# test for enrichment of ecdna in clusters
all_test <- list()
for(i in 1:num_cluster){
  print(i)
  
  cluster_with_ecdna      <- nrow(data[data$cluster_assignment == i & data$feature_type == "ecDNA", ])
  cluster_no_ecdna        <- nrow(data[data$cluster_assignment == i & data$feature_type == "no_ecdna", ])
  not_cluster_with_ecdna  <- nrow(data[data$cluster_assignment != i & data$feature_type == "ecDNA", ])
  not_cluster_no_ecdna    <- nrow(data[data$cluster_assignment != i & data$feature_type == "no_ecdna", ])
  
  dat <- data.frame(
    "ecdna_yes" = c(cluster_with_ecdna, not_cluster_with_ecdna),
    "ecdna_no" = c(cluster_no_ecdna, not_cluster_no_ecdna),
    row.names = c("in_cluster", "other_cluster"),
    stringsAsFactors = FALSE)
  
  test <- fisher.test(dat)
  
  test_df <- data.frame(cluster_id = i, 
                        in_cluster_with_ecdna = cluster_with_ecdna,
                        in_cluster_no_ecdna = cluster_no_ecdna,
                        other_cluster_with_ecdna = not_cluster_with_ecdna,
                        other_cluster_no_ecdna = not_cluster_no_ecdna,
                        fisher_p = test$p.value,
                        fisher_odd = test$estimate,
                        fisher_loCI = test$conf.int[1],
                        fisher_hiCI = test$conf.int[2])
  
  all_test[[i]] <- test_df
  
  
}

gene_test <- do.call(rbind, all_test)
gene_test$p_adjust <- p.adjust(gene_test$fisher_p, method = "fdr")

gene_test$significance <- NULL
gene_test[which(gene_test$p_adjust < 0.0001), "significance"] <- "****"
gene_test[which(gene_test$p_adjust < 0.001 & gene_test$p_adjust >= 0.0001), "significance"] <- "***"
gene_test[which(gene_test$p_adjust < 0.01 & gene_test$p_adjust >= 0.001), "significance"] <- "**"
gene_test[which(gene_test$p_adjust < 0.05 & gene_test$p_adjust >= 0.01), "significance"] <- "*"
gene_test[which(gene_test$p_adjust > 0.05), "significance"] <- ""

gene_test <- gene_test[order(gene_test$significance,
                             gene_test$fisher_odd,
                             decreasing = T), ]

gene_test$cluster_id <- factor(gene_test$cluster_id, levels = unique(gene_test$cluster_id))

gene_test[gene_test$fisher_odd < 1, "enriched_in"] <- "no_ecdna"
gene_test[gene_test$fisher_odd >= 1, "enriched_in"] <- "ecdna"
gene_test[gene_test$significance != "", "is_significant"] <- "significant"
gene_test[gene_test$significance == "", "is_significant"] <- "not_significant"

p <- ggplot(gene_test, aes(x = cluster_id, y = fisher_odd, color = enriched_in, alpha = is_significant)) + 
  geom_pointrange(aes(ymin = fisher_loCI, ymax = fisher_hiCI)) +
  scale_color_manual(values = c("#fc8d59","#91bfdb"), labels = c("tumours with ecDNA", "tumours without ecDNA"), name = "enriched with") +
  scale_y_log10() +
  scale_alpha_manual(values = c(0.3, 1), labels = c("not_significant", "significant"), name = "") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", linewidth = 0.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),) +
  xlab("signature cluster") +
  ylab("odds ratio") 

pdf(paste0(output_path, "cluster_ecDNA_enrichment.pdf"), width = 4.5, height = 2)
p
dev.off()
