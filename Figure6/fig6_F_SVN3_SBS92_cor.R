# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                       SV signature exploration                      # # #   
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
####################################################################   libraries

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
library(grid)
library(cowplot)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

################################################################################
####################################################################  file paths

signatures_path        <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
sample_table_path      <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
smoking_path           <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/smoking_data.txt"
output_path            <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/"
mutburden_path         <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/patient_mutburden.RDS"
shatterseek_path       <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/gridss_SV_shatterseek_summary_object.Rdata"
SV_sig_input_path      <- "~/re_gecip/cancer_lung/analysisResults/GEL/release_v8/4.MutationSignatures/HDP/Primary_Mets/SV/input_files/input_32matrix_patient.rds"

################################################################################
####################################################################        MAIN

# load signatures
load(signatures_path)

# make into one data frame
signatures_df           <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure, MDS_exposure)

# make into matrix
signatures_df_weight <- dcast(signatures_df, sample ~ label, value.var = "weight")
signatures_df_weight[is.na(signatures_df_weight)] <- 0

SV_weight <- signatures_df_weight[, c("sample", grep("SV", colnames(signatures_df_weight), value = T))]
SV_weight <- melt(SV_weight, id.vars = "sample")
colnames(SV_weight) <- c("sample", "label", "weight")

signatures_df_exposure <- dcast(signatures_df, sample ~ label, value.var = "exposure")
signatures_df_exposure[is.na(signatures_df_exposure)] <- 0

SV_exposure <- signatures_df_exposure[, c("sample", grep("SV", colnames(signatures_df_exposure), value = T))]
SV_exposure <- melt(SV_exposure, id.vars = "sample")
colnames(SV_exposure) <- c("sample", "label", "exposure")

SV_exposure <- left_join(SV_exposure, SV_weight)

# add histology
sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)

SV_exposure <- left_join(SV_exposure, sample_table[, c("participant_id", "histology")], by = c("sample" = "participant_id"))

smoking <- read.table(smoking_path, head = T, sep = "\t")
smoking$participant_id <- as.character(smoking$participant_id)

SV_exposure <- left_join(SV_exposure, smoking[, c("participant_id", "classifier")], by = c("sample" = "participant_id"))


# do this as well for all signatures
# and adjust for SV burden
mutburden <- readRDS(mutburden_path)
SV_exposure <- left_join(SV_exposure, mutburden[, c("sample", "SVcount")])

fish_df_all <- list()
for(sig in unique(SV_exposure$label)){
  print(sig)
  
  df <- SV_exposure[SV_exposure$label == sig, ]
  
  smoker_with_sig <- nrow(df[df$weight >= 0.05 & df$classifier == "smoker", ])
  non_smoker_with_sig <- nrow(df[df$weight >= 0.05 & df$classifier == "non_smoker", ])
  
  smoker_without_sig <- nrow(df[df$weight < 0.05 & df$classifier == "smoker", ])
  non_smoker_without_sig <- nrow(df[df$weight < 0.05 & df$classifier == "non_smoker", ])
  
  dat <- data.frame(has_sig = c(smoker_with_sig, non_smoker_with_sig),
                    no_sig = c(smoker_without_sig, non_smoker_without_sig))
  
  rownames(dat) <- c("smoker", "non_smoker")
  
  test <- fisher.test(dat)
  
  test_flip <- fisher.test(dat[, c(2, 1)])
  
  df[df$weight >= 0.05, "presence"] <- "present"
  df[df$weight < 0.05, "presence"] <- "absent"
  
  glmtest <- glm(data = df, as.factor(classifier) ~ as.factor(presence) + SVcount,
                 family = binomial)
  
  df_out <- data.frame(signature = sig,
                       odds = test$estimate,
                       odds_ratio_flip = test_flip$estimate,
                       loCI = test$conf.int[1],
                       loCI_flip = test_flip$conf.int[1],
                       hiCI = test$conf.int[2],
                       hiCI_flip = test_flip$conf.int[2],
                       pval = test$p.value,
                       glm_p_value = coef(summary(glmtest))[,'Pr(>|z|)'][2])
  
  fish_df_all[[sig]] <- df_out
}

fish_df_all <- do.call(rbind, fish_df_all)
fish_df_all$p_adjust_fish <- p.adjust(fish_df_all$pval, method = "fdr")
fish_df_all$p_adjust <- p.adjust(fish_df_all$glm_p_value, method = "fdr")

fish_df_all$significance <- NULL
fish_df_all[which(fish_df_all$p_adjust < 0.0001), "significance"] <- "****"
fish_df_all[which(fish_df_all$p_adjust < 0.001 & fish_df_all$p_adjust >= 0.0001), "significance"] <- "***"
fish_df_all[which(fish_df_all$p_adjust < 0.01 & fish_df_all$p_adjust >= 0.001), "significance"] <- "**"
fish_df_all[which(fish_df_all$p_adjust < 0.05 & fish_df_all$p_adjust >= 0.01), "significance"] <- "*"
fish_df_all[which(fish_df_all$p_adjust > 0.05), "significance"] <- ""

fish_df_all[fish_df_all$odds < 1, "enriched_in"] <- "non_smoker"
fish_df_all[fish_df_all$odds >= 1, "enriched_in"] <- "smoker"
fish_df_all[fish_df_all$significance != "", "is_significant"] <- "significant"
fish_df_all[fish_df_all$significance == "", "is_significant"] <- "not_significant"

fish_df_all[which(fish_df_all$odds > 1), "odds_ratio_plot"] <- fish_df_all[which(fish_df_all$odds > 1), "odds_ratio_flip"]
fish_df_all[which(fish_df_all$odds > 1), "lowCI_plot"] <- fish_df_all[which(fish_df_all$odds > 1), "loCI_flip"]
fish_df_all[which(fish_df_all$odds > 1), "hiCI_plot"] <- fish_df_all[which(fish_df_all$odds > 1), "hiCI_flip"]

fish_df_all[which(fish_df_all$odds < 1), "odds_ratio_plot"] <- 0-fish_df_all[which(fish_df_all$odds < 1), "odds"]
fish_df_all[which(fish_df_all$odds < 1), "lowCI_plot"] <- 0-fish_df_all[which(fish_df_all$odds < 1), "loCI"]
fish_df_all[which(fish_df_all$odds < 1), "hiCI_plot"] <- 0-fish_df_all[which(fish_df_all$odds < 1), "hiCI"]

fish_df_all$signature <- sub(":unknown", "", fish_df_all$signature)
fish_df_all <- fish_df_all[order(fish_df_all$significance,
                                 fish_df_all$odds_ratio_flip,
                                 decreasing = T), ]

fish_df_all$signature <- factor(fish_df_all$signature, levels = unique(fish_df_all$signature))

# plot
fish_df_all$signature <- sub("_", "", fish_df_all$signature)
fish_df_all$signature <- factor(fish_df_all$signature, levels = unique(fish_df_all$signature))

p_odd <- ggplot(fish_df_all, aes(x = signature, y = odds_ratio_flip, color = enriched_in, alpha = is_significant)) + 
  geom_pointrange(aes(ymin = loCI_flip, ymax = hiCI_flip)) +
  scale_color_manual(values = c("#7ad8fa","#46788a"), labels = c("inferred\nnon-smoker", "inferred\nsmoker"), name = "enriched in") +
  # scale_y_log10() +
  scale_alpha_manual(values = c(0.4, 1), labels = c("not significant", "significant"), name = "") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", linewidth = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),) +
  xlab("") +
  ylab("odds ratio") 

pdf(paste0(output_path, "SV_signature_smoker_enrichment.pdf"), width = 5, height = 3)
p_odd
dev.off()

# is SV N3 correlated with smoking signatures within LUSCS?
load(signatures_path)

all_signatures <- rbind(SBS_exposure, DBS_exposure, ID_exposure, SV_exposure)

sigs_interest <- c("SBS4", "SBS92", "ID3", "DBS2")

smoke_sig <- all_signatures[all_signatures$sig %in% sigs_interest, ]

SV_N3 <- SV_exposure[SV_exposure$signature == "SV_N3", c("sample", "weight", "exposure")]
colnames(SV_N3)[c(2,3)] <- c("SV_N3_weight", "SV_N3_exposure")

smoke_sig <- left_join(smoke_sig, SV_N3)
smoke_sig[is.na(smoke_sig$weight), "weight"] <- 0
smoke_sig[is.na(smoke_sig$exposure), "exposure"] <- 0
smoke_sig[is.na(smoke_sig$SV_N3_weight), "SV_N3_weight"] <- 0
smoke_sig[is.na(smoke_sig$SV_N3_exposure), "SV_N3_exposure"] <- 0

smoke_sig <- left_join(smoke_sig, sample_table[, c("participant_id", "histology")], by = c("sample" = "participant_id"))

p_cor_w <- ggscatter(smoke_sig[smoke_sig$weight > 0 & smoke_sig$SV_N3_weight > 0 & smoke_sig$label == "SBS92:smoking", ], x = "weight", y = "SV_N3_weight",
                     add = "reg.line",
                     conf.int = TRUE,
                     add.params = list(color = "blue",
                                       fill = "lightgray"))+
  stat_cor(method = "pearson") +
  xlab("% SBS92") +
  ylab("% SVN3")

pdf(paste0(output_path, "SVN3_SBS92_all_cor.pdf"), width = 3, height = 3)
p_cor_w
dev.off()
