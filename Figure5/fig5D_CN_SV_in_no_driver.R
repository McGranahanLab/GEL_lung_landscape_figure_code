# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #           CN and SV drivers in no SNV driver tumours                  # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#################################################################################
####################################################################### libraries

.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
.libPaths(c( .libPaths(), "/re_gecip/cancer_lung/R_packages_4_1/"))

library(ggpubr)
library(dplyr)
library(reshape2)

options(stringsAsFactors = F)
options(bitmapType = "cairo")


#################################################################################
#######################################################################file paths

driver_count_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/patient_driver_binary.txt"
SV_driver_mut_path      <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/SV/fishhook/50kb_bins_with_TRA/SV_driver_genes_in_samples.txt"
CN_driver_mut_path      <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/copynumber_peak_driver_per_patient.txt"

output_path             <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/CancerGenes/genomic_disruption/CN_SV_driver_in_samples_no_SNV_driver.pdf"

#################################################################################
#######################################################################     MAIN

driver_count <- read.table(driver_count_path, head = T, sep = "\t")
driver_count[driver_count$coding == TRUE | driver_count$non_coding == TRUE | driver_count$SVdriver == TRUE | driver_count$CNdriver == TRUE, "any_driver"] <- TRUE
driver_count[is.na(driver_count$any_driver), "any_driver"] <- FALSE

samples <- driver_count[driver_count$coding == FALSE & driver_count$non_coding == FALSE & driver_count$any_driver == TRUE, "participant_id"]

SV_driver <- read.table(SV_driver_mut_path, head = T, sep = "\t")
CN_driver <- read.table(CN_driver_mut_path, head = T, sep = "\t")

CN_driver_sample <- CN_driver[CN_driver$patient %in% samples, ]
SV_driver_sample <- SV_driver[SV_driver$participant_id %in% samples, ]

SVclass_count <- data.frame(table(SV_driver_sample$SV_CLASS))
SVclass_gene <- data.frame(table(SV_driver_sample$gene_name))
SVclass_element <- data.frame(table(SV_driver_sample$element_type))

SV_driver_sample <- unique(SV_driver_sample[,c("participant_id", "gene_name")])
SV_count <- data.frame(table(SV_driver_sample$gene_name))
SV_count$type <- "SV"

# for copy number need to be careful not to count the same peak that was found in adenos and in pancan twice
# maybe just use pancan peaks for this

CN_driver_sample <- CN_driver_sample[grep("PANCAN", CN_driver_sample$new_peak_name), ]
CN_driver_sample$name <- paste0(CN_driver_sample$new_peak_name, ":", CN_driver_sample$significant_CN_genes_on_peak)
CN_driver_sample$name <- sub("PANCAN_", "", CN_driver_sample$name)
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "TBL1XR1_PIK3CA_SOX2_MAP3K13_BCL6_TP63_MB21D2_MUC4_MECOM_STAG1_PIK3CB_ATR", "gene"] <- "PIK2CA, SOX2"
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "SDHA_TERT_CTNND2_CDH10_DROSHA_IL7R_NIPBL", "gene"] <- "TERT"
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "CSMD1_RPL23AP53_ZNF596_FAM87A_FBXO25_TDRP_ERICH1", "gene"] <- NA
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "SLC16A12_PANK1_FLJ37201_KIF20B_LINC00865_LINC01374", "gene"] <- NA
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "RECQL4_FAM135B_CDH17_PABPC1_UBR5_CSMD3_RAD21_EXT1_MYC", "gene"] <- "MYC"
CN_driver_sample[CN_driver_sample$significant_CN_genes_on_peak == "STK11_POLRMT_CD209_ZNF429_EEF2_CREB3L3_MAP2K2_KEAP1", "gene"] <- "STK11, KEAP1"
CN_driver_sample[grep("_", CN_driver_sample$significant_CN_genes_on_peak, invert = T), "gene"] <- CN_driver_sample[grep("_", CN_driver_sample$significant_CN_genes_on_peak, invert = T), "significant_CN_genes_on_peak"]
CN_driver_sample[nchar(CN_driver_sample$significant_CN_genes_on_peak) <= 20, "gene"] <- CN_driver_sample[nchar(CN_driver_sample$significant_CN_genes_on_peak) <= 20, "significant_CN_genes_on_peak"]
CN_driver_sample$short_peak_name <- sub("PANCAN_", "", CN_driver_sample$new_peak_name)
CN_driver_sample$short_peak_name <- sub("_p.*", "", CN_driver_sample$short_peak_name)
CN_driver_sample$short_peak_name <- sub("_q.*", "", CN_driver_sample$short_peak_name)

CN_driver_sample$short_name <- paste0(CN_driver_sample$short_peak_name, ":", CN_driver_sample$gene)
CN_driver_sample$short_name <- sub(":NA", "", CN_driver_sample$short_name)

CN_count <- data.frame(table(CN_driver_sample$short_name))
CN_count$type <- "CN"

count <- rbind(CN_count, SV_count)

# make as proportion
count$prop <- count$Freq / length(samples)
count$prop <- count$prop *100
count <- count[order(count$prop, decreasing = T), ]
count$Var1 <- factor(count$Var1, levels = unique(count$Var1))

p <- ggplot(count, aes(Var1, prop, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#756bb1", "#bcbddc")) +
  facet_grid(~type, scale = "free_x", space = "free_x") +
  xlab("") +
  ylab("% tumours without SNV driver") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_rect(fill="white"),
        legend.position = "none")

pdf(output_path, width = 6, height = 5)
p
dev.off()
