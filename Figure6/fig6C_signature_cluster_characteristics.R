# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                        signature cluster plot                         # # #    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
###################################################################   file paths

################################################################################
####################################################################   libraries

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
library(ComplexHeatmap)
library(circlize)
library(stringr)

options(stringsAsFactors = F)
options(bitmapType = "cairo")

################################################################################
####################################################################  file paths

signatures_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/HDP/HDP_primary_met_timepoint_patient_summary_object.RData"
sample_table_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/LungLandscape_sampleList_20221223.txt"
driver_path           <- "/re_gecip/cancer_lung/pipelines/Noncoding_drivers_Nextflow/completed_runs/2023_12_25/results/tables/driver_mutations/driverMutations-Panlung--hg19.csv"
TCRA_path             <- "/re_gecip/cancer_lung/shared/ImmuneLENS_pancan_data_to_use_20240202.RDS"
survival_path         <- "/re_gecip/cancer_lung/shared/pancan_survival_release_v16.RDS"
purity_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/CopyNumber/combined_ascatTable_segments_joined.rds"
smoking_path          <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/smoking_data.txt"
clone_table_path      <- '/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/clonal_subclonal_mutation_count.txt'
battenberg_path       <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/3.chromosomeAberrations/G.Battenberg/hg38_withSVs/combined_battenbergTable.rds"
ploidy_purity_path    <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/3.chromosomeAberrations/K.ASCAT/hg38_withSVs/combined_ascatTable.rds"
timing_path           <- "/re_gecip/cancer_lung/kthol/ComplexTiming/GRITIC/output_withSubclonalData_newBattenbergWithSV/"
WGD_path              <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/CopyNumber/patient_WGD.txt"
germline_path         <- "/re_gecip/cancer_lung/shared/germline_clinvar_variants.RDS"
HR_gene_path          <- "/re_gecip/cancer_lung/kthol/general/HR_genes.txt"
SV_path               <- "~/re_gecip/cancer_lung/kthol/SV/input/SV_list.Rdata"
clin_data_path        <- "/re_gecip/cancer_lung/analysisResults/GEL/release_v8/sampleTables/v14_lung_patient_clinical_data_2022-07-05_useful.tsv"
mut_burden_path       <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/patient_mutburden.RDS"
HLA_path              <- "~/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Immune/Lung_LOHHLA_20230425.RDS"
timed_mutation_path   <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/Timing/mutation_timing/timed_mutation_count.txt"
clonal_tumour_sig_path      <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/timed/clonal_tumour_sigs.txt"
clonal_all_sig_path         <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/timed/clonal_all_sigs.txt"
subclonal_tumour_sig_path   <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/timed/subclonal_tumour_sigs.txt"
subclonal_all_sig_path      <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/timed/subclonal_all_sigs.txt"
driver_count_path     <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/all_driver_count.txt"
SV_driver_mut_path    <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/SV/fishhook/50kb_bins_with_TRA/SV_driver_genes_in_samples.txt"
CN_driver_mut_path    <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/CancerGenes/CN/copynumber_peak_driver_per_patient.txt"
dnds_input_dir        <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/signature_clustering/results/"
CPI_path              <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Cohort/CPI_treated_patients_with_OS_list.csv"
wgii_path             <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/CopyNumber/wGII_scores.txt"
CPI_response_path     <- "/re_gecip/cancer_lung/bsimpson/CPI_cohort/Scripts/Analysis_version_lock_v13/Mastersheet_V31.csv"
chrom_mut_path        <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/MutBurden/chromatin_modifier_gene_mutburden.RData"
fusion_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Mutations/SV/fusions.txt"

output_path           <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/plots/Signatures/"
cluster_output_path   <- "/re_gecip/cancer_lung/analysisResults/landscape_paper/release_20230324/Signatures/"


num_clusters          <- 12

################################################################################
####################################################################        MAIN
signature_order <- c("SBS1:clock-like", "SBS2:APOBEC", "SBS3:HRd", "SBS4:smoking", "SBS5:clock-like", "SBS13:APOBEC", "SBS15:MMRd", "SBS17b:5FU-chemotherapy", "SBS18:ROS", "SBS21:MMRd", "SBS25:chemotherapy", "SBS44:MMRd", "SBS92:smoking",
                     "DBS2:smoking", "DBS4:unknown", "DBS11:unknown/APOBEC", "DBS_N1:unknown", "DBS_N3:unknown/platinum-chemo", "DBS_N4:unknown/MMRd",
                     "ID2:DNA-slippage", "ID3:smoking", "ID4:unknown", "ID6:HRd", "ID7:MMRd", "ID_N1:unknown", "ID_N2:unknown/NHEJ", "ID_N3:unknown",
                     "CN1:diploid", "CN2:tetraploid", "CN8:chromothripsis", "CN9:CIN", "CN17:HRd", "CN18:unknown/HRd", "CN20:unknown", "CN_N1:unknown", "CN_N2:unknown", "CN_N6:unknown", "CN_N8:unknown", "CN_N11:unknown", "CN_N12:unknown", "CN_N14:unknown", "CN_N15:unknown",
                     "SV2:unknown", "SV4:unknown", "SV6:unknown", "SV_N1:unknown", "SV_N2:unknown", "SV_N3:unknown", "SV_N4:unknown", "SV_N8:unknown", "SV_N9:unknown", "SV_N10:unknown", "SV_N11:unknown")

sample_table <- read.table(sample_table_path, head = T, sep = "\t")
sample_table$participant_id <- as.character(sample_table$participant_id)
sample_table[sample_table$histology == "MET_OTHER", "histology"] <- "OTHER"

# load signatures
load(signatures_path)

all_sigs <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure)
all_sigs$weight <- as.numeric(as.character(all_sigs$weight))
all_sigs$exposure <- as.numeric(as.character(all_sigs$exposure))

# make NAs into 0s
mat0 <- dcast(all_sigs, sample ~ label, value.var = "weight")
mat0 <- melt(mat0)
colnames(mat0) <- c("sample", "label", "weight")

mat0[is.na(mat0$weight), "weight"] <- 0
all_sigs <- mat0
all_sigs$label <- as.character(all_sigs$label)

mat_exposure <- rbind(SBS_exposure, DBS_exposure, ID_exposure, CN_exposure, SV_exposure)
mat_exposure$weight <- as.numeric(as.character(mat_exposure$weight))
mat_exposure$exposure <- as.numeric(as.character(mat_exposure$exposure))
mat_exposure <- dcast(mat_exposure, sample ~ label, value.var = "exposure")
mat_exposure <- melt(mat_exposure)
colnames(mat_exposure) <- c("sample", "label", "exposure")
mat_exposure[is.na(mat_exposure$weight), "exposure"] <- 0
mat_exposure$label <- as.character(mat_exposure$label)

# for each signature calculate the z scores

all_sigs$sig_type <- sapply(all_sigs$label, function(x){strsplit(x, ":")[[1]][1]})
all_sigs$sig_type <- sub("_N.", "", all_sigs$sig_type )
all_sigs$sig_type <- gsub("[[:digit:]]", "", all_sigs$sig_type)
all_sigs$sig_type <- gsub("b", "", all_sigs$sig_type)

all_sigs_z <- list()
for(sig in unique(all_sigs$label)){
  
  df <- all_sigs[all_sigs$label == sig, ]
  df$weight_z <- scale(df$weight)
  # df$exposure_z <- scale(log10(df$exposure +1))
  
  all_sigs_z[[sig]] <- df
  
}
all_sigs_z <- do.call(rbind, all_sigs_z)

# remove signatures present in less than 5% of tumours
signatures_to_keep <- c()
for(sig in unique(all_sigs_z$label)){
  
  df <- all_sigs_z[all_sigs_z$label == sig, ]
  num_sample_with_sig <- length(unique(df[df$weight >= 0.05, "sample"]))
  
  # in how many samples were we able to get signatures for this sig_type
  sig_type <- unique(df$sig_type)
  
  if(sig_type == "SBS"){
    num_sample_sig_type <- length(unique(SBS_exposure$sample))
  }
  
  if(sig_type == "DBS"){
    num_sample_sig_type <- length(unique(DBS_exposure$sample))
  }
  
  if(sig_type == "ID"){
    num_sample_sig_type <- length(unique(ID_exposure$sample))
  }
  
  if(sig_type == "CN"){
    num_sample_sig_type <- length(unique(CN_exposure$sample))
  }
  
  if(sig_type == "SV"){
    num_sample_sig_type <- length(unique(SV_exposure$sample))
  }
  
  # what is the porportion of samples that have this signature out of all of them 
  percentage <- (num_sample_with_sig/1011)*100
  if(percentage >= 3){
    signatures_to_keep[[sig]] <- sig
  }
  
}
signatures_to_keep <- as.character(signatures_to_keep)

all_sigs_a <- all_sigs_z[all_sigs_z$label %in% signatures_to_keep, ]

# z scaled weight
mat <- dcast(all_sigs_a, sample ~ label, value.var = "weight_z")
rownames(mat) <- mat$sample
mat$sample <- NULL

distance <- dist(mat, method = "euclidean")
cluster <- hclust(distance, method = "ward.D2")

cluster_assignment <- cutree(cluster, k = num_clusters)
order_cluster <- cluster$labels[cluster$order]

# now remake the tile plot splitting up the clusters
cluster_assignment <- data.frame(cluster_assignment)
cluster_assignment$variable <- rownames(cluster_assignment)

# now we can make a tiled plot with this prder
all_sigs_inf_clust_order <- mat[order_cluster, ] 
all_sigs_inf_clust_order$sample <- rownames(all_sigs_inf_clust_order)
all_sigs_inf_clust_order <- melt(all_sigs_inf_clust_order)
all_sigs_inf_clust_order$sample <- factor(all_sigs_inf_clust_order$sample, levels = order_cluster)
all_sigs_inf_clust_order$variable <- factor(all_sigs_inf_clust_order$variable, levels = signature_order)

# everything above 1 and below -1 gets made into 1.5
all_sigs_inf_clust_order[which(all_sigs_inf_clust_order$value < -1), "value"] <- -1.5
all_sigs_inf_clust_order[which(all_sigs_inf_clust_order$value > 1), "value"] <- 1.5

all_sigs_inf_clust_order <- left_join(all_sigs_inf_clust_order, cluster_assignment, by = c("sample" = "variable"))
all_sigs_inf_clust_order$cluster_assignment <- as.character(all_sigs_inf_clust_order$cluster_assignment)
all_sigs_inf_clust_order$cluster_assignment <- factor(all_sigs_inf_clust_order$cluster_assignment, unique(all_sigs_inf_clust_order$cluster_assignment))
all_sigs_inf_clust_order$sample <- factor(all_sigs_inf_clust_order$sample, levels = order_cluster)

# make a prettier plot
mat_t <- t(mat)
colnames(mat_t) <- NULL
mat_t <- as.matrix(mat_t)
colnames(mat_t) <- NULL

col_fun = colorRamp2(c(-1, 0, 1), c("#3f007d", "white", "#67000d"))

heat <- Heatmap(mat_t, clustering_method_columns = "ward.D2", clustering_distance_columns = "euclidean", col = col_fun, column_split = num_clusters, name = "signature\nweight")

dir.create(paste0(output_path, "clustering/ALL_new/"), recursive = T)

pdf(paste0(output_path, "clustering/ALL_new/signature_clustering_heatmap.pdf"), width = 12, height = 6)
heat
dev.off()

# the clusters here are the same but the order is different, so need to get the order from that
cluster_samples <- rownames(mat)
col_order <- column_order(heat)

cluster_name <- c(rep(1, length(col_order[[1]])), rep(2, length(col_order[[2]])),
                  rep(3, length(col_order[[3]])), rep(4, length(col_order[[4]])), 
                  rep(5, length(col_order[[5]])), rep(6, length(col_order[[6]])),
                  rep(7, length(col_order[[7]])), rep(8, length(col_order[[8]])),
                  rep(9, length(col_order[[9]])), rep(10, length(col_order[[10]])),
                  rep(11, length(col_order[[11]])), rep(12, length(col_order[[12]])))

complex_assignment <- data.frame(cluster_assignment = cluster_name,
                                 index = c(col_order[[1]], col_order[[2]], col_order[[3]], col_order[[4]], col_order[[5]], col_order[[6]], col_order[[7]], col_order[[8]], 
                                           col_order[[9]], col_order[[10]], col_order[[11]], col_order[[12]]))


cluster_samples <- data.frame(cluster_samples)
cluster_samples$index <- c(1:nrow(cluster_samples))

complex_assignment <- left_join(complex_assignment, cluster_samples)
complex_assignment$index <- NULL
colnames(complex_assignment) <- c("cluster_assignment", "variable")

cluster_assignment <- complex_assignment

# write.table(cluster_assignment, paste0(cluster_output_path, "new_signature_cluster_patient_assignment.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
cluster_assignment <- read.table(paste0(cluster_output_path, "new_signature_cluster_patient_assignment.txt"), head = T, sep = "\t")
################################################################################

# what are the predominant signatures in each group
cluster_sig_df <- mat
cluster_sig_df$sample <- rownames(cluster_sig_df)
cluster_sig_df$sample <- as.character(cluster_sig_df$sample)
cluster_assignment$variable <- as.character(cluster_assignment$variable)
cluster_sig_df <- left_join(cluster_sig_df, cluster_assignment, by = c("sample" = "variable"))

# for each cluster, calculate the median and sd of all signatures
sig_per_cluster <- list()
for(cluster in unique(cluster_sig_df$cluster_assignment)){
  df <- cluster_sig_df[cluster_sig_df$cluster_assignment == cluster, ]
  df$sample <- NULL
  df$cluster_assignment <- NULL
  
  all_median <- data.frame(apply(df, 2, median))
  colnames(all_median) <- "median"
  all_median$signature <- rownames(all_median)
  
  all_sd <- data.frame(apply(df, 2, sd))
  colnames(all_sd) <- "sd"
  all_sd$signature <- rownames(all_sd)
  
  all_out <- left_join(all_median, all_sd)
  all_out$cluster <- cluster
  
  sig_per_cluster[[cluster]] <- all_out
  
}

sig_per_cluster <- do.call(rbind, sig_per_cluster)

# # # # # # # # # # # # 
# test if there is an enrichment of any of the drivers in any of the groups

# driver mutations
driver_muts <- read.csv(driver_path, head = T, sep = "\t")
driver_muts <- driver_muts[driver_muts$var_class != "amp", ]
driver_muts <- driver_muts[driver_muts$var_class != "hom.del.", ]
driver_muts$gene_name <- paste0(driver_muts$gene_name, "_", driver_muts$gr_id)

# only keep those gene elements that occur in more than 10 samples
gene_count <- unique(driver_muts[, c("participant_id", "gene_name")])
gene_count <- data.frame(table(gene_count$gene_name))
gene_count <- gene_count[gene_count$Freq >= 10, ]

driver_muts <- driver_muts[driver_muts$gene_name %in% gene_count$Var1, ]

driver_muts <- unique(driver_muts[, c("participant_id", "gene_name")])
driver_muts$participant_id <- as.character(driver_muts$participant_id)

# add cluster 
driver_muts <- left_join(driver_muts, cluster_assignment, by = c("participant_id" = "variable"))

mutburden_df <- readRDS(mut_burden_path)
mutburden_df$SBS_DBS_ID_mutburden <- mutburden_df$SBScount + mutburden_df$DBScount + mutburden_df$IDcount

# for each of these genes do a fishers test of proportion of samples with and without mutation within or outside of the cluster
gene_test <- list()
for(gene in unique(driver_muts$gene_name)){
  
  print(gene)
  df <- driver_muts[driver_muts$gene_name == gene, ]
  
  cluster_test <- list()
  # iterate over the cluster
  for(i in 1:max(driver_muts$cluster_assignment)){
    cluster_df <- df[df$cluster_assignment == i, ]
    # how many of the samples with this gene mutation are in cluster of interest
    mutation_cluster <- nrow(cluster_df)
    if(length(mutation_cluster) == 0){mutation_cluster <- 0}
    
    # how many of the samples with this gene mutation are not in the cluster of interest
    no_cluster_df <- df[df$cluster_assignment != i, ]
    mutation_non_cluster <- nrow(no_cluster_df)
    if(length(mutation_non_cluster) == 0){mutation_non_cluster <- 0}
    
    # how many of all samples that don't have this mutation are in cluster of interest
    no_mutation_cluster_df <- cluster_assignment[-which(cluster_assignment$variable %in% df$participant_id), ]
    no_mutation_cluster <- nrow(no_mutation_cluster_df[no_mutation_cluster_df$cluster_assignment == i, ])
    if(length(no_mutation_cluster) == 0){no_mutation_cluster <- 0}
    
    # how many of all samples that don't have this mutation are not in the cluster
    no_mutation_no_cluster_df <- cluster_assignment[-which(cluster_assignment$variable %in% df$participant_id), ]
    no_mutation_no_cluster <- nrow(no_mutation_no_cluster_df[no_mutation_no_cluster_df$cluster_assignment != i, ])
    if(length(no_mutation_no_cluster) == 0){no_mutation_no_cluster <- 0}
    
    
    dat <- data.frame(cluster = c(mutation_cluster, no_mutation_cluster),
                      no_cluster = c(mutation_non_cluster, no_mutation_no_cluster))
    
    rownames(dat) <- c("mutation", "no_mutation")
    
    test <- fisher.test(dat)
    
    test_flip <- fisher.test(dat[, c(2, 1)])
    
    # also do a glm to adjust for mutation burden
    all_sample_info <- mutburden_df[, c("sample", "SBS_DBS_ID_mutburden")]
    all_sample_info[all_sample_info$sample %in% df$participant_id, "gene_mutated"] <- TRUE
    all_sample_info[is.na(all_sample_info$gene_mutated), "gene_mutated"] <- FALSE
    all_sample_info[all_sample_info$sample %in% cluster_assignment[cluster_assignment$cluster_assignment == i, "variable"], "in_cluster"] <- TRUE
    all_sample_info[is.na(all_sample_info$in_cluster), "in_cluster"] <- FALSE
    
    model <- glm(data = all_sample_info, as.factor(in_cluster) ~ as.factor(gene_mutated) + SBS_DBS_ID_mutburden,
                 family = binomial)
    
    
    data_out <- data.frame(gene_element = gene,
                           fisher_p_value = test$p.value,
                           glm_p_value = coef(summary(model))[,'Pr(>|z|)'][2],
                           odds_ratio = test$estimate,
                           odds_ratio_flip = test_flip$estimate,
                           cluster = mutation_cluster,
                           no_cluster = mutation_non_cluster,
                           cluster = i)
    
    cluster_test[[i]] <- data_out
    
  }
  cluster_test <- do.call(rbind, cluster_test)
  gene_test[[gene]] <- cluster_test
}

gene_test <- do.call(rbind, gene_test)

gene_test$p_adjust <- p.adjust(gene_test$glm_p_value, method = "fdr")

gene_test$significance <- NULL
gene_test[which(gene_test$p_adjust < 0.0001), "significance"] <- "****"
gene_test[which(gene_test$p_adjust < 0.001 & gene_test$p_adjust >= 0.0001), "significance"] <- "***"
gene_test[which(gene_test$p_adjust < 0.01 & gene_test$p_adjust >= 0.001), "significance"] <- "**"
gene_test[which(gene_test$p_adjust < 0.05 & gene_test$p_adjust >= 0.01), "significance"] <- "*"
gene_test[which(gene_test$p_adjust > 0.05), "significance"] <- ""

# stuff for plotting
gene_test$enrichment <- NULL
gene_test[gene_test$odds_ratio < 1, "enrichment"] <- "depleted"
gene_test[gene_test$odds_ratio >= 1, "enrichment"] <- "enriched"

gene_test[gene_test$odds_ratio == 0, "odds_ratio"] <- NA

# odds ratio plot
gene_test[which(gene_test$odds_ratio > 1), "odds_ratio_plot"] <- gene_test[which(gene_test$odds_ratio > 1), "odds_ratio_flip"]
gene_test[which(gene_test$odds_ratio < 1), "odds_ratio_plot"] <- 0-gene_test[which(gene_test$odds_ratio < 1), "odds_ratio"]

gene_test[gene_test$significance == "", "enrichment"] <- NA

gene_test$enrichment <- factor(gene_test$enrichment, levels = c("enriched", "depleted"))

# make a plot indicating the drivers
gene_test$cluster.1 <- factor(gene_test$cluster.1, levels = c(1:num_clusters))

plot_driver <- ggplot()+
  geom_tile(data = gene_test[gene_test$gene_element %in% unique(gene_test[gene_test$significance != "", "gene_element", ]), ], aes(cluster.1, gene_element, fill = odds_ratio_plot)) +
  geom_tile(data = gene_test[gene_test$gene_element %in% unique(gene_test[gene_test$significance != "", "gene_element", ]) & gene_test$significance == "", ], aes(cluster.1, gene_element), color = "transparent", fill = "white", linewidth = 1) +
  geom_tile(data = gene_test[gene_test$gene_element %in% unique(gene_test[gene_test$significance != "", "gene_element", ]), ], aes(cluster.1, gene_element), color = "transparent", fill = "transparent", linewidth = 1) +
  geom_tile(data = gene_test[gene_test$gene_element %in% unique(gene_test[gene_test$significance != "", "gene_element", ]), ], aes(cluster.1, gene_element, color = enrichment), fill = "transparent", linewidth = 1) +
  scale_color_manual(name = "significantly", values = c("#941107", "blue"), na.value = "transparent") +
  scale_fill_gradient2(name = "odds", low = "blue", mid = "white",
                       high = "#941107", midpoint = 0, na.value = "white") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, size = 8),
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "right") +
  ylab("SNV driver") +
  labs(x=NULL)

# reformat the dataframe so we can out it all together into one dataframe in the end
gene_test <- gene_test[, c("cluster.1", "gene_element", "odds_ratio", "odds_ratio_flip", "glm_p_value", "p_adjust", "cluster", "no_cluster")]
colnames(gene_test) <- c("cluster_assignment", "SNV_gene_element", "SNV_odds_ratio", "SNV_odds_ratio_flip", "SNV_glm_p_value", "SNV_p_adjust", "cluster_with_SNV", "cluster_without_SNV")


# test for enrrichment of SV drivers
SV_driver_muts <- read.table(SV_driver_mut_path, head = T, sep = "\t")
SV_driver_muts$gene_name <- paste0(SV_driver_muts$gene_name, "_", SV_driver_muts$element_type)

SV_driver_muts <- unique(SV_driver_muts[, c("participant_id", "gene_name")])
SV_driver_muts$participant_id <- as.character(SV_driver_muts$participant_id)

SV_driver_count <- data.frame(table(SV_driver_muts$gene_name))
SV_driver_count <- SV_driver_count[SV_driver_count$Freq >= 10, ]

SV_driver_muts <- SV_driver_muts[SV_driver_muts$gene_name %in% SV_driver_count$Var1, ]

# add cluster 
SV_driver_muts <- left_join(SV_driver_muts, cluster_assignment, by = c("participant_id" = "variable"))

# for each of these genes do a fishers test of proportion of samples with and without mutation within or outside of the cluster
SV_gene_test <- list()
for(gene in unique(SV_driver_muts$gene_name)){
  
  print(gene)
  df <- SV_driver_muts[SV_driver_muts$gene_name == gene, ]
  
  cluster_test <- list()
  # iterate over the cluster
  for(i in 1:max(SV_driver_muts$cluster_assignment)){
    cluster_df <- df[df$cluster_assignment == i, ]
    # how many of the samples with this gene mutation are in cluster of interest
    mutation_cluster <- nrow(cluster_df)
    if(length(mutation_cluster) == 0){mutation_cluster <- 0}
    
    # how many of the samples with this gene mutation are not in the cluster of interest
    no_cluster_df <- df[df$cluster_assignment != i, ]
    mutation_non_cluster <- nrow(no_cluster_df)
    if(length(mutation_non_cluster) == 0){mutation_non_cluster <- 0}
    
    # how many of all samples that don't have this mutation are in cluster of interest
    no_mutation_cluster_df <- cluster_assignment[-which(cluster_assignment$variable %in% df$participant_id), ]
    no_mutation_cluster <- nrow(no_mutation_cluster_df[no_mutation_cluster_df$cluster_assignment == i, ])
    if(length(no_mutation_cluster) == 0){no_mutation_cluster <- 0}
    
    # how many of all samples that don't have this mutation are not in the cluster
    no_mutation_no_cluster_df <- cluster_assignment[-which(cluster_assignment$variable %in% df$participant_id), ]
    no_mutation_no_cluster <- nrow(no_mutation_no_cluster_df[no_mutation_no_cluster_df$cluster_assignment != i, ])
    if(length(no_mutation_no_cluster) == 0){no_mutation_no_cluster <- 0}
    
    dat <- data.frame(cluster = c(mutation_cluster, no_mutation_cluster),
                      no_cluster = c(mutation_non_cluster, no_mutation_no_cluster))
    
    rownames(dat) <- c("mutation", "no_mutation")
    
    test <- fisher.test(dat)
    
    # also flip this around to get the odds ratio for the other way around
    test_flip <- fisher.test(dat[, c(2, 1)])
    
    
    # also do a glm to adjust for mutation burden
    all_sample_info <- mutburden_df[, c("sample", "SVcount")]
    all_sample_info[all_sample_info$sample %in% df$participant_id, "gene_mutated"] <- TRUE
    all_sample_info[is.na(all_sample_info$gene_mutated), "gene_mutated"] <- FALSE
    all_sample_info[all_sample_info$sample %in% cluster_assignment[cluster_assignment$cluster_assignment == i, "variable"], "in_cluster"] <- TRUE
    all_sample_info[is.na(all_sample_info$in_cluster), "in_cluster"] <- FALSE
    
    model <- glm(data = all_sample_info, as.factor(in_cluster) ~ as.factor(gene_mutated) + SVcount,
                 family = binomial)
    
    
    data_out <- data.frame(gene_element = gene,
                           fisher_p_value = test$p.value,
                           glm_p_value = coef(summary(model))[,'Pr(>|z|)'][2],
                           odds_ratio = test$estimate,
                           odds_ratio_flip = test_flip$estimate,
                           cluster = mutation_cluster,
                           no_cluster = mutation_non_cluster,
                           cluster = i)
    
    cluster_test[[i]] <- data_out
    
  }
  cluster_test <- do.call(rbind, cluster_test)
  SV_gene_test[[gene]] <- cluster_test
}

SV_gene_test <- do.call(rbind, SV_gene_test)

SV_gene_test$p_adjust <- p.adjust(SV_gene_test$glm_p_value, method = "fdr")

SV_gene_test$significance <- NULL
SV_gene_test[which(SV_gene_test$p_adjust < 0.0001), "significance"] <- "****"
SV_gene_test[which(SV_gene_test$p_adjust < 0.001 & SV_gene_test$p_adjust >= 0.0001), "significance"] <- "***"
SV_gene_test[which(SV_gene_test$p_adjust < 0.01 & SV_gene_test$p_adjust >= 0.001), "significance"] <- "**"
SV_gene_test[which(SV_gene_test$p_adjust < 0.05 & SV_gene_test$p_adjust >= 0.01), "significance"] <- "*"
SV_gene_test[which(SV_gene_test$p_adjust > 0.05), "significance"] <- ""

# plot SV driver
SV_gene_test$enrichment <- NULL
SV_gene_test[SV_gene_test$odds_ratio < 1, "enrichment"] <- "depleted"
SV_gene_test[SV_gene_test$odds_ratio >= 1, "enrichment"] <- "enriched"

SV_gene_test[SV_gene_test$odds_ratio == 0, "odds_ratio"] <- NA

# odds ratio plot
SV_gene_test[which(SV_gene_test$odds_ratio > 1), "odds_ratio_plot"] <- SV_gene_test[which(SV_gene_test$odds_ratio > 1), "odds_ratio_flip"]
SV_gene_test[which(SV_gene_test$odds_ratio < 1), "odds_ratio_plot"] <- 0-SV_gene_test[which(SV_gene_test$odds_ratio < 1), "odds_ratio"]

SV_gene_test[SV_gene_test$significance == "", "enrichment"] <- NA

SV_gene_test$enrichment <- factor(SV_gene_test$enrichment, levels = c("enriched", "depleted"))

SV_gene_test$cluster.1 <- factor(SV_gene_test$cluster.1, levels = c(1:num_clusters))

plot_SVdriver <- ggplot(SV_gene_test[SV_gene_test$gene_element %in% unique(SV_gene_test[SV_gene_test$significance != "", "gene_element", ]), ], aes(cluster.1, gene_element, fill = odds_ratio_plot,))+
  geom_tile() +
  geom_tile(data = SV_gene_test[SV_gene_test$gene_element %in% unique(SV_gene_test[SV_gene_test$significance != "", "gene_element", ]) & SV_gene_test$significance == "", ], color = "white", fill = "white", linewidth = 1) +
  geom_tile(data = SV_gene_test[SV_gene_test$gene_element %in% unique(SV_gene_test[SV_gene_test$significance != "", "gene_element", ]), ], color = "white", fill = NA, linewidth = 1) +
  geom_tile(data = SV_gene_test[!is.na(SV_gene_test$enrichment), ], aes(color = enrichment), fill = NA, linewidth = 1) +
  scale_fill_gradient2(name = "odds", low = "blue", mid = "white",
                       high = "#941107", midpoint = 0, na.value = "white") +
  scale_color_manual(name = "significantly", values = c("#941107", "blue")) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8),
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "right") +
  ylab("SV driver") +
  labs(x=NULL)

SV_gene_test <- SV_gene_test[, c("cluster.1", "gene_element", "odds_ratio", "odds_ratio_flip", "glm_p_value", "p_adjust", "cluster", "no_cluster")]
colnames(SV_gene_test) <- c("cluster_assignment", "SV_gene_element", "SV_odds_ratio", "SV_odds_ratio_flip", "SV_glm_p_value", "SV_p_adjust", "cluster_with_SV", "cluster_without_SV")

# and also copy number driver

CN_driver_muts <- read.table(CN_driver_mut_path, head = T, sep = "\t")
CN_driver_muts <- CN_driver_muts[grep("PANCAN", CN_driver_muts$new_peak_name), ]

# make as matrix
CN_driver_muts_mat <- dcast(CN_driver_muts, patient ~ new_peak_name, value.var = "event")

# add on all patient that are missing here
missing_patient <- sample_table[-which(sample_table$participant_id %in% CN_driver_muts_mat$patient), "participant_id"]

add_mat <- matrix(NA, nrow = length(missing_patient), ncol = ncol(CN_driver_muts_mat))
add_mat <- data.frame(add_mat)
colnames(add_mat) <- colnames(CN_driver_muts_mat)
add_mat$patient <- missing_patient

CN_driver_muts_mat <- rbind(CN_driver_muts_mat, add_mat)

CN_driver_muts_all <- melt(CN_driver_muts_mat, id.vars = "patient")
CN_driver_muts_long <- melt(CN_driver_muts_mat, id.vars = "patient")
colnames(CN_driver_muts_long) <- c("patient", "gene", "mutated")
CN_driver_muts_long <- CN_driver_muts_long[which(CN_driver_muts_long$mutated), ]

CN_driver_muts <- unique(CN_driver_muts_long[, c("patient", "gene")])
CN_driver_muts$patient <- as.character(CN_driver_muts$patient)
colnames(CN_driver_muts) <- c("participant_id", "gene_name")

CN_driver_count <- data.frame(table(CN_driver_muts$gene))
CN_driver_count <- CN_driver_count[CN_driver_count$Freq >= 10, ]

CN_driver_muts <- CN_driver_muts_all[CN_driver_muts_all$variable %in% CN_driver_count$Var1, ]

# add cluster 
CN_driver_muts <- left_join(CN_driver_muts, cluster_assignment, by = c("patient" = "variable"))
colnames(CN_driver_muts) <- c("participant_id", "gene_name", "mutation", "cluster_assignment")
CN_driver_muts[is.na(CN_driver_muts$mutation), "mutation"] <- FALSE

wgii <- read.table(wgii_path, head = T, sep = "\t")

# for each of these genes do a fishers test of proportion of samples with and without mutation within or outside of the cluste
CN_gene_test <- list()
for(gene in unique(CN_driver_muts$gene_name)){
  
  print(gene)
  df <- CN_driver_muts[CN_driver_muts$gene_name == gene, ]
  
  cluster_test <- list()
  # iterate over the cluster
  for(i in 1:max(CN_driver_muts$cluster_assignment)){
    cluster_df <- df[df$cluster_assignment == i, ]
    # how many of the samples with this gene mutation are in cluster of interest
    mutation_cluster <- nrow(cluster_df[cluster_df$mutation == TRUE, ])
    if(length(mutation_cluster) == 0){mutation_cluster <- 0}
    
    # how many of the samples with this gene mutation are not in the cluster of interest
    no_cluster_df <- df[df$cluster_assignment != i, ]
    mutation_non_cluster <- nrow(no_cluster_df[no_cluster_df$mutation == TRUE, ])
    if(length(mutation_non_cluster) == 0){mutation_non_cluster <- 0}
    
    # how many of all samples that don't have this mutation are in cluster of interest
    no_mutation_cluster_df <- df[df$cluster_assignment == i, ]
    no_mutation_cluster <- nrow(no_mutation_cluster_df[no_mutation_cluster_df$mutation == FALSE, ])
    if(length(no_mutation_cluster) == 0){no_mutation_cluster <- 0}
    
    # how many of all samples that don't have this mutation are not in the cluster
    no_mutation_no_cluster_df <- df[df$cluster_assignment != i, ]
    no_mutation_no_cluster <- nrow(no_mutation_no_cluster_df[no_mutation_no_cluster_df$mutation == FALSE, ])
    if(length(no_mutation_no_cluster) == 0){no_mutation_no_cluster <- 0}
    
    dat <- data.frame(cluster = c(mutation_cluster, no_mutation_cluster),
                      no_cluster = c(mutation_non_cluster, no_mutation_no_cluster))
    
    rownames(dat) <- c("mutation", "no_mutation")
    
    test <- fisher.test(dat)
    
    
    # also flip this around to get the odds ratio for the other way around
    test_flip <- fisher.test(dat[, c(2, 1)])
    
    # also do a glm to adjust for mutation burden
    all_sample_info <- wgii[, c("patient", "wGII")]
    all_sample_info[all_sample_info$patient %in% df[df$mutation == TRUE, "participant_id"], "gene_mutated"] <- TRUE
    all_sample_info[all_sample_info$patient %in% df[df$mutation == FALSE, "participant_id"], "gene_mutated"] <- FALSE
    all_sample_info[all_sample_info$patient %in% cluster_assignment[cluster_assignment$cluster_assignment == i, "variable"], "in_cluster"] <- TRUE
    all_sample_info[is.na(all_sample_info$in_cluster), "in_cluster"] <- FALSE
    
    model <- glm(data = all_sample_info, as.factor(in_cluster) ~ as.factor(gene_mutated) + wGII,
                 family = binomial)
    
    
    data_out <- data.frame(gene_element = gene,
                           fisher_p_value = test$p.value,
                           glm_p_value = coef(summary(model))[,'Pr(>|z|)'][2],
                           odds_ratio = test$estimate,
                           odds_ratio_flip = test_flip$estimate,
                           cluster = mutation_cluster,
                           no_cluster = mutation_non_cluster,
                           cluster = i)
    
    cluster_test[[i]] <- data_out
    
  }
  cluster_test <- do.call(rbind, cluster_test)
  CN_gene_test[[gene]] <- cluster_test
}

CN_gene_test <- do.call(rbind, CN_gene_test)

CN_gene_test$p_adjust <- p.adjust(CN_gene_test$glm_p_value, method = "fdr")

CN_gene_test$significance <- NULL
CN_gene_test[which(CN_gene_test$p_adjust < 0.0001), "significance"] <- "****"
CN_gene_test[which(CN_gene_test$p_adjust < 0.001 & CN_gene_test$p_adjust >= 0.0001), "significance"] <- "***"
CN_gene_test[which(CN_gene_test$p_adjust < 0.01 & CN_gene_test$p_adjust >= 0.001), "significance"] <- "**"
CN_gene_test[which(CN_gene_test$p_adjust < 0.05 & CN_gene_test$p_adjust >= 0.01), "significance"] <- "*"
CN_gene_test[which(CN_gene_test$p_adjust > 0.05), "significance"] <- ""

# plot CN driver
CN_gene_test$enrichment <- NULL
CN_gene_test[CN_gene_test$odds_ratio < 1, "enrichment"] <- "depleted"
CN_gene_test[CN_gene_test$odds_ratio >= 1, "enrichment"] <- "enriched"

CN_gene_test[CN_gene_test$odds_ratio == 0, "odds_ratio"] <- NA

# odds ratio plot
CN_gene_test[which(CN_gene_test$odds_ratio > 1), "odds_ratio_plot"] <- CN_gene_test[which(CN_gene_test$odds_ratio > 1), "odds_ratio_flip"]
CN_gene_test[which(CN_gene_test$odds_ratio < 1), "odds_ratio_plot"] <- 0-CN_gene_test[which(CN_gene_test$odds_ratio < 1), "odds_ratio"]

CN_gene_test[CN_gene_test$significance == "", "enrichment"] <- NA

CN_gene_test$enrichment <- factor(CN_gene_test$enrichment, levels = c("enriched", "depleted"))

CN_gene_test$gene_element <- sub("PANCAN_", "", CN_gene_test$gene_element)

CN_gene_test$cluster.1 <- factor(CN_gene_test$cluster.1, levels = c(1:num_clusters))

plot_CNdriver <- ggplot(CN_gene_test[CN_gene_test$gene_element %in% unique(CN_gene_test[CN_gene_test$significance != "", "gene_element", ]), ], aes(cluster.1, gene_element, fill = odds_ratio_plot))+
  geom_tile() +
  geom_tile(data = CN_gene_test[CN_gene_test$gene_element %in% unique(CN_gene_test[CN_gene_test$significance != "", "gene_element", ]) & CN_gene_test$significance == "", ], color = "white", fill = "white", linewidth = 1) +
  geom_tile(data = CN_gene_test[CN_gene_test$gene_element %in% unique(CN_gene_test[CN_gene_test$significance != "", "gene_element", ]), ], color = "white", fill = NA, linewidth = 1) +
  geom_tile(data = CN_gene_test[!is.na(CN_gene_test$enrichment), ], aes(color = enrichment), fill = NA, linewidth = 1) +
  scale_fill_gradient2(name = "odds", low = "blue", mid = "white",
                       high = "#941107", midpoint = 0, na.value = "white") +
  scale_color_manual(name = "significantly", values = c("#941107", "blue")) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5, size = 8),
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "right") +
  ylab("CN driver") +
  labs(x=NULL)

CN_gene_test <- CN_gene_test[, c("cluster.1", "gene_element", "odds_ratio", "odds_ratio_flip", "glm_p_value", "p_adjust", "cluster", "no_cluster")]
colnames(CN_gene_test) <- c("cluster_assignment", "CN_gene_element", "CN_odds_ratio", "CN_odds_ratio_flip", "CN_glm_p_value", "CN_p_adjust", "cluster_with_CN", "cluster_without_CN")

# also check for enrichment of germline mutations
HR_genes <- read.table(HR_gene_path, head = F, sep = "\t")
HR_genes <- c(HR_genes$V1, "CHEK2")

germline_all <- readRDS(germline_path)
germline_all <- germline_all[germline_all$gene %in% HR_genes, ]
germline_all <- germline_all[germline_all$patient %in% sample_table$participant_id, ]

# remove BRCA mutation thta's not trustworthy
germline_all <- germline_all[grep("rs80359307", germline_all$ID, invert = T), ]

germline_all <- germline_all[, c("patient", "gene")]
germline_all <- as.data.frame(germline_all)
germline_all <- unique(germline_all)

germline_muts <- left_join(germline_all, cluster_assignment, by = c("patient" = "variable"))

germ_gene_test <- list()

for(gene in "BRCA2"){
  
  print(gene)
  df <- germline_muts[germline_muts$gene == gene, ]
  
  cluster_test <- list()
  # iterate over the cluster
  for(i in 1:max(germline_muts$cluster_assignment)){
    cluster_df <- df[which(df$cluster_assignment == i), ]
    # how many of the samples with this gene mutation are in cluster of interest
    mutation_cluster <- nrow(cluster_df)
    if(length(mutation_cluster) == 0){mutation_cluster <- 0}
    
    # how many of the samples with this gene mutation are not in the cluster of interest
    no_cluster_df <- df[which(df$cluster_assignment != i), ]
    mutation_non_cluster <- nrow(no_cluster_df)
    if(length(mutation_non_cluster) == 0){mutation_non_cluster <- 0}
    
    # how many of all samples that don't have this mutation are in cluster of interest
    no_mutation_cluster_df <- cluster_assignment[-which(cluster_assignment$variable %in% df$patient), ]
    no_mutation_cluster <- nrow(no_mutation_cluster_df[no_mutation_cluster_df$cluster_assignment == i, ])
    if(length(no_mutation_cluster) == 0){no_mutation_cluster <- 0}
    
    # how many of all samples that don't have this mutation are not in the cluster
    no_mutation_no_cluster_df <- cluster_assignment[-which(cluster_assignment$variable %in% df$patient), ]
    no_mutation_no_cluster <- nrow(no_mutation_no_cluster_df[no_mutation_no_cluster_df$cluster_assignment != i, ])
    if(length(no_mutation_no_cluster) == 0){no_mutation_no_cluster <- 0}
    
    dat <- data.frame(cluster = c(mutation_cluster, no_mutation_cluster),
                      no_cluster = c(mutation_non_cluster, no_mutation_no_cluster))
    
    rownames(dat) <- c("mutation", "no_mutation")
    
    test <- fisher.test(dat)
    
    # also flip this around to get the odds ratio for the other way around
    test_flip <- fisher.test(dat[, c(2, 1)])
    
    
    data_out <- data.frame(gene_element = gene,
                           p_value = test$p.value,
                           odds_ratio = test$estimate,
                           odds_ratio_flip = test_flip$estimate,
                           cluster = mutation_cluster,
                           no_cluster = mutation_non_cluster,
                           cluster = i)
    
    cluster_test[[i]] <- data_out
    
  }
  cluster_test <- do.call(rbind, cluster_test)
  germ_gene_test[[gene]] <- cluster_test
}

germ_gene_test <- do.call(rbind, germ_gene_test)

germ_gene_test$p_adjust <- p.adjust(germ_gene_test$p_value, method = "fdr")

# plot this
germ_gene_test$significance <- NULL
germ_gene_test[which(germ_gene_test$p_adjust < 0.0001), "significance"] <- "****"
germ_gene_test[which(germ_gene_test$p_adjust < 0.001 & germ_gene_test$p_adjust >= 0.0001), "significance"] <- "***"
germ_gene_test[which(germ_gene_test$p_adjust < 0.01 & germ_gene_test$p_adjust >= 0.001), "significance"] <- "**"
germ_gene_test[which(germ_gene_test$p_adjust < 0.05 & germ_gene_test$p_adjust >= 0.01), "significance"] <- "*"
germ_gene_test[which(germ_gene_test$p_adjust > 0.05), "significance"] <- ""

germ_gene_test <- germ_gene_test[order(germ_gene_test$p_adjust, decreasing = F), ]
germ_gene_test_p <- germ_gene_test[, c("gene_element", "cluster", "no_cluster", "significance", "cluster.1")]
germ_gene_test_p <- melt(germ_gene_test_p, id.vars = c("gene_element", "significance", "cluster.1"))
germ_gene_test_p$gene_element <- factor(germ_gene_test_p$gene_element, levels = unique(germ_gene_test_p$gene_element))
germ_gene_test_p$value_log <- log10(germ_gene_test_p$value)

# only plot significant ones
germ_gene_test_p <- germ_gene_test_p[germ_gene_test_p$significance != "", ]
germ_gene_test_p$variable_cluster <- paste0(germ_gene_test_p$variable, "_", germ_gene_test_p$cluster.1)

germ_gene_test$enrichment <- NULL
germ_gene_test[germ_gene_test$odds_ratio < 1, "enrichment"] <- "depleted"
germ_gene_test[germ_gene_test$odds_ratio >= 1, "enrichment"] <- "enriched"

germ_gene_test[germ_gene_test$odds_ratio == 0, "odds_ratio"] <- NA

# odds ratio plot
germ_gene_test[which(germ_gene_test$odds_ratio > 1), "odds_ratio_plot"] <- germ_gene_test[which(germ_gene_test$odds_ratio > 1), "odds_ratio_flip"]
germ_gene_test[which(germ_gene_test$odds_ratio < 1), "odds_ratio_plot"] <- 0-germ_gene_test[which(germ_gene_test$odds_ratio < 1), "odds_ratio"]

germ_gene_test[germ_gene_test$significance == "", "enrichment"] <- NA

germ_gene_test$enrichment <- factor(germ_gene_test$enrichment, levels = c("enriched", "depleted"))

germ_gene_test$cluster.1 <- factor(germ_gene_test$cluster.1, levels = c(1:num_clusters))

plot_germ_driver <- ggplot(germ_gene_test[germ_gene_test$gene_element %in% unique(germ_gene_test[germ_gene_test$significance != "", "gene_element", ]), ], aes(cluster.1, gene_element, fill = odds_ratio))+
  geom_tile() +
  geom_tile(data = germ_gene_test[germ_gene_test$gene_element %in% unique(germ_gene_test[germ_gene_test$significance != "", "gene_element", ]) & germ_gene_test$significance == "", ], color = "white", fill = "white", linewidth = 1) +
  geom_tile(data = germ_gene_test[germ_gene_test$gene_element %in% unique(germ_gene_test[germ_gene_test$significance != "", "gene_element", ]), ], color = "white", fill = NA, linewidth = 1) +
  geom_tile(data = germ_gene_test[!is.na(germ_gene_test$enrichment), ], aes(color = enrichment), fill = NA, linewidth = 1) +
  scale_fill_gradient2(name = "odds", low = "blue", mid = "white",
                       high = "#941107", midpoint = 0, na.value = "white") +
  scale_color_manual(name = "significantly", values = c("#941107", "blue")) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8),
        plot.margin = unit(c(0.1,0,0.1,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "right") +
  ylab("germline driver") +
  labs(x=NULL)

germ_gene_test <- germ_gene_test[, c("cluster.1", "gene_element", "odds_ratio", "odds_ratio_flip", "p_value", "p_adjust", "cluster", "no_cluster")]
colnames(germ_gene_test) <- c("cluster_assignment", "germline_gene_element", "germline_odds_ratio", "germline_odds_ratio_flip", "germline_fisher_p_value", "germline_p_adjust", "cluster_with_germline_var", "cluster_without_germline_var")


############
# TCRA in different clusters
TCRA <- readRDS(TCRA_path)
TCRA <- TCRA[TCRA$type == "tumour", ]

immune <- left_join(cluster_assignment, TCRA[, c("patient", "TCRA.tcell.fraction.adj","IGH.bcell.fraction.adj", "igMD.frac.adj")], by = c("variable" = "patient"))
immune$class_switchedB <- immune$IGH.bcell.fraction.adj - immune$igMD.frac.adj
immune$igMD.frac.adj <- NULL

immune <- melt(immune, id.vars = c("cluster_assignment", "variable"))
colnames(immune) <- c("cluster_assignment", "variable", "immune_score", "value")
immune$cluster_assignment <- as.character(immune$cluster_assignment)
immune$cluster_assignment <- factor(immune$cluster_assignment,  levels = c(1:num_clusters))
tcra <- immune[immune$immune_score == "TCRA.tcell.fraction.adj", ]

ggplot(tcra, aes(cluster_assignment, value)) +
  geom_boxplot() +
  stat_compare_means()

tcra <- tcra %>% group_by(cluster_assignment) %>% mutate(mean = median(value, na.rm = T))
tcra <- unique(tcra[,c("cluster_assignment", "mean")])
tcra$cluster_assignment <- factor(tcra$cluster_assignment, levels = c(1:num_clusters))

plot_immunet <- ggplot(tcra, aes(cluster_assignment, 1, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(name = "T-cell fraction", low = "#e5f5e0", high = "#00441b") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("T-cell fraction") +
  labs(x=NULL)

colnames(tcra)[2] <- "mean_TCRAscore"

tcrb <- immune[immune$immune_score == "IGH.bcell.fraction.adj", ]
tcrb <- tcrb %>% group_by(cluster_assignment) %>% mutate(mean = median(value, na.rm = T))
tcrb <- unique(tcrb[,c("cluster_assignment", "mean")])
tcrb$cluster_assignment <- factor(tcrb$cluster_assignment, levels = c(1:num_clusters))

plot_immuneb <- ggplot(tcrb, aes(cluster_assignment, 1, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(name = "B-cell fraction", low = "#e5f5e0", high = "#00441b") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("naive B-cell fraction") +
  labs(x=NULL)


colnames(tcrb)[2] <- "mean_TCRBscore"

#############
# check histology

histology <- left_join(cluster_assignment, sample_table[, c("participant_id", "histology")], by = c("variable" = "participant_id"))

# how big is each cluster
clust_size <- data.frame(table(histology$cluster_assignment))
colnames(clust_size) <- c("cluster_assignment", "cluster_size")

hist_count <- data.frame(table(histology[, c("cluster_assignment", "histology")]))
hist_count <- left_join(hist_count, clust_size)
hist_count$proportion <- hist_count$Freq/hist_count$cluster_size 

hist_count$cluster_assignment <- as.character(hist_count$cluster_assignment)
hist_count$cluster_assignment <- factor(hist_count$cluster_assignment,  levels = c(1:num_clusters))

plot_hist <- ggplot(hist_count, aes(cluster_assignment, proportion, fill = histology))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("ADENOCARCINOMA" = "#67001f",
                               "MET_ADENOCARCINOMA" = "#67001f80",
                               "SQUAMOUS_CELL" = "#053061",
                               "MET_SQUAMOUS_CELL" = "#05306180",
                               "ADENOSQUAMOUS" = "#66c2a5",
                               "CARCINOID" = "#fc8d62",
                               "LARGE_CELL" = "#54278f",
                               "MESOTHELIOMA" = "#e3309e",
                               "SMALL_CELL" = "#b89704",
                               "MET_SMALL_CELL" = "#b8970480",
                               "NEUROENDOCRINE_CARCINOMA" = "#0b8ca3",
                               "OTHER" = "#7a7979")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("% tumours") +
  labs(x=NULL)

hist_count <- hist_count[, c('cluster_assignment', 'histology', 'proportion')]
colnames(hist_count)[3] <- "histology_proportion"

clust_size$cluster_assignment <- factor(clust_size$cluster_assignment, levels = c(1:num_clusters))
plot_size <- ggplot(clust_size, aes(cluster_assignment, cluster_size, fill = cluster_assignment))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#a6cee3",
                                        "#1f78b4",
                                        "#b2df8a",
                                        "#33a02c",
                                        "#fb9a99",
                                        "#e31a1c",
                                        "#fdbf6f",
                                        "#ff7f00",
                                        "#cab2d6",
                                        "#6a3d9a",
                                        "#b15928",
                                        "#ffff99",
                                        "#c51b8a")) +
                                          scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0,0.1,0), "cm"),
        legend.justification = c(0,1),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8),
        panel.spacing = unit(0,'lines'),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("# tumours") +
  xlab("cluster assignment") 

# ratio of primary to mets
time_point <- left_join(cluster_assignment, sample_table[, c("participant_id", "time_point")], by = c("variable" = "participant_id"))
time_point <- data.frame(table(time_point$cluster_assignment, time_point$time_point))
colnames(time_point) <- c("cluster_assignment", "time_point", "Freq")
time_point <- left_join(time_point, clust_size)

time_point$prop  <- time_point$Freq / time_point$cluster_size
time_point$time_point <- factor(time_point$time_point, levels = c("PRIMARY", "non_PRIMARY"))

time_point$cluster_assignment <- factor(time_point$cluster_assignment,  levels = c(1:num_clusters))

plot_timepoint <- ggplot(time_point, aes(cluster_assignment, prop, fill = time_point))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#bcbddc","#4a1486"), breaks = c("PRIMARY", "non_PRIMARY"), label = c("primary", "recurrence/\nmetastasis"), name = "time point") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("% tumours") +
  labs(x=NULL)

time_point <- time_point[, c("cluster_assignment", "time_point", "prop")]
time_point <- dcast(time_point, cluster_assignment ~ time_point)
colnames(time_point) <- c("cluster_assignment", "proportion_primary", "proportion_not_primary")


# and check smoking status
smoking <- read.table(smoking_path, head = T, sep = "\t")
smoking$participant_id <- as.character(smoking$participant_id)

smoking <- left_join(smoking[, c("participant_id", "classifier")], cluster_assignment, by = c("participant_id" = "variable"))

# how big is each cluster
clust_size <- data.frame(table(smoking$cluster_assignment))
colnames(clust_size) <- c("cluster_assignment", "cluster_size")

smok_count <- data.frame(table(smoking[, c("cluster_assignment", "classifier")]))
smok_count <- left_join(smok_count, clust_size)
smok_count$proportion <- smok_count$Freq/smok_count$cluster_size 

smok_count$cluster_assignment <- as.character(smok_count$cluster_assignment)
smok_count$cluster_assignment <- factor(smok_count$cluster_assignment,  levels = c(1:num_clusters))

plot_smok <- ggplot(smok_count, aes(cluster_assignment, proportion, fill = classifier))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#7ad8fa","#46788a"), breaks = c("non_smoker", "smoker"), label = c("inferred never-smoker", "inferred smoker"), name = "inferred smoking status") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("% tumours") +
  labs(x=NULL)

smok_count <- smok_count[, c("cluster_assignment", "classifier", "proportion")]
smok_count <- dcast(smok_count, cluster_assignment ~ classifier)
colnames(smok_count) <- c("cluster_assignment", "proportion_inferred_never_smoker", "proportion_inferred_smoker")

# also add heterogeneity measures
clone_table <- read.table(clone_table_path, head = T, sep = "\t")
clone_table$subclonal_clonal_fraction <- clone_table$number_subclonal_mutations/clone_table$number_clonal_mutationes
clone_table$sample <- as.character(clone_table$sample)

clone_cluster_assignment <- left_join(cluster_assignment, clone_table, by = c("variable" = "sample"))

clone_cluster_assignment$cluster_assignment <- as.character(clone_cluster_assignment$cluster_assignment)

clone_cluster_assignment$cluster_assignment <- factor(clone_cluster_assignment$cluster_assignment,  levels = c(1:num_clusters))

SNV_ITH <- clone_cluster_assignment[clone_cluster_assignment$number_subclonal_mutations > 0, c("cluster_assignment", "subclonal_clonal_fraction")]
SNV_ITH <- SNV_ITH %>% group_by(cluster_assignment) %>% mutate(mean = median(subclonal_clonal_fraction, na.rm = T))
SNV_ITH <- unique(SNV_ITH[,c("cluster_assignment", "mean")])
SNV_ITH$cluster_assignment <- factor(SNV_ITH$cluster_assignment, levels = c(1:num_clusters))

plot_SNV_ITH <- ggplot(SNV_ITH, aes(cluster_assignment, 1, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(name = "SNV ITH", low = "#fde0dd", high = "#ae017e") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("SNV ITH") +
  labs(x=NULL)

colnames(SNV_ITH) <- c("cluster_assignment", "median_SNV_ITH")

# make this back into a supplementary plot for supplementary
plot_SNV_ITH_b <- ggplot(clone_cluster_assignment[clone_cluster_assignment$number_subclonal_mutations > 0, ], aes(cluster_assignment, subclonal_clonal_fraction, fill = cluster_assignment))+
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(alpha = 0.8) +
  stat_compare_means(label.y = log10(30), label.x = 2) +
  scale_y_log10(expand = c(0,0)) +
  scale_fill_manual(values = c("#a6cee3",
                                        "#1f78b4",
                                        "#b2df8a",
                                        "#33a02c",
                                        "#fb9a99",
                                        "#e31a1c",
                                        "#fdbf6f",
                                        "#ff7f00",
                                        "#cab2d6",
                                        "#6a3d9a",
                                        "#b15928",
                                        "#ffff99",
                                        "#c51b8a")) +
                                          theme(legend.position = "none") +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("# SNV ITH") +
  labs(x=NULL)


plot_subclonal_SNV <- ggplot(clone_cluster_assignment[clone_cluster_assignment$number_subclonal_mutations > 0, ], aes(cluster_assignment, number_subclonal_mutations, fill = cluster_assignment))+
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(alpha = 0.8) +
  stat_compare_means(label.y = log10(15000), label.x = 2) +
  scale_y_log10(expand = c(0,0)) +
  scale_fill_manual(values = c("#a6cee3",
                                        "#1f78b4",
                                        "#b2df8a",
                                        "#33a02c",
                                        "#fb9a99",
                                        "#e31a1c",
                                        "#fdbf6f",
                                        "#ff7f00",
                                        "#cab2d6",
                                        "#6a3d9a",
                                        "#b15928",
                                        "#ffff99",
                                        "#c51b8a")) +
                                          theme(legend.position = "none") +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("# subclonal SNVs") +
  labs(x=NULL)


plot_clonal_SNV <- ggplot(clone_cluster_assignment[clone_cluster_assignment$number_subclonal_mutations > 0, ], aes(cluster_assignment, number_clonal_mutationes, fill = cluster_assignment))+
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(alpha = 0.8) +
  stat_compare_means(label.y = log10(55000), label.x = 2) +
  scale_y_log10(expand = c(0,0)) +
  scale_fill_manual(values = c("#a6cee3",
                                        "#1f78b4",
                                        "#b2df8a",
                                        "#33a02c",
                                        "#fb9a99",
                                        "#e31a1c",
                                        "#fdbf6f",
                                        "#ff7f00",
                                        "#cab2d6",
                                        "#6a3d9a",
                                        "#b15928",
                                        "#ffff99",
                                        "#c51b8a")) +
                                          theme(legend.position = "none") +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("# clonal SNVs") +
  labs(x=NULL)



p_ITH <- plot_grid(plot_SNV_ITH_b, plot_subclonal_SNV, plot_clonal_SNV, ncol = 1, align = "v")

pdf(paste0(output_path, "clustering/ALL_new/signature_cluster_SNV_ITH.pdf"), width = 6, height = 4)
p_ITH
dev.off()


WGD_df <- read.table(WGD_path, head = T, sep = "\t")
WGD_df <- WGD_df[, c("patient", "first_gd")]
WGD_df$patient <- as.character(WGD_df$patient)
WGD_df <- left_join(WGD_df, cluster_assignment, by = c("patient" = "variable"))

WGD_df <- WGD_cluster[, c("cluster_assignment", "WGD_Status")]
WGD_df <- data.frame(table(WGD_df[, c("first_gd", "cluster_assignment")]))
WGD_df <- WGD_df[WGD_df$first_gd == TRUE, ]
WGD_df <- left_join(WGD_df, clust_size)
WGD_df$WGD_prop <- WGD_df$Freq/WGD_df$cluster_size
WGD_df$WGD_prop <- WGD_df$WGD_prop *100

WGD_df$cluster_assignment <- factor(WGD_df$cluster_assignment, levels = c(1:num_clusters))

plot_WGD <- ggplot(WGD_df, aes(cluster_assignment, 1, fill = WGD_prop)) +
  geom_tile() +
  scale_fill_gradient(name = "% WGD", low = "#fee0d2", high = "#a50f15") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("% WGD") +
  labs(x=NULL)

WGD_df <- WGD_df[, c("cluster_assignment", "WGD_prop")]
colnames(WGD_df) <- c("cluster_assignment", "WGD_percent")

age_df <- sample_table[, c("participant_id", "patient_age")]
age_df <- left_join(cluster_assignment, age_df, by = c("variable" = "participant_id"))
age_df$cluster_assignment <- as.character(age_df$cluster_assignment)

# ggplot(age_df,  aes(cluster_assignment, patient_age)) +
#   geom_boxplot()+
#   stat_compare_means()

age_df <- age_df %>% group_by(cluster_assignment) %>% mutate(mean = median(patient_age, na.rm = T))
age_df <- unique(age_df[,c("cluster_assignment", "mean")])
age_df$cluster_assignment <- factor(age_df$cluster_assignment, levels = c(1:num_clusters))

plot_age <- ggplot(age_df, aes(cluster_assignment, 1, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(name = "age", low = "#bfd3e6", high = "#88419d", limits = c(50, 85)) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("age") +
  labs(x=NULL)

colnames(age_df)[2] <- "median_age"

# 
# add sex
sex_df <- sample_table[, c("participant_id", "sex")]
sex_df <- left_join(cluster_assignment, sex_df, by = c("variable" = "participant_id"))
sex_df <- data.frame(table(sex_df[, c("cluster_assignment", "sex")]))
sex_df <- left_join(sex_df, clust_size)
sex_df$prop <- sex_df$Freq/sex_df$cluster_size

sex_df$cluster_assignment <- factor(sex_df$cluster_assignment, levels = c(1:num_clusters))

plot_sex <- ggplot(sex_df, aes(cluster_assignment, prop, fill = sex)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("black", "white"), breaks = c("MALE", "FEMALE"), label = c("male", "female"), name = "sex") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        legend.position = "right") +
  ylab("% samples") +
  labs(x=NULL)

sex_df <- sex_df[, c("cluster_assignment", "sex", "prop")]
sex_df <- dcast(sex_df, cluster_assignment ~ sex)
colnames(sex_df) <- c("cluster_assignment", "proportion_female", "proportion_male")

# mutburden
mutburden_df <- readRDS(mut_burden_path)
mutburden_df <- left_join(cluster_assignment, mutburden_df, by = c("variable" = "sample"))
mutburden_df$cluster_assignment <- as.character(mutburden_df$cluster_assignment)

# make the mutation burdens into mutation per megabase, assuming that an average whole genome has 2800 Mb with sufficient coverage
mutburden_df$SBScount  <- mutburden_df$SBScount/2800
mutburden_df$IDcount   <- mutburden_df$IDcount/2800
mutburden_df$DBScount  <- mutburden_df$DBScount/2800
mutburden_df$SVcount   <- mutburden_df$SVcount/2800

mutburden_df[is.na(mutburden_df$SBScount), "SBScount"] <- 0
mutburden_df[is.na(mutburden_df$DBScount), "DBScount"] <- 0
mutburden_df[is.na(mutburden_df$IDcount), "IDcount"] <- 0
mutburden_df[is.na(mutburden_df$SVcount), "SVcount"] <- 0

mutburden_df <- mutburden_df %>% group_by(cluster_assignment) %>%
  mutate(mean_SBS = median(SBScount)) %>%
  mutate(mean_ID = median(IDcount)) %>%
  mutate(mean_DBS = median(DBScount)) %>%
  mutate(mean_SV = median(SVcount))

mutburden_df <- unique(mutburden_df[, c("cluster_assignment", "mean_SBS", "mean_DBS", "mean_ID", "mean_SV"), ])

mutburden_df$cluster_assignment <- factor(mutburden_df$cluster_assignment, levels = c(1:num_clusters))

plot_SBSburden <- ggplot(mutburden_df, aes(cluster_assignment, 1, fill = mean_SBS)) +
  geom_tile() +
  scale_fill_gradient(name = "SBS/megabase", low = "#fde0dd", high = "#ae017e") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("SBS/megabase") +
  labs(x=NULL)

plot_SVburden <- ggplot(mutburden_df, aes(cluster_assignment, 1, fill = mean_SV)) +
  geom_tile() +
  scale_fill_gradient(name = "SV/megabase", low = "#fde0dd", high = "#ae017e") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("SV/megabase") +
  labs(x=NULL)

# add stage
clin <- read.table(clin_data_path, head = T, sep = "\t")
clin$participant_id <- as.character(clin$participant_id)

stage <- clin[, c("participant_id", "stage_stage_best")]

stage <- left_join(cluster_assignment, stage, by = c("variable" = "participant_id"))
stage$cluster_assignment <- as.character(stage$cluster_assignment)
stage$stage <- stage$stage_stage_best
stage[which(stage$stage == ""), "stage"] <- NA
stage[which(stage$stage == "?"), "stage"] <- NA
stage[which(stage$stage == "U"), "stage"] <- NA
stage[which(stage$stage %in% c("1A", "1A_1", "1A3", "1A_3B_", "1A1", "1A2", "1A2", "1B")), "stage"] <- 1
stage[which(stage$stage %in% c("2", "2A", "2B")), "stage"] <- 2
stage[which(stage$stage %in% c("3", "3A", "3B", "3C")), "stage"] <- 3
stage[which(stage$stage %in% c("4", "4A", "4B")), "stage"] <- 4
stage$stage <- as.character(stage$stage)
stage$stage <- factor(stage$stage, levels = c("1", "2", "3", "4"))
stage <- data.frame(table(stage[, c("cluster_assignment", "stage")]))

# calculate proportion of all tumours in this cluster
stage <- left_join(stage, clust_size)
stage$prop <- stage$Freq / stage$cluster_size

stage$cluster_assignment <- factor(stage$cluster_assignment, levels = c(1:num_clusters))

plot_stage <- ggplot(stage, aes(cluster_assignment, prop, fill = stage))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ffffcc",
                                        "#a1dab4",
                                        "#41b6c4",
                                        "#225ea8")) +
                                          theme(legend.position = "none") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("% tumours") +
  labs(x=NULL)

stage <- stage[, c("cluster_assignment", "stage", "prop")]
stage <- dcast(stage, cluster_assignment ~ stage)

colnames(stage) <- c("cluster_assignment", "proportion_stage1", "proportion_stage2", "proportion_stage3", "proportion_stage4")

# LOHHLA
LOHHLA <- readRDS(HLA_path)

LOHHLA <- left_join(LOHHLA, sample_table[, c("participant_id", "tumour_sample_platekey")], by = c("BAM_ID" = "tumour_sample_platekey"))

LOHHLA_ct <- LOHHLA[, c("Gene", "HLA.loss", "participant_id")]
LOHHLA_ct <- LOHHLA_ct[!is.na(LOHHLA_ct$participant_id), ]

LOHHLA_ct <- LOHHLA_ct %>% group_by(participant_id) %>% mutate(any_LOHHLA = ifelse(any(HLA.loss), TRUE, FALSE))
LOHHLA_ct <- data.frame(LOHHLA_ct)
LOHHLA_ct <- LOHHLA_ct[, c("participant_id", "any_LOHHLA")]
LOHHLA_ct <- unique(LOHHLA_ct)

LOHHLA_ct <- left_join(cluster_assignment, LOHHLA_ct, by = c("variable" = "participant_id"))
LOHHLA_ct$cluster_assignment <- as.character(LOHHLA_ct$cluster_assignment)
LOHHLA_ct <- data.frame(table(LOHHLA_ct[, c("cluster_assignment", "any_LOHHLA")]))

LOHHLA_ct <- left_join(LOHHLA_ct, clust_size)
LOHHLA_ct$prop <- LOHHLA_ct$Freq / LOHHLA_ct$cluster_size
LOHHLA_ct$any_LOHHLA <- factor(LOHHLA_ct$any_LOHHLA, levels = c("TRUE", "FALSE"))
LOHHLA_ct <- LOHHLA_ct[LOHHLA_ct$any_LOHHLA == TRUE, ]

LOHHLA_ct$cluster_assignment <- factor(LOHHLA_ct$cluster_assignment, levels = c(1:num_clusters))

plot_lohhla <- ggplot(LOHHLA_ct, aes(cluster_assignment, 1, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(name = "% tumours\nwith HLA LOH", low = "#ffffcc", high = "#7fcdbb") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("% tumours\nwith HLA LOH") +
  labs(x=NULL)

LOHHLA_ct <- LOHHLA_ct[, c("cluster_assignment", "any_LOHHLA", "prop")]
LOHHLA_ct <- dcast(LOHHLA_ct, cluster_assignment ~ any_LOHHLA)
colnames(LOHHLA_ct) <- c("cluster_assignment", "proportion_LOHHLA")

# plot the number of drivers per group
driver_count <- read.table(driver_count_path, head = T, sep = "\t")
driver_count$sample <- as.character(driver_count$sample)
driver_count <- left_join(driver_count, cluster_assignment, by = c("sample" = "variable"))
driver_count$cluster_assignment <- as.character(driver_count$cluster_assignment)
driver_count$SNV <- driver_count$coding + driver_count$non_coding

driver_count <- driver_count %>% group_by(cluster_assignment) %>% 
  mutate(mean_coding = mean(coding)) %>%
  mutate(mean_noncoding = mean(non_coding)) %>%
  mutate(mean_SV = mean(SV_driver_count)) %>%
  mutate(mean_CN = mean(CN_driver_count)) %>%
  mutate(mean_total = mean(total_driver_count)) %>%
  mutate(mean_SNV = mean(SNV))

driver_count_mean <- unique(driver_count[, c("cluster_assignment", "mean_coding", "mean_noncoding", "mean_SV", "mean_CN", "mean_total", "mean_SNV")])

driver_count_mean$cluster_assignment <- factor(driver_count_mean$cluster_assignment, levels = c(1:num_clusters))

plot_cn_driver_count <- ggplot(driver_count_mean, aes(cluster_assignment, 1, fill = mean_CN)) +
  geom_tile() +
  scale_fill_gradient(name = "# CN cancer gene aberrations", low = "#ece7f2", high = "#023858") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("# CN cancer gene aberrations") +
  labs(x=NULL)

plot_sv_driver_count <- ggplot(driver_count_mean, aes(cluster_assignment, 1, fill = mean_SV)) +
  geom_tile() +
  scale_fill_gradient(name = "# SV cancer gene aberrations", low = "#ece7f2", high = "#023858") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("# SV cancer gene aberrations") +
  labs(x=NULL)

plot_snv_driver_count <- ggplot(driver_count_mean, aes(cluster_assignment, 1, fill = mean_SNV)) +
  geom_tile() +
  scale_fill_gradient(name = "# SNV cancer gene aberrations", low = "#ece7f2", high = "#023858") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("# SNV cancer gene aberrations") +
  labs(x=NULL)


colnames(driver_count_mean) <- c("cluster_assignment", "mean_coding_driver_mut", "mean_noncoding_driver_mut", "mean_SV_driver", "mean_CN_driver" , "mean_any_driver", "mean_coding_plus_non_coding")


# plot proportion of tumours that have zero drivers
prop_no_drivers <- data.frame(table(driver_count[driver_count$total_driver_count == 0, "cluster_assignment"]))

# cluster 7 and 9 are missing, add manually
add <- data.frame(Var1 = c("7", "9"), Freq = 0)
prop_no_drivers <- rbind(prop_no_drivers, add)

prop_no_drivers <- left_join(prop_no_drivers, clust_size, by = c("Var1" = "cluster_assignment"))
prop_no_drivers$prop <- prop_no_drivers$Freq / prop_no_drivers$cluster_size
prop_no_drivers$prop <- prop_no_drivers$prop * 100
prop_no_drivers$Var1 <- as.character(prop_no_drivers$Var1)

prop_no_drivers$Var1 <- factor(prop_no_drivers$Var1, levels = c(1:num_clusters))

plot_nodriver <- ggplot(prop_no_drivers, aes(Var1, 1, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(name = "% tumours without cancer gene aberration", low = "#fee6ce", high = "#d94801") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("% tumours without cancer gene aberration") +
  labs(x=NULL)

prop_no_drivers <- prop_no_drivers[, c("Var1", "prop")]
colnames(prop_no_drivers) <- c("cluster_assignment", "proportion_no_driver")

# chemo and radiotherapy
treat <- left_join(cluster_assignment, clin[, c("participant_id", "chemo_before_sample", "RT_before_sample")], by = c("variable" = "participant_id"))
treat[is.na(treat$chemo_before_sample), "chemo_before_sample"] <- FALSE
treat[is.na(treat$RT_before_sample), "RT_before_sample"] <- FALSE

# make a column indicating, whether only RT, only chemo, both or neither
treat[treat$chemo_before_sample == FALSE & treat$RT_before_sample == FALSE, "treatment"] <- "no_prior_treatment"
treat[treat$chemo_before_sample == TRUE & treat$RT_before_sample == FALSE, "treatment"] <- "chemo"
treat[treat$chemo_before_sample == FALSE & treat$RT_before_sample == TRUE, "treatment"] <- "RT"
treat[treat$chemo_before_sample == TRUE & treat$RT_before_sample == TRUE, "treatment"] <- "chemo_and_RT"

treat_count <- data.frame(table(treat[, c("cluster_assignment", "treatment")]))
treat_count <- left_join(treat_count, clust_size)
treat_count$prop <- treat_count$Freq / treat_count$cluster_size
treat_count$treatment <- factor(treat_count$treatment, levels = rev(c("chemo", "RT", "chemo_and_RT", "no_prior_treatment")))

treat_count$cluster_assignment <- factor(treat_count$cluster_assignment, levels = c(1:num_clusters))

plot_treatment <- ggplot(treat_count, aes(cluster_assignment, prop, fill = treatment))+
  # geom_quasirandom(alpha = 0.5) + 
  geom_bar(stat = "identity") +
  # scale_y_log10(expand = c(0,0)) +
  # stat_compare_means(label.y = 7, label.x = 2) +
  scale_fill_manual(breaks = rev(c("chemo", "RT", "chemo_and_RT", "no_prior_treatment")),
                    labels = rev(c("chemotherapy", "radiotherapy", "chemo- and\nradiotherapy", "no prior\ntreatment")),
                    values = rev(c("#ffeda0", "#feb24c", "#f03b20", "white"))) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8),
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("% samples\nwith treatment") +
  labs(x=NULL)

treat_count <- treat_count[, c("cluster_assignment", "treatment", "prop")]
treat_count <- dcast(treat_count, cluster_assignment ~ treatment)
colnames(treat_count) <- c("cluster_assignment", "proportion_no_prior_treatment", "proportion_chemo_and_RT", "proportion_RT", "proportion_chemo")


# number of smoking associated mutations

smoking_exposure <- mat_exposure[mat_exposure$label %in% c("SBS4:smoking", "SBS92:smoking", "DBS2:smoking", "ID3:smoking"), ]
smoking_exposure[is.na(smoking_exposure$exposure), "exposure"] <- 0
smoking_exposure <- smoking_exposure %>% group_by(sample) %>% summarize(sum_smoking = sum(exposure))
smoking_exposure <- left_join(smoking_exposure, cluster_assignment, by = c("sample" = "variable"))
smoking_exposure$cluster_assignment <- as.character(smoking_exposure$cluster_assignment)

smoking_exposure <- smoking_exposure %>% group_by(cluster_assignment) %>%  mutate(mean = median(sum_smoking))
smoking_exposure <- unique(smoking_exposure[, c("cluster_assignment", "mean")])
smoking_exposure$cluster_assignment <- factor(smoking_exposure$cluster_assignment, levels = c(1:num_clusters))

plot_smoking_exposure <- ggplot(smoking_exposure, aes(cluster_assignment, 1, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(name = "# smoking mutations", low = "#7ad8fa", high = "#46788a") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"),
        legend.justification = c(0,1),
        panel.spacing = unit(0,'lines'),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  ylab("# smoking mutations") +
  labs(x=NULL)

colnames(smoking_exposure) <- c("cluster_assignment", "median_smoking_exposure")

# sort out the legends
theme_margin <- theme(legend.box.margin = margin(1, 1, 2, 1))

legend_hist  <- cowplot::get_legend(plot_hist + theme_margin)
legend_stage  <- cowplot::get_legend(plot_stage + theme_margin)
legend_timepoint  <- cowplot::get_legend(plot_timepoint + theme_margin)
legend_smok  <- cowplot::get_legend(plot_smok + theme_margin)
legend_smok_exposure  <- cowplot::get_legend(plot_smoking_exposure + theme_margin)
legend_treatment  <- cowplot::get_legend(plot_treatment + theme_margin)
legend_age  <- cowplot::get_legend(plot_age + theme_margin)
legend_SBSburden  <- cowplot::get_legend(plot_SBSburden + theme_margin)
legend_SVburden  <- cowplot::get_legend(plot_SVburden + theme_margin)
legend_SNV_ITH <- cowplot::get_legend(plot_SNV_ITH + theme_margin)
legend_immunet  <- cowplot::get_legend(plot_immunet + theme_margin)
legend_immuneb  <- cowplot::get_legend(plot_immuneb + theme_margin)
legend_lohhla  <- cowplot::get_legend(plot_lohhla + theme_margin)
legend_snv_driver_count  <- cowplot::get_legend(plot_snv_driver_count + theme_margin)
legend_plot_cn_driver_count  <- cowplot::get_legend(plot_cn_driver_count + theme_margin)
legend_plot_sv_driver_count  <- cowplot::get_legend(plot_sv_driver_count + theme_margin)
legend_plot_driver  <- cowplot::get_legend(plot_CNdriver + theme_margin)
legend_plot_WGD  <- cowplot::get_legend(plot_WGD + theme_margin)
legend_plot_nodriver <- cowplot::get_legend(plot_nodriver + theme_margin)


combineLegend1 <- cowplot::plot_grid(legend_hist,
                                     legend_stage,
                                     legend_timepoint,
                                     legend_smok,
                                     legend_treatment,
                                     legend_age,
                                     ncol = 1, rel_heights = c(7, 2.5, 1.7, 1.6, 2.6, 3), align = "v", axis = "rlbt", scale = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5))

combine_Legend2 <- cowplot::plot_grid(legend_SBSburden, legend_SVburden, legend_SNV_ITH, ncol = 3)
combine_Legend3 <- cowplot::plot_grid(legend_immunet, legend_immuneb, legend_lohhla, ncol = 3)
combine_Legend4 <- cowplot::plot_grid(legend_snv_driver_count, legend_plot_cn_driver_count, legend_plot_sv_driver_count, ncol = 3)
combine_Legend5 <- cowplot::plot_grid(legend_plot_nodriver, legend_smok_exposure, legend_plot_WGD, ncol = 3)
combine_Legend6 <- cowplot::plot_grid(legend_plot_driver, ncol = 1)

combine_Legendother <- cowplot::plot_grid(combine_Legend2, combine_Legend3, combine_Legend4, combine_Legend5, combine_Legend6, 
                                          ncol = 1, align = "v", axis = "l", rel_heights = c(1, 1, 1, 1, 1.5), 
                                          scale =c(0.5, 0.5, 0.5, 0.5, 0.5))

all_legend <- cowplot::plot_grid(combineLegend1, combine_Legendother,
                                 ncol = 2)

# remove legends from plots
plot_hist <- plot_hist + theme(legend.position = "none")
plot_stage <- plot_stage + theme(legend.position = "none")
plot_timepoint <- plot_timepoint + theme(legend.position = "none")
plot_smok <- plot_smok + theme(legend.position = "none")
plot_smoking_exposure <- plot_smoking_exposure + theme(legend.position = "none")
plot_treatment <- plot_treatment + theme(legend.position = "none")
plot_age <- plot_age + theme(legend.position = "none")
plot_SBSburden <- plot_SBSburden + theme(legend.position = "none")
plot_SVburden <- plot_SVburden + theme(legend.position = "none")
plot_SNV_ITH <- plot_SNV_ITH + theme(legend.position = "none")
plot_immunet <- plot_immunet + theme(legend.position = "none")
plot_immuneb <- plot_immuneb + theme(legend.position = "none")
plot_lohhla <- plot_lohhla + theme(legend.position = "none")
plot_snv_driver_count <- plot_snv_driver_count + theme(legend.position = "none")
plot_cn_driver_count <- plot_cn_driver_count + theme(legend.position = "none")
plot_sv_driver_count <- plot_sv_driver_count + theme(legend.position = "none")
plot_nodriver <- plot_nodriver + theme(legend.position = "none")
plot_CNdriver <- plot_CNdriver + theme(legend.position = "none")
plot_driver <- plot_driver + theme(legend.position = "none")
plot_SVdriver <- plot_SVdriver + theme(legend.position = "none")
plot_germ_driver <- plot_germ_driver + theme(legend.position = "none")
plot_WGD <- plot_WGD + theme(legend.position = "none")

plot <- plot_grid(plot_size, plot_hist, plot_stage, plot_timepoint, plot_smok,
                  plot_treatment, plot_age,
                  plot_SBSburden, plot_SNV_ITH, plot_SVburden, plot_WGD, plot_smoking_exposure, plot_immunet, plot_immuneb, plot_lohhla,
                  plot_snv_driver_count, plot_cn_driver_count, plot_sv_driver_count, plot_nodriver,
                  plot_driver, plot_CNdriver, plot_SVdriver, plot_germ_driver,
                  rel_heights = c(1.05, 0.6, 0.6, 0.6, 0.6,
                                  0.6, 0.3,
                                  0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
                                  0.3, 0.3, 0.3, 0.3,
                                  1, 1, 0.3, 0.4), ncol = 1, align = "v", axis = "rlbt")

pdf(paste0(output_path, "clustering/ALL_new/signature_cluster_characteristics_heatmap.pdf"), width = 6, height = 9)
plot
dev.off()

pdf(paste0(output_path, "clustering/ALL_new/signature_cluster_characteristics_heatmap_LEGEND.pdf"), width = 12, height = 12)
all_legend
dev.off()


# make all the data going into this plot into one data frame

all_data <- list(gene_test,
                 SV_gene_test,
                 CN_gene_test,
                 germ_gene_test,
                 tcra,
                 tcrb,
                 hist_count,
                 clust_size,
                 time_point,
                 smok_count,
                 SNV_ITH,
                 WGD_df,
                 age_df,
                 sex_df,
                 mutburden_df,
                 stage,
                 LOHHLA_ct,
                 driver_count_mean,
                 prop_no_drivers,
                 treat_count,
                 smoking_exposure)
